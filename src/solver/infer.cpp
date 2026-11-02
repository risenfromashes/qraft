#include "internal.hpp"
#include "divide/bipartition_test.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <condition_variable>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string_view>
#include <thread>
#include <type_traits>

namespace qraft::internal {
namespace {

struct InferenceRunState {
  std::mutex mutex;
  std::condition_variable changed;
  int active_shared_runs = 0;
  int waiting_writers = 0;
  int configured_workers = 0;
  bool exclusive_run = false;
};

InferenceRunState& inference_run_state() {
  static InferenceRunState state;
  return state;
}

int requested_worker_count(const Options& options, int worker_cap) {
  const int requested =
      options.threads > 0
          ? options.threads
          : static_cast<int>(std::thread::hardware_concurrency());
  return std::min(std::max(1, worker_cap), std::max(1, requested));
}

bool requires_exclusive_inference_run(const Options& options) {
  return options.progress != ProgressMode::Off ||
         !options.stats_json.empty();
}

/// Admits one inference run to the shared executor. A run with progress not
/// Off, or with --stats-json, is exclusive, since both use process-wide
/// state. Quiet runs may overlap while they use the same worker count, and
/// the run statistics are disabled while any is active. The executor is
/// reconfigured only when no run is active. A waiting exclusive run, or a
/// quiet run that needs another worker count, stops further quiet runs from
/// joining (waiting_writers).
class ScopedInferenceRun {
 public:
  ScopedInferenceRun(const Options& options, int worker_cap)
      : requested_workers_(requested_worker_count(options, worker_cap)),
        exclusive_(requires_exclusive_inference_run(options)) {
    if (exclusive_) {
      acquire_exclusive(options);
    } else {
      acquire_shared(options);
    }
  }

  ScopedInferenceRun(const ScopedInferenceRun&) = delete;
  ScopedInferenceRun& operator=(const ScopedInferenceRun&) = delete;

  ~ScopedInferenceRun() noexcept {
    InferenceRunState& state = inference_run_state();
    {
      std::lock_guard<std::mutex> lock(state.mutex);
      if (exclusive_) {
        state.exclusive_run = false;
      } else {
        --state.active_shared_runs;
        if (state.active_shared_runs == 0) {
          set_run_stats_enabled(true);
        }
      }
    }
    state.changed.notify_all();
  }

  int workers() const noexcept { return workers_; }

  bool owns_process_telemetry() const noexcept { return exclusive_; }

 private:
  int configure_requested_workers(const Options& options) const {
    Options effective = options;
    effective.threads = requested_workers_;
    return configure_threads(effective);
  }

  void acquire_exclusive(const Options& options) {
    InferenceRunState& state = inference_run_state();
    std::unique_lock<std::mutex> lock(state.mutex);
    ++state.waiting_writers;
    state.changed.wait(lock, [&] {
      return !state.exclusive_run && state.active_shared_runs == 0;
    });
    --state.waiting_writers;
    state.exclusive_run = true;
    set_run_stats_enabled(true);
    try {
      workers_ = configure_requested_workers(options);
      state.configured_workers = workers_;
    } catch (...) {
      state.exclusive_run = false;
      lock.unlock();
      state.changed.notify_all();
      throw;
    }
  }

  void acquire_shared(const Options& options) {
    InferenceRunState& state = inference_run_state();
    std::unique_lock<std::mutex> lock(state.mutex);
    bool waiting_for_reconfiguration = false;
    for (;;) {
      if (!state.exclusive_run && state.active_shared_runs > 0 &&
          state.configured_workers == requested_workers_ &&
          state.waiting_writers == 0) {
        ++state.active_shared_runs;
        workers_ = state.configured_workers;
        return;
      }

      if (!state.exclusive_run && state.active_shared_runs == 0 &&
          (state.waiting_writers == 0 || waiting_for_reconfiguration)) {
        if (waiting_for_reconfiguration) {
          --state.waiting_writers;
          waiting_for_reconfiguration = false;
        }

        if (state.configured_workers == requested_workers_ &&
            get_num_workers() == requested_workers_) {
          set_run_stats_enabled(false);
          state.active_shared_runs = 1;
          workers_ = state.configured_workers;
          lock.unlock();
          state.changed.notify_all();
          return;
        }

        state.exclusive_run = true;
        try {
          workers_ = configure_requested_workers(options);
          state.configured_workers = workers_;
        } catch (...) {
          state.exclusive_run = false;
          lock.unlock();
          state.changed.notify_all();
          throw;
        }
        state.exclusive_run = false;
        set_run_stats_enabled(false);
        state.active_shared_runs = 1;
        lock.unlock();
        state.changed.notify_all();
        return;
      }

      if (!waiting_for_reconfiguration &&
          state.configured_workers != requested_workers_) {
        ++state.waiting_writers;
        waiting_for_reconfiguration = true;
      }
      state.changed.wait(lock);
    }
  }

  int requested_workers_ = 1;
  int workers_ = 1;
  bool exclusive_ = false;
};

/// When progress is on, finishes progress output on every exit. Declared
/// after the ScopedInferenceRun, it stops the progress heartbeat, which reads
/// process-wide state, before an exception releases the exclusive run.
class ScopedInferenceProgressFinish {
 public:
  explicit ScopedInferenceProgressFinish(bool active) : active_(active) {}
  ScopedInferenceProgressFinish(const ScopedInferenceProgressFinish&) = delete;
  ScopedInferenceProgressFinish& operator=(
      const ScopedInferenceProgressFinish&) = delete;

  ~ScopedInferenceProgressFinish() noexcept {
    if (!active_) return;
    try {
      progress_finish();
    } catch (...) {
    }
  }

  void finish() {
    if (!active_) return;
    progress_finish();
    active_ = false;
  }

 private:
  bool active_ = false;
};

const char* divide_score_name(DivideScore score) {
  return score == DivideScore::Quartet ? "quartet" : "balanced";
}

const char* regraft_score_name(RegraftScore score) {
  return score == RegraftScore::Normalized ? "normalized" : "raw";
}

const char* flag(bool value) { return value ? "true" : "false"; }

// Sorted distinct labels of every gene, as TaxonRegistry(labels) orders them,
// without copying each physical label: each gene's distinct labels are sorted
// in parallel, then merged pairwise in parallel rounds.
std::vector<std::string> sorted_unique_gene_labels(
    const std::vector<GeneTree>& genes,
    long long& physical_labels) {
  std::vector<std::vector<std::string_view>> lists(genes.size());
  std::vector<long long> counts(genes.size(), 0);
  parallel_for_inner(0, static_cast<int>(genes.size()), [&](int g) {
    std::vector<std::string_view>& list = lists[static_cast<std::size_t>(g)];
    for (const Tree::Node& node : genes[static_cast<std::size_t>(g)].nodes) {
      if (!node.label.empty()) list.push_back(node.label);
    }
    counts[static_cast<std::size_t>(g)] = static_cast<long long>(list.size());
    std::sort(list.begin(), list.end());
    list.erase(std::unique(list.begin(), list.end()), list.end());
  }, 1);
  physical_labels = std::accumulate(counts.begin(), counts.end(), 0LL);
  while (lists.size() > 1) {
    std::vector<std::vector<std::string_view>> merged((lists.size() + 1) / 2);
    parallel_for_inner(0, static_cast<int>(merged.size()), [&](int i) {
      const std::size_t left = 2u * static_cast<std::size_t>(i);
      std::vector<std::string_view>& out = merged[static_cast<std::size_t>(i)];
      if (left + 1u == lists.size()) {
        out = std::move(lists[left]);
        return;
      }
      out.reserve(std::max(lists[left].size(), lists[left + 1u].size()));
      std::set_union(lists[left].begin(), lists[left].end(),
                     lists[left + 1u].begin(), lists[left + 1u].end(),
                     std::back_inserter(out));
    }, 1);
    lists = std::move(merged);
  }
  std::vector<std::string> unique;
  if (!lists.empty()) {
    unique.reserve(lists.front().size());
    for (std::string_view label : lists.front()) unique.emplace_back(label);
  }
  return unique;
}

/// The --support setting as written in the statistics: none, localpp-gene or
/// localpp-quartet.
std::string support_name(const Options& options) {
  if (options.support == SupportMode::None) return "none";
  return options.support_frequencies == SupportFrequencies::Gene ? "localpp-gene"
                                                                  : "localpp-quartet";
}

std::string startup_settings_summary(const Options& options) {
  std::ostringstream out;
  out << "settings divide=" << flag(options.divide)
      << " divide_score=" << divide_score_name(options.divide_score)
      << " bipartition_test=" << flag(options.bipartition_test.enabled)
      << " regraft_score=" << regraft_score_name(options.regraft_score)
      << " nni=" << flag(options.nni)
      << " nni_move_cap=" << flag(options.nni_move_cap)
      << " nni_support_gate=" << flag(options.nni_support_gate)
      << " nni_batch=" << flag(options.nni_batch)
      << " spr_refine=" << flag(options.spr_refine)
      << " guide_tree=" << flag(!options.guide_tree_path.empty())
      << " guide_seed=" << flag(options.guide_seed)
      << " support=" << support_name(options);
  return out.str();
}

/// The public species tree of a solved tree: unlabelled vertices of degree
/// two suppressed, active vertices renumbered in order, leaves labelled.
/// With `support_ctx`, its internal edges carry the support it asks for.
SpeciesTree to_species(const SolverTree& solved, const TaxonRegistry& taxa,
                       const SolverContext* support_ctx = nullptr) {
  SolverTree compact = solved;
  suppress_solver_degree_two(compact);
  SpeciesTree species;
  std::vector<int> map(compact.nodes.size(), -1);
  for (int i = 0; i < static_cast<int>(compact.nodes.size()); ++i) {
    if (!compact.nodes[i].active) continue;
    const int taxon = compact.nodes[i].taxon;
    if (taxon >= taxa.size()) continue;
    map[i] = static_cast<int>(species.nodes.size());
    SpeciesTree::Node node;
    if (taxon >= 0) node.label = taxa.label_for(taxon);
    species.nodes.push_back(std::move(node));
  }
  for (int i = 0; i < static_cast<int>(compact.nodes.size()); ++i) {
    if (map[i] < 0) continue;
    for (int nb : compact.nodes[i].neighbors) {
      if (map[nb] >= 0) species.nodes[map[i]].neighbors.push_back(map[nb]);
    }
  }
  if (support_ctx != nullptr &&
      support_ctx->options.support == SupportMode::LocalPP) {
    annotate_branch_support(*support_ctx, compact, map, species);
  }
  return species;
}

/// Writes the run's options and counters as one flat JSON object.
void write_stats_json(const Options& options, const SolverContext& ctx) {
  if (options.stats_json.empty()) return;
  std::ofstream out(options.stats_json);
  if (!out) {
    throw std::runtime_error("failed to open stats file: " +
                             options.stats_json);
  }
  const BipartitionTestOptions& test = options.bipartition_test;
  const std::shared_ptr<const BipartitionTestIndex>& index =
      ctx.bipartition_test_index;
  out << "{\n";
  if (!options.program_version.empty()) {
    out << "  \"qraft_version\": \"" << options.program_version << "\",\n";
  }
  out << "  \"genes\": " << base_genes(ctx).size() << ",\n";
  out << "  \"taxa\": " << context_taxa(ctx).size() << ",\n";
  out << "  \"workers\": " << ctx.workers << ",\n";
  out << "  \"progress\": \"" << progress_mode_name(options.progress)
      << "\",\n";
  out << "  \"divide\": " << flag(options.divide) << ",\n";
  out << "  \"divide_score\": \"" << divide_score_name(options.divide_score)
      << "\",\n";
  out << "  \"bipartition_test\": " << flag(test.enabled) << ",\n";
  out << "  \"bipartition_test_available\": " << flag(index != nullptr)
      << ",\n";
  out << "  \"bipartition_test_complete_genes\": "
      << flag(index != nullptr && index->complete_single_copy) << ",\n";
  out << "  \"bipartition_test_coefficient\": " << test.coefficient << ",\n";
  out << "  \"bipartition_test_crossover\": " << test.crossover << ",\n";
  out << "  \"bipartition_test_proposals\": " << test.proposals << ",\n";
  out << "  \"bipartition_test_anchor_pairs\": " << test.anchor_pairs
      << ",\n";
  // On genes that miss taxa every pair of taxa has its own threshold.
  out << "  \"bipartition_test_threshold\": ";
  if (index != nullptr && index->complete_single_copy) {
    out << bipartition_test_threshold(static_cast<int>(index->genes.size()),
                                      index->taxon_count, test.coefficient,
                                      test.crossover);
  } else {
    out << "null";
  }
  out << ",\n";
  out << "  \"regraft_score\": \""
      << regraft_score_name(options.regraft_score) << "\",\n";
  out << "  \"nni\": " << flag(options.nni) << ",\n";
  out << "  \"nni_move_cap\": " << flag(options.nni_move_cap) << ",\n";
  out << "  \"nni_support_gate\": " << flag(options.nni_support_gate)
      << ",\n";
  out << "  \"nni_batch\": " << flag(options.nni_batch) << ",\n";
  out << "  \"spr_refine\": " << flag(options.spr_refine) << ",\n";
  out << "  \"spr_reads_per_pass\": " << options.spr_reads << ",\n";
  out << "  \"spr_min_clade\": " << options.spr_min_clade << ",\n";
  out << "  \"spr_candidates\": " << options.spr_candidates << ",\n";
  out << "  \"spr_max_rounds\": " << options.spr_rounds << ",\n";
  out << "  \"guide_tree\": " << flag(!options.guide_tree_path.empty()) << ",\n";
  out << "  \"guide_seed\": " << flag(options.guide_seed) << ",\n";
  out << "  \"support\": \"" << support_name(options) << "\",\n";
  const RunStats& stats = run_stats();
#define QRAFT_WRITE_STAT(name)                                   \
  out << "  \"" #name "\": "                                     \
      << stats.name.load(std::memory_order_relaxed) << ",\n";
  QRAFT_RUN_STATS(QRAFT_WRITE_STAT)
#undef QRAFT_WRITE_STAT
  out << "  \"nni_gain\": " << wide_decimal(stats.nni_gain) << ",\n";
  const bool negative_spr_gain = stats.spr_gain < 0;
  // Convert the magnitude without negating the smallest signed integer.
  const RegraftCount magnitude_base = negative_spr_gain
      ? -(stats.spr_gain + 1) : stats.spr_gain;
  WideCount magnitude{static_cast<std::uint64_t>(magnitude_base >> 64),
                      static_cast<std::uint64_t>(magnitude_base)};
  if (negative_spr_gain && ++magnitude.low == 0) ++magnitude.high;
  out << "  \"spr_gain\": " << (negative_spr_gain ? "-" : "")
      << wide_decimal(magnitude) << ",\n";
  out << "  \"spr_gain_exact\": " << flag(stats.spr_gain_exact) << ",\n";
  out << "  \"complexity_contract_satisfied\": "
      << flag(stats.complexity_contract_satisfied.load(
             std::memory_order_relaxed) != 0)
      << "\n}\n";
}

#undef QRAFT_ALWAYS_INLINE

}  // namespace

}  // namespace qraft::internal

namespace qraft {

/// Runs the QRAFT pipeline: indexes the genes, builds or reads the starting
/// tree, then runs the recursive divide and regraft and NNI refinement, or
/// under --no-divide refines the starting tree directly.
void validate_options(const Options& options) {
  if (!std::isfinite(options.bipartition_test.coefficient) ||
      options.bipartition_test.coefficient <= 0.0 ||
      !std::isfinite(options.bipartition_test.crossover) ||
      options.bipartition_test.crossover < 0.0 ||
      options.bipartition_test.proposals < 1 ||
      options.bipartition_test.proposals > 16 ||
      options.bipartition_test.anchor_pairs < 1 ||
      options.bipartition_test.anchor_pairs > 4) {
    throw std::invalid_argument(
        "bipartition test requires a finite positive coefficient, a finite "
        "nonnegative crossover, 1..16 proposals and 1..4 anchor pairs");
  }
  if (options.spr_reads < 1 || options.spr_min_clade < 2 ||
      options.spr_candidates < 1 || options.spr_rounds < 1) {
    throw std::invalid_argument(
        "SPR refinement requires positive reads, candidates and rounds and "
        "a minimum clade of at least 2 taxa");
  }
}

SpeciesTree infer_species_tree(std::vector<GeneTree> genes, const Options& options) {
  if (genes.empty()) {
    throw std::invalid_argument("at least one gene tree is required");
  }
  if (genes.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("gene count exceeds the integer-indexed domain");
  }
  validate_options(options);
  // The worker count is capped by the input's physical nodes, so executor
  // startup costs no more than indexing the input.
  std::size_t physical_nodes = 0;
  for (const GeneTree& gene : genes) {
    physical_nodes = std::min(
        static_cast<std::size_t>(std::numeric_limits<int>::max()),
        physical_nodes + std::min(
                             gene.nodes.size(),
                             static_cast<std::size_t>(
                                 std::numeric_limits<int>::max()) -
                                 physical_nodes));
  }
  const int worker_cap = std::max(1, static_cast<int>(physical_nodes));
  internal::ScopedInferenceRun inference_run(options, worker_cap);
  internal::ScopedInferenceProgressFinish progress_finish(
      options.progress == ProgressMode::On);
  const auto startup_started = std::chrono::steady_clock::now();
  if (inference_run.owns_process_telemetry()) {
    internal::reset_run_stats();
  }
  internal::SolverContext ctx;
  ctx.options = options;
  const Options& run_options = ctx.options;
  ctx.workers = inference_run.workers();
  const std::string startup_settings =
      internal::startup_settings_summary(run_options);
  internal::progress_startup_status(
      run_options, "taxa", startup_settings + " building registry");
  {
    const auto stage_started = std::chrono::steady_clock::now();
    long long physical_labels = 0;
    ctx.taxa = TaxonRegistry::from_ordered_unique(
        internal::sorted_unique_gene_labels(genes, physical_labels));
    internal::run_stats().taxon_registry_physical_labels.store(
        physical_labels, std::memory_order_relaxed);
    internal::run_stats().taxon_registry_unique_labels.store(
        ctx.taxa.size(), std::memory_order_relaxed);
    internal::run_stats().taxon_registry_sort_items.store(
        ctx.taxa.size(), std::memory_order_relaxed);
    internal::stats_add_time(
        internal::run_stats().taxon_registry_ns,
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - stage_started));
  }
  internal::progress_startup_status(
      run_options, "index", startup_settings + " indexing gene trees");
  ctx.genes.resize(genes.size());
  const auto index_started = std::chrono::steady_clock::now();
  internal::parallel_for_inner(
      0, static_cast<int>(genes.size()),
      [&](int i) {
        ctx.genes[static_cast<std::size_t>(i)] =
            internal::index_gene(genes[static_cast<std::size_t>(i)], ctx.taxa);
        // The parsed tree is not read again; release it as soon as it is
        // indexed rather than after every gene is.
        genes[static_cast<std::size_t>(i)] = GeneTree{};
      },
      1);
  std::vector<GeneTree>().swap(genes);
  internal::release_retained_allocator_pages();
  internal::stats_add_time(
      internal::run_stats().gene_index_ns,
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now() - index_started));

  internal::progress_startup_status(
      run_options, "starting-tree",
      startup_settings + " building MRE");
  if (run_options.guide_seed &&
      (run_options.guide_tree_path.empty() ||
       !run_options.starting_tree_path.empty())) {
    throw std::invalid_argument(
        "guide_seed needs a guide tree and builds the MRE, so not with a "
        "starting tree file");
  }
  std::optional<Tree> guide_tree;
  if (!run_options.guide_tree_path.empty()) {
    std::ifstream guide_in(run_options.guide_tree_path);
    if (!guide_in) {
      throw std::runtime_error("failed to open guide tree file: " +
                               run_options.guide_tree_path);
    }
    const std::string guide_text((std::istreambuf_iterator<char>(guide_in)),
                                 std::istreambuf_iterator<char>());
    guide_tree = parse_newick(guide_text);
  }
  if (run_options.starting_tree_path.empty()) {
    // --guide-seed: the MRE skips the clusters that conflict with the guide.
    std::unique_ptr<internal::GuideFilter> guide_filter;
    if (run_options.guide_seed && guide_tree) {
      guide_filter = std::make_unique<internal::GuideFilter>(*guide_tree, ctx.taxa);
    }
    ctx.starting_tree = internal::build_mre_tree(
        ctx.genes, ctx.taxa.size(), run_options, guide_filter.get());
  } else {
    std::ifstream starting_tree_in(run_options.starting_tree_path);
    if (!starting_tree_in) {
      throw std::runtime_error("failed to open starting tree file: " +
                               run_options.starting_tree_path);
    }
    const std::string starting_tree_text(
        (std::istreambuf_iterator<char>(starting_tree_in)),
        std::istreambuf_iterator<char>());
    ctx.starting_tree =
        internal::starting_tree_from(parse_newick(starting_tree_text), ctx.taxa);
  }
  if (guide_tree) {
    ctx.starting_tree =
        internal::apply_guide_tree(ctx.starting_tree, *guide_tree, ctx.taxa);
  }
  if (!run_options.write_starting_tree_path.empty()) {
    std::ofstream starting_tree_out(run_options.write_starting_tree_path);
    starting_tree_out << internal::starting_tree_newick(ctx.starting_tree, ctx.taxa);
    if (!starting_tree_out) {
      throw std::runtime_error("failed to write starting tree file: " +
                               run_options.write_starting_tree_path);
    }
  }
  // Built after MRE so its arrays are not live at the MRE peak; it depends
  // only on the indexed genes. The test needs every gene to hold each of its
  // taxa once; genes may miss taxa. Otherwise it is off for the run.
  if (run_options.divide && run_options.bipartition_test.enabled) {
    ctx.bipartition_test_index = internal::make_bipartition_test_index(
        ctx.genes, ctx.taxa.size());
    if (!ctx.bipartition_test_index->single_copy) {
      ctx.bipartition_test_index.reset();
    }
  }
  internal::stats_add_time(
      internal::run_stats().startup_ns,
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now() - startup_started));
  internal::release_retained_allocator_pages();
  internal::progress_start(run_options, std::max(0, ctx.taxa.size() - 3));

  internal::SolverTree solved;
  tf::Taskflow solve_flow;
  solve_flow.emplace([&](tf::Subflow& sf) {
    if (run_options.divide) {
      solved = internal::solve_qraft(ctx, &sf);
    } else {
      solved = internal::refine_starting_tree(ctx, run_options.nni);
    }
  });
  internal::get_executor().run(solve_flow).get();
  progress_finish.finish();
  // Support is counted on the final tree, before the statistics are written.
  SpeciesTree species =
      internal::to_species(solved, internal::context_taxa(ctx), &ctx);
  internal::write_stats_json(run_options, ctx);
  return species;
}

}  // namespace qraft
