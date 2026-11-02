#include "solver/solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>

namespace qraft::internal {

void rethrow_parallel_child_failure(
    const std::exception_ptr& first,
    const std::exception_ptr& second) {
  if (first) std::rethrow_exception(first);
  if (second) std::rethrow_exception(second);
}

namespace {

constexpr int kParallelChildMinTaxa = 64;
// Below this many scoring genes per worker, the last genes of a gene loop
// leave most workers idle, so siblings run in parallel to fill them.
constexpr int kParallelChildGenesPerWorker = 16;
// Progress output names a subproblem by a hash of its path from the root,
// seeded with this value.
constexpr uint64_t kRootPath = 0x6a09e667f3bcc909ULL;

// Payload accounting. When a run records statistics, every structure the
// recursion keeps alive across a descent holds a lease on one live_* counter,
// so peak_payload bounds the recursion's live space; the recursion tests
// check these peaks against closed forms.

/// Holds `amount` units of one live/peak payload category.
class PayloadLease {
 public:
  PayloadLease() = default;
  PayloadLease(const PayloadLease&) = delete;
  PayloadLease& operator=(const PayloadLease&) = delete;
  PayloadLease(PayloadLease&& other) noexcept { move_from(other); }
  PayloadLease& operator=(PayloadLease&& other) noexcept {
    if (this != &other) {
      release();
      move_from(other);
    }
    return *this;
  }
  ~PayloadLease() { release(); }

  void acquire(bool enabled,
               std::atomic<long long>& live,
               std::atomic<long long>& peak,
               long long amount) {
    release();
    if (!enabled || amount <= 0) return;
    live_ = &live;
    amount_ = amount;
    stats_acquire_payload(live, peak, amount);
  }

  void release() {
    if (live_ == nullptr) return;
    stats_release_payload(*live_, amount_);
    live_ = nullptr;
    amount_ = 0;
  }

  long long amount() const { return amount_; }

 private:
  void move_from(PayloadLease& other) noexcept {
    live_ = std::exchange(other.live_, nullptr);
    amount_ = std::exchange(other.amount_, 0);
  }

  std::atomic<long long>* live_ = nullptr;
  long long amount_ = 0;
};

/// The lease on one subproblem's taxon list. While a sibling is solved first,
/// the list is also counted as a pending continuation.
class TaxaLease {
 public:
  TaxaLease() = default;
  TaxaLease(bool enabled, long long taxa) {
    RunStats& stats = run_stats();
    taxa_.acquire(enabled, stats.live_taxa, stats.peak_taxa, taxa);
  }
  TaxaLease(const TaxaLease&) = delete;
  TaxaLease& operator=(const TaxaLease&) = delete;
  TaxaLease(TaxaLease&& other) noexcept
      : taxa_(std::move(other.taxa_)),
        continuation_(std::exchange(other.continuation_, false)) {}
  TaxaLease& operator=(TaxaLease&& other) noexcept {
    if (this != &other) {
      release();
      taxa_ = std::move(other.taxa_);
      continuation_ = std::exchange(other.continuation_, false);
    }
    return *this;
  }
  ~TaxaLease() { release(); }

  void mark_continuation() {
    if (continuation_ || taxa_.amount() == 0) return;
    RunStats& stats = run_stats();
    stats_acquire_gauge(stats.live_continuation_taxa,
                        stats.peak_continuation_taxa, taxa_.amount());
    continuation_ = true;
  }

  void unmark_continuation() {
    if (!continuation_) return;
    stats_release_gauge(run_stats().live_continuation_taxa, taxa_.amount());
    continuation_ = false;
  }

  void release() {
    unmark_continuation();
    taxa_.release();
  }

 private:
  PayloadLease taxa_;
  bool continuation_ = false;
};

struct TreePayload {
  long long nodes = 0;
  long long adjacencies = 0;
};

TreePayload tree_payload(const SolverTree& tree, bool accounting) {
  TreePayload payload;
  if (!accounting) return payload;
  payload.nodes = static_cast<long long>(tree.nodes.size());
  for (const SolverTree::Node& node : tree.nodes) {
    payload.adjacencies += static_cast<long long>(node.neighbors.size());
  }
  return payload;
}

/// The payload connect_at_edges adds to a copy of `copied_tree`.
TreePayload connect_payload(const SolverTree& copied_tree,
                            const RegraftEdge& left_edge,
                            const RegraftEdge& right_edge,
                            bool accounting) {
  TreePayload payload = tree_payload(copied_tree, accounting);
  if (!accounting) return payload;
  const auto note_endpoint = [&](const RegraftEdge& edge) {
    if (edge.a < 0) {
      // An empty tree gains one isolated connector.
      ++payload.nodes;
    } else if (edge.b >= 0) {
      // Subdividing an edge adds a node and replaces one edge by two.
      ++payload.nodes;
      payload.adjacencies += 2;
    }
    // A singleton reuses its only vertex.
  };
  note_endpoint(left_edge);
  note_endpoint(right_edge);
  // The edge joining the two connectors.
  payload.adjacencies += 2;
  return payload;
}

class TreeLease {
 public:
  TreeLease() = default;
  TreeLease(bool enabled, TreePayload payload) { acquire(enabled, payload); }

  void acquire(bool enabled, TreePayload payload) {
    RunStats& stats = run_stats();
    nodes_.acquire(enabled, stats.live_tree_nodes, stats.peak_tree_nodes,
                   payload.nodes);
    adjacencies_.acquire(enabled, stats.live_tree_adjacencies,
                         stats.peak_tree_adjacencies, payload.adjacencies);
  }

  void release() {
    adjacencies_.release();
    nodes_.release();
  }

 private:
  PayloadLease nodes_;
  PayloadLease adjacencies_;
};

struct GeneSetPayload {
  long long nodes = 0;
  // Per-gene headers plus every stored array other than the nodes.
  long long entries = 0;
};

GeneSetPayload gene_set_payload(const std::vector<GeneIndex>& genes,
                                bool accounting) {
  GeneSetPayload payload;
  if (!accounting) return payload;
  payload.entries = static_cast<long long>(genes.size());
  for (const GeneIndex& gene : genes) {
    payload.nodes += static_cast<long long>(gene.nodes.size());
    payload.entries += static_cast<long long>(
        gene.taxon_to_node.size() + gene.root_traversal.parent.size() +
        gene.root_traversal.preorder.size() + gene.adjacency_offsets.size() +
        gene.adjacency.size() + gene.leaf_mass.size());
  }
  return payload;
}

class GeneSetLease {
 public:
  GeneSetLease() = default;
  GeneSetLease(bool enabled, GeneSetPayload payload) {
    RunStats& stats = run_stats();
    nodes_.acquire(enabled, stats.live_gene_set_nodes,
                   stats.peak_gene_set_nodes, payload.nodes);
    entries_.acquire(enabled, stats.live_gene_set_entries,
                     stats.peak_gene_set_entries, payload.entries);
  }

  void release() {
    entries_.release();
    nodes_.release();
  }

 private:
  PayloadLease nodes_;
  PayloadLease entries_;
};

struct RecursionControl {
  bool accounting = false;
  bool force_shape = false;
  ForcedSplitShape shape = ForcedSplitShape::Balanced;
  // The recursion tests make every subproblem this deep throw; -1 for none.
  int fail_depth = -1;
};

/// Sorts taxa in four fixed byte passes.
void radix_sort_taxa(std::vector<int>& taxa, bool accounting) {
  if (taxa.size() <= 1) return;
  RunStats& stats = run_stats();
  PayloadLease scratch_lease;
  scratch_lease.acquire(accounting, stats.live_taxon_scratch,
                        stats.peak_taxon_scratch,
                        static_cast<long long>(taxa.size()) * 2);
  std::vector<int> scratch(taxa.size());
  std::vector<int>* source = &taxa;
  std::vector<int>* destination = &scratch;
  for (int shift = 0; shift < 32; shift += 8) {
    std::array<std::size_t, 256> counts{};
    for (int taxon : *source) {
      ++counts[(static_cast<unsigned int>(taxon) >> shift) & 0xffU];
    }
    std::size_t offset = 0;
    for (std::size_t& count : counts) {
      const std::size_t bucket_size = count;
      count = offset;
      offset += bucket_size;
    }
    for (int taxon : *source) {
      (*destination)[counts[(static_cast<unsigned int>(taxon) >> shift) &
                            0xffU]++] = taxon;
    }
    std::swap(source, destination);
  }
  if (source != &taxa) taxa.swap(*source);
  if (accounting) {
    stats_add(stats.recursion_radix_entries,
              static_cast<long long>(taxa.size()) * 8);
    stats_add(stats.recursion_radix_bucket_ops, 4 * 256);
  }
}

struct TrackedTaxa {
  std::vector<int> values;
  PayloadLease lease;
};

/// The ascending taxa of a solved subtree. Recomputing them after the descent
/// avoids keeping each child's taxon list alive while it is solved.
TrackedTaxa taxa_of_tree(const SolverTree& tree, bool accounting) {
  TrackedTaxa tracked;
  std::vector<int>& taxa = tracked.values;
  for (const SolverTree::Node& node : tree.nodes) {
    if (node.active && node.taxon >= 0) taxa.push_back(node.taxon);
  }
  radix_sort_taxa(taxa, accounting);
  taxa.erase(std::unique(taxa.begin(), taxa.end()), taxa.end());
  if (accounting) {
    RunStats& stats = run_stats();
    stats_add(stats.recursion_tree_node_visits,
              static_cast<long long>(tree.nodes.size()));
    stats_add(stats.recursion_unique_taxa,
              static_cast<long long>(taxa.size()));
    tracked.lease.acquire(true, stats.live_taxon_scratch,
                          stats.peak_taxon_scratch,
                          static_cast<long long>(taxa.size()));
  }
  return tracked;
}

uint64_t child_path(uint64_t path, int child_offset, int depth) {
  return stable_hash_u64(path, static_cast<uint64_t>(child_offset),
                         static_cast<uint64_t>(depth));
}

/// Number of taxa the scoring genes of `ctx` are restricted to.
int source_taxon_count(const SolverContext& ctx) {
  return ctx.local_taxa != nullptr
             ? static_cast<int>(ctx.local_taxa->size())
             : context_taxa(ctx).size();
}

/// Records the taxa of the scoring genes. More than twice the subproblem's
/// taxa, above the base case, counts as an inflation violation.
void note_scoring_source(const SolverContext& ctx, int active_taxa) {
  const int source_taxa = std::max(1, source_taxon_count(ctx));
  RunStats& stats = run_stats();
  stats_max(stats.max_scoring_source_taxa, source_taxa);
  if (active_taxa > 3 && static_cast<long long>(source_taxa) >
                             2LL * static_cast<long long>(active_taxa)) {
    stats_add(stats.scoring_source_inflation_violations);
  }
}

/// A shared copy of `taxa` whose entries count as gene-set payload.
std::shared_ptr<const std::vector<int>> tracked_local_taxa(
    const std::vector<int>& taxa, bool accounting) {
  const long long entries = static_cast<long long>(taxa.size());
  if (accounting) {
    RunStats& stats = run_stats();
    stats_acquire_payload(stats.live_gene_set_taxa, stats.peak_gene_set_taxa,
                          entries);
  }
  return std::shared_ptr<const std::vector<int>>(
      new std::vector<int>(taxa), [accounting, entries](const std::vector<int>* ptr) {
        delete ptr;
        if (accounting) {
          stats_release_payload(run_stats().live_gene_set_taxa, entries);
        }
      });
}

SolverContext make_child_context_base(const SolverContext& parent) {
  SolverContext child;
  child.taxa_ref = &context_taxa(parent);
  child.genes_ref = &base_genes(parent);
  child.bipartition_test_index = parent.bipartition_test_index;
  child.starting_tree_ref = &context_starting_tree(parent);
  child.options = parent.options;
  child.workers = parent.workers;
  child.child_budget = parent.child_budget;
  return child;
}

/// A child context scoring `gene_set` when it is non-null, and otherwise the
/// parent's scoring genes.
SolverContext make_child_context(
    const SolverContext& parent,
    std::shared_ptr<const std::vector<GeneIndex>> gene_set,
    const std::vector<int>& child_taxa,
    bool accounting) {
  SolverContext child = make_child_context_base(parent);
  if (gene_set != nullptr) {
    child.owned_active_genes = std::move(gene_set);
    child.local_taxa = tracked_local_taxa(child_taxa, accounting);
    return child;
  }
  child.owned_active_genes = parent.owned_active_genes;
  child.active_genes = parent.active_genes;
  if (child.owned_active_genes == nullptr && child.active_genes == nullptr) {
    child.active_genes = &base_genes(parent);
  }
  child.local_taxa = parent.local_taxa;
  stats_add(run_stats().child_gene_set_reuses);
  return child;
}

long long estimate_gene_set_bytes(const std::vector<GeneIndex>& genes) {
  long long bytes = static_cast<long long>(genes.capacity()) *
                    static_cast<long long>(sizeof(GeneIndex));
  auto vector_bytes = [](const auto& values) -> long long {
    return static_cast<long long>(values.capacity()) *
           static_cast<long long>(
               sizeof(typename std::decay_t<decltype(values)>::value_type));
  };
  for (const GeneIndex& gene : genes) {
    bytes += vector_bytes(gene.nodes) + vector_bytes(gene.taxon_to_node) +
             vector_bytes(gene.root_traversal.parent) +
             vector_bytes(gene.root_traversal.preorder) +
             vector_bytes(gene.adjacency_offsets) +
             vector_bytes(gene.adjacency) + vector_bytes(gene.leaf_mass);
  }
  return bytes;
}

long long gene_set_node_count(const std::vector<GeneIndex>& genes) {
  long long nodes = 0;
  for (const GeneIndex& gene : genes) {
    nodes += static_cast<long long>(gene.nodes.size());
  }
  return nodes;
}

/// Publishes a restricted gene set, holding `taxa` slots of the run's child
/// budget until the last context using it is destroyed.
std::shared_ptr<const std::vector<GeneIndex>> track_child_gene_set(
    std::unique_ptr<std::vector<GeneIndex>> genes,
    int taxa,
    bool accounting,
    std::shared_ptr<ChildSetBudget> budget) {
  RunStats& stats = run_stats();
  budget->live_sets.fetch_add(1, std::memory_order_relaxed);
  budget->live_slots.fetch_add(taxa, std::memory_order_acq_rel);
  const GeneSetPayload payload = gene_set_payload(*genes, accounting);
  if (accounting) {
    stats_acquire_payload(stats.live_gene_set_nodes, stats.peak_gene_set_nodes,
                          payload.nodes);
    stats_acquire_payload(stats.live_gene_set_entries,
                          stats.peak_gene_set_entries, payload.entries);
  }
  stats_max(stats.max_child_gene_sets,
            stats.live_child_gene_sets.fetch_add(1, std::memory_order_relaxed) +
                1);
  stats_max(stats.max_child_gene_set_slots,
            stats.live_child_gene_set_slots.fetch_add(
                taxa, std::memory_order_relaxed) +
                taxa);
  if (run_stats_enabled()) {
    stats_max(stats.max_child_gene_set_bytes, estimate_gene_set_bytes(*genes));
  }
  return std::shared_ptr<const std::vector<GeneIndex>>(
      genes.release(), [taxa, accounting, payload, budget = std::move(budget)](
                           const std::vector<GeneIndex>* ptr) {
        delete ptr;
        budget->live_sets.fetch_sub(1, std::memory_order_relaxed);
        budget->live_slots.fetch_sub(taxa, std::memory_order_acq_rel);
        RunStats& live_stats = run_stats();
        if (accounting) {
          stats_release_payload(live_stats.live_gene_set_entries,
                                payload.entries);
          stats_release_payload(live_stats.live_gene_set_nodes,
                                payload.nodes);
        }
        live_stats.live_child_gene_sets.fetch_sub(1,
                                                  std::memory_order_relaxed);
        live_stats.live_child_gene_set_slots.fetch_sub(
            taxa, std::memory_order_relaxed);
      });
}

bool should_restrict_genes(const std::vector<GeneIndex>& source,
                           const std::vector<int>& child_taxa) {
  return child_taxa.size() > 3 && !source.empty() &&
         gene_set_node_count(source) > 0;
}

long long projected_child_set_slots(const SolverContext& ctx,
                                    const std::vector<int>& child_taxa) {
  return should_restrict_genes(scoring_genes(ctx), child_taxa)
             ? static_cast<long long>(child_taxa.size())
             : 0;
}

/// Whether the scoring genes hold at most 8 nodes per gene and per taxon they
/// are restricted to, the physical-size premise for solving children in
/// parallel.
bool source_is_parallel_linear(const SolverContext& ctx) {
  constexpr long long kNodesPerGeneTaxon = 8;
  const std::vector<GeneIndex>& source = scoring_genes(ctx);
  const long long genes = static_cast<long long>(source.size());
  const long long source_taxa =
      static_cast<long long>(std::max(1, source_taxon_count(ctx)));
  if (genes == 0) return true;
  const long long nodes = gene_set_node_count(source);
  if (genes > std::numeric_limits<long long>::max() / kNodesPerGeneTaxon) {
    return false;
  }
  const long long per_taxon = kNodesPerGeneTaxon * genes;
  if (source_taxa > std::numeric_limits<long long>::max() / per_taxon) {
    return false;
  }
  return nodes <= per_taxon * source_taxa;
}

long long child_set_slot_budget(const SolverContext& ctx) {
  // Live restricted gene sets may hold up to four times the root taxa.
  return 4LL * static_cast<long long>(
                   std::max<std::size_t>(1, context_taxa(ctx).size()));
}

/// Reserves child-budget slots for gene sets not yet built.
class ChildSetReservation {
 public:
  ChildSetReservation() = default;
  ChildSetReservation(const ChildSetReservation&) = delete;
  ChildSetReservation& operator=(const ChildSetReservation&) = delete;
  ~ChildSetReservation() { release(); }

  bool acquire(const SolverContext& ctx, long long slots) {
    if (slots <= 0) return true;
    const long long budget = child_set_slot_budget(ctx);
    ChildSetBudget& state = *ctx.child_budget;
    long long observed = state.reserved_slots.load(std::memory_order_acquire);
    while (true) {
      const long long live = state.live_slots.load(std::memory_order_acquire);
      if (observed < 0 || live < 0) {
        throw std::logic_error("negative child gene-set slot accounting");
      }
      if (observed > budget || slots > budget - observed) return false;
      const long long next = observed + slots;
      if (live > budget || next > budget - live) return false;
      if (state.reserved_slots.compare_exchange_weak(
              observed, next, std::memory_order_acq_rel,
              std::memory_order_acquire)) {
        slots_ = slots;
        budget_ = ctx.child_budget;
        stats_max(run_stats().max_parallel_child_slots, next);
        return true;
      }
    }
  }

  void release() {
    if (slots_ <= 0) return;
    budget_->reserved_slots.fetch_sub(slots_, std::memory_order_acq_rel);
    slots_ = 0;
    budget_.reset();
  }

 private:
  long long slots_ = 0;
  std::shared_ptr<ChildSetBudget> budget_;
};

bool try_reserve_parallel_children(const SolverContext& ctx,
                                   const std::vector<int>& left,
                                   const std::vector<int>& right,
                                   ChildSetReservation& reservation) {
  const long long left_slots = projected_child_set_slots(ctx, left);
  const long long right_slots = projected_child_set_slots(ctx, right);
  if (right_slots > std::numeric_limits<long long>::max() - left_slots ||
      !reservation.acquire(ctx, left_slots + right_slots)) {
    stats_add(run_stats().parallel_child_budget_rejections);
    return false;
  }
  stats_add(run_stats().parallel_children);
  return true;
}

/// The position of `taxon` in the ascending `taxa`, or -1.
int position_in(const std::vector<int>& taxa, int taxon) {
  const auto it = std::lower_bound(taxa.begin(), taxa.end(), taxon);
  return it != taxa.end() && *it == taxon
             ? static_cast<int>(it - taxa.begin())
             : -1;
}

enum class ChildGeneSetDecision : std::uint8_t {
  NotNeeded,
  Materialized,
  HeavyReuse,
};

struct ChildGeneSet {
  std::shared_ptr<const std::vector<GeneIndex>> genes;
  ChildGeneSetDecision decision = ChildGeneSetDecision::NotNeeded;
  bool attempted = false;
  long long source_nodes = 0;
  long long candidate_nodes = 0;
};

/// Restricts the scoring genes of `ctx` to `child_taxa`. The result is kept
/// only if it has at most three quarters of the parent's nodes; otherwise the
/// child reuses the parent's genes (a heavy reuse).
ChildGeneSet build_child_gene_set(const SolverContext& ctx,
                                  const std::vector<int>& child_taxa,
                                  bool accounting) {
  const std::vector<GeneIndex>& source = scoring_genes(ctx);
  ChildGeneSet result;
  result.source_nodes = gene_set_node_count(source);
  result.attempted = should_restrict_genes(source, child_taxa);
  if (!result.attempted) return result;

  auto genes = std::make_unique<std::vector<GeneIndex>>(source.size());
  RunStats& stats = run_stats();
  stats_add(stats.child_gene_set_attempts);
  stats_add(stats.child_gene_set_input_nodes, result.source_nodes);
  const auto started = std::chrono::steady_clock::now();
  // The index maps and per-gene workspaces are freed before the set is
  // published; the loop bounds the live workspaces by the workers. Genes
  // differ widely in size, so they are handed out one at a time.
  if (ctx.local_taxa != nullptr) {
    const std::vector<int>& parent_taxa = *ctx.local_taxa;
    std::vector<int> local_to_child(parent_taxa.size());
    for (std::size_t i = 0; i < parent_taxa.size(); ++i) {
      local_to_child[i] = position_in(child_taxa, parent_taxa[i]);
    }
    parallel_for_dynamic(0, static_cast<int>(source.size()), [&](int g) {
      (*genes)[static_cast<std::size_t>(g)] = restrict_gene_by_local(
          source[static_cast<std::size_t>(g)], local_to_child);
    }, 1);
  } else {
    std::vector<int> taxon_to_child(
        static_cast<std::size_t>(context_taxa(ctx).size()), -1);
    for (int i = 0; i < static_cast<int>(child_taxa.size()); ++i) {
      taxon_to_child[static_cast<std::size_t>(
          child_taxa[static_cast<std::size_t>(i)])] = i;
    }
    parallel_for_dynamic(0, static_cast<int>(source.size()), [&](int g) {
      (*genes)[static_cast<std::size_t>(g)] = restrict_gene_by_taxon(
          source[static_cast<std::size_t>(g)], taxon_to_child);
    }, 1);
  }
  result.candidate_nodes = gene_set_node_count(*genes);
  stats_add(stats.child_gene_set_output_nodes, result.candidate_nodes);
  stats_add_time(stats.child_gene_set_ns,
                 std::chrono::steady_clock::now() - started);
  GeneSetLease candidate_lease(accounting,
                               gene_set_payload(*genes, accounting));
  if (run_stats_enabled()) {
    stats_max(stats.max_child_gene_set_bytes, estimate_gene_set_bytes(*genes));
  }

  // Keeping a set must shrink the physical nodes by a quarter; the pair
  // argument in validate_child_gene_sets shows at most one child can fail.
  const long long state_limit =
      result.source_nodes - (result.source_nodes + 3ll) / 4ll;
  if (result.candidate_nodes > state_limit) {
    stats_add(stats.child_gene_set_heavy_reuses);
    result.decision = ChildGeneSetDecision::HeavyReuse;
    return result;
  }

  stats_max(stats.child_gene_set_max_accepted_ratio_ppm,
            static_cast<long long>(std::llround(
                static_cast<long double>(result.candidate_nodes) * 1000000.0L /
                static_cast<long double>(result.source_nodes))));
  stats_add(stats.child_gene_set_accepts);
  result.decision = ChildGeneSetDecision::Materialized;
  candidate_lease.release();
  result.genes =
      track_child_gene_set(std::move(genes), static_cast<int>(child_taxa.size()),
                           accounting, ctx.child_budget);
  return result;
}

void validate_child_gene_sets(const ChildGeneSet& left,
                              const ChildGeneSet& right) {
  if (!left.attempted || !right.attempted || left.source_nodes == 0 ||
      right.source_nodes == 0 ||
      left.decision != ChildGeneSetDecision::HeavyReuse ||
      right.decision != ChildGeneSetDecision::HeavyReuse) {
    return;
  }
  // In each parent gene let A and B be the leaves kept by the two children,
  // U the other leaves and I the internal nodes, so L = A+B+U leaves. A
  // restricted child keeps its leaves and at most every internal node, so
  // |left| + |right| <= A+B+2I. Every unlabelled internal node has degree at
  // least three, except possibly one unary root, so I <= L and
  // 2(|left|+|right|) <= 2L+4I <= 3(L+I). Summed over the genes, the two
  // children cannot both exceed three quarters of the parent.
  RunStats& stats = run_stats();
  stats_add(stats.child_gene_set_pair_violations);
  stats.complexity_contract_satisfied.store(0, std::memory_order_relaxed);
  throw std::logic_error("both child gene sets failed the size gate");
}

Candidate forced_candidate(int taxa, ForcedSplitShape shape) {
  Candidate candidate;
  candidate.split = Bitset(taxa);
  const int left = shape == ForcedSplitShape::Balanced ? taxa / 2 : 1;
  for (int i = 0; i < left; ++i) candidate.split.set(i);
  return candidate;
}

/// One subproblem of the divide and regraft, on the ascending `taxa`. The
/// frames live on the heap and the recursion runs as a loop over them, so a
/// skewed recursion, such as a starting-tree polytomy peeled one arm per
/// level, uses no call stack per level.
struct SubproblemFrame {
  enum class Stage : std::uint8_t {
    Start,          // not divided yet
    AwaitLeft,      // solving the left child; the right one waits
    AwaitRight,     // solving the right child; the left tree is complete
    AwaitChildren,  // solving both children in parallel
  };

  // Owned by the parent frame, or by the caller for the root.
  SolverContext* ctx = nullptr;
  std::vector<int> taxa;
  uint64_t path = kRootPath;
  int depth = 0;
  // Whether the children may be solved in parallel. A serial subproblem
  // solves its whole subtree serially, so only a parallel subproblem's
  // children may be parallel themselves.
  bool parallel = false;
  TaxaLease taxa_lease;
  SubproblemFrame* parent = nullptr;
  // This frame's index among its parent's children: 0 left, 1 right.
  int slot = 0;
  Stage stage = Stage::Start;

  int active = 0;
  std::vector<int> left;
  std::vector<int> right;
  TaxaLease left_lease;
  TaxaLease right_lease;
  uint64_t left_path = 0;
  uint64_t right_path = 0;
  ChildGeneSet left_set;
  std::optional<ChildGeneSet> right_set;
  // Indexed by slot. A serial subproblem holds one child at a time.
  std::array<std::unique_ptr<SolverContext>, 2> child_ctx;
  std::array<std::unique_ptr<SubproblemFrame>, 2> child;
  std::array<SolverTree, 2> trees;
  std::array<TreeLease, 2> tree_leases;
  // A parallel subproblem's children that have not finished, and their
  // failures.
  std::atomic<int> pending{0};
  std::array<std::exception_ptr, 2> failures;
};

std::unique_ptr<SubproblemFrame> make_child_frame(SubproblemFrame& parent,
                                                  int slot,
                                                  SolverContext& ctx,
                                                  std::vector<int> taxa,
                                                  uint64_t path,
                                                  bool parallel,
                                                  TaxaLease taxa_lease) {
  auto child = std::make_unique<SubproblemFrame>();
  child->ctx = &ctx;
  child->taxa = std::move(taxa);
  child->path = path;
  child->depth = parent.depth + 1;
  child->parallel = parallel;
  child->taxa_lease = std::move(taxa_lease);
  child->parent = &parent;
  child->slot = slot;
  return child;
}

/// The divide and regraft of one run: each subproblem chooses a bipartition,
/// solves both sides and regrafts them. A thread works through the frames in
/// a loop: it descends into a child, and when a subproblem is solved hands
/// its tree to the parent and continues there. A parallel subproblem solves
/// its left child on the current thread and its right child as an
/// asynchronous job; whichever finishes last regrafts the two and continues
/// upwards, so no thread waits for a child and no call nests per level.
/// Serial subproblems perform the same steps in the same order as a
/// depth-first recursion, which the payload accounting relies on.
class DivideAndRegraft : public std::enable_shared_from_this<DivideAndRegraft> {
 public:
  explicit DivideAndRegraft(const RecursionControl& control)
      : control_(control) {}

  /// Solves the subproblem on the ascending `taxa` of `ctx`; `parallel` lets
  /// its children run in parallel. Returns once every job has handed back
  /// its tree, and rethrows the failure the recursion ends with.
  SolverTree solve(SolverContext& ctx,
                   std::vector<int> taxa,
                   bool parallel,
                   TaxaLease taxa_lease) {
    root_ = std::make_unique<SubproblemFrame>();
    SubproblemFrame* root = root_.get();
    root->ctx = &ctx;
    root->taxa = std::move(taxa);
    root->parallel = parallel;
    root->taxa_lease = std::move(taxa_lease);
    drive(root);
    if (!done_.load(std::memory_order_acquire)) {
      // Jobs on other workers still hold parts of the recursion. A worker
      // runs tasks while it waits; any other thread blocks.
      tf::Executor& executor = get_executor();
      if (executor.this_worker_id() >= 0) {
        executor.corun_until(
            [this] { return done_.load(std::memory_order_acquire); });
      } else {
        done_.wait(false, std::memory_order_acquire);
      }
    }
    if (root_failure_) std::rethrow_exception(root_failure_);
    return std::move(root_tree_);
  }

 private:
  /// What a stage asks for next: the child frame to solve, or, when `child`
  /// is null, the solved subproblem's tree.
  struct Step {
    SubproblemFrame* child = nullptr;
    SolverTree tree;
  };

  /// How a subproblem ended.
  struct Outcome {
    SolverTree tree;
    std::exception_ptr failure;
  };

  /// Divides the subproblem of a new frame and starts its children, or
  /// solves a base case.
  Step begin(SubproblemFrame& frame) {
    SolverContext& ctx = *frame.ctx;
    const RecursionControl& control = control_;
    const int active = static_cast<int>(frame.taxa.size());
    frame.active = active;
    RunStats& stats = run_stats();
    if (control.accounting) {
      stats_add(stats.recursion_calls);
      stats_add(stats.recursion_taxon_mass, active);
      stats_max(stats.recursion_max_depth, frame.depth);
    }
    if (frame.depth == control.fail_depth) {
      throw std::runtime_error("injected recursion failure at path " +
                               std::to_string(frame.path));
    }
    if (active <= 3) {
      Step step;
      step.tree = base_tree(frame.taxa);
      TreeLease result_lease(control.accounting,
                             tree_payload(step.tree, control.accounting));
      frame.taxa_lease.release();
      return step;
    }

    if (control.accounting) stats_add(stats.recursion_nonbase_calls);
    stats_add(stats.subproblems);
    stats_max(stats.max_depth, frame.depth);
    stats_max(stats.max_active_taxa, active);
    note_scoring_source(ctx, active);
    progress_note_subproblem(frame.path, frame.depth, active);

    progress_note_phase(frame.path, "divide");
    std::vector<int>& left = frame.left;
    std::vector<int>& right = frame.right;
    {
      BipartitionChoice split =
          choose_bipartition(ctx, frame.taxa, frame.path, frame.depth);
      if (control.force_shape) {
        split.candidate = forced_candidate(active, control.shape);
      }
      const Bitset& side = split.candidate.split;
      const int left_count = std::clamp(side.count(), 0, active);
      if (control.accounting) stats_add(stats.recursion_split_tests, active);
      if (left_count < 1 || left_count == active) {
        throw std::logic_error("bipartition has an empty side");
      }
      left.reserve(static_cast<std::size_t>(left_count));
      right.reserve(static_cast<std::size_t>(active - left_count));
      for (int i = 0; i < active; ++i) {
        (side.test(i) ? left : right)
            .push_back(frame.taxa[static_cast<std::size_t>(i)]);
      }
      std::vector<int>().swap(frame.taxa);
      frame.taxa_lease.release();
      frame.left_lease = TaxaLease(control.accounting,
                                   static_cast<long long>(left.size()));
      frame.right_lease = TaxaLease(control.accounting,
                                    static_cast<long long>(right.size()));
      if (control.accounting) {
        stats_add(stats.recursion_partition_moves, active);
        const int left_size = static_cast<int>(left.size());
        const int right_size = static_cast<int>(right.size());
        if (left_size == right_size) {
          stats_add(stats.recursion_balanced_splits);
        }
        if (std::min(left_size, right_size) == 1) {
          stats_add(stats.recursion_skewed_splits);
        }
      }
    }
    progress_note_split(frame.path, active, static_cast<int>(left.size()),
                        static_cast<int>(right.size()));

    frame.left_path = child_path(frame.path, 1, frame.depth);
    frame.right_path = child_path(frame.path, 2, frame.depth);

    // Children run in parallel only when the gene loops have too few genes
    // per worker to keep every worker busy, the scoring genes are linear in
    // their taxa, and the run's child budget can hold both restricted gene
    // sets.
    ChildSetReservation reservation;
    const long long scoring_gene_count =
        static_cast<long long>(scoring_genes(ctx).size());
    const bool wants_parallel =
        frame.parallel && get_num_workers() > 1 &&
        active >= kParallelChildMinTaxa &&
        scoring_gene_count <
            static_cast<long long>(kParallelChildGenesPerWorker) *
                get_num_workers();
    const bool size_ok = !wants_parallel || source_is_parallel_linear(ctx);
    if (wants_parallel && !size_ok) {
      stats_add(stats.parallel_child_size_rejections);
    }
    const bool parallel_children =
        wants_parallel && size_ok &&
        try_reserve_parallel_children(ctx, left, right, reservation);

    if (parallel_children) {
      progress_note_phase(frame.path, "child gene sets");
      ChildGeneSet left_set =
          build_child_gene_set(ctx, left, control.accounting);
      ChildGeneSet right_set =
          build_child_gene_set(ctx, right, control.accounting);
      validate_child_gene_sets(left_set, right_set);
      reservation.release();
      frame.child_ctx[0] = std::make_unique<SolverContext>(make_child_context(
          ctx, std::move(left_set.genes), left, control.accounting));
      frame.child_ctx[1] = std::make_unique<SolverContext>(make_child_context(
          ctx, std::move(right_set.genes), right, control.accounting));
      progress_note_phase(frame.path, "children");
      frame.child[0] = make_child_frame(
          frame, 0, *frame.child_ctx[0], std::move(left), frame.left_path,
          true, std::move(frame.left_lease));
      frame.child[1] = make_child_frame(
          frame, 1, *frame.child_ctx[1], std::move(right), frame.right_path,
          true, std::move(frame.right_lease));
      frame.stage = SubproblemFrame::Stage::AwaitChildren;
      frame.pending.store(2, std::memory_order_relaxed);
      SubproblemFrame* first = frame.child[0].get();
      launch(frame.child[1].get());
      return Step{first, {}};
    }

    stats_add(stats.serial_children);
    // The right gene set is built before the left child only when the left
    // one is a heavy reuse, so that the pair gate can be checked.
    progress_note_phase(frame.path, "left gene set");
    frame.left_set = build_child_gene_set(ctx, left, control.accounting);
    if (frame.left_set.decision == ChildGeneSetDecision::HeavyReuse) {
      progress_note_phase(frame.path, "right gene set");
      frame.right_set.emplace(
          build_child_gene_set(ctx, right, control.accounting));
      validate_child_gene_sets(frame.left_set, *frame.right_set);
    }

    frame.right_lease.mark_continuation();
    frame.child_ctx[0] = std::make_unique<SolverContext>(make_child_context(
        ctx, std::move(frame.left_set.genes), left, control.accounting));
    progress_note_stage("left child");
    frame.child[0] = make_child_frame(frame, 0, *frame.child_ctx[0],
                                      std::move(left), frame.left_path,
                                      false, std::move(frame.left_lease));
    frame.stage = SubproblemFrame::Stage::AwaitLeft;
    return Step{frame.child[0].get(), {}};
  }

  /// Takes a serial subproblem's solved child: after the left one starts the
  /// right one, after the right one regrafts the two.
  Step resume(SubproblemFrame& frame, SolverTree child_tree) {
    SolverContext& ctx = *frame.ctx;
    const RecursionControl& control = control_;
    if (frame.stage == SubproblemFrame::Stage::AwaitLeft) {
      frame.trees[0] = std::move(child_tree);
      frame.tree_leases[0].acquire(
          control.accounting, tree_payload(frame.trees[0], control.accounting));
      frame.child_ctx[0].reset();
      frame.right_lease.unmark_continuation();
      release_after_subproblem(frame.active);

      progress_note_phase(frame.path, "right gene set");
      if (!frame.right_set.has_value()) {
        frame.right_set.emplace(
            build_child_gene_set(ctx, frame.right, control.accounting));
        validate_child_gene_sets(frame.left_set, *frame.right_set);
      }
      frame.child_ctx[1] = std::make_unique<SolverContext>(make_child_context(
          ctx, std::move(frame.right_set->genes), frame.right,
          control.accounting));
      progress_note_stage("right child");
      frame.child[1] = make_child_frame(frame, 1, *frame.child_ctx[1],
                                        std::move(frame.right),
                                        frame.right_path, false,
                                        std::move(frame.right_lease));
      frame.stage = SubproblemFrame::Stage::AwaitRight;
      return Step{frame.child[1].get(), {}};
    }

    frame.trees[1] = std::move(child_tree);
    frame.tree_leases[1].acquire(
        control.accounting, tree_payload(frame.trees[1], control.accounting));
    frame.child_ctx[1].reset();
    release_after_subproblem(frame.active);

    TrackedTaxa left_taxa = taxa_of_tree(frame.trees[0], control.accounting);
    TrackedTaxa right_taxa = taxa_of_tree(frame.trees[1], control.accounting);
    progress_note_stage("regraft");
    progress_note_regraft_start(frame.active,
                                static_cast<int>(left_taxa.values.size()),
                                static_cast<int>(right_taxa.values.size()));
    const auto edges =
        choose_regraft_edges(ctx, frame.trees[0], frame.trees[1],
                             left_taxa.values, right_taxa.values);
    progress_note_regraft_done(frame.active);
    stats_add(run_stats().completed_subproblems);
    TreeLease connect_lease(
        control.accounting,
        connect_payload(frame.trees[1], edges.first, edges.second,
                        control.accounting));
    Step step;
    step.tree = connect_at_edges(std::move(frame.trees[0]),
                                 std::move(frame.trees[1]), edges.first,
                                 edges.second);
    frame.tree_leases[0].release();
    frame.tree_leases[1].release();
    connect_lease.release();
    TreeLease result_lease(control.accounting,
                           tree_payload(step.tree, control.accounting));
    return step;
  }

  /// Keeps a parallel subproblem's finished child and releases its context.
  void record(SubproblemFrame& frame, int slot, Outcome outcome) noexcept {
    const std::size_t index = static_cast<std::size_t>(slot);
    if (outcome.failure) {
      frame.failures[index] = std::move(outcome.failure);
    } else {
      frame.trees[index] = std::move(outcome.tree);
      frame.tree_leases[index].acquire(
          control_.accounting,
          tree_payload(frame.trees[index], control_.accounting));
    }
    frame.child_ctx[index].reset();
  }

  /// Regrafts a parallel subproblem's two children once both have finished.
  /// A failure of the left child takes precedence over one of the right, so
  /// the rethrown failure does not depend on which finished first.
  SolverTree join(SubproblemFrame& frame) {
    ScopedInnerParallelBudgetRestore full_budget;
    SolverContext& ctx = *frame.ctx;
    const RecursionControl& control = control_;
    rethrow_parallel_child_failure(frame.failures[0], frame.failures[1]);
    SolverTree out;
    {
      TrackedTaxa left_taxa = taxa_of_tree(frame.trees[0], control.accounting);
      TrackedTaxa right_taxa =
          taxa_of_tree(frame.trees[1], control.accounting);
      progress_note_stage("regraft");
      progress_note_regraft_start(frame.active,
                                  static_cast<int>(left_taxa.values.size()),
                                  static_cast<int>(right_taxa.values.size()));
      const auto edges =
          choose_regraft_edges(ctx, frame.trees[0], frame.trees[1],
                               left_taxa.values, right_taxa.values);
      progress_note_regraft_done(frame.active);
      TreeLease connect_lease(
          control.accounting,
          connect_payload(frame.trees[1], edges.first, edges.second,
                          control.accounting));
      out = connect_at_edges(std::move(frame.trees[0]),
                             std::move(frame.trees[1]), edges.first,
                             edges.second);
      frame.tree_leases[0].release();
      frame.tree_leases[1].release();
      connect_lease.release();
    }
    stats_add(run_stats().completed_subproblems);
    TreeLease result_lease(control.accounting,
                           tree_payload(out, control.accounting));
    return out;
  }

  /// Solves `frame` on another worker. The job keeps this object alive until
  /// it has handed its tree back.
  void launch(SubproblemFrame* frame) {
    get_executor().silent_async([self = shared_from_this(), frame]() {
      ScopedInnerParallelBudgetRestore full_budget;
      self->drive(frame);
    });
  }

  /// Works on the recursion from the new `frame` until this thread's branch
  /// waits for a sibling solved elsewhere or the root is solved. A failure
  /// ends every enclosing serial subproblem; a parallel one ends with it once
  /// its other child has finished.
  void drive(SubproblemFrame* frame) noexcept {
    bool entering = true;
    Outcome outcome;
    while (true) {
      if (entering) {
        outcome = Outcome{};
        try {
          Step step = begin(*frame);
          if (step.child != nullptr) {
            frame = step.child;
            continue;
          }
          outcome.tree = std::move(step.tree);
        } catch (...) {
          outcome.failure = std::current_exception();
        }
        entering = false;
      }

      // `frame` has ended with `outcome`; hand it to the parent.
      SubproblemFrame* parent = frame->parent;
      if (parent == nullptr) {
        finish(std::move(outcome));
        return;
      }
      const int slot = frame->slot;
      parent->child[static_cast<std::size_t>(slot)].reset();
      frame = parent;
      if (frame->stage == SubproblemFrame::Stage::AwaitChildren) {
        record(*frame, slot, std::move(outcome));
        outcome = Outcome{};
        if (frame->pending.fetch_sub(1, std::memory_order_acq_rel) != 1) {
          return;
        }
        try {
          outcome.tree = join(*frame);
        } catch (...) {
          outcome.failure = std::current_exception();
        }
        continue;
      }
      if (outcome.failure) continue;
      try {
        Step step = resume(*frame, std::move(outcome.tree));
        outcome = Outcome{};
        if (step.child != nullptr) {
          frame = step.child;
          entering = true;
          continue;
        }
        outcome.tree = std::move(step.tree);
      } catch (...) {
        outcome = Outcome{};
        outcome.failure = std::current_exception();
      }
    }
  }

  /// Keeps the root's outcome and wakes the thread in solve. Every other
  /// frame has been released by then.
  void finish(Outcome outcome) noexcept {
    root_tree_ = std::move(outcome.tree);
    root_failure_ = std::move(outcome.failure);
    root_.reset();
    done_.store(true, std::memory_order_release);
    done_.notify_all();
  }

  RecursionControl control_;
  std::unique_ptr<SubproblemFrame> root_;
  SolverTree root_tree_;
  std::exception_ptr root_failure_;
  std::atomic<bool> done_{false};
};

std::vector<int> all_taxa(const SolverContext& ctx) {
  std::vector<int> taxa(static_cast<std::size_t>(context_taxa(ctx).size()));
  std::iota(taxa.begin(), taxa.end(), 0);
  return taxa;
}

SolverTree solve_with_control(SolverContext& ctx,
                              tf::Subflow* subflow,
                              const RecursionControl& control) {
  std::vector<int> taxa = all_taxa(ctx);
  TaxaLease taxa_lease(control.accounting,
                       static_cast<long long>(taxa.size()));
  const GeneSetPayload root_payload =
      gene_set_payload(scoring_genes(ctx), control.accounting);
  GeneSetLease root_lease(control.accounting, root_payload);
  if (control.accounting) {
    RunStats& stats = run_stats();
    stats_add(stats.recursion_root_gene_set_nodes, root_payload.nodes);
    stats_add(stats.recursion_root_gene_set_entries, root_payload.entries);
  }
  SolverTree tree = std::make_shared<DivideAndRegraft>(control)->solve(
      ctx, std::move(taxa), subflow != nullptr, std::move(taxa_lease));
  TreeLease tree_lease(control.accounting,
                       tree_payload(tree, control.accounting));
  refine_nni(ctx, tree, all_taxa(ctx));
  refine_spr(ctx, tree, all_taxa(ctx));
  return tree;
}

}  // namespace

SolverTree solve_qraft(SolverContext& ctx, tf::Subflow* subflow) {
  const RecursionControl control{!ctx.options.stats_json.empty(), false,
                                 ForcedSplitShape::Balanced, -1};
  return solve_with_control(ctx, subflow, control);
}

SolverTree solve_qraft_for_test(SolverContext& ctx,
                                ForcedSplitShape shape,
                                int fail_depth) {
  const RecursionControl control{true, true, shape, fail_depth};
  if (ctx.workers > 1 && get_num_workers() > 1) {
    SolverTree tree;
    std::exception_ptr failure;
    tf::Taskflow flow;
    flow.emplace([&](tf::Subflow& sf) {
      try {
        tree = solve_with_control(ctx, &sf, control);
      } catch (...) {
        failure = std::current_exception();
      }
    });
    run_inner(flow);
    if (failure) std::rethrow_exception(failure);
    return tree;
  }
  return solve_with_control(ctx, nullptr, control);
}

SolverTree refine_starting_tree(SolverContext& ctx, bool run_nni) {
  const StartingTree& starting_tree = context_starting_tree(ctx);
  SolverTree tree;
  if (starting_tree.root < 0 || starting_tree.nodes.empty()) return tree;

  tree.nodes.resize(starting_tree.nodes.size());
  for (int i = 0; i < static_cast<int>(starting_tree.nodes.size()); ++i) {
    const StartingTree::Node& source =
        starting_tree.nodes[static_cast<std::size_t>(i)];
    SolverTree::Node& node = tree.nodes[static_cast<std::size_t>(i)];
    node.taxon = source.taxon;
    node.neighbors.reserve(source.children.size() +
                           (source.parent >= 0 ? 1 : 0));
    if (source.parent >= 0) node.neighbors.push_back(source.parent);
    for (int child : source.children) node.neighbors.push_back(child);
  }
  if (run_nni) {
    // NNI scores only edges whose two ends have three arms; a binary root
    // would leave the edges around it unscored and fixed.
    suppress_solver_degree_two(tree);
    refine_nni(ctx, tree, all_taxa(ctx));
  }
  refine_spr(ctx, tree, all_taxa(ctx));
  return tree;
}

}  // namespace qraft::internal
