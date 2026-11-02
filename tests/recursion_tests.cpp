#include "qraft/qraft.hpp"
#include "qraft/tree.hpp"

#include "solver/solver.hpp"
#include "counting/hdt.hpp"
#include "counting/query_taxa.hpp"
#include "support/small_stack.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <barrier>
#include <bit>
#include <charconv>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#if defined(__linux__)
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace {

[[noreturn]] void fail(const std::string& message) {
  throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
  if (!condition) fail(message);
}

void require_finite(double value, const std::string& message) {
  if (!std::isfinite(value)) {
    fail(message + ": non-finite value=" + std::to_string(value));
  }
}

void require_close(double actual,
                   double expected,
                   const std::string& message) {
  require_finite(actual, message + " actual");
  require_finite(expected, message + " expected");
  const double scale =
      std::max({1.0, std::fabs(actual), std::fabs(expected)});
  if (std::fabs(actual - expected) > 1e-10 * scale) {
    fail(message + ": actual=" + std::to_string(actual) +
         " expected=" + std::to_string(expected));
  }
}

template <typename Fn>
void require_invalid_argument(Fn&& fn,
                              const std::string& expected_message,
                              const std::string& context) {
  bool caught = false;
  try {
    std::forward<Fn>(fn)();
  } catch (const std::invalid_argument& error) {
    caught = true;
    require(error.what() == expected_message,
            context + ": unexpected message=" + error.what());
  } catch (const std::exception& error) {
    fail(context + ": wrong exception class: " + error.what());
  } catch (...) {
    fail(context + ": non-standard exception class");
  }
  require(caught, context + ": expected std::invalid_argument");
}

void atomic_min(std::atomic<int>& target, int value) {
  int current = target.load(std::memory_order_relaxed);
  while (value < current &&
         !target.compare_exchange_weak(current, value,
                                       std::memory_order_relaxed,
                                       std::memory_order_relaxed)) {
  }
}

class LogicalFailure : public std::runtime_error {
 public:
  explicit LogicalFailure(int logical_index)
      : std::runtime_error("logical failure at " +
                           std::to_string(logical_index)),
        logical_index(logical_index) {}

  int logical_index;
};

enum class FailureMode { Caller, Helper, Simultaneous };

const char* failure_mode_name(FailureMode mode) {
  switch (mode) {
    case FailureMode::Caller:
      return "caller";
    case FailureMode::Helper:
      return "helper";
    case FailureMode::Simultaneous:
      return "simultaneous";
  }
  return "unknown";
}

void exercise_parallel_failure_case(bool dynamic, FailureMode mode) {
  using namespace qraft::internal;
  constexpr int kWorkers = 4;
  constexpr int kIterations = 32;
  constexpr int kStaticChunk = kIterations / kWorkers;

  std::array<std::atomic<int>, kIterations> visits{};
  std::atomic<int> expected_index{std::numeric_limits<int>::max()};
  std::atomic<int> caller_failures{0};
  std::atomic<int> helper_failures{0};
  std::barrier start_gate(kWorkers);
  const std::thread::id caller = std::this_thread::get_id();

  auto body = [&](int index) {
    visits[static_cast<std::size_t>(index)].fetch_add(
        1, std::memory_order_relaxed);
    if ((dynamic && index < kWorkers) ||
        (!dynamic && index % kStaticChunk == 0)) {
      start_gate.arrive_and_wait();
    }

    const bool on_caller = std::this_thread::get_id() == caller;
    const bool should_fail =
        mode == FailureMode::Simultaneous ||
        (mode == FailureMode::Caller && on_caller) ||
        (mode == FailureMode::Helper && !on_caller);
    if (!should_fail) return;

    atomic_min(expected_index, index);
    (on_caller ? caller_failures : helper_failures)
        .fetch_add(1, std::memory_order_relaxed);
    throw LogicalFailure(index);
  };

  int rethrown_index = -1;
  try {
    if (dynamic) {
      parallel_for_dynamic(0, kIterations, body, 1);
    } else {
      parallel_for_inner(0, kIterations, body, 1);
    }
    fail(std::string(dynamic ? "dynamic" : "static") + " " +
         failure_mode_name(mode) + " failure was not rethrown");
  } catch (const LogicalFailure& ex) {
    rethrown_index = ex.logical_index;
  }

  const int expected = expected_index.load(std::memory_order_relaxed);
  require(expected != std::numeric_limits<int>::max(),
          "parallel failure fixture did not execute a failing callback");
  require(rethrown_index == expected,
          std::string(dynamic ? "dynamic" : "static") + " " +
              failure_mode_name(mode) +
              " loop did not rethrow the lowest logical failure");
  for (int index = 0; index < kIterations; ++index) {
    require(visits[static_cast<std::size_t>(index)].load(
                std::memory_order_relaxed) == 1,
            std::string(dynamic ? "dynamic" : "static") +
                " loop did not complete every logical iteration exactly once");
  }
  if (mode != FailureMode::Helper) {
    require(caller_failures.load(std::memory_order_relaxed) > 0,
            "caller failure path was not exercised");
  }
  if (mode != FailureMode::Caller) {
    require(helper_failures.load(std::memory_order_relaxed) > 0,
            "helper failure path was not exercised");
  }
}

void exercise_nested_worker_loop(bool dynamic) {
  using namespace qraft::internal;
  constexpr int kIterations = 24;
  auto completion = std::make_shared<std::promise<void>>();
  std::future<void> completed = completion->get_future();

  get_executor().silent_async([completion, dynamic]() {
    try {
      std::array<std::atomic<int>, kIterations> visits{};
      auto body = [&](int index) {
        visits[static_cast<std::size_t>(index)].fetch_add(
            1, std::memory_order_relaxed);
        if (index == 2 || index == 17) throw LogicalFailure(index);
      };
      int rethrown = -1;
      try {
        if (dynamic) {
          parallel_for_dynamic(0, kIterations, body, 1);
        } else {
          parallel_for_inner(0, kIterations, body, 1);
        }
      } catch (const LogicalFailure& ex) {
        rethrown = ex.logical_index;
      }
      require(rethrown == 2,
              "nested executor loop did not select the lowest failure");
      for (const std::atomic<int>& visit : visits) {
        require(visit.load(std::memory_order_relaxed) == 1,
                "nested executor loop left an iteration incomplete");
      }
      completion->set_value();
    } catch (...) {
      completion->set_exception(std::current_exception());
    }
  });

  require(completed.wait_for(std::chrono::seconds(10)) ==
              std::future_status::ready,
          std::string(dynamic ? "dynamic" : "static") +
              " nested executor loop timed out");
  completed.get();
}

void test_parallel_loop_exception_safety() {
  qraft::Options options;
  options.threads = 4;
  require(qraft::internal::configure_threads(options) == 4,
          "failed to configure four executor workers");

  for (int repetition = 0; repetition < 20; ++repetition) {
    for (bool dynamic : {false, true}) {
      for (FailureMode mode : {FailureMode::Caller, FailureMode::Helper,
                               FailureMode::Simultaneous}) {
        exercise_parallel_failure_case(dynamic, mode);
      }
    }
  }
  for (bool dynamic : {false, true}) {
    exercise_nested_worker_loop(dynamic);
  }
}

void test_parallel_child_exception_determinism() {
  using namespace qraft::internal;
  for (bool second_finishes_first : {false, true}) {
    struct FailureFrame {
      std::exception_ptr first;
      std::exception_ptr second;
      std::atomic<bool> first_captured{false};
      std::atomic<bool> second_captured{false};
    };
    auto frame = std::make_shared<FailureFrame>();
    tf::Executor executor(2);
    tf::Taskflow flow;
    tf::Task first = flow.emplace([frame, second_finishes_first] {
      if (second_finishes_first) {
        while (!frame->second_captured.load(std::memory_order_acquire)) {
          std::this_thread::yield();
        }
      }
      try {
        throw std::runtime_error("first structural child failed");
      } catch (...) {
        frame->first = std::current_exception();
      }
      frame->first_captured.store(true, std::memory_order_release);
    });
    tf::Task second = flow.emplace([frame, second_finishes_first] {
      if (!second_finishes_first) {
        while (!frame->first_captured.load(std::memory_order_acquire)) {
          std::this_thread::yield();
        }
      }
      try {
        throw std::runtime_error("second structural child failed");
      } catch (...) {
        frame->second = std::current_exception();
      }
      frame->second_captured.store(true, std::memory_order_release);
    });
    tf::Task join = flow.emplace([frame] {
      rethrow_parallel_child_failure(frame->first, frame->second);
    });
    first.precede(join);
    second.precede(join);

    bool selected_first = false;
    try {
      executor.run(flow).get();
    } catch (const std::runtime_error& error) {
      selected_first =
          std::string(error.what()) == "first structural child failed";
    }
    require(frame->first_captured.load(std::memory_order_acquire) &&
                frame->second_captured.load(std::memory_order_acquire),
            "structural child failure probe did not join both children");
    require(selected_first,
            "structural child failure selection depended on completion order");
  }
}

std::vector<qraft::GeneTree> inference_fixture_a() {
  return {
      qraft::parse_newick("((A,B),(C,D));"),
      qraft::parse_newick("((A,B),(C,D));"),
      qraft::parse_newick("((A,C),(B,D));"),
      qraft::parse_newick("((A,B),(E,F));"),
      qraft::parse_newick("((C,D),(E,F));"),
      qraft::parse_newick("((A,E),(B,F));"),
  };
}

std::vector<qraft::GeneTree> inference_fixture_b() {
  return {
      qraft::parse_newick("((P,Q),(R,(S,T)));"),
      qraft::parse_newick("((P,R),(Q,(S,T)));"),
      qraft::parse_newick("((P,Q),(S,(R,T)));"),
      qraft::parse_newick("((P,T),(Q,(R,S)));"),
  };
}

class TemporaryDirectory {
 public:
  TemporaryDirectory() {
    static std::atomic<unsigned long long> sequence{0};
    const auto tick = std::chrono::steady_clock::now().time_since_epoch().count();
    path_ = std::filesystem::temp_directory_path() /
            ("qraft-instrumentation-" + std::to_string(tick) + "-" +
             std::to_string(sequence.fetch_add(1, std::memory_order_relaxed)));
    require(std::filesystem::create_directories(path_),
            "failed to create temporary regression directory");
  }

  TemporaryDirectory(const TemporaryDirectory&) = delete;
  TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;

  ~TemporaryDirectory() {
    std::error_code ignored;
    std::filesystem::remove_all(path_, ignored);
  }

  const std::filesystem::path& path() const { return path_; }

 private:
  std::filesystem::path path_;
};

std::string read_file(const std::filesystem::path& path) {
  std::ifstream in(path);
  require(static_cast<bool>(in), "failed to open " + path.string());
  return std::string((std::istreambuf_iterator<char>(in)),
                     std::istreambuf_iterator<char>());
}

long long json_integer(const std::string& json, const std::string& key) {
  const std::string token = "\"" + key + "\"";
  std::size_t position = json.find(token);
  require(position != std::string::npos, "missing JSON key " + key);
  position = json.find(':', position + token.size());
  require(position != std::string::npos, "missing JSON separator for " + key);
  ++position;
  while (position < json.size() &&
         (json[position] == ' ' || json[position] == '\n' ||
          json[position] == '\r' || json[position] == '\t')) {
    ++position;
  }
  long long value = 0;
  const char* begin = json.data() + static_cast<std::ptrdiff_t>(position);
  const char* end = json.data() + static_cast<std::ptrdiff_t>(json.size());
  const auto parsed = std::from_chars(begin, end, value);
  require(parsed.ec == std::errc{}, "invalid JSON integer for " + key);
  return value;
}

qraft::Options inference_options(int workers,
                                std::filesystem::path stats_path) {
  qraft::Options options;
  options.threads = workers;
  options.progress = qraft::ProgressMode::Off;
  options.nni = false;
  options.stats_json = std::move(stats_path).string();
  return options;
}

struct InferenceResult {
  std::string tree;
  std::string telemetry;
};

InferenceResult run_inference(std::vector<qraft::GeneTree> genes,
                              const qraft::Options& options) {
  InferenceResult result;
  result.tree = qraft::to_newick(
      qraft::infer_species_tree(std::move(genes), options));
  result.telemetry = read_file(options.stats_json);
  return result;
}

void require_run_telemetry(const InferenceResult& actual,
                           const InferenceResult& expected,
                           int genes,
                           int taxa,
                           int workers,
                           const std::string& label) {
  require(json_integer(actual.telemetry, "genes") == genes,
          label + " telemetry used another run's gene count");
  require(json_integer(actual.telemetry, "taxa") == taxa,
          label + " telemetry used another run's taxon count");
  require(json_integer(actual.telemetry, "workers") == workers,
          label + " telemetry used another run's worker setting");

  constexpr std::array<const char*, 19> kStableCounters{
      "subproblems",
      "completed_subproblems",
      "candidate_sets",
      "starting_tree_index_builds",
      "starting_tree_index_nodes",
      "active_taxa_entries",
      "active_taxa_index_entries",
      "query_tree_selected_leaves",
      "query_tree_lca_queries",
      "query_tree_retained_nodes",
      "gene_set_scoring_nodes",
      "gene_set_scoring_calls",
      "recursion_calls",
      "recursion_taxon_mass",
      "recursion_split_tests",
      "recursion_tree_node_visits",
      "peak_payload",
      "peak_gene_set_nodes",
      "accounting_violations",
  };
  for (const char* counter : kStableCounters) {
    require(json_integer(actual.telemetry, counter) ==
                json_integer(expected.telemetry, counter),
            label + " counter contamination at " + counter);
  }
  const long long subproblems = json_integer(actual.telemetry, "subproblems");
  require(subproblems > 0 &&
              subproblems ==
                  json_integer(actual.telemetry, "completed_subproblems"),
          label + " solver-owned subproblem counters disagree");
  require(json_integer(actual.telemetry, "starting_tree_index_builds") == 1,
          label + " did not build exactly one run-owned starting-tree index");
  require(json_integer(actual.telemetry, "live_payload") == 0 &&
              json_integer(actual.telemetry,
                           "accounting_violations") == 0,
          label + " benchmark JSON reports unbalanced live ownership");
}

void test_concurrent_inference() {
  TemporaryDirectory temporary;
  const qraft::Options baseline_options_a = inference_options(
      1, temporary.path() / "baseline-a.json");
  const qraft::Options baseline_options_b = inference_options(
      3, temporary.path() / "baseline-b.json");
  const InferenceResult expected_a =
      run_inference(inference_fixture_a(), baseline_options_a);
  const InferenceResult expected_b =
      run_inference(inference_fixture_b(), baseline_options_b);
  require(expected_a.tree != expected_b.tree,
          "concurrent inference fixtures did not produce distinct outputs");

  for (int round = 0; round < 3; ++round) {
    const qraft::Options options_a = inference_options(
        1,
        temporary.path() / ("round-" + std::to_string(round) + "-a.json"));
    const qraft::Options options_b = inference_options(
        3,
        temporary.path() / ("round-" + std::to_string(round) + "-b.json"));
    std::barrier simultaneous_start(3);
    std::future<InferenceResult> future_a = std::async(
        std::launch::async, [&, options_a] {
          simultaneous_start.arrive_and_wait();
          return run_inference(inference_fixture_a(), options_a);
        });
    std::future<InferenceResult> future_b = std::async(
        std::launch::async, [&, options_b] {
          simultaneous_start.arrive_and_wait();
          return run_inference(inference_fixture_b(), options_b);
        });
    simultaneous_start.arrive_and_wait();

    require(future_a.wait_for(std::chrono::seconds(20)) ==
                    std::future_status::ready &&
                future_b.wait_for(std::chrono::seconds(20)) ==
                    std::future_status::ready,
            "concurrent public inference exceeded its bounded wait");
    const InferenceResult actual_a = future_a.get();
    const InferenceResult actual_b = future_b.get();
    require(actual_a.tree == expected_a.tree,
            "fixture A output changed during concurrent inference");
    require(actual_b.tree == expected_b.tree,
            "fixture B output changed during concurrent inference");
    require_run_telemetry(actual_a, expected_a, 6, 6, 1,
                          "fixture A round " + std::to_string(round));
    require_run_telemetry(actual_b, expected_b, 4, 5, 3,
                          "fixture B round " + std::to_string(round));
  }
}

qraft::internal::StartingTree balanced_starting_tree(int leaf_count) {
  using qraft::internal::StartingTree;
  StartingTree starting_tree;
  starting_tree.taxon_to_leaf.assign(static_cast<std::size_t>(leaf_count), -1);
  auto build = [&](auto&& self, int begin, int end) -> int {
    const int node = static_cast<int>(starting_tree.nodes.size());
    starting_tree.nodes.emplace_back();
    if (end - begin == 1) {
      starting_tree.nodes[static_cast<std::size_t>(node)].taxon = begin;
      starting_tree.taxon_to_leaf[static_cast<std::size_t>(begin)] = node;
      return node;
    }
    const int middle = begin + (end - begin) / 2;
    const int left = self(self, begin, middle);
    const int right = self(self, middle, end);
    starting_tree.nodes[static_cast<std::size_t>(node)].children = {left, right};
    starting_tree.nodes[static_cast<std::size_t>(left)].parent = node;
    starting_tree.nodes[static_cast<std::size_t>(right)].parent = node;
    return node;
  };
  starting_tree.root = leaf_count > 0 ? build(build, 0, leaf_count) : -1;
  return starting_tree;
}

qraft::internal::StartingTree random_polytomy_starting_tree(
    int leaf_count,
    std::mt19937& random) {
  using qraft::internal::StartingTree;
  StartingTree starting_tree;
  starting_tree.taxon_to_leaf.assign(static_cast<std::size_t>(leaf_count), -1);
  std::vector<int> roots;
  roots.reserve(static_cast<std::size_t>(leaf_count));
  for (int taxon = 0; taxon < leaf_count; ++taxon) {
    const int leaf = static_cast<int>(starting_tree.nodes.size());
    starting_tree.nodes.emplace_back();
    starting_tree.nodes.back().taxon = taxon;
    starting_tree.taxon_to_leaf[static_cast<std::size_t>(taxon)] = leaf;
    roots.push_back(leaf);
  }
  while (roots.size() > 1) {
    std::shuffle(roots.begin(), roots.end(), random);
    const int maximum_degree =
        std::min(5, static_cast<int>(roots.size()));
    std::uniform_int_distribution<int> degree_distribution(2, maximum_degree);
    const int degree = degree_distribution(random);
    const int parent = static_cast<int>(starting_tree.nodes.size());
    starting_tree.nodes.emplace_back();
    for (int i = 0; i < degree; ++i) {
      const int child = roots.back();
      roots.pop_back();
      starting_tree.nodes[static_cast<std::size_t>(parent)].children.push_back(child);
      starting_tree.nodes[static_cast<std::size_t>(child)].parent = parent;
    }
    roots.push_back(parent);
  }
  starting_tree.root = roots.empty() ? -1 : roots.front();
  return starting_tree;
}

void test_restriction_rejects_disconnected_forest() {
  using namespace qraft::internal;
  GeneIndex forest;
  forest.nodes.resize(8);
  std::vector<std::vector<int>> forest_adjacency(forest.nodes.size());
  auto add_edge = [&](int first, int second) {
    forest_adjacency[static_cast<std::size_t>(first)].push_back(second);
    forest_adjacency[static_cast<std::size_t>(second)].push_back(first);
  };
  for (int component = 0; component < 2; ++component) {
    const int center = component * 4;
    for (int offset = 1; offset <= 3; ++offset) {
      const int leaf = center + offset;
      forest.nodes[static_cast<std::size_t>(leaf)].set_taxon(
          component * 3 + offset - 1);
      add_edge(center, leaf);
    }
  }
  set_gene_adjacency(forest, forest_adjacency);
  forest.root = 0;
  forest.reference_root = 1;

  require_invalid_argument(
      [&] {
        static_cast<void>(
            restrict_gene_by_taxon(forest, std::vector<int>(6, 0)));
      },
      "gene adjacency must be connected",
      "restriction by taxon accepted a disconnected forest");

  GeneIndex local_forest = forest;
  for (GeneIndex::Node& node : local_forest.nodes) {
    node.set_taxon(-1);
  }
  for (int component = 0; component < 2; ++component) {
    const int center = component * 4;
    for (int offset = 1; offset <= 3; ++offset) {
      GeneIndex::Node& leaf = local_forest.nodes[
          static_cast<std::size_t>(center + offset)];
      leaf.set_local(0);
    }
  }
  require_invalid_argument(
      [&] {
        static_cast<void>(
            restrict_gene_by_local(local_forest, std::vector<int>{0}));
      },
      "gene adjacency must be connected",
      "restriction by local index accepted a disconnected forest");
}

void test_gene_restriction_duplicates() {
  using namespace qraft::internal;
  test_restriction_rejects_disconnected_forest();
  constexpr int kDuplicateLeaves = 32;
  constexpr int kUnmappedLeaves = 3;
  constexpr int kOwnerCount = 4;
  constexpr int kInternalNodes = 3;
  constexpr int kMappedSingletons = 3;
  constexpr int kTaxa =
      kDuplicateLeaves + kUnmappedLeaves + kMappedSingletons;

  GeneIndex original;
  original.nodes.resize(static_cast<std::size_t>(kInternalNodes + kTaxa));
  original.taxon_to_node.assign(static_cast<std::size_t>(kTaxa), -1);
  std::vector<std::vector<int>> original_adjacency(original.nodes.size());
  auto add_edge = [&](int first, int second) {
    original_adjacency[static_cast<std::size_t>(first)].push_back(second);
    original_adjacency[static_cast<std::size_t>(second)].push_back(first);
  };
  auto add_taxon_leaf = [&](int taxon, int parent) {
    const int node = kInternalNodes + taxon;
    original.nodes[static_cast<std::size_t>(node)].set_taxon(taxon);
    original.taxon_to_node[static_cast<std::size_t>(taxon)] = node;
    add_edge(parent, node);
  };

  add_edge(0, 1);
  add_edge(1, 2);
  for (int taxon = 0; taxon < kDuplicateLeaves; ++taxon) {
    add_taxon_leaf(taxon, 0);
  }
  for (int taxon = kDuplicateLeaves;
       taxon < kDuplicateLeaves + kUnmappedLeaves; ++taxon) {
    add_taxon_leaf(taxon, 0);
  }
  const int owner_one_taxon = kDuplicateLeaves + kUnmappedLeaves;
  const int owner_two_taxon = owner_one_taxon + 1;
  const int owner_three_taxon = owner_two_taxon + 1;
  add_taxon_leaf(owner_one_taxon, 1);
  add_taxon_leaf(owner_two_taxon, 2);
  add_taxon_leaf(owner_three_taxon, 2);
  set_gene_adjacency(original, original_adjacency);
  original.root = 1;
  original.reference_root = original.taxon_to_node.front();
  original.root_traversal = traversal(original, original.root);

  // Owner 0 has every duplicate leaf; unmapped leaves carry a large mass that
  // must not survive the restriction.
  std::vector<int> taxon_to_local(static_cast<std::size_t>(kTaxa), -1);
  std::vector<double> taxon_mass(static_cast<std::size_t>(kTaxa), 128.0);
  std::array<double, kOwnerCount> expected_mass{};
  for (int taxon = 0; taxon < kDuplicateLeaves; ++taxon) {
    const double mass = 1.0 + 0.25 * static_cast<double>(taxon % 4);
    taxon_to_local[static_cast<std::size_t>(taxon)] = 0;
    taxon_mass[static_cast<std::size_t>(taxon)] = mass;
    expected_mass[0] += mass;
  }
  for (const std::pair<int, int> mapped :
       {std::pair{owner_one_taxon, 1}, std::pair{owner_two_taxon, 2},
        std::pair{owner_three_taxon, 3}}) {
    const double mass = static_cast<double>(mapped.second + 1);
    taxon_to_local[static_cast<std::size_t>(mapped.first)] = mapped.second;
    taxon_mass[static_cast<std::size_t>(mapped.first)] = mass;
    expected_mass[static_cast<std::size_t>(mapped.second)] = mass;
  }
  for (int taxon = 0; taxon < kTaxa; ++taxon) {
    set_gene_leaf_mass(original, kInternalNodes + taxon,
                       taxon_mass[static_cast<std::size_t>(taxon)]);
  }

  const GeneIndex restricted = restrict_gene_by_taxon(original, taxon_to_local);
  constexpr int kMappedLeaves = kDuplicateLeaves + kMappedSingletons;
  require(restricted.nodes.size() == kInternalNodes + kMappedLeaves,
          "exact pruning coalesced or lost a mapped physical leaf");
  require(restricted.taxon_to_node.empty(),
          "restricted star retained original taxon slots");
  require(restricted.root >= 0 &&
              restricted.root < static_cast<int>(restricted.nodes.size()) &&
              restricted.reference_root >= 0 &&
              restricted.reference_root <
                  static_cast<int>(restricted.nodes.size()),
          "restricted star has an invalid root");
  GeneTraversal reference_scratch;
  require(restricted.root_traversal.preorder.size() == restricted.nodes.size() &&
              gene_reference_traversal(restricted, reference_scratch)
                      .preorder.size() == restricted.nodes.size(),
          "restricted star traversal does not cover the canonical tree");

  std::array<int, kOwnerCount> owner_leaf_count{};
  std::vector<double> owner_zero_masses;
  double total_mass = 0.0;
  int internal_count = 0;
  std::vector<int> internal_degrees;
  for (int node = 0; node < static_cast<int>(restricted.nodes.size()); ++node) {
    const GeneIndex::Node& current =
        restricted.nodes[static_cast<std::size_t>(node)];
    require(current.taxon() < 0,
            "restricted star retained an unmapped taxon");
    for (int neighbor : gene_neighbors(restricted, node)) {
      require(neighbor >= 0 &&
                  neighbor < static_cast<int>(restricted.nodes.size()),
              "restricted star has an out-of-range edge");
      const std::span<const int> reverse =
          gene_neighbors(restricted, neighbor);
      require(std::count(reverse.begin(), reverse.end(), node) == 1,
              "restricted star adjacency is not symmetric");
    }
    if (current.local() < 0) {
      ++internal_count;
      require(gene_degree(restricted, node) >= 3,
              "exact pruned tree retained an uninformative internal node");
      internal_degrees.push_back(gene_degree(restricted, node));
      continue;
    }
    require(current.local() < kOwnerCount,
            "restricted star produced an out-of-range owner");
    require(gene_degree(restricted, node) == 1,
            "restricted owner occurrence is not a leaf");
    ++owner_leaf_count[static_cast<std::size_t>(current.local())];
    const double mass = gene_leaf_mass(restricted, node);
    if (current.local() == 0) owner_zero_masses.push_back(mass);
    total_mass += mass;
  }
  require(internal_count == kInternalNodes,
          "exact pruned tree changed its informative internal topology");
  std::sort(internal_degrees.begin(), internal_degrees.end());
  require(internal_degrees == std::vector<int>({3, 3, kDuplicateLeaves + 1}),
          "exact pruned tree changed the high-degree input topology");
  require(owner_leaf_count ==
              std::array<int, kOwnerCount>{kDuplicateLeaves, 1, 1, 1},
          "exact pruning changed duplicate-owner physical multiplicity");
  std::vector<double> expected_owner_zero_masses;
  for (int taxon = 0; taxon < kDuplicateLeaves; ++taxon) {
    expected_owner_zero_masses.push_back(
        taxon_mass[static_cast<std::size_t>(taxon)]);
  }
  std::sort(owner_zero_masses.begin(), owner_zero_masses.end());
  std::sort(expected_owner_zero_masses.begin(),
            expected_owner_zero_masses.end());
  require(owner_zero_masses == expected_owner_zero_masses,
          "exact pruning regrouped duplicate-owner masses");
  require_close(total_mass,
                std::accumulate(expected_mass.begin(), expected_mass.end(),
                                0.0),
                "mapped mass changed while pruning unmapped leaves");

  GeneIndex expected;
  expected.nodes.resize(static_cast<std::size_t>(kInternalNodes + kMappedLeaves));
  std::vector<std::vector<int>> expected_adjacency(expected.nodes.size());
  auto add_expected_edge = [&](int first, int second) {
    expected_adjacency[static_cast<std::size_t>(first)].push_back(second);
    expected_adjacency[static_cast<std::size_t>(second)].push_back(first);
  };
  add_expected_edge(0, 1);
  add_expected_edge(1, 2);
  for (int occurrence = 0; occurrence < kDuplicateLeaves; ++occurrence) {
    const int leaf_node = kInternalNodes + occurrence;
    add_expected_edge(0, leaf_node);
    expected.nodes[static_cast<std::size_t>(leaf_node)].set_local(0);
    set_gene_leaf_mass(expected, leaf_node,
                       taxon_mass[static_cast<std::size_t>(occurrence)]);
  }
  for (int owner = 1; owner < kOwnerCount; ++owner) {
    const int leaf_node = kInternalNodes + kDuplicateLeaves + owner - 1;
    add_expected_edge(owner == 1 ? 1 : 2, leaf_node);
    expected.nodes[static_cast<std::size_t>(leaf_node)].set_local(owner);
    set_gene_leaf_mass(expected, leaf_node,
                       expected_mass[static_cast<std::size_t>(owner)]);
  }
  set_gene_adjacency(expected, expected_adjacency);
  expected.root = 1;
  expected.reference_root = 1;
  expected.root_traversal = traversal(expected, expected.root);

  // Negative control: a mass-only tree that sums each owner's leaves into one
  // atom. It cannot preserve terms that select two distinct physical leaves of
  // the same owner, so its two-color score must differ.
  GeneIndex collapsed;
  collapsed.nodes.resize(6);
  std::vector<std::vector<int>> collapsed_adjacency(collapsed.nodes.size());
  auto add_collapsed_edge = [&](int first, int second) {
    collapsed_adjacency[static_cast<std::size_t>(first)].push_back(second);
    collapsed_adjacency[static_cast<std::size_t>(second)].push_back(first);
  };
  add_collapsed_edge(0, 1);
  add_collapsed_edge(0, 2);
  add_collapsed_edge(0, 3);
  add_collapsed_edge(1, 4);
  add_collapsed_edge(1, 5);
  for (int owner = 0; owner < kOwnerCount; ++owner) {
    collapsed.nodes[static_cast<std::size_t>(2 + owner)].set_local(owner);
    set_gene_leaf_mass(collapsed, 2 + owner,
                       expected_mass[static_cast<std::size_t>(owner)]);
  }
  set_gene_adjacency(collapsed, collapsed_adjacency);
  collapsed.root = 0;
  collapsed.reference_root = 0;
  collapsed.root_traversal = traversal(collapsed, collapsed.root);

  OwnerProjection owners;
  owners.local_keys.mode = OwnerKeyMode::Dense;
  owners.local_keys.dense.resize(static_cast<std::size_t>(kOwnerCount));
  std::iota(owners.local_keys.dense.begin(), owners.local_keys.dense.end(), 0);
  owners.owner_count = kOwnerCount;
  owners.live_owner.assign(static_cast<std::size_t>(kOwnerCount), 1);

  auto four_color_counts = [&](const GeneIndex& gene) {
    HdtColorPolicy colors;
    colors.current_by_owner = {0, 1, 2, 3};
    colors.background_by_owner.assign(kOwnerCount, 0);
    const HdtState state = build_hdt_state(
        build_base_component_graph(gene, owners, 4, colors));
    const std::array<double, 3> counts = hdt_root_counts(state);
    return ChannelCounts{counts[0], counts[1], counts[2]};
  };
  const ChannelCounts restricted_counts = four_color_counts(restricted);
  const ChannelCounts expected_counts = four_color_counts(expected);
  require_close(restricted_counts.q01_23, expected_counts.q01_23,
                "restricted/exact q01|23 score mismatch");
  require_close(restricted_counts.q02_13, expected_counts.q02_13,
                "restricted/exact q02|13 score mismatch");
  require_close(restricted_counts.q03_12, expected_counts.q03_12,
                "restricted/exact q03|12 score mismatch");
  const double expected_score =
      std::accumulate(expected_mass.begin(), expected_mass.end(), 1.0,
                      std::multiplies<double>());
  require_close(restricted_counts.q01_23, expected_score,
                "restricted quartet score differs from mapped mass product");
  require_close(restricted_counts.q02_13, 0.0,
                "restricted quartet gained q02|13 mass");
  require_close(restricted_counts.q03_12, 0.0,
                "restricted quartet gained q03|12 mass");

  auto two_color_scores = [&](const GeneIndex& gene) {
    HdtColorPolicy colors;
    colors.current_by_owner = {1, 1, 0, 0};
    colors.background_by_owner.assign(kOwnerCount, 0);
    const HdtState state = build_hdt_state(
        build_base_component_graph(gene, owners, 2, colors));
    const std::array<double, 3> counts = hdt_root_counts(state);
    const double score = counts[0] - counts[1];
    require_finite(score, "two-color fixed-vector score");
    return score;
  };
  const double restricted_score = two_color_scores(restricted);
  const double expected_two_color_score = two_color_scores(expected);
  const double collapsed_score = two_color_scores(collapsed);
  require_close(restricted_score, expected_two_color_score,
                "restricted/exact two-color score mismatch");
  require_finite(collapsed_score, "collapsed two-color negative control");
  require(std::fabs(restricted_score - collapsed_score) > 1e-9,
          "mass-only owner coalescing unexpectedly preserved pair terms");
}

qraft::internal::GeneIndex synthetic_local_gene(int active_count) {
  qraft::internal::GeneIndex gene;
  const int internal_count = std::max(0, active_count - 1);
  gene.nodes.resize(static_cast<std::size_t>(2 * active_count - 1));
  for (int local = 0; local < active_count; ++local) {
    gene.nodes[static_cast<std::size_t>(internal_count + local)]
        .set_local(local);
  }
  return gene;
}

struct SetupCounts {
  long long index_builds = 0;
  long long index_nodes = 0;
  long long active_taxon_entries = 0;
  long long active_index_entries = 0;
  long long selected_leaves = 0;
  long long lca_queries = 0;
  long long retained_nodes = 0;
  long long localized_nodes = 0;
  long long localized_calls = 0;

  long long aggregate_work() const {
    return index_nodes + active_taxon_entries + active_index_entries +
           selected_leaves + lca_queries + retained_nodes + localized_nodes +
           localized_calls;
  }
};

SetupCounts setup_counts() {
  const qraft::internal::RunStats& stats =
      qraft::internal::run_stats();
  auto load = [](const std::atomic<long long>& value) {
    return value.load(std::memory_order_relaxed);
  };
  return SetupCounts{
      load(stats.starting_tree_index_builds),
      load(stats.starting_tree_index_nodes),
      load(stats.active_taxa_entries),
      load(stats.active_taxa_index_entries),
      load(stats.query_tree_selected_leaves),
      load(stats.query_tree_lca_queries),
      load(stats.query_tree_retained_nodes),
      load(stats.gene_set_scoring_nodes),
      load(stats.gene_set_scoring_calls),
  };
}

void exercise_active_subset(const qraft::internal::StartingTree& starting_tree,
                            int root_taxa,
                            int begin,
                            int end) {
  using namespace qraft::internal;
  std::vector<int> taxa(static_cast<std::size_t>(end - begin));
  std::iota(taxa.begin(), taxa.end(), begin);
  const ActiveTaxa active = make_active_taxa(taxa, root_taxa);
  const QueryTree query =
      induce_query_tree(starting_tree, active);
  require(static_cast<int>(query.leaf_active_order.size()) == end - begin,
          "virtual query lost an active leaf");

  SolverContext context;
  context.genes.push_back(synthetic_local_gene(end - begin));
  context.local_taxa = std::make_shared<const std::vector<int>>(taxa);
  const OwnerProjection owners =
      make_owner_projection(context, active.active_by_taxon);
  require(owners.owner_count == end - begin,
          "restricted gene-set owner map changed active cardinality");
  for (int local = 0; local < end - begin; ++local) {
    require(owners.local_keys.find(local) == local,
            "restricted gene-set owner map changed local ownership");
  }
  require(context.genes.front().nodes.size() ==
              static_cast<std::size_t>(2 * (end - begin) - 1),
          "immutable gene-set fixture changed node cardinality");
  note_gene_set_use(context.genes);
}

struct ScalingSample {
  SetupCounts counts;
  long long calls = 0;
  long long active_mass = 0;
  double normalized = 0.0;
};

ScalingSample run_subset_family(int leaf_count, bool skew) {
  using namespace qraft::internal;
  StartingTree starting_tree = balanced_starting_tree(leaf_count);
  reset_run_stats();
  long long calls = 0;
  long long active_mass = 0;
  auto visit = [&](auto&& self, int begin, int end) -> void {
    ++calls;
    active_mass += end - begin;
    exercise_active_subset(starting_tree, leaf_count, begin, end);
    if (end - begin <= 1) return;
    const int middle = skew ? begin + 1 : begin + (end - begin) / 2;
    self(self, begin, middle);
    self(self, middle, end);
  };
  visit(visit, 0, leaf_count);

  const SetupCounts counts = setup_counts();
  const long long expected_calls = 2LL * leaf_count - 1;
  int log2_leaves = 0;
  for (int value = leaf_count; value > 1; value /= 2) ++log2_leaves;
  const long long expected_mass =
      skew ? static_cast<long long>(leaf_count) * (leaf_count + 1) / 2 +
                 leaf_count - 1
           : static_cast<long long>(leaf_count) * (log2_leaves + 1);
  require(calls == expected_calls && active_mass == expected_mass,
          "recursive subset enumeration does not match its closed form");
  require(counts.index_builds == 1 &&
              counts.index_nodes == 2LL * leaf_count - 1,
          "starting-tree index was not one-time linear work");
  require(counts.active_taxon_entries == active_mass &&
              counts.active_index_entries == active_mass,
          "active-local setup initialized entries outside active mass");
  require(counts.selected_leaves == active_mass &&
              counts.lca_queries == active_mass - calls,
          "virtual query primitive counts do not match selected leaves");
  require(counts.retained_nodes >= active_mass &&
              counts.retained_nodes <= 2 * active_mass - calls,
          "virtual query retained more than 2m-1 nodes");
  if (!skew) {
    require(counts.retained_nodes == 2 * active_mass - calls,
            "balanced clade queries did not retain exactly 2m-1 nodes");
  }
  require(counts.localized_calls == calls &&
              counts.localized_nodes == 2 * active_mass - calls,
          "localized scoring counters do not match copied gene nodes");
  require(counts.active_index_entries <
              static_cast<long long>(leaf_count) * calls,
          "active setup behaves like per-query root-slot initialization");

  const double envelope =
      skew ? static_cast<double>(leaf_count) * leaf_count
           : static_cast<double>(leaf_count) * (log2_leaves + 1);
  const double normalized =
      static_cast<double>(counts.aggregate_work()) / envelope;
  require_finite(normalized, "active-local normalized operation count");
  return ScalingSample{counts, calls, active_mass, normalized};
}

using Clusters = std::vector<std::vector<int>>;

Clusters full_scan_induced_clusters(
    const qraft::internal::StartingTree& starting_tree,
    const qraft::internal::ActiveTaxa& active) {
  std::vector<int> taxon_to_active(starting_tree.taxon_to_leaf.size(), -1);
  for (int active_id = 0; active_id < static_cast<int>(active.taxa.size());
       ++active_id) {
    const int taxon = active.taxa[static_cast<std::size_t>(active_id)];
    if (taxon >= 0 &&
        taxon < static_cast<int>(taxon_to_active.size())) {
      taxon_to_active[static_cast<std::size_t>(taxon)] = active_id;
    }
  }

  Clusters clusters;
  auto build = [&](auto&& self, int node) -> std::vector<int> {
    const auto& source = starting_tree.nodes[static_cast<std::size_t>(node)];
    if (source.taxon >= 0 &&
        source.taxon < static_cast<int>(taxon_to_active.size())) {
      const int active_id =
          taxon_to_active[static_cast<std::size_t>(source.taxon)];
      if (active_id >= 0) {
        clusters.push_back({active_id});
        return {active_id};
      }
    }
    std::vector<std::vector<int>> children;
    for (int child : source.children) {
      std::vector<int> induced = self(self, child);
      if (!induced.empty()) children.push_back(std::move(induced));
    }
    if (children.empty()) return {};
    if (children.size() == 1) return std::move(children.front());
    std::vector<int> combined;
    for (const std::vector<int>& child : children) {
      combined.insert(combined.end(), child.begin(), child.end());
    }
    std::sort(combined.begin(), combined.end());
    clusters.push_back(combined);
    return combined;
  };
  if (starting_tree.root >= 0) (void)build(build, starting_tree.root);
  std::sort(clusters.begin(), clusters.end());
  return clusters;
}

Clusters virtual_query_clusters(const qraft::internal::QueryTree& query) {
  Clusters clusters;
  clusters.reserve(query.nodes.size());
  for (int node = 0; node < static_cast<int>(query.nodes.size()); ++node) {
    const int begin =
        query.subtree_leaf_first[static_cast<std::size_t>(node)];
    const int end = begin +
                    query.subtree_leaf_count[static_cast<std::size_t>(node)];
    std::vector<int> cluster(
        query.leaf_active_order.begin() + static_cast<std::ptrdiff_t>(begin),
        query.leaf_active_order.begin() + static_cast<std::ptrdiff_t>(end));
    std::sort(cluster.begin(), cluster.end());
    clusters.push_back(std::move(cluster));
  }
  std::sort(clusters.begin(), clusters.end());
  return clusters;
}

void test_virtual_tree_parity_and_publication() {
  using namespace qraft::internal;
  std::mt19937 random(20260713u);
  int polytomy_trials = 0;
  for (int trial = 0; trial < 200; ++trial) {
    std::uniform_int_distribution<int> leaf_distribution(2, 64);
    const int leaves = leaf_distribution(random);
    StartingTree starting_tree = random_polytomy_starting_tree(leaves, random);
    if (std::any_of(starting_tree.nodes.begin(), starting_tree.nodes.end(),
                    [](const StartingTree::Node& node) {
                      return node.children.size() > 2;
                    })) {
      ++polytomy_trials;
    }
    std::vector<int> taxa(static_cast<std::size_t>(leaves));
    std::iota(taxa.begin(), taxa.end(), 0);
    std::shuffle(taxa.begin(), taxa.end(), random);
    std::uniform_int_distribution<int> active_distribution(1, leaves);
    taxa.resize(static_cast<std::size_t>(active_distribution(random)));
    const ActiveTaxa active =
        make_active_taxa(taxa, leaves);
    const QueryTree query =
        induce_query_tree(starting_tree, active);
    require(virtual_query_clusters(query) ==
                full_scan_induced_clusters(starting_tree, active),
            "virtual/full-scan induced-tree parity failed at seed 20260713 "
            "trial=" +
                std::to_string(trial));
  }
  require(polytomy_trials > 100,
          "random parity corpus did not exercise arbitrary polytomies");

  StartingTree retry = balanced_starting_tree(16);
  std::vector<int> retry_taxa(16);
  std::iota(retry_taxa.begin(), retry_taxa.end(), 0);
  const ActiveTaxa retry_active =
      make_active_taxa(retry_taxa, 16);
  retry.nodes.emplace_back();
  bool rejected = false;
  try {
    (void)induce_query_tree(retry, retry_active);
  } catch (const std::logic_error&) {
    rejected = true;
  }
  require(rejected, "malformed starting-tree index build did not fail");
  retry.nodes.pop_back();
  reset_run_stats();
  const QueryTree retried =
      induce_query_tree(retry, retry_active);
  require(!retried.nodes.empty() && setup_counts().index_builds == 1,
          "starting-tree index did not retry after a build exception");

  constexpr int kPublicationLeaves = 16384;
  constexpr int kReaders = 6;
  StartingTree shared = balanced_starting_tree(kPublicationLeaves);
  std::vector<int> all_taxa(kPublicationLeaves);
  std::iota(all_taxa.begin(), all_taxa.end(), 0);
  const ActiveTaxa all_active = make_active_taxa(
      all_taxa, kPublicationLeaves);
  reset_run_stats();
  std::barrier publication_start(kReaders + 1);
  std::vector<std::future<int>> readers;
  readers.reserve(kReaders);
  for (int reader = 0; reader < kReaders; ++reader) {
    readers.push_back(std::async(std::launch::async, [&] {
      publication_start.arrive_and_wait();
      const QueryTree query =
          induce_query_tree(shared, all_active);
      return static_cast<int>(query.nodes.size());
    }));
  }
  publication_start.arrive_and_wait();
  for (std::future<int>& reader : readers) {
    require(reader.wait_for(std::chrono::seconds(10)) ==
                std::future_status::ready,
            "starting-tree index publication reader timed out");
    require(reader.get() == 2 * kPublicationLeaves - 1,
            "concurrent starting-tree query retained an invalid tree");
  }
  const SetupCounts publication = setup_counts();
  require(publication.index_builds == 1 &&
              publication.index_nodes == 2LL * kPublicationLeaves - 1,
          "concurrent starting-tree queries did not share one published index");
}

std::string outer_taxon_label(int taxon) {
  std::string digits = std::to_string(taxon);
  return "T" + std::string(6 - digits.size(), '0') + digits;
}

qraft::GeneTree balanced_gene_tree(int leaf_count) {
  qraft::GeneTree tree;
  auto add_edge = [&](int first, int second) {
    tree.nodes[static_cast<std::size_t>(first)].neighbors.push_back(second);
    tree.nodes[static_cast<std::size_t>(second)].neighbors.push_back(first);
  };
  auto build = [&](auto&& self, int begin, int end) -> int {
    const int node = static_cast<int>(tree.nodes.size());
    tree.nodes.emplace_back();
    if (end - begin == 1) {
      tree.nodes[static_cast<std::size_t>(node)].label =
          outer_taxon_label(begin);
      return node;
    }
    const int middle = begin + (end - begin) / 2;
    const int left = self(self, begin, middle);
    const int right = self(self, middle, end);
    add_edge(node, left);
    add_edge(node, right);
    return node;
  };
  tree.root = build(build, 0, leaf_count);
  return tree;
}

qraft::GeneTree duplicate_heavy_gene_tree(int leaf_count,
                                         int repeated_pair_count) {
  require(leaf_count >= 4 && repeated_pair_count >= 2,
          "duplicate-heavy gene fixture has invalid dimensions");
  qraft::GeneTree tree;
  tree.nodes.resize(static_cast<std::size_t>(repeated_pair_count));
  auto add_edge = [&](int first, int second) {
    tree.nodes[static_cast<std::size_t>(first)].neighbors.push_back(second);
    tree.nodes[static_cast<std::size_t>(second)].neighbors.push_back(first);
  };
  auto add_leaf = [&](int center, int taxon) {
    const int leaf = static_cast<int>(tree.nodes.size());
    tree.nodes.emplace_back();
    tree.nodes.back().label = outer_taxon_label(taxon);
    add_edge(center, leaf);
  };

  for (int center = 1; center < repeated_pair_count; ++center) {
    add_edge(center - 1, center);
  }
  for (int center = 0; center < repeated_pair_count; ++center) {
    // Every corridor node sees both owners, so neither duplicate family is a
    // same-owner pendant component. Keeping T0/T1 together therefore keeps a
    // topology as large as the input gene even after most taxa are dropped.
    add_leaf(center, 0);
    add_leaf(center, 1);
  }
  for (int taxon = 2; taxon < leaf_count; ++taxon) {
    add_leaf((taxon - 2) % repeated_pair_count, taxon);
  }
  tree.root = 0;
  return tree;
}

long long estimated_gene_set_bytes(
    const std::vector<qraft::internal::GeneIndex>& genes) {
  using qraft::internal::GeneIndex;
  using qraft::internal::GeneTraversal;
  auto vector_bytes = [](const auto& values) {
    using Value = typename std::decay_t<decltype(values)>::value_type;
    return static_cast<long long>(values.capacity()) * sizeof(Value);
  };
  auto traversal_bytes = [&](const GeneTraversal& traversal) {
    return vector_bytes(traversal.parent) +
           vector_bytes(traversal.preorder);
  };

  long long bytes = vector_bytes(genes);
  for (const GeneIndex& gene : genes) {
    bytes += vector_bytes(gene.nodes) + vector_bytes(gene.taxon_to_node) +
             traversal_bytes(gene.root_traversal) +
             vector_bytes(gene.adjacency_offsets) +
             vector_bytes(gene.adjacency) + vector_bytes(gene.leaf_mass);
  }
  return bytes;
}

qraft::internal::SolverContext outer_solver_context(int leaf_count,
                                                    int gene_count,
                                                    int workers = 1) {
  using namespace qraft::internal;
  SolverContext context;
  std::vector<std::string> labels(static_cast<std::size_t>(leaf_count));
  for (int taxon = 0; taxon < leaf_count; ++taxon) {
    labels[static_cast<std::size_t>(taxon)] = outer_taxon_label(taxon);
  }
  context.taxa = qraft::TaxonRegistry::from_ordered_unique(std::move(labels));
  context.genes.reserve(static_cast<std::size_t>(gene_count));
  for (int gene = 0; gene < gene_count; ++gene) {
    context.genes.push_back(
        index_gene(balanced_gene_tree(leaf_count), context.taxa));
  }
  context.starting_tree = balanced_starting_tree(leaf_count);
  context.options.threads = workers;
  context.options.progress = qraft::ProgressMode::Off;
  context.options.nni = false;
  context.workers = workers;
  return context;
}

qraft::internal::SolverContext duplicate_heavy_outer_context(
    int leaf_count,
    int repeated_pair_count,
    int workers) {
  using namespace qraft::internal;
  SolverContext context;
  std::vector<std::string> labels(static_cast<std::size_t>(leaf_count));
  for (int taxon = 0; taxon < leaf_count; ++taxon) {
    labels[static_cast<std::size_t>(taxon)] = outer_taxon_label(taxon);
  }
  context.taxa = qraft::TaxonRegistry::from_ordered_unique(std::move(labels));
  context.genes.push_back(index_gene(
      duplicate_heavy_gene_tree(leaf_count, repeated_pair_count),
      context.taxa));
  context.starting_tree = balanced_starting_tree(leaf_count);
  context.options.threads = workers;
  context.options.progress = qraft::ProgressMode::Off;
  context.options.nni = false;
  context.workers = workers;
  return context;
}

// Every run counter, copied after a run.
struct OuterSnapshot {
#define QRAFT_SNAPSHOT_FIELD(name) long long name = 0;
  QRAFT_RUN_STATS(QRAFT_SNAPSHOT_FIELD)
#undef QRAFT_SNAPSHOT_FIELD
  long long complexity_contract_satisfied = 0;
};

OuterSnapshot outer_snapshot() {
  const qraft::internal::RunStats& stats = qraft::internal::run_stats();
  OuterSnapshot out;
#define QRAFT_SNAPSHOT_COPY(name) \
  out.name = stats.name.load(std::memory_order_relaxed);
  QRAFT_RUN_STATS(QRAFT_SNAPSHOT_COPY)
#undef QRAFT_SNAPSHOT_COPY
  out.complexity_contract_satisfied =
      stats.complexity_contract_satisfied.load(std::memory_order_relaxed);
  return out;
}

void require_outer_payload_released(const std::string& label) {
  const qraft::internal::RunStats& stats = qraft::internal::run_stats();
  auto load = [](const std::atomic<long long>& value) {
    return value.load(std::memory_order_relaxed);
  };
  const std::array<const std::atomic<long long>*, 9> live{
      &stats.live_payload,       &stats.live_taxa,
      &stats.live_gene_set_nodes, &stats.live_gene_set_entries,
      &stats.live_gene_set_taxa, &stats.live_continuation_taxa,
      &stats.live_tree_nodes,    &stats.live_tree_adjacencies,
      &stats.live_taxon_scratch,
  };
  for (const std::atomic<long long>* counter : live) {
    require(load(*counter) == 0, label + " left a payload live");
  }
  require(load(stats.live_child_gene_sets) == 0 &&
              load(stats.live_child_gene_set_slots) == 0,
          label + " left a child gene set live");
}

struct OuterRun {
  OuterSnapshot counts;
  int genes = 0;
  long long root_gene_set_bytes = 0;
  long long output_nodes = 0;
  long long output_adjacencies = 0;
};

OuterRun run_outer_context(
    qraft::internal::SolverContext context,
    qraft::internal::ForcedSplitShape shape) {
  using namespace qraft::internal;
  const int leaf_count = static_cast<int>(context.taxa.size());
  const int gene_count = static_cast<int>(context.genes.size());
  const long long root_gene_set_bytes =
      estimated_gene_set_bytes(context.genes);
  reset_run_stats();
  SolverTree tree = solve_qraft_for_test(context, shape);
  long long output_adjacencies = 0;
  for (const SolverTree::Node& node : tree.nodes) {
    output_adjacencies += static_cast<long long>(node.neighbors.size());
  }
  OuterRun run;
  run.counts = outer_snapshot();
  run.genes = gene_count;
  run.root_gene_set_bytes = root_gene_set_bytes;
  run.output_nodes = static_cast<long long>(tree.nodes.size());
  run.output_adjacencies = output_adjacencies;
  require_outer_payload_released("outer recursion n=" +
                                 std::to_string(leaf_count));
  return run;
}

OuterRun run_outer_shape(int leaf_count,
                         qraft::internal::ForcedSplitShape shape,
                         int gene_count,
                         int workers = 1) {
  return run_outer_context(
      outer_solver_context(leaf_count, gene_count, workers), shape);
}

struct OuterExpected {
  long long calls = 1;
  long long nonbase = 0;
  long long active_mass = 0;
  long long nonbase_mass = 0;
  long long tree_visits = 0;
  long long radix_mass = 0;
  long long radix_calls = 0;
  long long balanced_splits = 0;
  long long skew_splits = 0;
  long long max_depth = 0;
  long long output_nodes = 0;
  long long output_adjacencies = 0;
  long long continuation = 0;
  long long completed_nodes = 0;
  long long completed_adjacencies = 0;
  long long taxon_scratch = 0;
  long long regraft_scoring_calls = 0;
};

OuterExpected expected_outer_shape(int n, bool skew) {
  OuterExpected out;
  out.active_mass = n;
  out.output_nodes = n == 1 ? 1 : 2LL * n - 2;
  out.output_adjacencies = n == 1 ? 0 : 2 * out.output_nodes - 2;
  out.completed_nodes = out.output_nodes;
  out.completed_adjacencies = out.output_adjacencies;
  if (n <= 3) return out;

  const int left_size = skew ? 1 : n / 2;
  const int right_size = n - left_size;
  const OuterExpected left = expected_outer_shape(left_size, skew);
  const OuterExpected right = expected_outer_shape(right_size, skew);

  out.calls += left.calls + right.calls;
  out.nonbase = 1 + left.nonbase + right.nonbase;
  out.active_mass += left.active_mass + right.active_mass;
  out.nonbase_mass = n + left.nonbase_mass + right.nonbase_mass;
  out.tree_visits = left.tree_visits + right.tree_visits +
                    left.output_nodes + right.output_nodes;
  out.radix_mass = left.radix_mass + right.radix_mass +
                   (left_size > 1 ? left_size : 0) +
                   (right_size > 1 ? right_size : 0);
  out.radix_calls = left.radix_calls + right.radix_calls +
                    (left_size > 1 ? 1 : 0) +
                    (right_size > 1 ? 1 : 0);
  out.balanced_splits = left.balanced_splits + right.balanced_splits +
                        (left_size == right_size ? 1 : 0);
  out.skew_splits = left.skew_splits + right.skew_splits +
                    (left_size == 1 && right_size == n - 1 ? 1 : 0);
  out.max_depth = 1 + std::max(left.max_depth, right.max_depth);
  const long long connector_nodes =
      (left_size > 1 ? 1LL : 0LL) + (right_size > 1 ? 1LL : 0LL);
  const long long connector_adjacencies = 2 + 2 * connector_nodes;
  out.output_nodes =
      left.output_nodes + right.output_nodes + connector_nodes;
  out.output_adjacencies = left.output_adjacencies +
                           right.output_adjacencies +
                           connector_adjacencies;

  // Solved serially, the right side's taxa stay live while the left side is
  // solved, then the completed left tree stays live while the right side is.
  out.continuation =
      std::max<long long>(right_size + left.continuation,
                          right.continuation);
  out.completed_nodes = std::max(
      {left.completed_nodes,
       left.output_nodes + right.completed_nodes,
       left.output_nodes + 2 * right.output_nodes + connector_nodes,
       out.output_nodes});
  out.completed_adjacencies = std::max(
      {left.completed_adjacencies,
       left.output_adjacencies + right.completed_adjacencies,
       left.output_adjacencies + 2 * right.output_adjacencies +
           connector_adjacencies,
       out.output_adjacencies});

  const long long left_local_scratch =
      left_size > 1 ? 2LL * left_size : left_size;
  const long long right_local_scratch =
      left_size + (right_size > 1 ? 2LL * right_size : right_size);
  out.taxon_scratch =
      std::max({left.taxon_scratch, right.taxon_scratch,
                left_local_scratch, right_local_scratch,
                static_cast<long long>(n)});
  out.regraft_scoring_calls =
      left.regraft_scoring_calls + right.regraft_scoring_calls +
      (left_size >= 3 || right_size >= 3 ? 1 : 0);
  return out;
}

void require_outer_formula(const OuterRun& run, int n, bool skew) {
  const OuterExpected expected = expected_outer_shape(n, skew);
  const std::string label = std::string(skew ? "skew" : "balanced") +
                            " n=" + std::to_string(n) + " k=" +
                            std::to_string(run.genes);
  const OuterSnapshot& counts = run.counts;

  require(counts.starting_tree_index_builds == 1 &&
              counts.starting_tree_index_nodes == 2LL * n - 1,
          label + " did not use one starting-tree index");
  require(counts.query_tree_selected_leaves == expected.nonbase_mass &&
              counts.query_tree_lca_queries ==
                  expected.nonbase_mass - expected.nonbase,
          label + " virtual-tree setup differs from active mass");
  require(counts.query_tree_retained_nodes >= expected.nonbase_mass &&
              counts.query_tree_retained_nodes <=
                  2 * expected.nonbase_mass - expected.nonbase &&
              (skew || counts.query_tree_retained_nodes ==
                           2 * expected.nonbase_mass - expected.nonbase),
          label + " virtual-tree retention violates the 2m-1 invariant");
  require(counts.recursion_calls == expected.calls &&
              counts.recursion_nonbase_calls == expected.nonbase &&
              counts.recursion_max_depth == expected.max_depth,
          label + " recursive shape differs from its closed form");
  require(counts.recursion_taxon_mass == expected.active_mass,
          label + " active mass differs from its closed form");
  require(counts.recursion_split_tests == expected.nonbase_mass &&
              counts.recursion_partition_moves == expected.nonbase_mass,
          label + " split/partition work differs from nonbase mass");
  require(counts.recursion_balanced_splits == expected.balanced_splits &&
              counts.recursion_skewed_splits == expected.skew_splits,
          label + " did not force the requested bipartition shape");
  require(counts.recursion_tree_node_visits == expected.tree_visits,
          label + " completed-tree reconstruction work is not exact");
  require(counts.recursion_radix_entries == 8 * expected.radix_mass &&
              counts.recursion_radix_bucket_ops == 1024 * expected.radix_calls &&
              counts.recursion_unique_taxa == expected.nonbase_mass,
          label + " taxon ordering work differs from fixed-pass formulas");
  require(run.output_nodes == expected.output_nodes &&
              run.output_adjacencies == expected.output_adjacencies,
          label + " output tree payload differs from its recurrence");

  const long long expected_scoring_calls =
      expected.nonbase + expected.regraft_scoring_calls;
  require(counts.candidate_sets == expected.nonbase &&
              counts.bipartitions_scored >= expected.nonbase,
          label + " did not execute the bipartition scorer");
  require(counts.regraft_calls == 2 * expected.nonbase &&
              counts.regraft_vertices > 0,
          label + " did not execute regraft and vertex scoring");
  require(counts.completed_subproblems == expected.nonbase,
          label + " completed-subproblem count differs from nonbase recursion");
  require(counts.hdt_builds == run.genes * expected_scoring_calls,
          label + " HDT builds do not scale exactly with k");
  if (run.genes == 0) {
    require(counts.hdt_queries == 0,
            label + " gene-free scoring performed per-gene HDT work");
  } else {
    require(counts.hdt_queries > 0,
            label + " gene-bearing run bypassed HDT work");
    require(counts.worker_total_peak_bytes <=
                    counts.extract_scratch_peak_bytes +
                        counts.contract_scratch_peak_bytes +
                        counts.state_arena_peak_bytes,
            label + " scorer workspace is not fully bounded");
    require(counts.max_scoring_source_taxa <= n &&
                counts.scoring_source_inflation_violations == 0,
            label + " scoring genes exceeded their taxon bound");
  }

  require(counts.peak_taxa == n,
          label + " taxon lists are not handed down exactly once");
  require(counts.peak_continuation_taxa == expected.continuation,
          label + " pending-sibling payload differs from the DFS recurrence");
  require(counts.peak_tree_nodes == expected.completed_nodes &&
              counts.peak_tree_adjacencies ==
                  expected.completed_adjacencies,
          label + " completed-tree payload differs from exact recurrence");
  require(counts.peak_taxon_scratch == expected.taxon_scratch,
          label + " taxon reconstruction scratch differs from exact peak");
  require(counts.accounting_violations == 0,
          label + " live accounting underflowed");

  if (run.genes == 0) {
    require(counts.recursion_root_gene_set_nodes == 0 &&
                counts.recursion_root_gene_set_entries == 0 &&
                counts.peak_gene_set_nodes == 0 &&
                counts.peak_gene_set_entries == 0 &&
                counts.peak_gene_set_taxa == 0,
            label + " empty-gene run acquired gene-set payload");
  }
  const long long sum_of_category_peaks =
      counts.peak_taxa + counts.peak_gene_set_nodes +
      counts.peak_gene_set_entries + counts.peak_gene_set_taxa +
      counts.peak_tree_nodes + counts.peak_tree_adjacencies +
      counts.peak_taxon_scratch;
  require(counts.peak_payload <= sum_of_category_peaks,
          label + " total peak exceeds the sum of disjoint category peaks");
  require(counts.peak_payload >=
              n + counts.recursion_root_gene_set_nodes +
                  counts.recursion_root_gene_set_entries,
          label + " total peak missed the root-owned payload");
}

void require_gene_set_bounds(const OuterRun& run,
                                    int leaves,
                                    bool skew) {
  const std::string label = std::string(skew ? "skew" : "balanced") +
                            " n=" + std::to_string(leaves) + " k=" +
                            std::to_string(run.genes);
  const OuterExpected expected = expected_outer_shape(leaves, skew);
  const OuterSnapshot& counts = run.counts;
  const long long root_nodes =
      static_cast<long long>(run.genes) * (2LL * leaves - 2);
  // Per gene: one header entry, n taxon slots, parent/preorder of the root
  // traversal (2(2n-2); the reference traversal is rebuilt on demand),
  // 2(2n-3) CSR neighbor entries and 2n-1 CSR offsets: 11(n-1)+1.
  const long long root_other =
      static_cast<long long>(run.genes) * (11LL * (leaves - 1) + 1LL);

  require(counts.recursion_root_gene_set_nodes == root_nodes,
          label + " root gene-set nodes do not equal k(2n-2)");
  require(counts.recursion_root_gene_set_entries == root_other,
          label + " root gene-set payload does not equal k(11(n-1)+1)");
  require(counts.child_gene_set_accepts > 0 && counts.child_gene_set_reuses > 0 &&
              counts.child_gene_set_accepts + counts.child_gene_set_reuses ==
                  2 * expected.nonbase,
          label + " does not cover every child gene-set keep/reuse choice");
  require(counts.child_gene_set_attempts ==
                  counts.child_gene_set_accepts +
                      counts.child_gene_set_heavy_reuses &&
              counts.child_gene_set_pair_violations == 0 &&
              counts.child_gene_set_max_accepted_ratio_ppm <= 750000,
          label + " violated the physical child gene-set gate");
  require(counts.max_child_gene_sets > 0 &&
              counts.max_child_gene_set_slots > 0 &&
              counts.max_child_gene_set_slots <= 4LL * leaves,
          label + " serial child gene-set retention exceeds its budget");
  require(counts.peak_gene_set_nodes > root_nodes &&
              counts.peak_gene_set_nodes <=
                  2 * root_nodes +
                      2LL * run.genes * counts.max_child_gene_set_slots,
          label + " child construction/retention exceeds its physical bound");
  require(counts.peak_gene_set_entries > root_other &&
              counts.peak_gene_set_entries <=
                  2 * root_other +
                      24LL * run.genes * counts.max_child_gene_set_slots,
          label + " child construction/retention structural payload is unbounded");
  require(counts.peak_gene_set_taxa <= counts.max_child_gene_set_slots,
          label + " child gene-set taxa exceed one entry per slot");

  // The capacity counter covers the vectors of a retained gene set after
  // construction. The temporary index maps and restriction workspaces stay
  // inside build_child_gene_set and are checked only by the RSS smoke below.
  require(counts.max_child_gene_set_bytes > 0 &&
              counts.max_child_gene_set_bytes <= run.root_gene_set_bytes,
          label + " retained child gene-set capacity exceeds the root set");
}

void require_exact_gene_scaling(const OuterRun& one,
                                const OuterRun& two,
                                int leaves,
                                bool skew) {
  const std::string label = std::string(skew ? "skew" : "balanced") +
                            " n=" + std::to_string(leaves);
  const OuterSnapshot& k1 = one.counts;
  const OuterSnapshot& k2 = two.counts;
  require(one.genes == 1 && two.genes == 2,
          label + " malformed k-scaling fixture");
  require(k2.recursion_root_gene_set_nodes == 2 * k1.recursion_root_gene_set_nodes &&
              k2.recursion_root_gene_set_entries ==
                  2 * k1.recursion_root_gene_set_entries &&
              k2.peak_gene_set_nodes == 2 * k1.peak_gene_set_nodes &&
              k2.peak_gene_set_entries ==
                  2 * k1.peak_gene_set_entries,
          label + " root/child gene-set payload is not exactly linear in k");
  require(k2.peak_gene_set_nodes - k2.recursion_root_gene_set_nodes ==
                  2 * (k1.peak_gene_set_nodes - k1.recursion_root_gene_set_nodes) &&
              k2.peak_gene_set_entries - k2.recursion_root_gene_set_entries ==
                  2 * (k1.peak_gene_set_entries -
                       k1.recursion_root_gene_set_entries),
          label + " materialized child payload is not exactly linear in k");
  require(k2.hdt_builds == 2 * k1.hdt_builds &&
              k2.hdt_queries == 2 * k1.hdt_queries &&
              k2.hdt_recolor_ops == 2 * k1.hdt_recolor_ops &&
              k2.hdt_nodes_recomputed == 2 * k1.hdt_nodes_recomputed,
          label + " bipartition/regraft work is not exactly linear in k");
  require(k2.child_gene_set_accepts == k1.child_gene_set_accepts &&
              k2.child_gene_set_reuses == k1.child_gene_set_reuses &&
              k2.max_child_gene_set_slots == k1.max_child_gene_set_slots &&
              k2.max_child_gene_sets == k1.max_child_gene_sets &&
              k2.peak_gene_set_taxa ==
                  k1.peak_gene_set_taxa &&
              k2.worker_total_peak_bytes == k1.worker_total_peak_bytes,
          label + " gene-set policy changed with duplicate gene count");
  require(two.root_gene_set_bytes == 2 * one.root_gene_set_bytes,
          label + " root gene-set retained capacity is not linear in k");
}

void test_gene_set_accounting() {
  using namespace qraft::internal;
  for (int leaves : {9, 17}) {
    for (ForcedSplitShape shape :
         {ForcedSplitShape::Balanced,
          ForcedSplitShape::MaximallySkewed}) {
      const bool skew = shape == ForcedSplitShape::MaximallySkewed;
      const OuterRun one = run_outer_shape(leaves, shape, 1);
      const OuterRun two = run_outer_shape(leaves, shape, 2);
      require_outer_formula(one, leaves, skew);
      require_outer_formula(two, leaves, skew);
      require_gene_set_bounds(one, leaves, skew);
      require_gene_set_bounds(two, leaves, skew);
      require_exact_gene_scaling(one, two, leaves, skew);
      std::cout << "outer-gene-exact n=" << leaves
                << " shape=" << (skew ? "skew" : "balanced")
                << " k1_hdt=" << one.counts.hdt_queries
                << " k2_hdt=" << two.counts.hdt_queries
                << " k1_child_nodes="
                << one.counts.peak_gene_set_nodes -
                       one.counts.recursion_root_gene_set_nodes
                << " k2_child_nodes="
                << two.counts.peak_gene_set_nodes -
                       two.counts.recursion_root_gene_set_nodes
                << '\n';
    }
  }
}

void test_duplicate_heavy_gene_set_gate() {
  using namespace qraft::internal;
  constexpr int kLeaves = 65;
  constexpr int kRepeatedPairs = 8 * kLeaves;
  constexpr int kWorkers = 4;
  qraft::Options parallel_options;
  parallel_options.threads = kWorkers;
  require(configure_threads(parallel_options) == kWorkers,
          "failed to configure duplicate-heavy parallel gate fixture");
  const OuterRun run = run_outer_context(
      duplicate_heavy_outer_context(kLeaves, kRepeatedPairs, kWorkers),
      ForcedSplitShape::Balanced);
  const OuterSnapshot& counts = run.counts;

  require(counts.child_gene_set_attempts > 0 &&
              counts.child_gene_set_accepts > 0 &&
              counts.child_gene_set_heavy_reuses > 0,
          "duplicate-heavy fixture did not exercise both physical gate paths");
  require(counts.child_gene_set_attempts ==
                  counts.child_gene_set_accepts +
                      counts.child_gene_set_heavy_reuses,
          "duplicate-heavy fixture lost a child gene-set gate decision");
  require(counts.child_gene_set_pair_violations == 0 &&
              counts.child_gene_set_max_accepted_ratio_ppm <= 750000 &&
              counts.complexity_contract_satisfied == 1,
          "duplicate-heavy fixture violated the physical state gate");
  require(counts.scoring_source_inflation_violations > 0,
          "duplicate-heavy fixture did not keep scoring genes on more than "
          "twice the active taxa");
  require(counts.parallel_child_size_rejections > 0 &&
              counts.max_parallel_child_slots <= 4LL * kLeaves,
          "duplicate-heavy genes bypassed the structural-parallel gate");
  require(counts.max_child_gene_set_slots <= 4LL * kLeaves,
          "duplicate-heavy gene-set metadata exceeded the slot budget");
  require(counts.peak_gene_set_nodes <= 5 * counts.recursion_root_gene_set_nodes,
          "duplicate-heavy live gene-set state exceeded a linear envelope");
  require(counts.accounting_violations == 0,
          "duplicate-heavy physical gene-set accounting underflowed");

  std::cout << "outer-duplicate-heavy n=" << kLeaves
            << " root_nodes=" << counts.recursion_root_gene_set_nodes
            << " peak_nodes=" << counts.peak_gene_set_nodes
            << " attempts=" << counts.child_gene_set_attempts
            << " accepts=" << counts.child_gene_set_accepts
            << " heavy=" << counts.child_gene_set_heavy_reuses
            << " max_ratio_ppm="
            << counts.child_gene_set_max_accepted_ratio_ppm
            << " physical_parallel_rejections="
            << counts.parallel_child_size_rejections << '\n';

  qraft::Options serial_options;
  serial_options.threads = 1;
  require(configure_threads(serial_options) == 1,
          "failed to restore serial duplicate-heavy gate fixture");
}

void test_outer_structural_parallel_recursion() {
  using namespace qraft::internal;
  qraft::Options parallel_options;
  parallel_options.threads = 4;
  require(configure_threads(parallel_options) == 4,
          "failed to configure structural-parallel outer driver");

  constexpr int kLeaves = 65;
  constexpr int kGenes = 2;
  const OuterExpected expected = expected_outer_shape(kLeaves, false);
  const OuterRun run = run_outer_shape(
      kLeaves, ForcedSplitShape::Balanced, kGenes, 4);
  const OuterSnapshot& counts = run.counts;
  require(counts.parallel_children > 0 &&
              counts.serial_children +
                      counts.parallel_children ==
                  expected.nonbase,
          "controlled driver did not enter structural parallel recursion");
  require(counts.max_parallel_child_slots > 0 &&
              counts.max_parallel_child_slots <= 4LL * kLeaves &&
              counts.max_child_gene_set_slots <= 4LL * kLeaves,
          "structural parallel slot ownership exceeded the configured budget");
  require(counts.recursion_calls == expected.calls &&
              counts.recursion_nonbase_calls == expected.nonbase &&
              counts.recursion_taxon_mass == expected.active_mass &&
              counts.recursion_partition_moves == expected.nonbase_mass &&
              counts.completed_subproblems == expected.nonbase &&
              run.output_nodes == expected.output_nodes &&
              run.output_adjacencies == expected.output_adjacencies,
          "parallel structural recursion changed outer work or connection");
  require(counts.regraft_calls == 2 * expected.nonbase &&
              counts.hdt_builds ==
                  kGenes *
                      (expected.nonbase + expected.regraft_scoring_calls) &&
              counts.child_gene_set_accepts + counts.child_gene_set_reuses ==
                  2 * expected.nonbase,
          "parallel recursion bypassed scorer or gene-set ownership");
  require(counts.accounting_violations == 0,
          "parallel structural recursion violated logical live accounting");

  qraft::Options serial_options;
  serial_options.threads = 1;
  require(configure_threads(serial_options) == 1,
          "failed to restore serial outer-recursion driver");
}

// Makes every subproblem at one depth throw, serially and on four workers.
// The failure of the leftmost one must reach the caller whatever order the
// parallel children finish in, and unwinding the frames must release every
// payload and the child budget.
void test_recursion_failure_unwinding() {
  using namespace qraft::internal;
  struct FailureCase {
    int leaves;
    ForcedSplitShape shape;
    int genes;
    int fail_depth;
  };
  const std::array<FailureCase, 4> cases{{
      {65, ForcedSplitShape::Balanced, 2, 0},
      {65, ForcedSplitShape::Balanced, 2, 2},
      {600, ForcedSplitShape::MaximallySkewed, 0, 500},
      {600, ForcedSplitShape::MaximallySkewed, 2, 64},
  }};
  const auto run = [](const FailureCase& failure, int workers) {
    const std::string label =
        "failure at depth " + std::to_string(failure.fail_depth) + " n=" +
        std::to_string(failure.leaves) + " k=" +
        std::to_string(failure.genes) + " workers=" + std::to_string(workers);
    SolverContext context =
        outer_solver_context(failure.leaves, failure.genes, workers);
    const std::shared_ptr<ChildSetBudget> budget = context.child_budget;
    reset_run_stats();
    std::string message;
    try {
      (void)solve_qraft_for_test(context, failure.shape, failure.fail_depth);
    } catch (const std::runtime_error& error) {
      message = error.what();
    }
    require(message.rfind("injected recursion failure at path ", 0) == 0,
            label + " ended with '" + message + "'");
    require_outer_payload_released(label);
    require(budget->reserved_slots.load() == 0 &&
                budget->live_slots.load() == 0 &&
                budget->live_sets.load() == 0,
            label + " left the child budget held");
    return message;
  };

  std::vector<std::string> serial;
  for (const FailureCase& failure : cases) serial.push_back(run(failure, 1));

  qraft::Options parallel_options;
  parallel_options.threads = 4;
  require(configure_threads(parallel_options) == 4,
          "failed to configure four workers");
  for (int round = 0; round < 20; ++round) {
    for (std::size_t i = 0; i < cases.size(); ++i) {
      require(run(cases[i], 4) == serial[i],
              "a parallel recursion rethrew a different subproblem's "
              "failure from the serial one");
    }
  }
  qraft::Options serial_options;
  serial_options.threads = 1;
  require(configure_threads(serial_options) == 1,
          "failed to restore one worker");
}

// A maximally skewed recursion is as deep as it has taxa. The divide and
// regraft keeps its subproblems on the heap, so the serial recursion runs here
// on a thread with a 512 KiB stack, and the parallel one on executor workers
// with such stacks; a call per level overflowed both.
void test_deep_skewed_recursion() {
  using namespace qraft::internal;
  constexpr std::size_t kStack = 512 * 1024;
  const std::array<std::pair<int, int>, 2> shapes{{{2000, 0}, {600, 2}}};
  for (const auto& [leaves, genes] : shapes) {
    OuterRun serial;
    qraft::test::run_with_stack(kStack, [&] {
      serial = run_outer_shape(leaves, ForcedSplitShape::MaximallySkewed,
                               genes, 1);
    });
    require_outer_formula(serial, leaves, true);
    if (genes > 0) require_gene_set_bounds(serial, leaves, true);
    std::cout << "deep-skew-serial n=" << leaves << " k=" << genes
              << " depth=" << serial.counts.recursion_max_depth << '\n';
  }

#if defined(__GLIBC__)
  require(qraft::test::set_default_thread_stack(kStack),
          "failed to set the default thread stack");
  qraft::Options parallel_options;
  parallel_options.threads = 4;
  require(configure_threads(parallel_options) == 4,
          "failed to configure four small-stack workers");
  for (const auto& [leaves, genes] : shapes) {
    const OuterExpected expected = expected_outer_shape(leaves, true);
    OuterRun parallel;
    qraft::test::run_with_stack(kStack, [&] {
      parallel = run_outer_shape(leaves, ForcedSplitShape::MaximallySkewed,
                                 genes, 4);
    });
    const OuterSnapshot& counts = parallel.counts;
    const std::string label = "deep parallel skew n=" +
                              std::to_string(leaves) + " k=" +
                              std::to_string(genes);
    require(counts.parallel_children > 0 &&
                counts.serial_children + counts.parallel_children ==
                    expected.nonbase,
            label + " did not solve its children in parallel");
    require(counts.recursion_calls == expected.calls &&
                counts.recursion_nonbase_calls == expected.nonbase &&
                counts.recursion_max_depth == expected.max_depth &&
                counts.recursion_taxon_mass == expected.active_mass &&
                counts.recursion_partition_moves == expected.nonbase_mass &&
                counts.completed_subproblems == expected.nonbase &&
                parallel.output_nodes == expected.output_nodes &&
                parallel.output_adjacencies == expected.output_adjacencies,
            label + " changed the recursion's work or its output");
    require(counts.regraft_calls == 2 * expected.nonbase &&
                counts.child_gene_set_accepts +
                        counts.child_gene_set_reuses ==
                    2 * expected.nonbase,
            label + " bypassed regraft or gene-set ownership");
    require(counts.accounting_violations == 0,
            label + " violated the live accounting");
    std::cout << "deep-skew-parallel n=" << leaves << " k=" << genes
              << " depth=" << counts.recursion_max_depth
              << " parallel_children=" << counts.parallel_children << '\n';
  }
  qraft::Options serial_options;
  serial_options.threads = 1;
  require(configure_threads(serial_options) == 1,
          "failed to restore one worker");
#endif
}

#if defined(__linux__)
struct RssProbe {
  long long rss_kib = 0;
  long long logical_peak = 0;
  long long scorer_peak_bytes = 0;
  long long child_gene_set_peak_bytes = 0;
};

RssProbe run_outer_rss_child(int leaves,
                             qraft::internal::ForcedSplitShape shape,
                             bool solve,
                             int gene_count = 0) {
  struct ChildReport {
    int ok = 0;
    long long logical_peak = 0;
    long long scorer_peak_bytes = 0;
    long long child_gene_set_peak_bytes = 0;
  };
  int pipe_fds[2];
  require(pipe(pipe_fds) == 0, "failed to create RSS probe pipe");
  const pid_t child = fork();
  require(child >= 0, "failed to fork RSS probe");
  if (child == 0) {
    close(pipe_fds[0]);
    ChildReport report;
    try {
      qraft::internal::SolverContext context =
          outer_solver_context(leaves, gene_count);
      if (solve) {
        qraft::Options options;
        options.threads = 1;
        (void)qraft::internal::configure_threads(options);
        qraft::internal::reset_run_stats();
        qraft::internal::SolverTree tree =
            qraft::internal::solve_qraft_for_test(context, shape);
        report.logical_peak =
            qraft::internal::run_stats()
                .peak_payload.load(std::memory_order_relaxed);
        report.scorer_peak_bytes =
            qraft::internal::run_stats()
                .worker_total_peak_bytes.load(std::memory_order_relaxed);
        report.child_gene_set_peak_bytes =
            qraft::internal::run_stats()
                .max_child_gene_set_bytes.load(
                    std::memory_order_relaxed);
        report.ok =
            !tree.nodes.empty() &&
            qraft::internal::run_stats()
                    .live_payload.load(
                        std::memory_order_relaxed) == 0 &&
            qraft::internal::run_stats()
                    .accounting_violations.load(
                        std::memory_order_relaxed) == 0;
      } else {
        std::vector<int> taxa(static_cast<std::size_t>(context.taxa.size()));
        std::iota(taxa.begin(), taxa.end(), 0);
        long long touched = static_cast<long long>(context.starting_tree.nodes.size());
        for (int taxon : taxa) touched += taxon;
        report.logical_peak = touched;
        report.ok = touched > 0;
      }
    } catch (...) {
      report.ok = 0;
    }
    const char* data = reinterpret_cast<const char*>(&report);
    std::size_t written = 0;
    while (written < sizeof(report)) {
      const ssize_t count =
          write(pipe_fds[1], data + written, sizeof(report) - written);
      if (count <= 0) break;
      written += static_cast<std::size_t>(count);
    }
    close(pipe_fds[1]);
    _exit(report.ok ? 0 : 1);
  }

  close(pipe_fds[1]);
  ChildReport report;
  char* data = reinterpret_cast<char*>(&report);
  std::size_t received = 0;
  while (received < sizeof(report)) {
    const ssize_t count =
        read(pipe_fds[0], data + received, sizeof(report) - received);
    if (count <= 0) break;
    received += static_cast<std::size_t>(count);
  }
  close(pipe_fds[0]);
  int status = 0;
  struct rusage usage {};
  require(wait4(child, &status, 0, &usage) == child,
          "failed to reap RSS probe child");
  require(received == sizeof(report) && WIFEXITED(status) &&
              WEXITSTATUS(status) == 0 && report.ok == 1,
          "RSS probe child failed");
  return RssProbe{static_cast<long long>(usage.ru_maxrss),
                  report.logical_peak, report.scorer_peak_bytes,
                  report.child_gene_set_peak_bytes};
}

void test_outer_physical_high_water_smoke() {
  using qraft::internal::ForcedSplitShape;
#if defined(__SANITIZE_ADDRESS__) || defined(__SANITIZE_THREAD__)
  // Sanitizer runtimes add process-level allocator and shadow memory that is
  // unrelated to the logical payload. Keep the smoke bounded but use enough
  // fixed slack that sanitizer metadata is not treated as payload.
  constexpr int kSmallLeaves = 32;
  constexpr int kLargeLeaves = 64;
  constexpr int kGeneSetLeaves = 17;
  constexpr bool kProbeGeneSets = false;
  constexpr long long kRuntimeSlackKib = 20LL * 1024;
  constexpr long long kScalingNoiseKib = 4LL * 1024;
#else
  constexpr int kSmallLeaves = 128;
  constexpr int kLargeLeaves = 256;
  constexpr int kGeneSetLeaves = 64;
  constexpr bool kProbeGeneSets = true;
  constexpr long long kRuntimeSlackKib = 4LL * 1024;
  constexpr long long kScalingNoiseKib = 1LL * 1024;
#endif
  const RssProbe baseline_small = run_outer_rss_child(
      kSmallLeaves, ForcedSplitShape::Balanced, false);
  const RssProbe baseline_large = run_outer_rss_child(
      kLargeLeaves, ForcedSplitShape::Balanced, false);
  const RssProbe balanced_large = run_outer_rss_child(
      kLargeLeaves, ForcedSplitShape::Balanced, true);
  const RssProbe skew_small = run_outer_rss_child(
      kSmallLeaves, ForcedSplitShape::MaximallySkewed, true);
  const RssProbe skew_large = run_outer_rss_child(
      kLargeLeaves, ForcedSplitShape::MaximallySkewed, true);
  RssProbe gene_set_baseline;
  RssProbe gene_set_balanced;
  RssProbe gene_set_skew;
  if constexpr (kProbeGeneSets) {
    gene_set_baseline = run_outer_rss_child(
        kGeneSetLeaves, ForcedSplitShape::Balanced, false, 2);
    gene_set_balanced = run_outer_rss_child(
        kGeneSetLeaves, ForcedSplitShape::Balanced, true, 2);
    gene_set_skew = run_outer_rss_child(
        kGeneSetLeaves, ForcedSplitShape::MaximallySkewed, true, 2);
  }

  auto require_empirical_bound = [=](const RssProbe& baseline,
                                     const RssProbe& workload,
                                     const std::string& label) {
    const long long logical_allowance_kib =
        (128 * workload.logical_peak + workload.scorer_peak_bytes +
         workload.child_gene_set_peak_bytes + 1023) /
        1024;
    const long long allowed =
        baseline.rss_kib + kRuntimeSlackKib + logical_allowance_kib;
    require(workload.rss_kib > 0 &&
                workload.rss_kib <= allowed,
            label + " allocator high-water smoke exceeded its bounded "
                    "runtime/payload allowance: baseline_kib=" +
                std::to_string(baseline.rss_kib) + " measured_kib=" +
                std::to_string(workload.rss_kib) + " allowed_kib=" +
                std::to_string(allowed));
  };
  require_empirical_bound(baseline_large, balanced_large,
                          "balanced n=" + std::to_string(kLargeLeaves));
  require_empirical_bound(baseline_small, skew_small,
                          "skew n=" + std::to_string(kSmallLeaves));
  require_empirical_bound(baseline_large, skew_large,
                          "skew n=" + std::to_string(kLargeLeaves));
  if constexpr (kProbeGeneSets) {
    require_empirical_bound(gene_set_baseline, gene_set_balanced,
                            "balanced two-gene set n=" +
                                std::to_string(kGeneSetLeaves));
    require_empirical_bound(gene_set_baseline, gene_set_skew,
                            "skew two-gene set n=" +
                                std::to_string(kGeneSetLeaves));
  }
  const long long skew_delta_small =
      std::max(0LL, skew_small.rss_kib - baseline_small.rss_kib);
  const long long skew_delta_large =
      std::max(0LL, skew_large.rss_kib - baseline_large.rss_kib);
  require(skew_delta_large <= 2 * skew_delta_small + kScalingNoiseKib,
          "allocator high-water smoke shows gross skew retention");
  std::cout << "outer-rss-smoke small_n=" << kSmallLeaves
            << " baseline_small_kib=" << baseline_small.rss_kib
            << " skew_small_kib=" << skew_small.rss_kib
            << " large_n=" << kLargeLeaves
            << " baseline_large_kib=" << baseline_large.rss_kib
            << " balanced_large_kib=" << balanced_large.rss_kib
            << " skew_large_kib=" << skew_large.rss_kib
            << " gene_set_n=" << kGeneSetLeaves
            << " gene_set_balanced_kib=" << gene_set_balanced.rss_kib
            << " gene_set_skew_kib=" << gene_set_skew.rss_kib << '\n';
}
#else
void test_outer_physical_high_water_smoke() {}
#endif

void test_recursion_work_and_live_space() {
  using namespace qraft::internal;
  test_outer_physical_high_water_smoke();

  qraft::Options thread_options;
  thread_options.threads = 1;
  require(configure_threads(thread_options) == 1,
          "failed to configure serial outer-recursion driver");
  for (int leaves : {8, 16, 32, 64, 128, 256}) {
    const OuterRun balanced =
        run_outer_shape(leaves, ForcedSplitShape::Balanced, 0);
    const OuterRun skew = run_outer_shape(
        leaves, ForcedSplitShape::MaximallySkewed, 0);
    require_outer_formula(balanced, leaves, false);
    require_outer_formula(skew, leaves, true);
    std::cout << "outer-exact n=" << leaves
              << " balanced_mass=" << balanced.counts.recursion_taxon_mass
              << " balanced_peak=" << balanced.counts.peak_payload
              << " skew_mass=" << skew.counts.recursion_taxon_mass
              << " skew_peak=" << skew.counts.peak_payload << '\n';
  }
  test_gene_set_accounting();
  test_duplicate_heavy_gene_set_gate();
  test_outer_structural_parallel_recursion();
}

void test_recursion_complexity() {
  std::vector<ScalingSample> balanced;
  std::vector<ScalingSample> skew;
  for (int leaves : {16, 32, 64, 128, 256}) {
    balanced.push_back(run_subset_family(leaves, false));
    skew.push_back(run_subset_family(leaves, true));
    const double balanced_ratio =
        balanced.size() == 1
            ? 0.0
            : static_cast<double>(balanced.back().counts.aggregate_work()) /
                  balanced[balanced.size() - 2].counts.aggregate_work();
    const double skew_ratio =
        skew.size() == 1
            ? 0.0
            : static_cast<double>(skew.back().counts.aggregate_work()) /
                  skew[skew.size() - 2].counts.aggregate_work();
    require_finite(balanced.back().normalized,
                   "balanced active-local normalized work");
    require_finite(skew.back().normalized,
                   "skew active-local normalized work");
    if (balanced.size() > 1) {
      require_finite(balanced_ratio, "balanced active-local scaling ratio");
      require_finite(skew_ratio, "skew active-local scaling ratio");
      require(balanced_ratio <= 2.75,
              "balanced aggregate grows faster than n log n envelope");
      require(skew_ratio <= 4.5,
              "maximally skewed aggregate grows faster than n squared envelope");
    }
    std::cout << "active-local n=" << leaves
              << " balanced_work="
              << balanced.back().counts.aggregate_work()
              << " balanced_norm=" << balanced.back().normalized
              << " balanced_ratio=" << balanced_ratio
              << " skew_work=" << skew.back().counts.aggregate_work()
              << " skew_norm=" << skew.back().normalized
              << " skew_ratio=" << skew_ratio << '\n';
  }
  require(balanced.back().normalized <=
              1.25 * balanced.front().normalized,
          "balanced normalized active-local work is not stable");
  require(skew.back().normalized <= 1.25 * skew.front().normalized,
          "skew normalized active-local work is not stable");

  test_virtual_tree_parity_and_publication();
  test_recursion_work_and_live_space();
}

using Test = void (*)();

const std::map<std::string, Test>& tests() {
  static const std::map<std::string, Test> kTests{
      {"parallel_loop_exception_safety",
       test_parallel_loop_exception_safety},
      {"parallel_child_exception_determinism",
       test_parallel_child_exception_determinism},
      {"concurrent_inference", test_concurrent_inference},
      {"gene_restriction_duplicates",
       test_gene_restriction_duplicates},
      {"recursion_complexity", test_recursion_complexity},
      {"deep_skewed_recursion", test_deep_skewed_recursion},
      {"recursion_failure_unwinding", test_recursion_failure_unwinding},
  };
  return kTests;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: " << argv[0] << " <test-name>\n";
    return 2;
  }
  const auto test = tests().find(argv[1]);
  if (test == tests().end()) {
    std::cerr << "unknown test: " << argv[1] << '\n';
    return 2;
  }
  try {
    test->second();
  } catch (const std::exception& ex) {
    std::cerr << argv[1] << " failed: " << ex.what() << '\n';
    return 1;
  }
  return 0;
}
