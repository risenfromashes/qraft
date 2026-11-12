#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <exception>
#include <iosfwd>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "qraft/qraft.hpp"
#include "qraft/taxa.hpp"

#include <taskflow/algorithm/for_each.hpp>
#include <taskflow/algorithm/reduce.hpp>
#include <taskflow/taskflow.hpp>

#ifndef QRAFT_TASKFLOW_DETERMINISTIC_WORKER_LOOKUP
#error "the deterministic Taskflow worker-lookup overlay is required"
#endif

namespace qraft::internal {

/// Side of a bipartition.
enum class Side : uint8_t { Left = 0, Right = 1 };
/// Fixed-size bitset used for taxon sets and bipartitions.
class Bitset {
 public:
  Bitset() = default;
  /// A bitset of `size` bits, all clear.
  explicit Bitset(int size);

  int size() const { return size_; }
  bool empty() const;
  /// Throws std::out_of_range for a position outside [0, size()).
  void set(int pos);
  /// Throws std::out_of_range for a position outside [0, size()).
  void reset(int pos);
  /// False for a position outside [0, size()).
  bool test(int pos) const;
  int count() const;
  void clear();
  /// Sets the bits [0, size()); padding bits stay clear.
  void fill();
  /// Bitsets of different sizes are unequal.
  bool operator==(const Bitset& other) const;
  bool operator!=(const Bitset& other) const { return !(*this == other); }

 private:
  int size_ = 0;
  std::vector<uint64_t> words_;

  /// Clears the bits of the last word beyond size(), which count() and
  /// operator== read.
  void mask_unused_bits();
};

/// Ancestry and LCA index of a starting tree, in linear space. An LCA query
/// takes O(log n): it finds the shallowest node between two first occurrences
/// in the Euler tour.
struct StartingTreeIndex {
  /// The root is its own parent.
  std::vector<int> parent;
  std::vector<int> depth;
  /// Preorder; node v's subtree is dfs_order[tin[v], tout[v]).
  std::vector<int> dfs_order;
  std::vector<int> tin;
  std::vector<int> tout;
  /// The leaf of each taxon, or -1.
  std::vector<int> taxon_to_leaf;
  /// The first position of each node in `euler`.
  std::vector<int> first_occurrence;
  std::vector<int> euler;
  /// Segment tree over `euler` holding the shallowest node of each range;
  /// its leaves start at rmq_base, a power of two.
  std::vector<int> rmq_nodes;
  int rmq_base = 0;
};

/// Builds a StartingTreeIndex once and publishes it to every reader. The
/// first reader builds it while the others wait; a failed build returns the
/// state to empty so a later reader retries.
struct StartingTreeIndexCache {
  /// 0 = empty, 1 = building, 2 = published.
  std::atomic<int> state{0};
  std::shared_ptr<const StartingTreeIndex> value;
};

/// Rooted starting tree whose edges propose the bipartitions of the divide.
struct StartingTree {
  struct Node {
    int parent = -1;
    int taxon = -1;
    std::vector<int> children;
  };

  std::vector<Node> nodes;
  int root = -1;
  std::vector<int> taxon_to_leaf;
  std::vector<int> preorder;
  std::vector<int> tin;
  std::vector<int> tout;
  // A starting tree is immutable once built. Copies share this cache, and
  // child contexts read the root's tree through starting_tree_ref, so a run
  // builds the index once.
  std::shared_ptr<StartingTreeIndexCache> index_cache =
      std::make_shared<StartingTreeIndexCache>();
};

/// Unrooted tree built by the recursion; a leaf carries its taxon.
struct SolverTree {
  struct Node {
    std::vector<int> neighbors;
    int taxon = -1;
    // An inactive node has been deleted; readers skip it.
    bool active = true;
  };

  std::vector<Node> nodes;
};

/// Rooted traversal of a gene tree; the root is its own parent. Readers that
/// need children before parents iterate preorder in reverse, since no
/// postorder is stored.
struct GeneTraversal {
  std::vector<int> parent;
  std::vector<int> preorder;
};

/// Indexed gene tree. A leaf of an input gene carries its taxon; a leaf of a
/// gene restricted to a subproblem carries its local index, the position of
/// its taxon in that subproblem's taxa (SolverContext::local_taxa).
struct GeneIndex {
  struct Node {
    // Taxon t is stored as t, local index i as -2 - i, and neither as -1.
    int key = -1;

    int taxon() const noexcept { return key >= 0 ? key : -1; }
    int local() const noexcept { return key <= -2 ? -2 - key : -1; }
    void set_taxon(int taxon) {
      if (taxon < 0) {
        if (key >= 0) key = -1;
        return;
      }
      if (key <= -2) {
        throw std::logic_error("gene node cannot hold a taxon and a local index");
      }
      key = taxon;
    }
    void set_local(int local) {
      if (local < 0) {
        if (key <= -2) key = -1;
        return;
      }
      if (key >= 0) {
        throw std::logic_error("gene node cannot hold a taxon and a local index");
      }
      if (local > std::numeric_limits<int>::max() - 1) {
        throw std::length_error("gene local index exceeds its key domain");
      }
      key = -2 - local;
    }
  };

  std::vector<Node> nodes;
  /// Per-node leaf masses, or empty when every mass is exactly 1.0. Inference
  /// never stores masses other than 1.0; tests do. Read through
  /// gene_leaf_mass.
  std::vector<double> leaf_mass;
  /// Immutable CSR adjacency in exact incidence order: the neighbors of node
  /// v are adjacency[adjacency_offsets[v], adjacency_offsets[v + 1]). The
  /// offsets hold nodes.size() + 1 entries, or none for an empty gene.
  std::vector<int> adjacency_offsets;
  std::vector<int> adjacency;
  /// The node root_traversal starts at, or -1 for an empty gene.
  int root = 0;
  /// For an input gene, the leaf of its smallest taxon, so that the side of
  /// a split that a subtree represents does not depend on the Newick root.
  /// A restricted gene uses `root`.
  int reference_root = 0;
  /// For an input gene, one leaf of each taxon (the last in node order), or
  /// -1. Duplicate leaves are not listed. Empty for a restricted gene.
  std::vector<int> taxon_to_node;
  GeneTraversal root_traversal;
  // The traversal rooted at reference_root is not stored. Its readers, MRE
  // counting and the NNI split support, rebuild it per gene with
  // gene_reference_traversal.
};

static_assert(sizeof(GeneIndex::Node) == sizeof(int),
              "gene nodes hold one packed leaf key");

/// Neighbors of one gene node, in incidence order.
inline std::span<const int> gene_neighbors(const GeneIndex& gene, int node) {
  const int* base = gene.adjacency.data();
  const int begin = gene.adjacency_offsets[static_cast<std::size_t>(node)];
  const int end = gene.adjacency_offsets[static_cast<std::size_t>(node) + 1u];
  return {base + begin, static_cast<std::size_t>(end - begin)};
}

inline int gene_degree(const GeneIndex& gene, int node) {
  return gene.adjacency_offsets[static_cast<std::size_t>(node) + 1u] -
         gene.adjacency_offsets[static_cast<std::size_t>(node)];
}

/// Mass of one gene node; 1.0 unless the gene stores explicit masses.
inline double gene_leaf_mass(const GeneIndex& gene, int node) {
  return gene.leaf_mass.empty()
             ? 1.0
             : gene.leaf_mass[static_cast<std::size_t>(node)];
}

/// Sets one node's mass, materializing explicit masses only for a value other
/// than 1.0. Call after gene.nodes has its final size.
inline void set_gene_leaf_mass(GeneIndex& gene, int node, double mass) {
  if (gene.leaf_mass.empty()) {
    if (mass == 1.0) return;
    gene.leaf_mass.assign(gene.nodes.size(), 1.0);
  }
  gene.leaf_mass[static_cast<std::size_t>(node)] = mass;
}

/// Replaces a gene's adjacency with per-node lists, preserving list order.
void set_gene_adjacency(GeneIndex& gene,
                        const std::vector<std::vector<int>>& lists);

/// The traversal of `gene` rooted at `root`, or at gene.root when `root` is
/// negative.
GeneTraversal traversal(const GeneIndex& gene, int root);

/// The traversal rooted at gene.reference_root: root_traversal when the roots
/// agree, otherwise rebuilt into `scratch`.
inline const GeneTraversal& gene_reference_traversal(const GeneIndex& gene,
                                                     GeneTraversal& scratch) {
  if (gene.reference_root == gene.root) return gene.root_traversal;
  scratch = traversal(gene, gene.reference_root);
  return scratch;
}

/// A proposed bipartition of a subproblem's taxa: `split` holds the active
/// indices of one side. representative_order is the preorder position, in the
/// query tree, of the starting-tree edge that proposes it; it breaks ties
/// between equal scores, and the trivial bipartition uses INT_MAX.
struct Candidate {
  Bitset split;
  int representative_order = 0;
};

struct BipartitionTestIndex;

/// Per-run scheduling state for parallel child subproblems. It bounds the
/// taxa held by reserved and live restricted child gene sets, so the
/// decision never depends on another concurrent run.
struct ChildSetBudget {
  /// Taxa of gene sets reserved but not yet built.
  std::atomic<long long> reserved_slots{0};
  /// Taxa of the live restricted gene sets.
  std::atomic<long long> live_slots{0};
  std::atomic<long long> live_sets{0};
};

/// The state a subproblem reads. The root context owns the taxa, the indexed
/// input genes and the starting tree; a child context reaches the root's
/// through the *_ref pointers. Its scoring genes are owned_active_genes when
/// set, else *active_genes, else the input genes (see scoring_genes).
struct SolverContext {
  TaxonRegistry taxa;
  const TaxonRegistry* taxa_ref = nullptr;
  std::vector<GeneIndex> genes;
  const std::vector<GeneIndex>* genes_ref = nullptr;
  const std::vector<GeneIndex>* active_genes = nullptr;
  std::shared_ptr<const std::vector<GeneIndex>> owned_active_genes;
  // The taxon of each local index of a restricted gene set; null while the
  // scoring genes are the input genes.
  std::shared_ptr<const std::vector<int>> local_taxa;
  // Null when the bipartition test is disabled or unavailable for the input.
  std::shared_ptr<const BipartitionTestIndex> bipartition_test_index;
  StartingTree starting_tree;
  const StartingTree* starting_tree_ref = nullptr;
  Options options;
  int workers = 1;
  // Shared by every context of one run; child contexts copy the pointer.
  std::shared_ptr<ChildSetBudget> child_budget =
      std::make_shared<ChildSetBudget>();
};

// Every counter of a run, in JSON order. Values are integers; *_ns are
// nanoseconds and *_bytes are bytes. The live_*/peak_* pairs of the payload
// group count entries; live_payload and peak_payload sum every live_*
// category except live_continuation_taxa, whose taxa live_taxa already holds.
#define QRAFT_RUN_STATS(X) \
  /* Run shape */ \
  X(subproblems) \
  X(completed_subproblems) \
  X(max_depth) \
  X(max_active_taxa) \
  /* Input */ \
  X(taxon_registry_ns) \
  X(taxon_registry_physical_labels) \
  X(taxon_registry_unique_labels) \
  X(taxon_registry_sort_items) \
  X(gene_index_ns) \
  X(startup_ns) \
  /* MRE starting tree */ \
  X(mre_distinct_splits) \
  X(mre_accepted_splits) \
  X(mre_count_ns) \
  X(mre_ids_ns) \
  X(mre_ids_forest_ns) \
  X(mre_assembly_ns) \
  X(starting_tree_index_builds) \
  X(starting_tree_index_nodes) \
  /* Quartet-guided divide and the bipartition test */ \
  X(candidate_sets) \
  X(candidates) \
  X(bipartitions_scored) \
  X(bipartition_score_ns) \
  X(bipartition_test_index_ns) \
  X(bipartition_test_index_bytes) \
  X(bipartition_test_ns) \
  X(bipartition_tests) \
  X(bipartition_test_accepts) \
  X(bipartition_test_rejects) \
  X(bipartition_test_negative_accepts) \
  X(bipartition_test_max_negative_gap_twice_k) \
  X(bipartition_test_singletons) \
  X(bipartition_test_distance_queries) \
  X(bipartition_test_incomplete) \
  X(bipartition_test_uncovered_pairs) \
  X(bipartition_test_untestable) \
  X(bipartition_test_shared_genes) \
  /* Contextual regrafting */ \
  X(regraft_calls) \
  X(regraft_vertices) \
  X(regraft_edge_ties) \
  X(regraft_ns) \
  /* NNI refinement */ \
  X(nni_edges_scored) \
  X(nni_full_score_passes) \
  X(nni_local_score_passes) \
  X(nni_local_edges_scored) \
  X(nni_moves_accepted) \
  X(nni_batch_passes) \
  X(nni_support_stops) \
  X(nni_split_support_ns) \
  X(nni_ns) \
  /* SPR refinement: rounds run, core vertices read, two-sided re-joins */ \
  /* scored and applied, pass and total time */ \
  X(spr_rounds) \
  X(spr_reads) \
  X(spr_candidates_scored) \
  X(spr_moves) \
  X(spr_pass_ns) \
  X(spr_ns) \
  X(guide_shared_taxa) \
  X(guide_contracted_edges) \
  X(guide_added_groups) \
  X(guide_ns) \
  X(guide_seed_rejected) \
  X(support_edges) \
  X(support_ns) \
  /* Quartet counting */ \
  X(active_taxa_entries) \
  X(active_taxa_index_entries) \
  X(query_tree_selected_leaves) \
  X(query_tree_lca_queries) \
  X(query_tree_retained_nodes) \
  X(gene_set_scoring_calls) \
  X(gene_set_scoring_nodes) \
  X(hdt_builds) \
  X(hdt_build_ns) \
  X(hdt_queries) \
  X(hdt_recolor_ops) \
  X(hdt_nodes_recomputed) \
  X(extract_scratch_peak_bytes) \
  X(contract_scratch_peak_bytes) \
  X(state_arena_peak_bytes) \
  X(worker_total_peak_bytes) \
  X(max_scoring_source_taxa) \
  X(scoring_source_inflation_violations) \
  /* Restricted child gene sets and parallel children */ \
  X(child_gene_set_reuses) \
  X(child_gene_set_attempts) \
  X(child_gene_set_accepts) \
  X(child_gene_set_heavy_reuses) \
  X(child_gene_set_pair_violations) \
  X(child_gene_set_input_nodes) \
  X(child_gene_set_output_nodes) \
  X(child_gene_set_max_accepted_ratio_ppm) \
  X(child_gene_set_ns) \
  X(live_child_gene_sets) \
  X(max_child_gene_sets) \
  X(live_child_gene_set_slots) \
  X(max_child_gene_set_slots) \
  X(max_child_gene_set_bytes) \
  X(serial_children) \
  X(parallel_children) \
  X(parallel_child_budget_rejections) \
  X(parallel_child_size_rejections) \
  X(max_parallel_child_slots) \
  /* Recursion work, checked by the complexity tests */ \
  X(recursion_calls) \
  X(recursion_nonbase_calls) \
  X(recursion_max_depth) \
  X(recursion_taxon_mass) \
  X(recursion_split_tests) \
  X(recursion_partition_moves) \
  X(recursion_balanced_splits) \
  X(recursion_skewed_splits) \
  X(recursion_tree_node_visits) \
  X(recursion_radix_entries) \
  X(recursion_radix_bucket_ops) \
  X(recursion_unique_taxa) \
  X(recursion_root_gene_set_nodes) \
  X(recursion_root_gene_set_entries) \
  /* Live and peak payload of the recursion */ \
  X(live_payload) \
  X(peak_payload) \
  X(live_taxa) \
  X(peak_taxa) \
  X(live_continuation_taxa) \
  X(peak_continuation_taxa) \
  X(live_gene_set_nodes) \
  X(peak_gene_set_nodes) \
  X(live_gene_set_entries) \
  X(peak_gene_set_entries) \
  X(live_gene_set_taxa) \
  X(peak_gene_set_taxa) \
  X(live_tree_nodes) \
  X(peak_tree_nodes) \
  X(live_tree_adjacencies) \
  X(peak_tree_adjacencies) \
  X(live_taxon_scratch) \
  X(peak_taxon_scratch) \
  X(accounting_violations)

/// An exact nonnegative integer below 2^128, for sums that pass 2^63.
struct WideCount {
  std::uint64_t high = 0;
  std::uint64_t low = 0;
};

/// Adds a nonnegative integer-valued double exactly, after rounding it to
/// the nearest integer; throws std::overflow_error for NaN, a negative
/// value or a sum of 2^128 or more.
void wide_add(WideCount& sum, double value);

/// The decimal digits of `value`.
std::string wide_decimal(const WideCount& value);

/// Signed integer regraft scores and SPR gains.
__extension__ typedef __int128 RegraftCount;

/// Process-wide counters of one run, written by --stats-json. stats_add,
/// stats_max and stats_add_time are no-ops unless run_stats_enabled().
struct RunStats {
#define QRAFT_DECLARE_STAT(name) std::atomic<long long> name{0};
  QRAFT_RUN_STATS(QRAFT_DECLARE_STAT)
#undef QRAFT_DECLARE_STAT
  // Cleared when a checked complexity invariant fails.
  std::atomic<long long> complexity_contract_satisfied{1};
  // The summed gain of the accepted NNI moves over the run's NNI passes.
  WideCount nni_gain;
  // The summed integer gains of accepted SPR re-joins, in gene quartets.
  RegraftCount spr_gain = 0;
  // False if an applied move could not contribute an integer gain.
  bool spr_gain_exact = true;
};

/// The edge {a, b} of a solved side where contextual regrafting attaches the
/// other side. b = -1 names the only vertex a of a one-vertex tree, and
/// a = b = -1 an empty tree.
struct RegraftEdge {
  int a = -1;
  int b = -1;
  /// raw_score or normalized_score, as --regraft-score selects.
  double rank_score = 0.0;
  /// The best rank_score outside the tie band of the chosen edge.
  double second_rank_score = 0.0;
  /// rank_score - second_rank_score; infinite without a competing edge.
  double rank_margin = 0.0;
  double raw_score = 0.0;
  /// Integer reduction of per-gene scores, when requested by the SPR oracle.
  RegraftCount exact_raw_score = 0;
  bool exact_raw_score_valid = false;
  double normalized_score = 0.0;
  /// The edges within the tie band of the best score.
  int ties = 0;
  int rank = 0;
};


/// The genes that score this context's queries: its restricted gene set when
/// it has one, otherwise the input genes.
const std::vector<GeneIndex>& scoring_genes(const SolverContext& ctx);
/// The indexed input genes of the run.
const std::vector<GeneIndex>& base_genes(const SolverContext& ctx);
/// The run's taxon registry.
const TaxonRegistry& context_taxa(const SolverContext& ctx);
/// The run's starting tree.
const StartingTree& context_starting_tree(const SolverContext& ctx);
/// The counters of the current run.
RunStats& run_stats();
/// Whether counter updates are recorded; set only between runs.
bool run_stats_enabled() noexcept;
void set_run_stats_enabled(bool enabled) noexcept;
void reset_run_stats();
void stats_add(std::atomic<long long>& counter, long long value = 1);
void stats_max(std::atomic<long long>& counter, long long value);
void stats_add_time(std::atomic<long long>& counter,
                    std::chrono::nanoseconds elapsed);
/// Records the summed gain of the accepted NNI moves.
void stats_set_nni_gain(const WideCount& gain);
/// Live payload accounting: a category's live/peak pair, whose live values
/// also sum into live_payload/peak_payload. Live values change even while
/// the statistics are disabled, which keeps acquires and releases balanced;
/// peaks are recorded only while enabled.
void stats_acquire_payload(std::atomic<long long>& live,
                           std::atomic<long long>& peak, long long value);
void stats_release_payload(std::atomic<long long>& live, long long value);
/// A live/peak pair outside live_payload (it overlaps a payload category).
void stats_acquire_gauge(std::atomic<long long>& live,
                         std::atomic<long long>& peak, long long value);
void stats_release_gauge(std::atomic<long long>& live, long long value);
/// The name of a progress mode in --stats-json output.
const char* progress_mode_name(ProgressMode mode);
/// Initializes progress reporting for one inference run.
void progress_start(const Options& options, int total_splits);
/// Emits a startup status line before the recursion begins.
void progress_startup_status(const Options& options,
                             std::string stage,
                             std::string detail = {});
/// Starts the progress bar of the MRE starting tree.
void progress_mre_start(const Options& options,
                              int total_genes,
                              int taxon_count);
/// Records one gene whose splits MRE counting has processed.
void progress_mre_count_gene();
/// Updates the progress of MRE tree assembly.
void progress_mre_assembly(long long checked_entries,
                                 long long total_entries,
                                 int accepted,
                                 int max_accepted);
/// Finishes the MRE progress bar.
void progress_mre_finish(int accepted, int max_accepted);
/// Records the start of one recursive subproblem on `active_taxa` taxa.
void progress_note_subproblem(uint64_t path, int depth, int active_taxa);
/// Updates the live phase label for a recursive subproblem.
void progress_note_phase(uint64_t path, std::string phase);
/// Updates the global phase label shown when no recursive subproblem is active.
void progress_note_stage(std::string stage);
/// Records the start and end of one regraft of two child trees.
void progress_note_regraft_start(int parent_active_size,
                                 int left_taxa,
                                 int right_taxa);
void progress_note_regraft_done(int parent_active_size);
/// Updates NNI refinement progress.
void progress_note_nni(std::string mode,
                           int moves,
                           int max_moves,
                           long long edges_scored,
                           double gain,
                           double min_gain);
/// Records the bipartition chosen for a subproblem.
void progress_note_split(uint64_t path,
                         int active_taxa,
                         int left_size,
                         int right_size);
/// Finishes progress output for one inference run.
void progress_finish();

GeneTraversal traversal(const GeneIndex& gene, int root);
/// Indexes a parsed gene tree by taxon, suppressing unlabelled vertices of
/// degree two.
GeneIndex index_gene(const GeneTree& gene, const TaxonRegistry& taxa);
/// The gene induced on a subproblem's taxa: every leaf whose taxon has a
/// local index (taxon_to_local[taxon] >= 0) is kept, relabelled by that index,
/// and suppressed internal nodes are removed. Duplicate leaves of one taxon
/// are all kept.
GeneIndex restrict_gene_by_taxon(const GeneIndex& gene,
                                 const std::vector<int>& taxon_to_local);
/// The same for a gene that is already restricted: a leaf with local index i
/// is kept as local index local_to_local[i] when that is nonnegative.
GeneIndex restrict_gene_by_local(const GeneIndex& gene,
                                 const std::vector<int>& local_to_local);

/// --guide-tree: while the MRE is built, rejects the clusters that conflict
/// with the guide tree on the taxa both hold.
class GuideFilter {
 public:
  GuideFilter(const Tree& guide, const TaxonRegistry& taxa);
  ~GuideFilter();
  /// True when the cluster of these taxa conflicts with the guide.
  bool conflicts(const std::vector<int>& cluster_taxa) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

/// Builds the MRE starting tree, the majority-rule extended consensus of the
/// genes, skipping the clusters that `guide` rejects when it is set.
StartingTree build_mre_tree(const std::vector<GeneIndex>& genes,
                              int taxon_count,
                              const Options& options,
                              const GuideFilter* guide = nullptr);

/// Checks of the MRE builder's radix aggregation, split ids, AVL forest and
/// flat topology on synthetic split occurrences, for the complexity tests.
struct MreComplexityProbeResult {
  bool radix_matches_oracle = false;
  bool dense_split_ids = false;
  bool avl_invariants_hold = false;
  bool range_aggregates_match = false;
  int radix_passes = 0;
  long long radix_items = 0;
  int occurrence_radix_passes = 0;
  int priority_radix_passes = 0;
  long long occurrence_radix_items = 0;
  long long priority_radix_items = 0;
  long long emitted_occurrences = 0;
  long long deduplicated_gene_splits = 0;
  long long distinct_splits = 0;
  int max_avl_height = 0;
  int avl_height_limit = 0;
  long long empty_complement_visits = 0;
  long long empty_complement_visit_limit = 0;
  bool flat_topology_matches = false;
  int topology_peak_slots = 0;
  int topology_slot_limit = 0;
};

MreComplexityProbeResult mre_complexity_probe_for_tests(int item_count);
/// Gene-tree quartet frequencies around one internal edge {edge_a, edge_b}
/// of a solver tree whose two ends are binary: the quartets with one leaf in
/// each of the edge's four arms, summed over genes, that a gene resolves as
/// the tree's own topology (main) or as either alternative.
struct BranchQuartetFrequencies {
  int edge_a = -1;
  int edge_b = -1;
  double main = 0.0;
  double alternative1 = 0.0;
  double alternative2 = 0.0;
  /// The product of the four arms' taxon counts.
  double quartets_around = 0.0;
};
/// The frequencies of every internal edge of `tree` on `active_taxa` whose
/// ends are both binary, in one full K=4 counting pass. With
/// `per_gene_normalized`, each gene's three counts are divided by their sum
/// before the sum over genes (ASTRAL, ASTER); otherwise they are raw counts.
std::vector<BranchQuartetFrequencies> branch_quartet_frequencies(
    const SolverContext& ctx,
    const SolverTree& tree,
    const std::vector<int>& active_taxa,
    bool per_gene_normalized);

/// log I_x(a, b), the regularized incomplete beta function, accurate when
/// I_x(a, b) is tiny and when it is close to 1.
double log_regularized_incomplete_beta(double a, double b, double x);
/// The local posterior probabilities of an edge's three topologies from
/// their quartet frequencies f (the tree's own topology first), under the
/// Yule prior rate 0.5 of ASTRAL and ASTER.
std::array<double, 3> local_posterior(const std::array<double, 3>& f);
/// The edge length in coalescent units from the main topology's frequency f1
/// of n = f1 + f2 + f3: -ln(3/2 - 3/2 f1 / (n + 1)) when f1 > (n + 1) / 3,
/// else 0;
/// NaN when n is 0.
double coalescent_length(double f1, double n);
/// Writes support:length (or the detailed label) on every internal edge of
/// `species` whose ends are binary in `tree`; species_node maps the nodes of
/// `tree` to those of `species` (-1 for none).
void annotate_branch_support(const SolverContext& ctx,
                             const SolverTree& tree,
                             const std::vector<int>& species_node,
                             SpeciesTree& species);
/// A user starting tree (--starting-tree) as a StartingTree.
StartingTree starting_tree_from(const Tree& tree, const TaxonRegistry& taxa);
/// --guide-tree: `tree` with every edge that conflicts with the guide on the
/// shared taxa contracted, then its polytomies refined by the guide
/// bipartitions that keep their arms whole.
StartingTree apply_guide_tree(const StartingTree& tree, const Tree& guide, const TaxonRegistry& taxa);
/// A starting tree as Newick, children in their stored order.
std::string starting_tree_newick(const StartingTree& tree, const TaxonRegistry& taxa);
/// The side of each of the first `m` indices: Left where `split` is set.
std::vector<Side> assignment_from_split(const Bitset& split, int m);
/// QRAFT on every taxon: the divide and regraft, run as a loop over heap
/// frames, then NNI refinement when enabled. When `subflow` is set, the
/// recursion runs on an executor worker and independent children may be
/// solved in parallel as executor jobs.
SolverTree solve_qraft(SolverContext& ctx, tf::Subflow* subflow);
/// A split shape forced by the recursion tests in place of the chosen one.
enum class ForcedSplitShape : uint8_t { Balanced, MaximallySkewed };
/// solve_qraft with every bipartition replaced by `shape` and with the
/// recursion's payload accounting enabled. Every subproblem `fail_depth`
/// deep throws a std::runtime_error naming its path.
SolverTree solve_qraft_for_test(SolverContext& ctx,
                                ForcedSplitShape shape,
                                int fail_depth = -1);
/// The starting tree as a solver tree, refined by NNI when `run_nni`.
SolverTree refine_starting_tree(SolverContext& ctx, bool run_nni);

/// Deletes every active unlabelled vertex of degree two, joining its two
/// neighbours: a starting tree's binary root, for example.
void suppress_solver_degree_two(SolverTree& tree);
/// Replaces the shared Taskflow executor with one of options.threads workers
/// (every hardware thread when not positive), tunes the allocator for them,
/// and returns the worker count.
int configure_threads(const Options& options);
/// The shared Taskflow executor, created on first use.
tf::Executor& get_executor();
int get_num_workers();
/// Tunes glibc's malloc for `workers` threads; a no-op with the linked
/// jemalloc or another C library.
void configure_allocator(int workers);
/// Returns pages of already freed memory to the operating system. The
/// linked jemalloc otherwise keeps them resident for up to ten seconds;
/// callers use this after large one-time releases. A no-op without jemalloc.
void release_retained_allocator_pages();
/// Returns freed pages to the operating system after a child of a subproblem
/// on `active_taxa` taxa is solved; a no-op for small subproblems.
void release_after_subproblem(int active_taxa);
/// Runs `taskflow` to completion. A worker thread coruns it rather than
/// blocking, so nested flows cannot deadlock the executor.
void run_inner(tf::Taskflow& taskflow);
/// Mixes three integers into a deterministic 64-bit hash.
uint64_t stable_hash_u64(uint64_t a, uint64_t b, uint64_t c);

/// The worker limit of inner parallel loops started on this thread; zero
/// means every executor worker.
inline thread_local int tl_inner_parallel_worker_budget = 0;

/// Makes inner loops started on this thread serial for its lifetime. A worker
/// holds one while it coruns other tasks as it waits for a parallel loop, so
/// the loops those tasks start submit no further helpers.
class ScopedInnerParallelBudgetSuspend {
 public:
  ScopedInnerParallelBudgetSuspend()
      : saved_budget_(tl_inner_parallel_worker_budget) {
    tl_inner_parallel_worker_budget = 1;
  }

  ScopedInnerParallelBudgetSuspend(
      const ScopedInnerParallelBudgetSuspend&) = delete;
  ScopedInnerParallelBudgetSuspend& operator=(
      const ScopedInnerParallelBudgetSuspend&) = delete;

  ~ScopedInnerParallelBudgetSuspend() {
    tl_inner_parallel_worker_budget = saved_budget_;
  }

 private:
  int saved_budget_;
};

/// Lets inner loops started on this thread use every worker for its lifetime.
/// A child subproblem runs under one, so that its loops stay parallel when a
/// worker picks the child up while it coruns inside a suspended loop wait.
class ScopedInnerParallelBudgetRestore {
 public:
  ScopedInnerParallelBudgetRestore()
      : saved_budget_(tl_inner_parallel_worker_budget) {
    tl_inner_parallel_worker_budget = 0;
  }

  ScopedInnerParallelBudgetRestore(
      const ScopedInnerParallelBudgetRestore&) = delete;
  ScopedInnerParallelBudgetRestore& operator=(
      const ScopedInnerParallelBudgetRestore&) = delete;

  ~ScopedInnerParallelBudgetRestore() {
    tl_inner_parallel_worker_budget = saved_budget_;
  }

 private:
  int saved_budget_;
};

/// The number of workers an inner parallel loop started here may use.
inline int current_inner_parallel_worker_limit() {
  const int workers = std::max(1, get_num_workers());
  if (tl_inner_parallel_worker_budget <= 0) {
    return workers;
  }
  return std::max(1, std::min(workers, tl_inner_parallel_worker_budget));
}

/// The first exception one loop participant caught, with its index.
struct ParallelLoopFailure {
  int index = std::numeric_limits<int>::max();
  std::exception_ptr exception;

  void capture(int logical_index) noexcept {
    // Each slot belongs to one ascending participant, so its first failure is
    // also that participant's lowest-index failure.
    if (!exception) {
      index = logical_index;
      exception = std::current_exception();
    }
  }
};

/// Counts the helpers of a parallel loop that have not finished.
class ParallelLoopJoinState {
 public:
  explicit ParallelLoopJoinState(int helpers) : remaining_(helpers) {}

  void complete(int count = 1) noexcept {
    if (remaining_.fetch_sub(count, std::memory_order_acq_rel) == count) {
      remaining_.notify_one();
    }
  }

  int remaining() const noexcept {
    return remaining_.load(std::memory_order_acquire);
  }

  void wait(int observed) const noexcept {
    remaining_.wait(observed, std::memory_order_acquire);
  }

 private:
  std::atomic<int> remaining_;
};

/// Marks one helper finished when it leaves scope, even on failure.
class ParallelLoopCompletion {
 public:
  explicit ParallelLoopCompletion(ParallelLoopJoinState& join) noexcept
      : join_(join) {}
  ParallelLoopCompletion(const ParallelLoopCompletion&) = delete;
  ParallelLoopCompletion& operator=(const ParallelLoopCompletion&) = delete;
  ~ParallelLoopCompletion() { join_.complete(); }

 private:
  ParallelLoopJoinState& join_;
};

template <typename F>
void run_parallel_loop_range(int begin,
                             int end,
                             F& body,
                             ParallelLoopFailure& failure) noexcept {
  // Parallel loops finish all iterations after a callback failure. This keeps
  // completion bounded and makes the lowest failing logical index observable.
  for (int i = begin; i < end; ++i) {
    try {
      body(i);
    } catch (...) {
      failure.capture(i);
    }
  }
}

/// Waits until every helper has finished. A worker thread runs other tasks
/// meanwhile; any other thread blocks.
inline void wait_for_parallel_loop(
    tf::Executor& exec,
    const std::shared_ptr<ParallelLoopJoinState>& join) {
  if (exec.this_worker_id() >= 0) {
    ScopedInnerParallelBudgetSuspend suspend_budget;
    exec.corun_until([&] { return join->remaining() <= 0; });
  } else {
    int observed = join->remaining();
    while (observed > 0) {
      join->wait(observed);
      observed = join->remaining();
    }
  }
}

/// The failure with the lowest loop index, or null.
inline std::exception_ptr select_parallel_loop_failure(
    const std::vector<ParallelLoopFailure>& failures) noexcept {
  const ParallelLoopFailure* selected = nullptr;
  for (const ParallelLoopFailure& failure : failures) {
    if (failure.exception &&
        (selected == nullptr || failure.index < selected->index)) {
      selected = &failure;
    }
  }
  return selected == nullptr ? std::exception_ptr{} : selected->exception;
}

/// Runs body(i) for i in [begin, end) on the shared executor, one contiguous
/// range per worker; serially below `min_tasks` indices or with one worker.
/// In parallel every index runs even after a failure. The call rethrows the
/// failure with the lowest index, after every helper has joined.
template <typename F>
void parallel_for_inner(int begin, int end, F&& body, int min_tasks = 4) {
  const int n = end - begin;
  if (n <= 0) return;
  auto& exec = get_executor();
  const int nw = std::min(current_inner_parallel_worker_limit(), n);
  if (nw <= 1 || n < min_tasks) {
    for (int i = begin; i < end; ++i) body(i);
    return;
  }
  std::vector<ParallelLoopFailure> failures(static_cast<std::size_t>(nw));
  auto join = std::make_shared<ParallelLoopJoinState>(nw - 1);
  std::exception_ptr submission_failure;
  int first_unsubmitted = nw;
  for (int t = 1; t < nw; ++t) {
    int lo = begin + static_cast<int>(static_cast<int64_t>(n) * t / nw);
    int hi = begin + static_cast<int>(static_cast<int64_t>(n) * (t + 1) / nw);
    try {
      exec.silent_async(
          [lo, hi, &body, &failure = failures[t], join]() noexcept {
            ParallelLoopCompletion completion(*join);
            run_parallel_loop_range(lo, hi, body, failure);
          });
    } catch (...) {
      submission_failure = std::current_exception();
      first_unsubmitted = t;
      join->complete(nw - t);
      break;
    }
  }
  int hi0 = begin + static_cast<int>(static_cast<int64_t>(n) / nw);
  run_parallel_loop_range(begin, hi0, body, failures[0]);
  for (int t = first_unsubmitted; t < nw; ++t) {
    int lo = begin + static_cast<int>(static_cast<int64_t>(n) * t / nw);
    int hi = begin + static_cast<int>(static_cast<int64_t>(n) * (t + 1) / nw);
    run_parallel_loop_range(lo, hi, body, failures[t]);
  }
  wait_for_parallel_loop(exec, join);
  if (std::exception_ptr failure = select_parallel_loop_failure(failures)) {
    std::rethrow_exception(failure);
  }
  if (submission_failure) {
    std::rethrow_exception(submission_failure);
  }
}

/// Like parallel_for_inner, but workers take indices one at a time from a
/// shared atomic cursor, which suits indices of uneven cost. Each index should
/// be a coarse chunk of work that amortises the atomic. Which worker runs an
/// index is not deterministic, so callers must reduce results in index
/// order, never in completion order.
template <typename F>
void parallel_for_dynamic(int begin, int end, F&& body, int min_tasks = 2) {
  const int n = end - begin;
  if (n <= 0) return;
  auto& exec = get_executor();
  const int nw = std::min(current_inner_parallel_worker_limit(), n);
  if (nw <= 1 || n < min_tasks) {
    for (int i = begin; i < end; ++i) body(i);
    return;
  }
  std::atomic<int> cursor{begin};
  std::vector<ParallelLoopFailure> failures(static_cast<std::size_t>(nw));
  auto join = std::make_shared<ParallelLoopJoinState>(nw - 1);
  std::exception_ptr submission_failure;
  auto drain = [&](ParallelLoopFailure& failure) noexcept {
    int i;
    while ((i = cursor.fetch_add(1, std::memory_order_relaxed)) < end) {
      try {
        body(i);
      } catch (...) {
        failure.capture(i);
      }
    }
  };
  for (int t = 1; t < nw; ++t) {
    try {
      exec.silent_async([&drain, &failure = failures[t], join]() noexcept {
        ParallelLoopCompletion completion(*join);
        drain(failure);
      });
    } catch (...) {
      submission_failure = std::current_exception();
      join->complete(nw - t);
      break;
    }
  }
  drain(failures[0]);
  wait_for_parallel_loop(exec, join);
  if (std::exception_ptr failure = select_parallel_loop_failure(failures)) {
    std::rethrow_exception(failure);
  }
  if (submission_failure) {
    std::rethrow_exception(submission_failure);
  }
}

}  // namespace qraft::internal
