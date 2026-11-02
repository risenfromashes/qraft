#pragma once

#include "solver/solver.hpp"
#include "counting/extract_contract.hpp"

#include <array>
#include <cstdint>
#include <functional>
#include <vector>

namespace qraft::internal {

// A query's leaf order against the owners of one state. owner_query_position
// maps an owner to its query leaf position, or -1. represented_owner_order
// lists, in query order, the owners with at least one live atom.
// represented_prefix[p] and live_atom_prefix[p] count those owners and their
// live atoms among the first p query positions, so the selector of a query
// subtree follows from its leaf interval in O(1).
struct HdtQueryOwnerIndex {
  const std::vector<int>* query_owner_order = nullptr;
  std::vector<int> owner_query_position;
  std::vector<int> represented_owner_order;
  std::vector<int> represented_prefix;
  std::vector<int> live_atom_prefix;
};

// Work, size and memory counters of one worker's traversals.
// complexity_contract_satisfied is cleared when a checked size bound, shrink
// bound or HDT invariant fails.
struct HdtTraversalStats {
  long long states_visited = 0;
  long long score_reads = 0;

  long long extract_calls = 0;
  long long extract_source_hdt_nodes = 0;
  long long extract_marked_hdt_nodes = 0;
  long long extract_hanging_components = 0;
  long long extract_intermediate_units = 0;
  long long extract_live_leaves = 0;

  long long contract_calls = 0;
  long long contract_light_calls = 0;
  long long contract_heavy_calls = 0;
  long long contract_input_units = 0;
  long long contract_output_units = 0;
  long long contract_live_leaves = 0;
  long long contract_corridors_emitted = 0;
  long long contract_background_units_absorbed = 0;
  long long contract_polytomy_arms_folded = 0;
  long long contract_cc_ops = 0;
  long long contract_gg_ops = 0;
  long long contract_ig_ops = 0;
  long long contract_c_to_g_ops = 0;
  long long contract_bound_violations = 0;

  long long heavy_input_units = 0;
  long long heavy_output_units = 0;
  long long heavy_calls = 0;
  long long heavy_shrink_violations = 0;

  long long hdt_initial_atoms = 0;
  long long hdt_rebuild_nodes = 0;
  int hdt_max_recompute_degree = 0;
  long long hdt_local_balance_violations = 0;
  int hdt_max_update_height = 0;

  long long recolor_changed_leaves = 0;
  long long recolor_marked_nodes = 0;
  long long dp_transition_evals = 0;
  long long dp_scalar_fmas = 0;

  int max_state_units = 0;
  int max_state_hdt_nodes = 0;
  long long pending_state_units_peak = 0;
  long long pending_state_hdt_nodes_peak = 0;
  std::size_t extract_scratch_peak_bytes = 0;
  std::size_t contract_scratch_peak_bytes = 0;
  std::size_t state_arena_peak_bytes = 0;
  std::size_t worker_total_peak_bytes = 0;
  long long contract_mark_ns = 0;
  long long contract_background_fold_ns = 0;
  long long contract_corridor_fold_ns = 0;
  long long contract_emit_ns = 0;
  long long hdt_rebuild_ns = 0;
  bool complexity_contract_satisfied = true;
};

// Publishes the memory peaks and the complexity flag to the run statistics.
// Callers merge the remaining counters they report.
void merge_count_stats(const HdtTraversalStats& stats);

// The counting telemetry of one gene of a gene loop: HDT reads, recoloured
// leaves, recomputed nodes and the traversal statistics.
struct GeneCountTelemetry {
  long long queries = 0;
  long long recolors = 0;
  long long recomputed = 0;
  HdtTraversalStats traversal;
};

// One telemetry slot per gene while the run records statistics and none
// otherwise, so that quiet runs neither allocate nor merge them.
std::vector<GeneCountTelemetry> gene_count_telemetry_slots(int genes);

// Adds the genes' telemetry, in gene order, to the run statistics.
void merge_gene_count_telemetry(
    const std::vector<GeneCountTelemetry>& telemetry);

// Charges the InitialCompact Contract of `initial`, if any, to `stats` and
// marks the state so that the traversal does not charge its HDT again.
void hdt_account_initial_state(
    HdtInitialStateResult& initial,
    HdtTraversalStats& stats);

// Throws when an owner occupies two query leaves or when the root's leaf
// interval does not cover query_owner_order. A deferred state is indexed
// from its component graph.
HdtQueryOwnerIndex build_hdt_query_owner_index(
    const HdtState& base,
    const QueryTree& query,
    const std::vector<int>& query_owner_order);

// The O(1) selector of the owners under query_node. It borrows from `index`.
HdtLiveSelector hdt_query_node_selector(
    const HdtQueryOwnerIndex& index,
    const QueryTree& query,
    int query_node);

// Both traversals walk the query tree depth first with an explicit stack.
// A light child reuses its parent's state in place when that state has at
// most state_reuse_bound_factor times the child's Contract bound
// (hdt_contract_size_bound of its live atoms) in component units, and gets
// its own Extract/Contract state otherwise. The heavy child reuses the
// parent's state under the same test; otherwise that state is replaced by
// the heavy child's Extract/Contract, or, while the parent still needs it,
// the heavy child gets its own. So every non-root node is entered on a state
// of at most state_reuse_bound_factor times its Contract bound. A visitor
// may recolour temporarily but must restore every colour before it returns.

// K=2 counting for bipartitions. `state` enters with every query owner at
// colour 1 over background 0. The visitor runs at every non-root node that
// score_node marks (every one when score_node is empty) while exactly that
// node's taxa have colour 1, so the root channels count the quartets of the
// bipartition at the node's parent edge.
using BipartitionVisitor =
    std::function<void(int query_node, HdtState& state)>;

void count_bipartitions(
    HdtState state,
    const QueryTree& query,
    const std::vector<int>& query_owner_order,
    const std::vector<unsigned char>& score_node,
    HdtExtractWorkspace& workspace,
    HdtTraversalStats& stats,
    const BipartitionVisitor& visitor,
    int state_reuse_bound_factor = 4);

// K=4 counting at the vertices of the query tree, for regrafting and NNI.
// `state` enters with every query owner at colour 1; other taxa, such as a
// regrafting context, stay at their background colour. The entry contracts
// the state to the query owners when it holds other live leaves or exceeds
// the root's Contract bound, and builds a deferred state otherwise. At the
// root with three children, or a non-root node with two, the arms are
// coloured by vertex_arm_coloring, the visitor runs, and the colours are
// restored. The traversal never descends into a query leaf.
// requested_vertices, when given, has one byte per query node: only marked
// vertices are read, and subtrees without one are skipped.
using VertexVisitor = std::function<void(
    int vertex_node,
    const std::array<int, 3>& branch_nodes,
    const std::array<std::uint8_t, 3>& branch_colors,
    HdtState& state)>;

void count_vertices(
    HdtState state,
    const QueryTree& query,
    const std::vector<int>& query_owner_order,
    HdtExtractWorkspace& workspace,
    HdtTraversalStats& stats,
    const VertexVisitor& visitor,
    int state_reuse_bound_factor = 4,
    const std::vector<std::uint8_t>* requested_vertices = nullptr);

// K=4 counting of the clades of a query tree rooted at a vertex whose three
// arms carry colours 0, 1 and 2 (the SPR refinement's skeleton pass).
// `state` enters with every query owner at colour 3 over a background equal
// to its arm colour; rest_color[node] is the arm colour of every query node
// (the root's entry is unused). The visitor runs at every non-root node that
// score_node marks while exactly that node's taxa have colour 3 and every
// other taxon has its arm colour, so the root channels count the quartets
// with one leaf in the node's clade and one in each of the three arms of the
// root, less the clade. Subtrees without a marked node are not entered.
using ColoredCladeVisitor =
    std::function<void(int query_node, HdtState& state)>;

void count_clades_colored(
    HdtState state,
    const QueryTree& query,
    const std::vector<int>& query_owner_order,
    const std::vector<unsigned char>& score_node,
    const std::vector<std::uint8_t>& rest_color,
    HdtExtractWorkspace& workspace,
    HdtTraversalStats& stats,
    const ColoredCladeVisitor& visitor,
    int state_reuse_bound_factor = 4);

}  // namespace qraft::internal
