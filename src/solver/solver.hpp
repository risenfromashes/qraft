#pragma once

#include "internal.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace qraft::internal {

/// Local indices of a list of distinct taxa. taxa[i] is the taxon of local
/// index i; by_taxon holds the (taxon, local index) pairs sorted by taxon.
struct TaxonIndex {
  std::vector<int> taxa;
  std::vector<std::pair<int, int>> by_taxon;

  /// The local index of `taxon`, or -1; O(log m).
  int find(int taxon) const {
    const auto it = std::lower_bound(
        by_taxon.begin(), by_taxon.end(), taxon,
        [](const std::pair<int, int>& entry, int value) {
          return entry.first < value;
        });
    return it != by_taxon.end() && it->first == taxon ? it->second : -1;
  }
};

/// A subproblem's taxa; the active index of a taxon is its position in
/// `taxa`.
struct ActiveTaxa {
  std::vector<int> taxa;
  TaxonIndex active_by_taxon;
};

/// The bipartition the divide chose for a subproblem.
struct BipartitionChoice {
  Candidate candidate;
  /// The index of the chosen proposal, or -1 for the trivial bipartition.
  int candidate_index = -1;
  /// The K=2 score W(00|11) - W(01|01); zero for the balanced score and the
  /// trivial bipartition.
  double score = 0.0;
  int representative_order = 0;
  /// The number of distinct bipartitions the starting-tree edges propose.
  int candidate_count = 0;
  int left_count = 0;
  int right_count = 0;
};

/// A scored edge for select_bipartition_for_tests. Its side is the leaves at
/// query preorder leaf positions [leaf_intervals[0], leaf_intervals[1]) and
/// [leaf_intervals[2], leaf_intervals[3]).
struct BipartitionProbe {
  double score = 0.0;
  int left_count = 0;
  int representative_order = 0;
  int edge_child_node = 0;
  std::array<int, 4> leaf_intervals{};
};

/// The three K=4 channels at the root of a query: q01_23 is W(01|23),
/// q02_13 is W(02|13) and q03_12 is W(03|12).
struct ChannelCounts {
  double q01_23 = 0.0;
  double q02_13 = 0.0;
  double q03_12 = 0.0;
};

/// The chunk count for a gene loop that holds one counting workspace per
/// chunk: at most one chunk per worker the loop may use.
inline int hdt_workspace_worker_limit(const Options& options) {
  (void)options;
  return current_inner_parallel_worker_limit();
}

// One chunk, and so one score row, per gene. Rows are filled in any order and
// reduced in gene order, so scores do not depend on the worker count, and no
// fixed chunk count limits the workers.
inline int gene_partition_count(int gene_count) {
  return std::max(0, gene_count);
}

/// A table of `rows` rows of `width` values, stored one row after another and
/// left uninitialized until init_row(r) value-initializes row r. A gene loop
/// gives each gene one row, so the task that fills a row initializes it, and
/// no single thread touches the whole table first.
template <class T>
class RowTable {
  static_assert(std::is_trivially_copyable_v<T> &&
                std::is_trivially_destructible_v<T>);

 public:
  RowTable(std::size_t rows, std::size_t width)
      : width_(width), data_(allocate(rows, width)) {}

  T* init_row(std::size_t r) {
    T* row = data_.get() + r * width_;
    std::uninitialized_value_construct_n(row, width_);
    return row;
  }

  T* row(std::size_t r) { return data_.get() + r * width_; }
  const T* data() const { return data_.get(); }

 private:
  struct Free {
    void operator()(T* ptr) const {
      ::operator delete(static_cast<void*>(ptr), std::align_val_t{alignof(T)});
    }
  };

  static T* allocate(std::size_t rows, std::size_t width) {
    if (width != 0 &&
        rows > std::numeric_limits<std::size_t>::max() / width / sizeof(T)) {
      throw std::length_error("row table exceeds the addressable size domain");
    }
    const std::size_t count = std::max<std::size_t>(1, rows * width);
    return static_cast<T*>(
        ::operator new(count * sizeof(T), std::align_val_t{alignof(T)}));
  }

  std::size_t width_;
  std::unique_ptr<T, Free> data_;
};

/// Adds `row_count` rows of `width` values, stored one row after another, into
/// out[0, width) with add(total, part). Blocks of columns run in parallel, and
/// each column adds the rows in row order, so the sums equal a serial
/// row-by-row reduction for every worker count.
template <class T, class Add>
void add_rows_in_order(const T* rows,
                       int row_count,
                       std::size_t width,
                       T* out,
                       Add&& add) {
  if (row_count <= 0 || width == 0) return;
  // A block adds at least this many values, and has at least one column.
  constexpr std::size_t kAddsPerBlock = std::size_t{1} << 16;
  const std::size_t adds = width * static_cast<std::size_t>(row_count);
  const int blocks = static_cast<int>(std::clamp<std::size_t>(
      adds / kAddsPerBlock, 1,
      std::min<std::size_t>(
          width,
          4 * static_cast<std::size_t>(current_inner_parallel_worker_limit()))));
  parallel_for_inner(0, blocks, [&](int block) {
    const std::size_t first = width * static_cast<std::size_t>(block) /
                              static_cast<std::size_t>(blocks);
    const std::size_t last = width * static_cast<std::size_t>(block + 1) /
                             static_cast<std::size_t>(blocks);
    for (int r = 0; r < row_count; ++r) {
      const T* row = rows + static_cast<std::size_t>(r) * width;
      for (std::size_t c = first; c < last; ++c) add(out[c], row[c]);
    }
  }, 2);
}

/// A rooted tree on the active taxa of a query, induced from the starting
/// tree or from a solver tree. Leaves carry active indices, and every array is
/// indexed by query node. tin/tout bound each subtree in `preorder`, and the
/// leaves of node v's subtree are leaf_active_order[subtree_leaf_first[v],
/// subtree_leaf_first[v] + subtree_leaf_count[v]), in preorder.
struct QueryTree {
  struct Node {
    /// This node's neighbours are adjacency[first_adj, first_adj + degree).
    int first_adj = -1;
    int degree = 0;
    /// The active index of a leaf, or -1.
    int taxon_active = -1;
    /// The solver-tree node this node comes from, or -1.
    int source_node = -1;
  };

  std::vector<Node> nodes;
  std::vector<int> adjacency;
  std::vector<int> parent;
  std::vector<int> child_first;
  std::vector<int> child_count;
  std::vector<int> child_nodes;
  std::vector<int> preorder;
  std::vector<int> postorder;
  std::vector<int> tin;
  std::vector<int> tout;
  std::vector<int> subtree_size;
  std::vector<int> subtree_leaf_first;
  std::vector<int> subtree_leaf_count;
  std::vector<int> leaf_active_order;
  std::vector<int> heavy_child;
  std::vector<int> active_to_node;
  int root = -1;
  int leaf_count = 0;
};

/// The colours of the three arms of a query vertex for K=4 counting: the
/// parent, the larger child and the smaller child of an internal vertex, or
/// the three children of a root of degree three.
struct VertexArmColoring {
  std::array<int, 3> branch_nodes{-1, -1, -1};
  std::array<std::uint8_t, 3> branch_colors{0, 1, 2};
  int temporary_color_two_node = -1;
};

inline int query_leaf_count(const QueryTree& query, int node) {
  return node >= 0 && node < static_cast<int>(query.subtree_leaf_count.size())
             ? query.subtree_leaf_count[static_cast<std::size_t>(node)]
             : 0;
}

inline VertexArmColoring vertex_arm_coloring(
    const QueryTree& query,
    int node) {
  VertexArmColoring coloring;
  if (node < 0 || node >= static_cast<int>(query.nodes.size()) ||
      node >= static_cast<int>(query.child_first.size())) {
    return coloring;
  }
  const int child_begin = query.child_first[static_cast<std::size_t>(node)];
  if (node == query.root) {
    if (query.child_count[static_cast<std::size_t>(node)] != 3) {
      return coloring;
    }
    coloring.branch_nodes = std::array<int, 3>{
        query.child_nodes[static_cast<std::size_t>(child_begin)],
        query.child_nodes[static_cast<std::size_t>(child_begin + 1)],
        query.child_nodes[static_cast<std::size_t>(child_begin + 2)]};
    // The largest branch keeps the query colour 1, so a traversal recolours
    // only the two smaller branches. Readers map counts through
    // branch_colors, so this relabels exact values.
    const std::array<int, 3>& branches = coloring.branch_nodes;
    int keep = 1;
    for (int i = 0; i < 3; ++i) {
      if (query_leaf_count(query, branches[static_cast<std::size_t>(i)]) >
          query_leaf_count(query,
                               branches[static_cast<std::size_t>(keep)])) {
        keep = i;
      }
    }
    const int color_zero = keep == 0 ? 1 : 0;
    const int color_two = keep == 2 ? 1 : 2;
    coloring.branch_colors = std::array<std::uint8_t, 3>{1, 1, 1};
    coloring.branch_colors[static_cast<std::size_t>(color_zero)] = 0;
    coloring.branch_colors[static_cast<std::size_t>(color_two)] = 2;
    coloring.temporary_color_two_node =
        coloring.branch_nodes[static_cast<std::size_t>(color_two)];
    return coloring;
  }
  if (query.child_count[static_cast<std::size_t>(node)] != 2) {
    return coloring;
  }
  const int child0 = query.child_nodes[static_cast<std::size_t>(child_begin)];
  const int child1 =
      query.child_nodes[static_cast<std::size_t>(child_begin + 1)];
  int color_one_child = child0;
  int color_two_child = child1;
  if (query_leaf_count(query, child0) <
      query_leaf_count(query, child1)) {
    color_one_child = child1;
    color_two_child = child0;
  }
  coloring.branch_nodes = std::array<int, 3>{
      query.parent[static_cast<std::size_t>(node)], color_one_child,
      color_two_child};
  coloring.branch_colors = std::array<std::uint8_t, 3>{0, 1, 2};
  coloring.temporary_color_two_node = color_two_child;
  return coloring;
}

/// An NNI across the edge {edge_a, edge_b}: swap_a_neighbor, a neighbour of
/// edge_a, trades places with swap_b_neighbor, a neighbour of edge_b.
/// `alternative` (1 or 2) names which of the edge's two NNIs this is, and
/// `gain` is the change in quartet score.
struct NniMove {
  int edge_a = -1;
  int edge_b = -1;
  int swap_a_neighbor = -1;
  int swap_b_neighbor = -1;
  int alternative = 0;
  double gain = 0.0;
};

ActiveTaxa make_active_taxa(const std::vector<int>& taxa, int taxon_count);

/// Indexes the distinct taxa of `taxa` that lie in [0, taxon_count), in order
/// of first occurrence.
TaxonIndex make_taxon_index(const std::vector<int>& taxa,
                                   int taxon_count);

/// The active index of `taxon`, or -1.
int active_index(const ActiveTaxa& active, int taxon);

/// Counts one scoring pass over `genes` and their nodes in the run statistics.
void note_gene_set_use(const std::vector<GeneIndex>& genes);

/// The starting tree induced on `active`: the leaves of the active taxa and
/// the LCAs of leaves adjacent in preorder, rooted at the LCA of them all.
QueryTree induce_query_tree(
    const StartingTree& starting_tree,
    const ActiveTaxa& active);

/// The solver tree induced on `active_taxa`, with unlabelled vertices of
/// degree two suppressed.
QueryTree query_tree_from_solver_tree(
    const SolverTree& tree,
    const std::vector<int>& active_taxa,
    int taxon_count);

/// The query tree of a binary `tree` on all of `active_taxa`, rooted at its
/// internal vertex `vertex`; the root's children, the vertex's three arms,
/// are its neighbours' subtrees in node order.
QueryTree query_tree_rooted_at_vertex(const SolverTree& tree,
                                      const std::vector<int>& active_taxa,
                                      int taxon_count, int vertex);

/// The owner of each query leaf, in preorder leaf order.
std::vector<int> query_leaf_owner_order(
    const QueryTree& query,
    const std::vector<int>& query_active_to_owner);

/// The tree on at most three taxa: a leaf, an edge, or a star.
SolverTree base_tree(const std::vector<int>& taxa);

std::vector<int> active_neighbors(const SolverTree& tree, int node);

bool has_solver_edge(const SolverTree& tree, int a, int b);

void add_solver_edge(SolverTree& tree, int a, int b);

void remove_solver_edge(SolverTree& tree, int a, int b);

/// Joins two solved sides at their regraft edges. Each edge is subdivided by
/// a new connector vertex (a one-vertex tree connects at its vertex, an empty
/// tree through a new vertex), and the two connectors are joined. The right
/// tree's nodes follow the left's. Throws std::logic_error for an edge that
/// does not match its tree.
SolverTree connect_at_edges(SolverTree left,
                    SolverTree right,
                    RegraftEdge left_edge,
                    RegraftEdge right_edge);

// Rethrows the left child's failure, else the right child's, once both
// parallel children have joined.
void rethrow_parallel_child_failure(
    const std::exception_ptr& first,
    const std::exception_ptr& second);

/// The quartet-guided divide: offers the best-ranked bipartitions proposed by
/// the starting tree's edges on `taxa` (by K=2 score, or by balance) to the
/// bipartition test in turn, at most options.bipartition_test.proposals of
/// them, and takes the first accepted. Without one it takes the trivial
/// bipartition.
BipartitionChoice choose_bipartition(const SolverContext& ctx,
                                     const std::vector<int>& taxa,
                                     uint64_t path,
                                     int depth);

/// The index of the probe the divide would select, or -1.
int select_bipartition_for_tests(
    const std::vector<BipartitionProbe>& probes,
    int active_count,
    int min_side);

/// The edge of `child_tree`, on `active_taxa`, where the taxa
/// `context_taxa` attach.
RegraftEdge choose_regraft_edge(const SolverContext& ctx,
                                const SolverTree& child_tree,
                                const std::vector<int>& active_taxa,
                                const std::vector<int>& context_taxa);

/// Contextual regrafting of two solved sides: the edge of each side where the
/// other attaches, both scored in one pass over the genes. The optional
/// arrays receive every edge's score (a side without a scored edge reports
/// its trivial edges); force_raw ranks by raw scores whatever the run's
/// regraft score. A side's optional `allowed` mask, indexed by node, limits
/// its choice to the edges with both ends marked (a side without a scored
/// edge ignores it).
std::pair<RegraftEdge, RegraftEdge> choose_regraft_edges(
    const SolverContext& ctx,
    const SolverTree& left_tree,
    const SolverTree& right_tree,
    const std::vector<int>& left_taxa,
    const std::vector<int>& right_taxa,
    std::vector<RegraftEdge>* left_scores = nullptr,
    std::vector<RegraftEdge>* right_scores = nullptr,
    bool force_raw = false,
    bool exact_scores = false,
    const std::vector<char>* left_allowed = nullptr,
    const std::vector<char>* right_allowed = nullptr);

/// The K=2 score W(00|11) - W(01|01), summed over the scoring genes, of the
/// bipartition cut by the edge above each of `edge_child_nodes`.
std::vector<double> score_bipartitions(
    const SolverContext& ctx,
    const ActiveTaxa& active,
    const QueryTree& query,
    const std::vector<int>& edge_child_nodes);

/// NNI refinement of `tree` on `active_taxa` when options.nni is set: applies
/// the best-gaining move, or under options.nni_batch the best set of moves
/// whose gains add, while gains reach the threshold and the support gate
/// allows them, up to a fixed number of moves unless options.nni_move_cap is
/// false.
void refine_nni(SolverContext& ctx,
                        SolverTree& tree,
                        const std::vector<int>& active_taxa);

/// A quartet count summed over the genes: one skeleton channel holds up to
/// |X| |A| |B| |C| quartets per gene, beyond 64 bits at 200,000 taxa.
__extension__ typedef __int128 SkeletonCount;

/// The skeleton pass of the SPR refinement. Its core is the set of
/// endpoints of the heavy edges, those with at least `resolution` taxa on
/// both sides, at the smallest resolution (at least 2) whose endpoints
/// number at most the read budget; the core is connected. For every core
/// vertex u, its three arms coloured 0, 1 and 2, and every scored edge e
/// (at least `min_clade` taxa on its smaller side), the channels count over
/// the scoring genes the quartets with one leaf x in the side of e away from
/// u and one leaf in each arm of u outside that side, by the arm whose leaf
/// the gene pairs with x. These are the terms of the regraft sigma identity
/// at u, so sums of their differences along a core path are exact quartet
/// gains of re-attaching that side.
struct SkeletonPass {
  SolverTree tree;   ///< the input with degree-two vertices suppressed
  QueryTree query;   ///< reference query tree; scored edges are its nodes
  int resolution = 0;
  int min_clade = 0;
  std::vector<int> edge_nodes;          ///< per scored edge, the node below it
  std::vector<int> edge_child_source;   ///< that node's solver vertex
  std::vector<int> edge_parent_source;  ///< its parent's solver vertex
  std::vector<int> edge_min_side;       ///< taxa on the edge's smaller side
  std::vector<std::uint8_t> edge_heavy; ///< both ends are core vertices
  std::vector<int> core;                ///< core vertices, in node order
  /// Core vertices traversed: all of them, or none when no edge is scored.
  int reads = 0;
  std::vector<int> core_index_of;       ///< solver vertex -> core index or -1
  /// Solver vertex -> the core vertex its group hangs from (itself for a
  /// core vertex).
  std::vector<int> attach;
  /// Per core vertex and arm c (the root's child order): the neighbour the
  /// arm starts with, and that neighbour's core index or -1.
  std::vector<std::array<int, 3>> arm_neighbor;
  std::vector<std::array<int, 3>> arm_core;
  /// channels[(i * edges + e) * 3 + c]: W(x y_c | y_a y_b) at core vertex i
  /// for the side of e away from it; each gene's counts (the counter's
  /// doubles, exact below 2^53) rounded and summed in gene order.
  std::vector<SkeletonCount> channels;
  std::vector<std::uint8_t> arm_of;         ///< [i * edges + e]: e's arm at i
  /// [i * edges + e]: whether the side away from core vertex i is the side
  /// below e's reference node (1) or its complement (0).
  std::vector<std::uint8_t> away_is_child;
};

/// The skeleton pass of `tree` on `taxa` with at most `read_budget` core
/// vertices (more only when more tie at the largest edge), scoring the edges
/// with at least `min_clade` taxa on their smaller side. Tests pass a
/// smaller `state_reuse_bound_factor` to force Extract/Contract.
SkeletonPass skeleton_pass(const SolverContext& ctx, const SolverTree& tree,
                           const std::vector<int>& taxa, int read_budget,
                           int min_clade, int state_reuse_bound_factor = 4);

/// One destination of a moved side: the edge from the core vertex `vertex`
/// to its neighbour `neighbor`.
struct SkeletonDestination {
  int vertex = -1;
  int neighbor = -1;
  /// The quartet gain of re-attaching the side there, over the quartets
  /// whose outside triple is centred at a core vertex.
  SkeletonCount gain = 0;
};

/// The skeleton gains of one scored edge. The child side is the side below
/// the edge's reference node; each side's best destination has the largest
/// gain (none, with gain zero, when no destination gains), and a side has
/// no data when the other side holds no core vertex.
struct SkeletonEdgeGains {
  bool child_has_data = false;
  bool parent_has_data = false;
  SkeletonDestination child_best, parent_best;
  SkeletonCount gain() const { return child_best.gain + parent_best.gain; }
};

/// The skeleton gains of scored edge `edge`; the profiles, when given,
/// receive every destination of each side.
SkeletonEdgeGains skeleton_edge_gains(
    const SkeletonPass& pass, int edge,
    std::vector<SkeletonDestination>* child_profile = nullptr,
    std::vector<SkeletonDestination>* parent_profile = nullptr);

/// The SPR refinement when options.spr_refine is set: up to
/// options.spr_rounds rounds, each a skeleton pass over the current tree
/// scoring the clades of at least options.spr_min_clade taxa, the scored
/// edges with a positive skeleton gain ranked by it, and for the first
/// options.spr_candidates of them, in order, the exact two-sided re-join of
/// the clade (the edge's smaller side, found again by its taxa; skipped when
/// no longer a clade), applied when it improves the quartet score; NNI
/// follows a round that moved. A round without a move ends the refinement.
/// Returns whether any move was applied. Requires a binary tree.
bool refine_spr(SolverContext& ctx, SolverTree& tree,
                const std::vector<int>& taxa);

// Number of genes displaying the bipartition `side_taxa | rest`, as read by
// the NNI support gate; -1 when the bipartition is trivial.
int nni_split_support_for_tests(const SolverContext& ctx,
                                    const std::vector<int>& side_taxa);
// Runs `moves` in order on `tree` through the NNI support gate: for each, the
// support delta the gate reads before applying it, or INT_MIN when it no
// longer applies.
std::vector<int> nni_support_gate_trace_for_tests(
    const SolverContext& ctx,
    SolverTree& tree,
    const std::vector<NniMove>& moves);

// Every scored NNI alternative of `tree`, with its quartet-score gain.
std::vector<NniMove> nni_scored_moves_for_tests(
    const SolverContext& ctx,
    const SolverTree& tree,
    const std::vector<int>& active_taxa);

// The same for the alternatives of the given central edges only, scored as a
// local pass after a move scores them.
std::vector<NniMove> nni_local_scored_moves_for_tests(
    const SolverContext& ctx,
    const SolverTree& tree,
    const std::vector<int>& active_taxa,
    const std::vector<std::pair<int, int>>& central_edges);

// The set of moves one --nni-batch pass applies to `tree`, chosen from every
// alternative of positive gain, without the support gate.
std::vector<NniMove> nni_batch_for_tests(
    const SolverContext& ctx,
    const SolverTree& tree,
    const std::vector<int>& active_taxa);

}  // namespace qraft::internal
