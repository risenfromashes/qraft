// The skeleton pass of the SPR refinement.
//
// The gain of re-attaching a clade X at an edge d of the rest of the tree
// is, by the regraft sigma identity, a sum over the vertices u on the path
// from X's edge to d: at each such u, the quartets with one leaf in X and one
// in each arm of u gain when the gene pairs the X leaf with the arm of d and
// lose when it pairs it with the arm X hangs in. Vertices off the path see
// X's edge and d in one arm and contribute nothing. Scoring every clade at
// every vertex would be the oracle for every clade, quadratic. The pass reads
// the vertices of a core only: the endpoints of the heavy edges, those with
// at least G taxa on both sides, at the smallest G whose endpoints fit the
// read budget. For each core vertex one K=4 Count traversal over the clades
// of at least M taxa (the clade coloured 3, the arms 0, 1 and 2) reads the
// three channels of every clade at once.
//
// The heavy edges form a subtree, so the core is connected and the path
// between two core edges lies in it: summed along a core path, the channel
// differences give the exact gain of the move over the quartets whose
// outside triple is centred at a core vertex, which is every changed quartet
// when the clade's own edge is heavy. A clade inside a group (a subtree
// hanging from a core vertex by a light edge) misses the quartets centred at
// the group's vertices on its way out. The best destination's gain ranks the
// clade; the exact two-sided oracle verifies it (spr.cpp).
//
// A vertex of more than three arms (a polytomy kept from the starting tree)
// is never a core vertex and is not counted against the read budget, so core
// paths, and the moves ranked along them, stay within three-arm vertices. A
// clade is scored only when both ends of its edge have three arms, as the
// oracle's cut needs.

#include "solver/solver.hpp"
#include "counting/count.hpp"
#include "counting/query_taxa.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace qraft::internal {
namespace {

// A scored or heavy edge keyed by its unordered solver endpoints.
struct EdgeKey {
  int a, b, edge;
  bool operator<(const EdgeKey& o) const {
    return std::tie(a, b) < std::tie(o.a, o.b);
  }
};

int lookup_edge(const std::vector<EdgeKey>& keys, int a, int b) {
  if (a > b) std::swap(a, b);
  const auto it =
      std::lower_bound(keys.begin(), keys.end(), EdgeKey{a, b, -1});
  return it != keys.end() && it->a == a && it->b == b ? it->edge : -1;
}

// The smallest resolution G >= 2 whose heavy-edge endpoints number at most
// `read_budget`. With h(u) the largest smaller side among the edges at an
// internal vertex u, u is an endpoint at G exactly when h(u) >= G, so
// G = h_(S+1) + 1, the (S+1)-th largest h plus one (vertices tied with it
// leave together; every vertex stays when there are at most S). When more
// than S vertices tie at the largest h, the top edges are read, over the
// budget by the tie.
int core_resolution(const QueryTree& query, int read_budget) {
  const int n = query.leaf_count;
  const int nodes = static_cast<int>(query.nodes.size());
  const auto leaves = [&](int v) {
    return query.subtree_leaf_count[static_cast<std::size_t>(v)];
  };
  std::vector<int> h(static_cast<std::size_t>(nodes), 0);
  for (int v = 0; v < nodes; ++v) {
    if (v == query.root) continue;
    const int side = std::min(leaves(v), n - leaves(v));
    if (query.nodes[static_cast<std::size_t>(v)].taxon_active < 0) {
      h[static_cast<std::size_t>(v)] =
          std::max(h[static_cast<std::size_t>(v)], side);
    }
    const int p = query.parent[static_cast<std::size_t>(v)];
    h[static_cast<std::size_t>(p)] = std::max(h[static_cast<std::size_t>(p)], side);
  }
  std::vector<int> values;
  for (int v = 0; v < nodes; ++v) {
    const int arms = query.child_count[static_cast<std::size_t>(v)] +
                     (v == query.root ? 0 : 1);
    if (query.nodes[static_cast<std::size_t>(v)].taxon_active < 0 && arms == 3) {
      values.push_back(h[static_cast<std::size_t>(v)]);
    }
  }
  std::sort(values.begin(), values.end(), std::greater<int>());
  int resolution = static_cast<int>(values.size()) > read_budget
                       ? values[static_cast<std::size_t>(read_budget)] + 1
                       : 2;
  resolution = std::max(resolution, 2);
  if (!values.empty() && values.front() < resolution) {
    resolution = std::max(2, values.front());
  }
  return resolution;
}

}  // namespace

SkeletonPass skeleton_pass(const SolverContext& ctx, const SolverTree& input,
                           const std::vector<int>& taxa, int read_budget,
                           int min_clade, int state_reuse_bound_factor) {
  if (read_budget < 1 || min_clade < 1 || state_reuse_bound_factor < 1) {
    throw std::invalid_argument(
        "skeleton pass needs a positive read budget, minimum clade and "
        "state reuse bound factor");
  }
  SkeletonPass pass;
  pass.tree = input;
  suppress_solver_degree_two(pass.tree);
  const SolverTree& tree = pass.tree;
  const int taxon_count = static_cast<int>(context_taxa(ctx).size());
  pass.min_clade = min_clade;
  pass.query = query_tree_from_solver_tree(tree, taxa, taxon_count);
  const QueryTree& query = pass.query;
  if (query.root < 0) return pass;
  const int n = query.leaf_count;
  const int nodes = static_cast<int>(query.nodes.size());
  // The solver vertices with more than three arms.
  std::vector<std::uint8_t> polytomy(tree.nodes.size(), 0);
  for (int v = 0; v < nodes; ++v) {
    if (query.nodes[static_cast<std::size_t>(v)].taxon_active >= 0) continue;
    const int arms = query.child_count[static_cast<std::size_t>(v)] +
                     (v == query.root ? 0 : 1);
    const int at = query.nodes[static_cast<std::size_t>(v)].source_node;
    if (arms > 3 && at >= 0) polytomy[static_cast<std::size_t>(at)] = 1;
  }
  const auto leaves = [&](int v) {
    return query.subtree_leaf_count[static_cast<std::size_t>(v)];
  };
  const auto source = [&](int v) {
    return query.nodes[static_cast<std::size_t>(v)].source_node;
  };
  pass.resolution = core_resolution(query, read_budget);

  // Scored edges, heavy edges and the core.
  const int tree_nodes = static_cast<int>(tree.nodes.size());
  std::vector<std::uint8_t> is_core(static_cast<std::size_t>(tree_nodes), 0);
  std::vector<EdgeKey> keys;
  for (int v = 0; v < nodes; ++v) {
    if (v == query.root) continue;
    const int below = leaves(v);
    const int side = std::min(below, n - below);
    const int a = source(v);
    const int b = source(query.parent[static_cast<std::size_t>(v)]);
    if (a < 0 || b < 0) {
      throw std::logic_error("skeleton pass: query edge without source nodes");
    }
    if (side >= pass.resolution) {
      if (polytomy[static_cast<std::size_t>(a)] == 0) is_core[static_cast<std::size_t>(a)] = 1;
      if (polytomy[static_cast<std::size_t>(b)] == 0) is_core[static_cast<std::size_t>(b)] = 1;
    }
    if (side < min_clade) continue;
    if (polytomy[static_cast<std::size_t>(a)] != 0 ||
        polytomy[static_cast<std::size_t>(b)] != 0) {
      continue;
    }
    const int edge = static_cast<int>(pass.edge_nodes.size());
    pass.edge_nodes.push_back(v);
    pass.edge_child_source.push_back(a);
    pass.edge_parent_source.push_back(b);
    pass.edge_min_side.push_back(side);
    keys.push_back(EdgeKey{std::min(a, b), std::max(a, b), edge});
  }
  std::sort(keys.begin(), keys.end());
  const int edge_count = static_cast<int>(pass.edge_nodes.size());
  pass.core_index_of.assign(static_cast<std::size_t>(tree_nodes), -1);
  for (int u = 0; u < tree_nodes; ++u) {
    if (is_core[static_cast<std::size_t>(u)] == 0) continue;
    if (tree.nodes[static_cast<std::size_t>(u)].taxon >= 0) {
      throw std::logic_error("skeleton pass: a leaf is a core vertex");
    }
    pass.core_index_of[static_cast<std::size_t>(u)] =
        static_cast<int>(pass.core.size());
    pass.core.push_back(u);
  }
  const int core_count = static_cast<int>(pass.core.size());
  pass.edge_heavy.assign(static_cast<std::size_t>(edge_count), 0);
  for (int e = 0; e < edge_count; ++e) {
    pass.edge_heavy[static_cast<std::size_t>(e)] =
        is_core[static_cast<std::size_t>(
            pass.edge_child_source[static_cast<std::size_t>(e)])] != 0 &&
                is_core[static_cast<std::size_t>(
                    pass.edge_parent_source[static_cast<std::size_t>(e)])] != 0
            ? 1
            : 0;
  }
  // Each group's attachment: a breadth-first search from the core vertices
  // over the other vertices.
  pass.attach.assign(static_cast<std::size_t>(tree_nodes), -1);
  {
    std::vector<int> frontier;
    for (int u : pass.core) {
      pass.attach[static_cast<std::size_t>(u)] = u;
      frontier.push_back(u);
    }
    for (std::size_t i = 0; i < frontier.size(); ++i) {
      const int u = frontier[i];
      // A polytomy is attached but not crossed: past it, no core move reaches.
      if (polytomy[static_cast<std::size_t>(u)] != 0) continue;
      for (int v : tree.nodes[static_cast<std::size_t>(u)].neighbors) {
        if (!tree.nodes[static_cast<std::size_t>(v)].active ||
            pass.attach[static_cast<std::size_t>(v)] >= 0) {
          continue;
        }
        pass.attach[static_cast<std::size_t>(v)] =
            pass.attach[static_cast<std::size_t>(u)];
        frontier.push_back(v);
      }
    }
  }
  pass.arm_neighbor.assign(static_cast<std::size_t>(core_count),
                           std::array<int, 3>{-1, -1, -1});
  pass.arm_core.assign(static_cast<std::size_t>(core_count),
                       std::array<int, 3>{-1, -1, -1});
  const std::size_t cells =
      static_cast<std::size_t>(core_count) * static_cast<std::size_t>(edge_count);
  pass.channels.assign(cells * 3, 0);
  pass.arm_of.assign(cells, 0);
  pass.away_is_child.assign(cells, 0);
  if (core_count == 0 || edge_count == 0 || n < 4) return pass;
  pass.reads = core_count;

  const ActiveTaxa active = make_active_taxa(taxa, taxon_count);
  const OwnerProjection owners = make_owner_projection(ctx, active.active_by_taxon);
  const std::vector<GeneIndex>& genes = scoring_genes(ctx);
  const int gene_count = static_cast<int>(genes.size());
  const std::size_t row_width = static_cast<std::size_t>(edge_count) * 3;
  // One row per gene; a cell is written once per traversal (the counter's
  // double), and the rows are summed exactly in gene order.
  std::vector<double> partial(static_cast<std::size_t>(gene_count) * row_width,
                              0.0);
  std::vector<std::uint8_t> score_u, rest;
  std::vector<int> node_to_edge_u;
  std::vector<std::uint8_t> background(
      static_cast<std::size_t>(owners.owner_count), 0);
  for (int i = 0; i < core_count; ++i) {
    const int u = pass.core[static_cast<std::size_t>(i)];
    const QueryTree qu = query_tree_rooted_at_vertex(tree, taxa, taxon_count, u);
    const int qn = static_cast<int>(qu.nodes.size());
    if (qu.root < 0 || qu.child_count[static_cast<std::size_t>(qu.root)] != 3) {
      throw std::logic_error("skeleton pass: could not root at a core vertex");
    }
    const int root_first = qu.child_first[static_cast<std::size_t>(qu.root)];
    rest.assign(static_cast<std::size_t>(qn), 0);
    for (int c = 0; c < 3; ++c) {
      const int child = qu.child_nodes[static_cast<std::size_t>(root_first + c)];
      rest[static_cast<std::size_t>(child)] = static_cast<std::uint8_t>(c);
      const int nb = qu.nodes[static_cast<std::size_t>(child)].source_node;
      pass.arm_neighbor[static_cast<std::size_t>(i)][static_cast<std::size_t>(c)] = nb;
      pass.arm_core[static_cast<std::size_t>(i)][static_cast<std::size_t>(c)] =
          nb >= 0 ? pass.core_index_of[static_cast<std::size_t>(nb)] : -1;
    }
    for (int v : qu.preorder) {
      const int p = qu.parent[static_cast<std::size_t>(v)];
      if (v == qu.root || p == qu.root || p < 0) continue;
      rest[static_cast<std::size_t>(v)] = rest[static_cast<std::size_t>(p)];
    }
    score_u.assign(static_cast<std::size_t>(qn), 0);
    node_to_edge_u.assign(static_cast<std::size_t>(qn), -1);
    for (int v = 0; v < qn; ++v) {
      if (v == qu.root) continue;
      const int p = qu.parent[static_cast<std::size_t>(v)];
      const int a = qu.nodes[static_cast<std::size_t>(v)].source_node;
      if (a < 0) throw std::logic_error("skeleton pass: query node without a source");
      const int b = p == qu.root ? u : qu.nodes[static_cast<std::size_t>(p)].source_node;
      const int e = lookup_edge(keys, a, b);
      if (e < 0) continue;
      score_u[static_cast<std::size_t>(v)] = 1;
      node_to_edge_u[static_cast<std::size_t>(v)] = e;
      const std::size_t cell = static_cast<std::size_t>(i) *
                                   static_cast<std::size_t>(edge_count) +
                               static_cast<std::size_t>(e);
      pass.arm_of[cell] = rest[static_cast<std::size_t>(v)];
      pass.away_is_child[cell] =
          a == pass.edge_child_source[static_cast<std::size_t>(e)] ? 1 : 0;
    }
    for (int v = 0; v < qn; ++v) {
      const int owner = qu.nodes[static_cast<std::size_t>(v)].taxon_active;
      if (owner >= 0 && owner < owners.owner_count) {
        background[static_cast<std::size_t>(owner)] =
            rest[static_cast<std::size_t>(v)];
      }
    }
    std::fill(partial.begin(), partial.end(), 0.0);
    parallel_for_dynamic(0, gene_count, [&](int g) {
      double* row = partial.data() + static_cast<std::size_t>(g) * row_width;
      const GeneIndex& gene = genes[static_cast<std::size_t>(g)];
      if (!gene_has_owner_quartet(gene, owners)) return;
      HdtColorPolicy colors;
      colors.default_current = 3;
      colors.default_background = 0;
      colors.current_by_owner.assign(
          static_cast<std::size_t>(owners.owner_count), 3);
      colors.background_by_owner = background;
      HdtTraversalStats hdt_stats;
      HdtExtractWorkspace hdt_workspace;
      HdtInitialStateResult initial = hdt_build_initial_state(
          build_base_component_graph(gene, owners, 4, colors));
      hdt_account_initial_state(initial, hdt_stats);
      count_clades_colored(
          std::move(initial.state), qu, qu.leaf_active_order, score_u, rest,
          hdt_workspace, hdt_stats,
          [&](int v, HdtState& state) {
            const int e = node_to_edge_u[static_cast<std::size_t>(v)];
            if (e < 0) return;
            // The root counts are W(01|23), W(02|13), W(03|12); the clade
            // is colour 3, so channel c pairs it with arm c.
            const std::array<double, 3> counts = hdt_root_counts(state);
            double* cell = row + static_cast<std::size_t>(e) * 3;
            cell[0] += counts[2];
            cell[1] += counts[1];
            cell[2] += counts[0];
          },
          state_reuse_bound_factor);
    }, 1);
    SkeletonCount* out = pass.channels.data() + static_cast<std::size_t>(i) * row_width;
    for (int g = 0; g < gene_count; ++g) {
      const double* row = partial.data() + static_cast<std::size_t>(g) * row_width;
      for (std::size_t j = 0; j < row_width; ++j) {
        out[j] += static_cast<SkeletonCount>(std::nearbyint(row[j]));
      }
    }
  }
  return pass;
}

SkeletonEdgeGains skeleton_edge_gains(const SkeletonPass& pass, int edge,
                                      std::vector<SkeletonDestination>* child_profile,
                                      std::vector<SkeletonDestination>* parent_profile) {
  SkeletonEdgeGains out;
  const std::size_t edge_count = pass.edge_nodes.size();
  if (edge < 0 || static_cast<std::size_t>(edge) >= edge_count) return out;
  const auto cell = [&](int i) {
    return static_cast<std::size_t>(i) * edge_count + static_cast<std::size_t>(edge);
  };
  const auto q = [&](int i, int c) {
    return pass.channels[cell(i) * 3 + static_cast<std::size_t>(c)];
  };
  // The best move of the side away from the core vertices whose
  // away_is_child equals `away`, into the side holding them: a search over
  // the core from the vertex the side enters it at, summing the channel
  // differences along each path.
  const auto run = [&](std::uint8_t away, SkeletonDestination& best,
                       std::vector<SkeletonDestination>* profile) {
    best = SkeletonDestination{};
    const int end = away ? pass.edge_parent_source[static_cast<std::size_t>(edge)]
                         : pass.edge_child_source[static_cast<std::size_t>(edge)];
    int entry = pass.core_index_of[static_cast<std::size_t>(end)];
    if (entry < 0) {
      const int at = pass.attach[static_cast<std::size_t>(end)];
      entry = at >= 0 ? pass.core_index_of[static_cast<std::size_t>(at)] : -1;
    }
    if (entry < 0 || pass.away_is_child[cell(entry)] != away) return false;
    std::vector<std::tuple<int, int, SkeletonCount>> stack{
        {entry, pass.arm_of[cell(entry)], SkeletonCount{0}}};
    while (!stack.empty()) {
      const auto [i, from, acc] = stack.back();
      stack.pop_back();
      for (int c = 0; c < 3; ++c) {
        if (c == from) continue;
        const SkeletonCount gain = acc + q(i, c) - q(i, from);
        SkeletonDestination dest;
        dest.vertex = pass.core[static_cast<std::size_t>(i)];
        dest.neighbor =
            pass.arm_neighbor[static_cast<std::size_t>(i)][static_cast<std::size_t>(c)];
        dest.gain = gain;
        if (profile != nullptr) profile->push_back(dest);
        if (gain > best.gain) best = dest;
        const int j = pass.arm_core[static_cast<std::size_t>(i)][static_cast<std::size_t>(c)];
        if (j >= 0 && pass.away_is_child[cell(j)] == away) {
          stack.emplace_back(j, pass.arm_of[cell(j)], gain);
        }
      }
    }
    return true;
  };
  out.child_has_data = run(1, out.child_best, child_profile);
  out.parent_has_data = run(0, out.parent_best, parent_profile);
  return out;
}

}  // namespace qraft::internal
