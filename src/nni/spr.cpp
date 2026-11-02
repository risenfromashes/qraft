// The SPR refinement: rounds of a skeleton pass (skeleton.cpp) that ranks the
// clades by their skeleton gain, then the exact two-sided re-join of the
// first ranked clades.
//
// The re-join cuts the edge above a clade into sides A and B and scores both
// in one call of the K=4 regraft scorer with raw weights. The raw score of
// attaching B on an edge e of A is 3 x (gene quartets with one taxon in B and
// three in A that the joined tree shows) minus a constant, so the best
// re-join and its exact gain over the current join follow from that call.
// The gain is the certificate: a move is applied only when it improves the
// quartet score, whatever the skeleton gain predicted.
//
// On a tree with polytomies (a starting tree kept under --no-divide) the
// scorer counts only vertices with three arms, so a polytomy's share of an
// edge's score is missing. It is the same for every edge in one arm of the
// polytomy, so differences stay exact between edges joined by a path of
// three-arm vertices. Moves are limited to those: the cut edge has two
// three-arm ends, and each side re-joins at an edge reached from its old join
// without passing through a polytomy. The polytomies themselves stay.

#include "solver/solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace qraft::internal {
namespace {

struct CutSide {
  SolverTree tree;
  std::vector<int> taxa;
  RegraftEdge original;  ///< the side's edge where the other side attached
};

// The regraft scorer requires a canonical unrooted component: suppresses its
// degree-two vertices, carrying the old attachment edge through each
// suppression, then compacts the indices.
void suppress_cut_roots(CutSide& side) {
  auto& tree = side.tree;
  bool suppressed = false;
  for (int u = 0; u < static_cast<int>(tree.nodes.size()); ++u) {
    auto& node = tree.nodes[static_cast<std::size_t>(u)];
    if (node.taxon >= 0 || node.neighbors.size() != 2) continue;
    const int a = node.neighbors[0], b = node.neighbors[1];
    if (side.original.a == u || side.original.b == u) {
      side.original.a = a;
      side.original.b = b;
    }
    remove_solver_edge(tree, u, a);
    remove_solver_edge(tree, u, b);
    add_solver_edge(tree, a, b);
    node.active = false;
    suppressed = true;
  }
  if (!suppressed) return;
  std::vector<int> local(tree.nodes.size(), -1);
  SolverTree compact;
  for (int u = 0; u < static_cast<int>(tree.nodes.size()); ++u) {
    if (!tree.nodes[static_cast<std::size_t>(u)].active) continue;
    local[static_cast<std::size_t>(u)] = static_cast<int>(compact.nodes.size());
    compact.nodes.push_back(std::move(tree.nodes[static_cast<std::size_t>(u)]));
  }
  for (auto& node : compact.nodes) {
    for (int& v : node.neighbors) v = local[static_cast<std::size_t>(v)];
  }
  side.original.a = local[static_cast<std::size_t>(side.original.a)];
  side.original.b = local[static_cast<std::size_t>(side.original.b)];
  tree = std::move(compact);
}

// The component of `endpoint` once the edge {endpoint, opposite} is cut, the
// degree-two endpoint omitted, with compact indices so repeated re-joins do
// not accumulate inactive vertices. Traversal and neighbour orders follow
// the input tree.
CutSide cut_side(const SolverTree& tree, int endpoint, int opposite) {
  auto arms = active_neighbors(tree, endpoint);
  arms.erase(std::remove(arms.begin(), arms.end(), opposite), arms.end());
  if (arms.size() != 2) throw std::logic_error("SPR cut endpoint is not binary");
  std::vector<int> order;
  std::vector<int> parent(tree.nodes.size(), -1);
  std::vector<int> local(tree.nodes.size(), -1);
  parent[static_cast<std::size_t>(endpoint)] = opposite;
  std::vector<int> stack{endpoint};
  while (!stack.empty()) {
    const int u = stack.back();
    stack.pop_back();
    if (u != endpoint) {
      local[static_cast<std::size_t>(u)] = static_cast<int>(order.size());
      order.push_back(u);
    }
    for (int v : tree.nodes[static_cast<std::size_t>(u)].neighbors) {
      if (!tree.nodes[static_cast<std::size_t>(v)].active ||
          v == parent[static_cast<std::size_t>(u)]) {
        continue;
      }
      parent[static_cast<std::size_t>(v)] = u;
      stack.push_back(v);
    }
  }
  CutSide side;
  side.tree.nodes.reserve(order.size());
  for (int u : order) {
    SolverTree::Node node;
    node.taxon = tree.nodes[static_cast<std::size_t>(u)].taxon;
    if (node.taxon >= 0) side.taxa.push_back(node.taxon);
    for (int v : tree.nodes[static_cast<std::size_t>(u)].neighbors) {
      if (!tree.nodes[static_cast<std::size_t>(v)].active) continue;
      if (v == endpoint) v = u == arms[0] ? arms[1] : arms[0];
      const int mapped = local[static_cast<std::size_t>(v)];
      if (mapped < 0) throw std::logic_error("SPR component has an external edge");
      node.neighbors.push_back(mapped);
    }
    side.tree.nodes.push_back(std::move(node));
  }
  std::sort(side.taxa.begin(), side.taxa.end());
  side.original.a = local[static_cast<std::size_t>(arms[0])];
  side.original.b = local[static_cast<std::size_t>(arms[1])];
  suppress_cut_roots(side);
  return side;
}

// The edge {a, b} (a < b) of `tree` one of whose sides holds exactly the
// `clade_size` taxa marked in `in_clade` (indexed by taxon id), if the tree
// has one.
bool find_clade_edge(const SolverTree& tree, const std::vector<char>& in_clade,
                     int clade_size, int& a, int& b) {
  const int n = static_cast<int>(tree.nodes.size());
  int root = -1;
  for (int u = 0; u < n && root < 0; ++u) {
    if (tree.nodes[static_cast<std::size_t>(u)].active) root = u;
  }
  if (root < 0) return false;
  std::vector<int> parent(static_cast<std::size_t>(n), -2);
  std::vector<int> taxa_below(static_cast<std::size_t>(n), 0);
  std::vector<int> clade_below(static_cast<std::size_t>(n), 0);
  std::vector<int> order, stack{root};
  parent[static_cast<std::size_t>(root)] = -1;
  while (!stack.empty()) {
    const int u = stack.back();
    stack.pop_back();
    order.push_back(u);
    for (int v : tree.nodes[static_cast<std::size_t>(u)].neighbors) {
      if (!tree.nodes[static_cast<std::size_t>(v)].active ||
          v == parent[static_cast<std::size_t>(u)]) {
        continue;
      }
      parent[static_cast<std::size_t>(v)] = u;
      stack.push_back(v);
    }
  }
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    const int u = *it;
    const int taxon = tree.nodes[static_cast<std::size_t>(u)].taxon;
    if (taxon >= 0) {
      ++taxa_below[static_cast<std::size_t>(u)];
      if (in_clade[static_cast<std::size_t>(taxon)] != 0) {
        ++clade_below[static_cast<std::size_t>(u)];
      }
    }
    const int p = parent[static_cast<std::size_t>(u)];
    if (p >= 0) {
      taxa_below[static_cast<std::size_t>(p)] += taxa_below[static_cast<std::size_t>(u)];
      clade_below[static_cast<std::size_t>(p)] += clade_below[static_cast<std::size_t>(u)];
    }
  }
  const int total = taxa_below[static_cast<std::size_t>(root)];
  for (int u = 0; u < n; ++u) {
    const int p = parent[static_cast<std::size_t>(u)];
    if (p < 0) continue;
    const int below = taxa_below[static_cast<std::size_t>(u)];
    const int in = clade_below[static_cast<std::size_t>(u)];
    if ((in == clade_size && below == clade_size) ||
        (in == 0 && total - below == clade_size)) {
      a = std::min(p, u);
      b = std::max(p, u);
      return true;
    }
  }
  return false;
}

// Marks the nodes of `side` reached from both ends of its old join without
// passing through a vertex of more than three arms (such a vertex is marked
// but not crossed). An edge with both ends marked lies in the old join's
// region.
std::vector<char> join_region(const CutSide& side) {
  const SolverTree& tree = side.tree;
  std::vector<char> marked(tree.nodes.size(), 0);
  std::vector<int> stack{side.original.a, side.original.b};
  for (int u : stack) marked[static_cast<std::size_t>(u)] = 1;
  while (!stack.empty()) {
    const int u = stack.back();
    stack.pop_back();
    if (active_neighbors(tree, u).size() > 3) continue;
    for (int v : tree.nodes[static_cast<std::size_t>(u)].neighbors) {
      if (!tree.nodes[static_cast<std::size_t>(v)].active ||
          marked[static_cast<std::size_t>(v)] != 0) {
        continue;
      }
      marked[static_cast<std::size_t>(v)] = 1;
      stack.push_back(v);
    }
  }
  return marked;
}

bool in_region(const std::vector<char>& region, const RegraftEdge& edge) {
  return edge.a >= 0 && edge.b >= 0 &&
         region[static_cast<std::size_t>(edge.a)] != 0 &&
         region[static_cast<std::size_t>(edge.b)] != 0;
}

bool has_polytomy(const SolverTree& tree) {
  for (int u = 0; u < static_cast<int>(tree.nodes.size()); ++u) {
    const auto& node = tree.nodes[static_cast<std::size_t>(u)];
    if (node.active && node.taxon < 0 && active_neighbors(tree, u).size() > 3) {
      return true;
    }
  }
  return false;
}

const RegraftEdge& original_score(const std::vector<RegraftEdge>& scores,
                                  const RegraftEdge& original) {
  for (const auto& edge : scores) {
    if ((edge.a == original.a && edge.b == original.b) ||
        (edge.a == original.b && edge.b == original.a)) {
      return edge;
    }
  }
  throw std::logic_error("SPR original join absent from regraft scores");
}

// The candidates of one pass: the first `limit` scored edges with a positive
// skeleton gain, largest first, ties by edge order, kept as their query
// nodes on a snapshot of the scored tree. Each candidate's clade is listed
// only when its turn comes, so the candidates take O(n + limit) space.
struct RankedClades {
  QueryTree query;
  ActiveTaxa active;
  std::vector<int> nodes;  ///< per candidate, the query node below its edge
};

RankedClades rank_clades(const SolverContext& ctx, const SolverTree& tree,
                         const std::vector<int>& taxa, int read_budget,
                         int min_clade, int limit) {
  SkeletonPass pass = skeleton_pass(ctx, tree, taxa, read_budget, min_clade);
  stats_add(run_stats().spr_reads, static_cast<long long>(pass.reads));
  std::vector<std::pair<SkeletonCount, int>> ranked;
  for (int e = 0; e < static_cast<int>(pass.edge_nodes.size()); ++e) {
    const SkeletonCount gain = skeleton_edge_gains(pass, e).gain();
    if (gain > 0) ranked.emplace_back(gain, e);
  }
  std::sort(ranked.begin(), ranked.end(), [](const auto& x, const auto& y) {
    if (x.first != y.first) return x.first > y.first;
    return x.second < y.second;
  });
  if (ranked.size() > static_cast<std::size_t>(limit)) {
    ranked.resize(static_cast<std::size_t>(limit));
  }
  RankedClades out;
  out.nodes.reserve(ranked.size());
  for (const auto& [gain, e] : ranked) {
    out.nodes.push_back(pass.edge_nodes[static_cast<std::size_t>(e)]);
  }
  out.query = std::move(pass.query);
  out.active = make_active_taxa(taxa, static_cast<int>(context_taxa(ctx).size()));
  return out;
}

// Marks in `in_clade` the clade of the edge above query node `v`: its
// smaller side, or the side below the edge when the halves are equal.
// Returns the clade's size.
int mark_clade(const RankedClades& ranked, int v, std::vector<char>& in_clade) {
  const QueryTree& query = ranked.query;
  const int n = query.leaf_count;
  const int first = query.subtree_leaf_first[static_cast<std::size_t>(v)];
  const int below = query.subtree_leaf_count[static_cast<std::size_t>(v)];
  const bool clade_below = below <= n - below;
  if (!clade_below) {
    for (int t : ranked.active.taxa) in_clade[static_cast<std::size_t>(t)] = 1;
  }
  for (int i = first; i < first + below; ++i) {
    in_clade[static_cast<std::size_t>(ranked.active.taxa[static_cast<std::size_t>(
        query.leaf_active_order[static_cast<std::size_t>(i)])])] = clade_below ? 1 : 0;
  }
  return clade_below ? below : n - below;
}

}  // namespace

bool refine_spr(SolverContext& ctx, SolverTree& tree,
                const std::vector<int>& taxa) {
  if (!ctx.options.spr_refine) return false;
  const Options& options = ctx.options;
  if (options.spr_reads < 1 || options.spr_min_clade < 2 ||
      options.spr_candidates < 1 || options.spr_rounds < 1) {
    throw std::invalid_argument(
        "SPR refinement requires positive reads, candidates and rounds and "
        "a minimum clade of at least 2 taxa");
  }
  if (taxa.size() < 4) return false;
  const auto started = std::chrono::steady_clock::now();
  progress_note_stage("spr");
  suppress_solver_degree_two(tree);
  // Moves neither create nor remove polytomies, so this holds for the run.
  const bool polytomies = has_polytomy(tree);
  bool changed = false;
  std::vector<char> in_clade(context_taxa(ctx).size(), 0);
  for (int round = 0; round < options.spr_rounds; ++round) {
    stats_add(run_stats().spr_rounds);
    const auto pass_started = std::chrono::steady_clock::now();
    const RankedClades ranked =
        rank_clades(ctx, tree, taxa, options.spr_reads, options.spr_min_clade,
                    options.spr_candidates);
    stats_add_time(run_stats().spr_pass_ns,
                   std::chrono::duration_cast<std::chrono::nanoseconds>(
                       std::chrono::steady_clock::now() - pass_started));
    int moves = 0;
    for (int v : ranked.nodes) {
      // An earlier move of the round may have broken the clade up.
      const int clade_size = mark_clade(ranked, v, in_clade);
      int edge_a = -1, edge_b = -1;
      const bool found = find_clade_edge(tree, in_clade, clade_size, edge_a, edge_b);
      for (int t : taxa) in_clade[static_cast<std::size_t>(t)] = 0;
      if (!found) continue;
      if (polytomies && (active_neighbors(tree, edge_a).size() != 3 ||
                         active_neighbors(tree, edge_b).size() != 3)) {
        continue;
      }
      CutSide left = cut_side(tree, edge_a, edge_b);
      CutSide right = cut_side(tree, edge_b, edge_a);
      std::vector<char> left_region, right_region;
      if (polytomies) {
        left_region = join_region(left);
        right_region = join_region(right);
      }
      std::vector<RegraftEdge> left_scores, right_scores;
      auto best = choose_regraft_edges(ctx, left.tree, right.tree, left.taxa,
          right.taxa, &left_scores, &right_scores, true, true,
          polytomies ? &left_region : nullptr,
          polytomies ? &right_region : nullptr);
      stats_add(run_stats().spr_candidates_scored);
      const auto& original_left = original_score(left_scores, left.original);
      const auto& original_right = original_score(right_scores, right.original);
      // A side without a three-arm vertex reports its first edge, which can
      // lie past a polytomy: such a side keeps its old join.
      if (polytomies && !in_region(left_region, best.first)) best.first = original_left;
      if (polytomies && !in_region(right_region, best.second)) best.second = original_right;
      const double old_left = original_left.raw_score;
      const double old_right = original_right.raw_score;
      const double gain = ((best.first.raw_score - old_left) +
                           (best.second.raw_score - old_right)) / 3.0;
      // Ties and rounding stay below this relative tolerance.
      const double threshold =
          1e-9 * (std::fabs(old_left) + std::fabs(old_right) + 1.0);
      if (!(gain > threshold)) continue;
      tree = connect_at_edges(std::move(left.tree), std::move(right.tree),
                              best.first, best.second);
      stats_add(run_stats().spr_moves);
      if (run_stats_enabled()) {
        const RegraftCount exact_gain =
            (best.first.exact_raw_score - original_left.exact_raw_score) +
            (best.second.exact_raw_score - original_right.exact_raw_score);
        if (best.first.exact_raw_score_valid && original_left.exact_raw_score_valid &&
            best.second.exact_raw_score_valid && original_right.exact_raw_score_valid &&
            exact_gain % 3 == 0) {
          run_stats().spr_gain += exact_gain / 3;
        } else {
          run_stats().spr_gain_exact = false;
        }
      }
      ++moves;
    }
    if (moves == 0) break;
    changed = true;
    refine_nni(ctx, tree, taxa);
    suppress_solver_degree_two(tree);
  }
  stats_add_time(run_stats().spr_ns,
                 std::chrono::duration_cast<std::chrono::nanoseconds>(
                     std::chrono::steady_clock::now() - started));
  return changed;
}

}  // namespace qraft::internal
