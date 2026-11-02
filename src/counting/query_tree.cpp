#include "solver/solver.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace qraft::internal {

// The owner of each query leaf, in leaf_active_order; -1 where the leaf's
// active taxon has no owner.
std::vector<int> query_leaf_owner_order(
    const QueryTree& query,
    const std::vector<int>& query_active_to_owner) {
  std::vector<int> owners(query.leaf_active_order.size(), -1);
  for (int pos = 0; pos < static_cast<int>(query.leaf_active_order.size());
       ++pos) {
    const int active = query.leaf_active_order[static_cast<std::size_t>(pos)];
    if (active >= 0 &&
        active < static_cast<int>(query_active_to_owner.size())) {
      owners[static_cast<std::size_t>(pos)] =
          query_active_to_owner[static_cast<std::size_t>(active)];
    }
  }
  return owners;
}

namespace {

// The shallower of two nodes, the smaller id on equal depth; -1 is empty.
int shallower_starting_tree_node(const StartingTreeIndex& index, int a, int b) {
  if (a < 0) return b;
  if (b < 0) return a;
  const int a_depth = index.depth[static_cast<std::size_t>(a)];
  const int b_depth = index.depth[static_cast<std::size_t>(b)];
  if (a_depth != b_depth) return a_depth < b_depth ? a : b;
  return std::min(a, b);
}

// Preorder intervals, an Euler tour and a segment tree over it that returns
// the shallowest node of a tour range: O(n) space, O(n) construction and
// O(log n) per LCA query.
std::shared_ptr<const StartingTreeIndex> build_starting_tree_index(
    const StartingTree& starting_tree) {
  constexpr std::size_t kMaxEulerNodes =
      (static_cast<std::size_t>(std::numeric_limits<int>::max()) + 1u) / 2u;
  if (starting_tree.nodes.size() > kMaxEulerNodes ||
      starting_tree.taxon_to_leaf.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("starting-tree index exceeds its integer domain");
  }
  const int node_count = static_cast<int>(starting_tree.nodes.size());
  if (starting_tree.root < 0 || starting_tree.root >= node_count) {
    throw std::logic_error("cannot index an empty starting tree");
  }

  auto index = std::make_shared<StartingTreeIndex>();
  index->parent.assign(static_cast<std::size_t>(node_count), -1);
  index->depth.assign(static_cast<std::size_t>(node_count), 0);
  index->tin.assign(static_cast<std::size_t>(node_count), -1);
  index->tout.assign(static_cast<std::size_t>(node_count), -1);
  index->first_occurrence.assign(static_cast<std::size_t>(node_count), -1);
  index->dfs_order.reserve(static_cast<std::size_t>(node_count));
  index->euler.reserve(static_cast<std::size_t>(node_count) * 2u - 1u);

  int taxon_slots = static_cast<int>(starting_tree.taxon_to_leaf.size());
  for (const StartingTree::Node& node : starting_tree.nodes) {
    if (node.taxon == std::numeric_limits<int>::max()) {
      throw std::length_error("starting-tree taxon exceeds its integer domain");
    }
    if (node.taxon >= 0) taxon_slots = std::max(taxon_slots, node.taxon + 1);
  }
  index->taxon_to_leaf.assign(static_cast<std::size_t>(taxon_slots), -1);

  struct DfsFrame {
    int node = -1;
    int next_child = 0;
  };
  std::vector<DfsFrame> stack{{starting_tree.root, 0}};
  index->parent[static_cast<std::size_t>(starting_tree.root)] = starting_tree.root;
  while (!stack.empty()) {
    DfsFrame& frame = stack.back();
    const int node = frame.node;
    if (frame.next_child == 0) {
      index->tin[static_cast<std::size_t>(node)] =
          static_cast<int>(index->dfs_order.size());
      index->dfs_order.push_back(node);
      index->first_occurrence[static_cast<std::size_t>(node)] =
          static_cast<int>(index->euler.size());
      index->euler.push_back(node);
      const int taxon = starting_tree.nodes[static_cast<std::size_t>(node)].taxon;
      if (taxon >= 0) {
        int& leaf = index->taxon_to_leaf[static_cast<std::size_t>(taxon)];
        if (leaf >= 0 && leaf != node) {
          throw std::logic_error("duplicate taxon in starting-tree index");
        }
        leaf = node;
      }
    }

    const auto& children =
        starting_tree.nodes[static_cast<std::size_t>(node)].children;
    if (frame.next_child < static_cast<int>(children.size())) {
      const int child = children[static_cast<std::size_t>(frame.next_child++)];
      if (child < 0 || child >= node_count ||
          index->parent[static_cast<std::size_t>(child)] >= 0) {
        throw std::logic_error("starting-tree children do not form a tree");
      }
      index->parent[static_cast<std::size_t>(child)] = node;
      index->depth[static_cast<std::size_t>(child)] =
          index->depth[static_cast<std::size_t>(node)] + 1;
      stack.push_back({child, 0});
      continue;
    }

    index->tout[static_cast<std::size_t>(node)] =
        static_cast<int>(index->dfs_order.size());
    stack.pop_back();
    if (!stack.empty()) index->euler.push_back(stack.back().node);
  }
  if (static_cast<int>(index->dfs_order.size()) != node_count) {
    throw std::logic_error("starting tree is disconnected");
  }

  index->rmq_base = 1;
  while (static_cast<std::size_t>(index->rmq_base) < index->euler.size()) {
    if (index->rmq_base > std::numeric_limits<int>::max() / 2) {
      throw std::length_error("starting-tree RMQ exceeds its integer domain");
    }
    index->rmq_base *= 2;
  }
  index->rmq_nodes.assign(static_cast<std::size_t>(index->rmq_base) * 2u, -1);
  for (int i = 0; i < static_cast<int>(index->euler.size()); ++i) {
    index->rmq_nodes[static_cast<std::size_t>(index->rmq_base + i)] =
        index->euler[static_cast<std::size_t>(i)];
  }
  for (int node = index->rmq_base - 1; node > 0; --node) {
    index->rmq_nodes[static_cast<std::size_t>(node)] =
        shallower_starting_tree_node(
            *index,
            index->rmq_nodes[static_cast<std::size_t>(2 * node)],
            index->rmq_nodes[static_cast<std::size_t>(2 * node + 1)]);
  }
  stats_add(run_stats().starting_tree_index_builds);
  stats_add(run_stats().starting_tree_index_nodes,
                        node_count);
  return index;
}

// Builds the index once per starting tree and shares it with every copy of
// the tree. The first caller moves the state from 0 (empty) to 1 (building)
// and publishes at 2; concurrent callers wait. A failed build returns the
// state to 0, so a later caller retries.
const StartingTreeIndex& shared_starting_tree_index(
    const StartingTree& starting_tree) {
  const std::shared_ptr<StartingTreeIndexCache> cache = starting_tree.index_cache;
  if (cache == nullptr) {
    throw std::logic_error("starting-tree index cache is missing");
  }
  for (;;) {
    const int state = cache->state.load(std::memory_order_acquire);
    if (state == 2) {
      if (cache->value == nullptr) {
        throw std::logic_error("published starting-tree index is missing");
      }
      return *cache->value;
    }
    if (state == 1) {
      cache->state.wait(1, std::memory_order_acquire);
      continue;
    }
    if (state != 0) {
      throw std::logic_error("invalid starting-tree index publication state");
    }
    int expected = 0;
    if (!cache->state.compare_exchange_strong(
            expected, 1, std::memory_order_acq_rel,
            std::memory_order_acquire)) {
      continue;
    }
    try {
      cache->value = build_starting_tree_index(starting_tree);
      cache->state.store(2, std::memory_order_release);
      cache->state.notify_all();
      return *cache->value;
    } catch (...) {
      cache->state.store(0, std::memory_order_release);
      cache->state.notify_all();
      throw;
    }
  }
}

bool starting_tree_is_ancestor(const StartingTreeIndex& index,
                           int ancestor,
                           int node) {
  return ancestor >= 0 && node >= 0 &&
         index.tin[static_cast<std::size_t>(ancestor)] <=
             index.tin[static_cast<std::size_t>(node)] &&
         index.tout[static_cast<std::size_t>(node)] <=
             index.tout[static_cast<std::size_t>(ancestor)];
}

int starting_tree_lca(const StartingTreeIndex& index, int a, int b) {
  if (a < 0 || b < 0 ||
      a >= static_cast<int>(index.first_occurrence.size()) ||
      b >= static_cast<int>(index.first_occurrence.size())) {
    throw std::logic_error("invalid starting-tree LCA endpoint");
  }
  int left = index.first_occurrence[static_cast<std::size_t>(a)];
  int right = index.first_occurrence[static_cast<std::size_t>(b)];
  if (left > right) std::swap(left, right);
  left += index.rmq_base;
  right += index.rmq_base + 1;
  int best = -1;
  while (left < right) {
    if ((left & 1) != 0) {
      best = shallower_starting_tree_node(
          index, best, index.rmq_nodes[static_cast<std::size_t>(left++)]);
    }
    if ((right & 1) != 0) {
      best = shallower_starting_tree_node(
          index, best,
          index.rmq_nodes[static_cast<std::size_t>(--right)]);
    }
    left /= 2;
    right /= 2;
  }
  if (best < 0) throw std::logic_error("empty starting-tree LCA range");
  return best;
}

struct QueryBuilder {
  std::vector<QueryTree::Node> nodes;
  std::vector<std::vector<int>> adj;

  int add_node(int active_taxon, int source_node = -1) {
    const int node = static_cast<int>(nodes.size());
    nodes.push_back(QueryTree::Node{-1, 0, active_taxon, source_node});
    adj.emplace_back();
    return node;
  }

  void add_edge(int a, int b) {
    if (a < 0 || b < 0 || a == b) return;
    adj[static_cast<std::size_t>(a)].push_back(b);
    adj[static_cast<std::size_t>(b)].push_back(a);
  }
};

// Stores the builder's adjacency as CSR with each neighbour list sorted and
// deduplicated, roots the tree at `root`, and fills the traversal orders,
// child lists and per-subtree data. Leaves are numbered in preorder, so the
// leaves of every subtree are one contiguous range of leaf_active_order.
// subtree_size counts active leaves, and the heavy child is the first child
// with the most of them.
void finalize_query_tree(QueryBuilder&& builder, int root, int active_count,
                         QueryTree& out) {
  out = QueryTree{};
  out.root = root;
  out.leaf_count = active_count;
  out.nodes = std::move(builder.nodes);
  out.active_to_node.assign(static_cast<std::size_t>(std::max(0, active_count)),
                            -1);
  for (int node = 0; node < static_cast<int>(out.nodes.size()); ++node) {
    const int active = out.nodes[static_cast<std::size_t>(node)].taxon_active;
    if (active >= 0 && active < active_count) {
      out.active_to_node[static_cast<std::size_t>(active)] = node;
    }
  }

  std::vector<int> offset(out.nodes.size() + 1, 0);
  for (int node = 0; node < static_cast<int>(out.nodes.size()); ++node) {
    std::vector<int>& nbs = builder.adj[static_cast<std::size_t>(node)];
    std::sort(nbs.begin(), nbs.end());
    nbs.erase(std::unique(nbs.begin(), nbs.end()), nbs.end());
    offset[static_cast<std::size_t>(node + 1)] =
        offset[static_cast<std::size_t>(node)] + static_cast<int>(nbs.size());
    out.nodes[static_cast<std::size_t>(node)].first_adj =
        offset[static_cast<std::size_t>(node)];
    out.nodes[static_cast<std::size_t>(node)].degree =
        static_cast<int>(nbs.size());
  }
  out.adjacency.assign(static_cast<std::size_t>(offset.back()), -1);
  for (int node = 0; node < static_cast<int>(out.nodes.size()); ++node) {
    const int begin = offset[static_cast<std::size_t>(node)];
    const std::vector<int>& nbs = builder.adj[static_cast<std::size_t>(node)];
    std::copy(nbs.begin(), nbs.end(),
              out.adjacency.begin() + static_cast<std::ptrdiff_t>(begin));
  }

  const int n = static_cast<int>(out.nodes.size());
  out.parent.assign(static_cast<std::size_t>(n), -1);
  out.child_first.assign(static_cast<std::size_t>(n), 0);
  out.child_count.assign(static_cast<std::size_t>(n), 0);
  out.tin.assign(static_cast<std::size_t>(n), 0);
  out.tout.assign(static_cast<std::size_t>(n), 0);
  out.subtree_size.assign(static_cast<std::size_t>(n), 0);
  out.subtree_leaf_first.assign(static_cast<std::size_t>(n), 0);
  out.subtree_leaf_count.assign(static_cast<std::size_t>(n), 0);
  out.heavy_child.assign(static_cast<std::size_t>(n), -1);
  if (root < 0 || root >= n) return;

  std::vector<std::pair<int, int>> stack{{root, 0}};
  out.parent[static_cast<std::size_t>(root)] = root;
  while (!stack.empty()) {
    const int node = stack.back().first;
    int& next = stack.back().second;
    if (next == 0) {
      out.tin[static_cast<std::size_t>(node)] =
          static_cast<int>(out.preorder.size());
      out.preorder.push_back(node);
    }
    const int first = out.nodes[static_cast<std::size_t>(node)].first_adj;
    const int degree = out.nodes[static_cast<std::size_t>(node)].degree;
    if (next < degree) {
      const int nb = out.adjacency[static_cast<std::size_t>(first + next++)];
      if (nb == out.parent[static_cast<std::size_t>(node)]) continue;
      out.parent[static_cast<std::size_t>(nb)] = node;
      stack.push_back({nb, 0});
      continue;
    }
    out.tout[static_cast<std::size_t>(node)] =
        static_cast<int>(out.preorder.size());
    out.postorder.push_back(node);
    stack.pop_back();
  }

  std::vector<int> child_offset(static_cast<std::size_t>(n + 1), 0);
  for (int node = 0; node < n; ++node) {
    const int first = out.nodes[static_cast<std::size_t>(node)].first_adj;
    const int degree = out.nodes[static_cast<std::size_t>(node)].degree;
    for (int i = 0; i < degree; ++i) {
      const int nb = out.adjacency[static_cast<std::size_t>(first + i)];
      if (out.parent[static_cast<std::size_t>(nb)] == node) {
        ++child_offset[static_cast<std::size_t>(node + 1)];
      }
    }
  }
  for (int node = 0; node < n; ++node) {
    child_offset[static_cast<std::size_t>(node + 1)] +=
        child_offset[static_cast<std::size_t>(node)];
    out.child_first[static_cast<std::size_t>(node)] =
        child_offset[static_cast<std::size_t>(node)];
    out.child_count[static_cast<std::size_t>(node)] =
        child_offset[static_cast<std::size_t>(node + 1)] -
        child_offset[static_cast<std::size_t>(node)];
  }
  out.child_nodes.assign(static_cast<std::size_t>(child_offset.back()), -1);
  std::vector<int> child_cursor = child_offset;
  for (int node = 0; node < n; ++node) {
    const int first = out.nodes[static_cast<std::size_t>(node)].first_adj;
    const int degree = out.nodes[static_cast<std::size_t>(node)].degree;
    for (int i = 0; i < degree; ++i) {
      const int nb = out.adjacency[static_cast<std::size_t>(first + i)];
      if (out.parent[static_cast<std::size_t>(nb)] != node) continue;
      out.child_nodes[static_cast<std::size_t>(
          child_cursor[static_cast<std::size_t>(node)]++)] = nb;
    }
  }

  std::vector<int> leaf_position(static_cast<std::size_t>(n), -1);
  out.leaf_active_order.reserve(static_cast<std::size_t>(active_count));
  for (int node : out.preorder) {
    const int active = out.nodes[static_cast<std::size_t>(node)].taxon_active;
    if (active < 0) continue;
    leaf_position[static_cast<std::size_t>(node)] =
        static_cast<int>(out.leaf_active_order.size());
    out.leaf_active_order.push_back(active);
  }

  for (int node : out.postorder) {
    int size = out.nodes[static_cast<std::size_t>(node)].taxon_active >= 0 ? 1 : 0;
    int heavy = -1;
    int heavy_size = -1;
    int leaf_first = leaf_position[static_cast<std::size_t>(node)];
    int leaf_count = leaf_first >= 0 ? 1 : 0;
    const int child_begin = out.child_first[static_cast<std::size_t>(node)];
    const int child_end =
        child_begin + out.child_count[static_cast<std::size_t>(node)];
    for (int child_pos = child_begin; child_pos < child_end; ++child_pos) {
      const int nb = out.child_nodes[static_cast<std::size_t>(child_pos)];
      const int child_size = out.subtree_size[static_cast<std::size_t>(nb)];
      size += child_size;
      if (child_size > heavy_size) {
        heavy_size = child_size;
        heavy = nb;
      }
      const int child_leaf_count =
          out.subtree_leaf_count[static_cast<std::size_t>(nb)];
      if (child_leaf_count > 0) {
        const int child_leaf_first =
            out.subtree_leaf_first[static_cast<std::size_t>(nb)];
        leaf_first = leaf_count == 0
                         ? child_leaf_first
                         : std::min(leaf_first, child_leaf_first);
        leaf_count += child_leaf_count;
      }
    }
    out.subtree_size[static_cast<std::size_t>(node)] = size;
    out.subtree_leaf_first[static_cast<std::size_t>(node)] =
        leaf_count > 0 ? leaf_first : 0;
    out.subtree_leaf_count[static_cast<std::size_t>(node)] = leaf_count;
    out.heavy_child[static_cast<std::size_t>(node)] = heavy;
  }
}

}  // namespace

// The starting tree induced on the active taxa, as a virtual tree: the active
// leaves in DFS order plus the LCAs of DFS-adjacent pairs, rooted at their
// common ancestor. It has O(m) nodes and costs O(m log n) without touching
// the rest of the starting tree. Active taxa absent from the starting tree
// are skipped.
QueryTree induce_query_tree(
    const StartingTree& starting_tree,
    const ActiveTaxa& active) {
  QueryTree out;
  if (starting_tree.root < 0 || starting_tree.nodes.empty() || active.taxa.empty()) {
    return out;
  }
  const StartingTreeIndex& index = shared_starting_tree_index(starting_tree);
  std::vector<int> selected_leaves;
  selected_leaves.reserve(active.taxa.size());
  for (int taxon : active.taxa) {
    if (taxon < 0 ||
        taxon >= static_cast<int>(index.taxon_to_leaf.size())) {
      continue;
    }
    const int leaf = index.taxon_to_leaf[static_cast<std::size_t>(taxon)];
    if (leaf >= 0) selected_leaves.push_back(leaf);
  }
  if (selected_leaves.empty()) return out;
  const auto dfs_less = [&](int a, int b) {
    const int a_tin = index.tin[static_cast<std::size_t>(a)];
    const int b_tin = index.tin[static_cast<std::size_t>(b)];
    return a_tin != b_tin ? a_tin < b_tin : a < b;
  };
  std::sort(selected_leaves.begin(), selected_leaves.end(), dfs_less);
  selected_leaves.erase(
      std::unique(selected_leaves.begin(), selected_leaves.end()),
      selected_leaves.end());
  stats_add(run_stats().query_tree_selected_leaves,
                        static_cast<long long>(selected_leaves.size()));

  std::vector<int> virtual_nodes = selected_leaves;
  virtual_nodes.reserve(selected_leaves.size() * 2u);
  stats_add(
      run_stats().query_tree_lca_queries,
      static_cast<long long>(selected_leaves.size() - 1u));
  for (int i = 1; i < static_cast<int>(selected_leaves.size()); ++i) {
    virtual_nodes.push_back(starting_tree_lca(
        index, selected_leaves[static_cast<std::size_t>(i - 1)],
        selected_leaves[static_cast<std::size_t>(i)]));
  }
  std::sort(virtual_nodes.begin(), virtual_nodes.end(), dfs_less);
  virtual_nodes.erase(std::unique(virtual_nodes.begin(), virtual_nodes.end()),
                      virtual_nodes.end());
  stats_add(run_stats().query_tree_retained_nodes,
                        static_cast<long long>(virtual_nodes.size()));

  QueryBuilder builder;
  for (int cnode : virtual_nodes) {
    const int taxon = starting_tree.nodes[static_cast<std::size_t>(cnode)].taxon;
    builder.add_node(taxon >= 0 ? active_index(active, taxon) : -1);
  }

  // In DFS order, each node's parent is the nearest stacked ancestor.
  std::vector<int> virtual_stack;
  virtual_stack.reserve(virtual_nodes.size());
  for (int node = 0; node < static_cast<int>(virtual_nodes.size()); ++node) {
    const int cnode = virtual_nodes[static_cast<std::size_t>(node)];
    while (!virtual_stack.empty() &&
           !starting_tree_is_ancestor(
               index,
               virtual_nodes[static_cast<std::size_t>(virtual_stack.back())],
               cnode)) {
      virtual_stack.pop_back();
    }
    if (!virtual_stack.empty()) {
      builder.add_edge(virtual_stack.back(), node);
    } else if (node != 0) {
      throw std::logic_error("starting-tree virtual tree has multiple roots");
    }
    virtual_stack.push_back(node);
  }

  const int root = virtual_nodes.empty() ? -1 : 0;
  finalize_query_tree(std::move(builder), root,
                      static_cast<int>(active.taxa.size()), out);
  return out;
}

// The subtree of a solver tree spanning the active taxa, over active nodes
// only, with unlabelled degree-two vertices suppressed. It is rooted at its
// first internal vertex in node order, or at node 0 when it has none.
QueryTree query_tree_from_solver_tree(
    const SolverTree& tree,
    const std::vector<int>& active_taxa,
    int taxon_count) {
  QueryTree out;
  if (tree.nodes.empty() || active_taxa.empty()) return out;
  const TaxonIndex active_by_taxon =
      make_taxon_index(active_taxa, taxon_count);
  const int n = static_cast<int>(tree.nodes.size());
  std::vector<unsigned char> active_leaf(static_cast<std::size_t>(n), 0);
  std::vector<int> active_leaf_node;
  for (int node = 0; node < n; ++node) {
    if (!tree.nodes[static_cast<std::size_t>(node)].active) continue;
    const int taxon = tree.nodes[static_cast<std::size_t>(node)].taxon;
    if (taxon >= 0 && active_by_taxon.find(taxon) >= 0) {
      active_leaf[static_cast<std::size_t>(node)] = 1;
      active_leaf_node.push_back(node);
    }
  }
  if (active_leaf_node.empty()) return out;

  std::vector<int> parent(static_cast<std::size_t>(n), -1);
  std::vector<int> order;
  order.reserve(static_cast<std::size_t>(n));
  const int traversal_root = active_leaf_node.front();
  parent[static_cast<std::size_t>(traversal_root)] = traversal_root;
  order.push_back(traversal_root);
  for (int i = 0; i < static_cast<int>(order.size()); ++i) {
    const int node = order[static_cast<std::size_t>(i)];
    for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
      if (nb < 0 || nb >= n || !tree.nodes[static_cast<std::size_t>(nb)].active ||
          parent[static_cast<std::size_t>(nb)] >= 0) {
        continue;
      }
      parent[static_cast<std::size_t>(nb)] = node;
      order.push_back(nb);
    }
  }
  std::vector<int> active_subtree_count(static_cast<std::size_t>(n), 0);
  std::vector<unsigned char> keep(static_cast<std::size_t>(n), 0);
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    const int node = *it;
    int count = active_leaf[static_cast<std::size_t>(node)] != 0 ? 1 : 0;
    for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
      if (nb >= 0 && nb < n && parent[static_cast<std::size_t>(nb)] == node) {
        count += active_subtree_count[static_cast<std::size_t>(nb)];
      }
    }
    active_subtree_count[static_cast<std::size_t>(node)] = count;
    if (count > 0 && count < static_cast<int>(active_leaf_node.size())) {
      keep[static_cast<std::size_t>(node)] = 1;
    }
  }
  keep[static_cast<std::size_t>(traversal_root)] = 1;
  for (int leaf : active_leaf_node) {
    int node = leaf;
    while (node >= 0 && keep[static_cast<std::size_t>(node)] == 0) {
      keep[static_cast<std::size_t>(node)] = 1;
      if (parent[static_cast<std::size_t>(node)] == node) break;
      node = parent[static_cast<std::size_t>(node)];
    }
  }

  QueryBuilder builder;
  std::vector<int> old_to_new(static_cast<std::size_t>(n), -1);
  for (int node = 0; node < n; ++node) {
    if (keep[static_cast<std::size_t>(node)] == 0) continue;
    const int taxon = tree.nodes[static_cast<std::size_t>(node)].taxon;
    const int active_id = taxon >= 0 ? active_by_taxon.find(taxon) : -1;
    const int kept_degree = static_cast<int>(std::count_if(
        tree.nodes[static_cast<std::size_t>(node)].neighbors.begin(),
        tree.nodes[static_cast<std::size_t>(node)].neighbors.end(), [&](int nb) {
          return nb >= 0 && nb < n && keep[static_cast<std::size_t>(nb)] != 0;
        }));
    if (active_id < 0 && kept_degree == 2) continue;
    old_to_new[static_cast<std::size_t>(node)] =
        builder.add_node(active_id, node);
  }
  auto endpoint = [&](int node, int from) {
    int prev = from;
    int cur = node;
    while (cur >= 0 && old_to_new[static_cast<std::size_t>(cur)] < 0) {
      int next = -1;
      for (int nb : tree.nodes[static_cast<std::size_t>(cur)].neighbors) {
        if (nb != prev && nb >= 0 && nb < n &&
            keep[static_cast<std::size_t>(nb)] != 0) {
          next = nb;
          break;
        }
      }
      prev = cur;
      cur = next;
    }
    return cur;
  };
  for (int node = 0; node < n; ++node) {
    if (old_to_new[static_cast<std::size_t>(node)] < 0) continue;
    for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
      if (nb < 0 || nb >= n || keep[static_cast<std::size_t>(nb)] == 0) continue;
      const int end = endpoint(nb, node);
      if (end >= 0 && old_to_new[static_cast<std::size_t>(end)] >= 0) {
        builder.add_edge(old_to_new[static_cast<std::size_t>(node)],
                         old_to_new[static_cast<std::size_t>(end)]);
      }
    }
  }
  int root = -1;
  for (int node = 0; node < static_cast<int>(builder.nodes.size()); ++node) {
    if (builder.nodes[static_cast<std::size_t>(node)].taxon_active < 0) {
      root = node;
      break;
    }
  }
  if (root < 0 && !builder.nodes.empty()) root = 0;
  finalize_query_tree(std::move(builder), root,
                      static_cast<int>(active_taxa.size()), out);
  return out;
}


// The query tree of a binary solver tree on all of `active_taxa`, rooted at
// its internal vertex `vertex`. The root is added after every other active
// node, so its children are the subtrees of the vertex's neighbours in node
// order; the SPR refinement colours them as the vertex's three arms.
QueryTree query_tree_rooted_at_vertex(const SolverTree& tree,
                                      const std::vector<int>& active_taxa,
                                      int taxon_count, int vertex) {
  QueryTree out;
  const int n = static_cast<int>(tree.nodes.size());
  if (n == 0 || active_taxa.empty()) return out;
  if (vertex < 0 || vertex >= n ||
      !tree.nodes[static_cast<std::size_t>(vertex)].active) {
    throw std::invalid_argument("query tree root is not an active vertex");
  }
  const TaxonIndex active_by_taxon =
      make_taxon_index(active_taxa, taxon_count);
  QueryBuilder builder;
  std::vector<int> map(static_cast<std::size_t>(n), -1);
  for (int node = 0; node < n; ++node) {
    if (!tree.nodes[static_cast<std::size_t>(node)].active || node == vertex) {
      continue;
    }
    const int taxon = tree.nodes[static_cast<std::size_t>(node)].taxon;
    const int active_id = taxon >= 0 ? active_by_taxon.find(taxon) : -1;
    map[static_cast<std::size_t>(node)] = builder.add_node(active_id, node);
  }
  const int root = builder.add_node(-1, vertex);
  for (int node = 0; node < n; ++node) {
    if (map[static_cast<std::size_t>(node)] < 0) continue;
    for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
      if (nb < 0 || nb >= n || !tree.nodes[static_cast<std::size_t>(nb)].active) {
        continue;
      }
      if (nb == vertex) {
        builder.add_edge(root, map[static_cast<std::size_t>(node)]);
      } else if (node < nb) {
        builder.add_edge(map[static_cast<std::size_t>(node)],
                         map[static_cast<std::size_t>(nb)]);
      }
    }
  }
  finalize_query_tree(std::move(builder), root,
                      static_cast<int>(active_taxa.size()), out);
  return out;
}

}  // namespace qraft::internal
