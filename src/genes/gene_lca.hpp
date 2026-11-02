#pragma once

#include <cstddef>
#include <utility>
#include <vector>

namespace qraft::internal {

struct GeneIndex;

// Per-gene ancestry index rooted at the gene's reference root. tin/tout give
// O(1) subtree tests and the child-range keys; LCA queries climb heavy paths,
// O(log n) each and O(n) words in total. An LCA is unique, so results do not
// depend on the choice of LCA structure.
struct StaticLcaIndex {
  std::vector<int> tin;
  std::vector<int> tout;
  std::vector<int> depth;
  std::vector<int> parent;
  std::vector<int> head;
  // The indexed gene's own taxon map; entries past its end read as absent.
  const std::vector<int>* taxon_to_node = nullptr;

  int node_for_taxon(int taxon) const {
    if (taxon_to_node == nullptr || taxon < 0 ||
        static_cast<std::size_t>(taxon) >= taxon_to_node->size()) {
      return -1;
    }
    return (*taxon_to_node)[static_cast<std::size_t>(taxon)];
  }

  bool in_subtree(int node, int root) const {
    if (node < 0 || root < 0) return false;
    return tin[static_cast<std::size_t>(root)] <=
               tin[static_cast<std::size_t>(node)] &&
           tin[static_cast<std::size_t>(node)] <
               tout[static_cast<std::size_t>(root)];
  }

  int lca(int a, int b) const {
    if (a < 0) return b;
    if (b < 0) return a;
    while (head[static_cast<std::size_t>(a)] !=
           head[static_cast<std::size_t>(b)]) {
      if (depth[static_cast<std::size_t>(head[static_cast<std::size_t>(a)])] <
          depth[static_cast<std::size_t>(head[static_cast<std::size_t>(b)])]) {
        std::swap(a, b);
      }
      a = parent[static_cast<std::size_t>(head[static_cast<std::size_t>(a)])];
    }
    return depth[static_cast<std::size_t>(a)] <=
                   depth[static_cast<std::size_t>(b)]
               ? a
               : b;
  }
};

// Builds the index of `gene` rooted at its reference root; taxon_count bounds
// the taxa the gene may hold.
StaticLcaIndex build_static_lca_index(const GeneIndex& gene, int taxon_count);

}  // namespace qraft::internal
