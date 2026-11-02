#include "genes/gene_lca.hpp"

#include "internal.hpp"

#include <limits>
#include <stdexcept>
#include <tuple>

namespace qraft::internal {

StaticLcaIndex build_static_lca_index(const GeneIndex& gene, int taxon_count) {
  StaticLcaIndex out;
  if (taxon_count < 0 ||
      gene.nodes.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      gene.taxon_to_node.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("gene LCA index exceeds its integer domain");
  }
  const int node_count = static_cast<int>(gene.nodes.size());
  out.tin.assign(static_cast<std::size_t>(node_count), 0);
  out.tout.assign(static_cast<std::size_t>(node_count), 0);
  out.depth.assign(static_cast<std::size_t>(node_count), 0);
  out.parent.assign(static_cast<std::size_t>(node_count), -1);
  out.head.assign(static_cast<std::size_t>(node_count), -1);
  out.taxon_to_node = &gene.taxon_to_node;
  if (node_count == 0) return out;
  if (gene.reference_root < 0 || gene.reference_root >= node_count) {
    throw std::invalid_argument("gene has an invalid reference root");
  }

  std::vector<int> order;
  order.reserve(static_cast<std::size_t>(node_count));
  int timer = 0;
  std::vector<std::tuple<int, int, int>> stack;
  stack.push_back({gene.reference_root, gene.reference_root, 0});
  while (!stack.empty()) {
    auto& [node, parent, next] = stack.back();
    if (next == 0) {
      out.tin[static_cast<std::size_t>(node)] = timer++;
      order.push_back(node);
    }
    if (next < gene_degree(gene, node)) {
      const int nb =
          gene_neighbors(gene, node)[static_cast<std::size_t>(next++)];
      if (nb == parent) continue;
      out.depth[static_cast<std::size_t>(nb)] =
          out.depth[static_cast<std::size_t>(node)] + 1;
      out.parent[static_cast<std::size_t>(nb)] = node;
      stack.push_back({nb, node, 0});
      continue;
    }
    out.tout[static_cast<std::size_t>(node)] = timer;
    stack.pop_back();
  }

  std::vector<int> size(static_cast<std::size_t>(node_count), 1);
  std::vector<int> heavy(static_cast<std::size_t>(node_count), -1);
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    const int node = *it;
    const int up = out.parent[static_cast<std::size_t>(node)];
    if (up < 0) continue;
    size[static_cast<std::size_t>(up)] += size[static_cast<std::size_t>(node)];
    const int current = heavy[static_cast<std::size_t>(up)];
    if (current < 0 || size[static_cast<std::size_t>(node)] >=
                           size[static_cast<std::size_t>(current)]) {
      heavy[static_cast<std::size_t>(up)] = node;
    }
  }
  for (int node : order) {
    const int up = out.parent[static_cast<std::size_t>(node)];
    out.head[static_cast<std::size_t>(node)] =
        up >= 0 && heavy[static_cast<std::size_t>(up)] == node
            ? out.head[static_cast<std::size_t>(up)]
            : node;
  }
  return out;
}

}  // namespace qraft::internal
