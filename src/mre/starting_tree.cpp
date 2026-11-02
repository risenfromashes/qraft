#include "internal.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace qraft::internal {
namespace {

/// Fills the preorder and, for each node, its preorder interval [tin, tout).
void fill_starting_tree_traversal(StartingTree& starting_tree) {
  starting_tree.preorder.clear();
  starting_tree.tin.assign(starting_tree.nodes.size(), 0);
  starting_tree.tout.assign(starting_tree.nodes.size(), 0);
  if (starting_tree.root < 0) return;
  std::vector<std::pair<int, int>> stack{{starting_tree.root, 0}};
  while (!stack.empty()) {
    const int node = stack.back().first;
    int& next = stack.back().second;
    if (next == 0) {
      starting_tree.tin[static_cast<std::size_t>(node)] =
          static_cast<int>(starting_tree.preorder.size());
      starting_tree.preorder.push_back(node);
    }
    if (next < static_cast<int>(
                   starting_tree.nodes[static_cast<std::size_t>(node)]
                       .children.size())) {
      stack.push_back(
          {starting_tree.nodes[static_cast<std::size_t>(node)]
               .children[static_cast<std::size_t>(next++)],
           0});
      continue;
    }
    starting_tree.tout[static_cast<std::size_t>(node)] =
        static_cast<int>(starting_tree.preorder.size());
    stack.pop_back();
  }
}

}  // namespace

/// Converts a user starting tree (--starting-tree) into a StartingTree rooted
/// at its parse root. Every registry taxon must label exactly one leaf; an
/// unknown label throws from id_for, and labels of internal nodes are ignored.
StartingTree starting_tree_from(const Tree& tree, const TaxonRegistry& taxa) {
  if (tree.nodes.empty()) {
    throw std::invalid_argument("starting tree is empty");
  }

  StartingTree out;
  out.nodes.resize(tree.nodes.size());
  out.root = tree.root >= 0 ? tree.root : 0;
  out.taxon_to_leaf.assign(taxa.size(), -1);

  std::vector<int> parent(tree.nodes.size(), -1);
  std::vector<int> stack{out.root};
  parent[static_cast<std::size_t>(out.root)] = out.root;
  while (!stack.empty()) {
    const int node = stack.back();
    stack.pop_back();
    out.nodes[static_cast<std::size_t>(node)].parent =
        node == out.root ? -1 : parent[static_cast<std::size_t>(node)];

    int child_count = 0;
    for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
      if (nb == parent[static_cast<std::size_t>(node)]) continue;
      if (parent[static_cast<std::size_t>(nb)] != -1) {
        throw std::invalid_argument("starting tree is not acyclic");
      }
      parent[static_cast<std::size_t>(nb)] = node;
      out.nodes[static_cast<std::size_t>(node)].children.push_back(nb);
      stack.push_back(nb);
      ++child_count;
    }

    if (!tree.nodes[static_cast<std::size_t>(node)].label.empty() &&
        child_count == 0) {
      const int taxon =
          taxa.id_for(tree.nodes[static_cast<std::size_t>(node)].label);
      if (out.taxon_to_leaf[static_cast<std::size_t>(taxon)] != -1) {
        throw std::invalid_argument(
            "duplicate taxon in starting tree: " +
            tree.nodes[static_cast<std::size_t>(node)].label);
      }
      out.nodes[static_cast<std::size_t>(node)].taxon = taxon;
      out.taxon_to_leaf[static_cast<std::size_t>(taxon)] = node;
    }
  }
  for (int i = 0; i < static_cast<int>(parent.size()); ++i) {
    if (parent[static_cast<std::size_t>(i)] == -1) {
      throw std::invalid_argument("starting tree is disconnected");
    }
  }
  for (int taxon = 0; taxon < taxa.size(); ++taxon) {
    if (out.taxon_to_leaf[static_cast<std::size_t>(taxon)] < 0) {
      throw std::invalid_argument("starting tree is missing taxon: " +
                                  taxa.label_for(taxon));
    }
  }

  fill_starting_tree_traversal(out);
  return out;
}

}  // namespace qraft::internal
