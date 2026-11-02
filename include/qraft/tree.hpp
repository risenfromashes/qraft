#pragma once

#include <string>
#include <string_view>
#include <vector>

namespace qraft {

/// An unrooted tree stored as adjacency lists, used for both gene trees and
/// species trees.
struct Tree {
  struct Node {
    std::vector<int> neighbors;  ///< indices into `nodes`; edges at both ends
    std::string label;           ///< taxon of a leaf; empty on internal nodes
  };

  /// The node the Newick text is rooted at, or -1. infer_species_tree needs
  /// a valid root on every non-empty gene tree; quartets do not depend on
  /// where it is.
  int root = -1;
  std::vector<Node> nodes;
  /// Text written after the lower end of an edge, such as "0.98:1.25"
  /// (support and length), keyed by the edge's node indices, smaller first,
  /// and sorted by them. Empty unless support was requested.
  struct EdgeLabel {
    int a = -1;
    int b = -1;
    std::string text;
  };
  std::vector<EdgeLabel> edge_labels;
  /// Every non-empty node label, in node order.
  std::vector<std::string> leaf_labels() const;
};

using GeneTree = Tree;
using SpeciesTree = Tree;

/// Parses one Newick tree; the final ';' is optional. Branch lengths,
/// internal-node labels and bracketed comments are skipped. A quoted label
/// may contain any character, with '' standing for a quote. Throws
/// std::invalid_argument on malformed text or an unlabeled leaf.
Tree parse_newick(std::string_view text);
/// Writes a tree as Newick, with the edge_labels after the lower end of
/// their edges. Children are sorted by their text, so the output does not
/// depend on the order of adjacency lists. Labels with whitespace or Newick
/// punctuation are quoted.
std::string to_newick(const Tree &tree);

} // namespace qraft
