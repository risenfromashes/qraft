#include "internal.hpp"

#include <cctype>
#include <sstream>

namespace qraft {

/// Every non-empty node label, in node order.
std::vector<std::string> Tree::leaf_labels() const {
  std::vector<std::string> labels;
  for (const auto &node : nodes) {
    if (!node.label.empty()) {
      labels.push_back(node.label);
    }
  }
  return labels;
}

namespace {

/// Quotes a label that is empty or contains whitespace or Newick punctuation,
/// doubling embedded single quotes.
std::string escape_newick_label(const std::string &label) {
  bool needs_quote = label.empty();
  for (char ch : label) {
    if (std::isspace(static_cast<unsigned char>(ch)) || ch == '(' ||
        ch == ')' || ch == ',' || ch == ':' || ch == ';' || ch == '\'' ||
        ch == '[' || ch == ']') {
      needs_quote = true;
      break;
    }
  }
  if (!needs_quote) {
    return label;
  }
  std::string out = "'";
  for (char ch : label) {
    if (ch == '\'') {
      out += "''";
    } else {
      out += ch;
    }
  }
  out += "'";
  return out;
}

/// Writes the subtree at root away from parent. Child strings are sorted, so
/// the output does not depend on neighbor order. The nodes whose children
/// are still being written are kept on an explicit stack, so the depth of
/// the tree costs no call stack.
std::string write_subtree(const SpeciesTree &tree, int root, int parent) {
  // The text of edge {u, v}, or nothing.
  const auto edge_text = [&](int u, int v) -> std::string {
    if (tree.edge_labels.empty() || v < 0) return {};
    const int a = std::min(u, v), b = std::max(u, v);
    const auto it = std::lower_bound(
        tree.edge_labels.begin(), tree.edge_labels.end(), std::pair<int, int>{a, b},
        [](const SpeciesTree::EdgeLabel &e, const std::pair<int, int> &key) {
          return std::pair<int, int>{e.a, e.b} < key;
        });
    if (it == tree.edge_labels.end() || it->a != a || it->b != b) return {};
    return it->text;
  };
  struct Open {
    int node;
    std::vector<int> children;
    std::size_t next = 0;
    std::vector<std::string> parts;
  };
  std::vector<Open> open;
  std::string text;
  // Writes a labelled leaf into `text` and returns false, or opens any other
  // node and returns true.
  const auto start = [&](int node, int from) {
    const auto &current = tree.nodes[node];
    std::vector<int> children;
    for (int neighbor : current.neighbors) {
      if (neighbor != from) {
        children.push_back(neighbor);
      }
    }
    if (!current.label.empty() && children.empty()) {
      text = escape_newick_label(current.label) + edge_text(node, from);
      return false;
    }
    Open entry{node, std::move(children), 0, {}};
    entry.parts.reserve(entry.children.size());
    open.push_back(std::move(entry));
    return true;
  };

  if (!start(root, parent)) {
    return text;
  }
  while (true) {
    Open &top = open.back();
    if (top.next < top.children.size()) {
      const int child = top.children[top.next++];
      if (!start(child, top.node)) {
        open.back().parts.push_back(std::move(text));
      }
      continue;
    }
    std::sort(top.parts.begin(), top.parts.end());
    std::ostringstream out;
    out << "(";
    for (std::size_t i = 0; i < top.parts.size(); ++i) {
      if (i != 0) {
        out << ",";
      }
      out << top.parts[i];
    }
    out << ")";
    const std::string &label = tree.nodes[top.node].label;
    if (!label.empty()) {
      out << escape_newick_label(label);
    }
    out << edge_text(top.node, open.size() > 1 ? open[open.size() - 2].node : parent);
    open.pop_back();
    if (open.empty()) {
      return out.str();
    }
    open.back().parts.push_back(out.str());
  }
}

} // namespace

/// Writes the tree rooted at its first node that has several neighbors or no
/// label, with the children of every node in sorted order.
std::string to_newick(const SpeciesTree &tree) {
  if (tree.nodes.empty()) {
    return ";";
  }
  if (tree.nodes.size() == 1) {
    return escape_newick_label(tree.nodes[0].label) + ";";
  }
  if (tree.nodes.size() == 2 && !tree.nodes[0].label.empty() &&
      !tree.nodes[1].label.empty()) {
    std::vector<std::string> labels = {
        escape_newick_label(tree.nodes[0].label),
        escape_newick_label(tree.nodes[1].label)};
    std::sort(labels.begin(), labels.end());
    return "(" + labels[0] + "," + labels[1] + ");";
  }
  int root = -1;
  for (int i = 0; i < static_cast<int>(tree.nodes.size()); ++i) {
    if (tree.nodes[i].neighbors.size() > 1 || tree.nodes[i].label.empty()) {
      root = i;
      break;
    }
  }
  if (root == -1) {
    root = 0;
  }
  return write_subtree(tree, root, -1) + ";";
}

} // namespace qraft

namespace qraft::internal {

/// FNV-style mix of three integers. Fixed-width arithmetic makes it identical
/// on every platform and run.
uint64_t stable_hash_u64(uint64_t a, uint64_t b, uint64_t c) {
  uint64_t h = 1469598103934665603ULL;
  auto mix = [&](uint64_t x) {
    h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h *= 1099511628211ULL;
  };
  mix(a);
  mix(b);
  mix(c);
  return h;
}

} // namespace qraft::internal
