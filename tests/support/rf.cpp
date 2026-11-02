#include "support/rf.hpp"

#include <algorithm>
#include <set>
#include <string>
#include <vector>

namespace qraft {
namespace {

using LabelSet = std::set<std::string>;
using Split = std::vector<std::string>;
using SplitSet = std::set<Split>;

/// Labels of the leaves of `tree`: labelled nodes of degree at most one.
LabelSet leaf_label_set(const Tree& tree) {
  LabelSet labels;
  for (const Tree::Node& node : tree.nodes) {
    if (!node.label.empty() && node.neighbors.size() <= 1) {
      labels.insert(node.label);
    }
  }
  return labels;
}

LabelSet set_intersection(LabelSet lhs, const LabelSet& rhs) {
  for (auto it = lhs.begin(); it != lhs.end();) {
    if (rhs.find(*it) == rhs.end()) {
      it = lhs.erase(it);
    } else {
      ++it;
    }
  }
  return lhs;
}

/// Canonical form of the split of `keep` with one side `side`: the smaller
/// side, sorted, or the lexicographically smaller side when both are equal.
Split canonical_split(Split side, const LabelSet& keep) {
  std::sort(side.begin(), side.end());
  Split other;
  other.reserve(keep.size() - side.size());
  for (const std::string& label : keep) {
    if (!std::binary_search(side.begin(), side.end(), label)) {
      other.push_back(label);
    }
  }
  if (side.size() > other.size() ||
      (side.size() == other.size() && other < side)) {
    return other;
  }
  return side;
}

/// Leaf labels in `keep` reachable from `start` without passing `blocked`.
Split edge_side(const Tree& tree, int start, int blocked, const LabelSet& keep) {
  Split side;
  if (start < 0 || start >= static_cast<int>(tree.nodes.size())) {
    return side;
  }

  std::vector<char> seen(tree.nodes.size(), false);
  std::vector<int> stack{start};
  if (blocked >= 0 && blocked < static_cast<int>(seen.size())) {
    seen[blocked] = true;
  }

  while (!stack.empty()) {
    const int node = stack.back();
    stack.pop_back();
    if (node < 0 || node >= static_cast<int>(tree.nodes.size()) || seen[node]) {
      continue;
    }
    seen[node] = true;
    const Tree::Node& current = tree.nodes[node];
    if (!current.label.empty() && current.neighbors.size() <= 1 &&
        keep.find(current.label) != keep.end()) {
      side.push_back(current.label);
    }
    for (int next : current.neighbors) {
      if (next >= 0 && next < static_cast<int>(seen.size()) && !seen[next]) {
        stack.push_back(next);
      }
    }
  }
  return side;
}

/// Canonical non-trivial splits that the edges of `tree` induce on `keep`.
SplitSet unrooted_splits(const Tree& tree, const LabelSet& keep) {
  SplitSet splits;
  if (keep.size() < 4) {
    return splits;
  }

  for (int u = 0; u < static_cast<int>(tree.nodes.size()); ++u) {
    for (int v : tree.nodes[u].neighbors) {
      if (u > v) {
        continue;
      }
      Split side = edge_side(tree, u, v, keep);
      const std::size_t other_size = keep.size() - side.size();
      if (side.size() <= 1 || other_size <= 1) {
        continue;
      }
      splits.insert(canonical_split(std::move(side), keep));
    }
  }
  return splits;
}

}  // namespace

std::pair<int, double> rf_distance(const Tree& lhs, const Tree& rhs) {
  const LabelSet common = set_intersection(leaf_label_set(lhs), leaf_label_set(rhs));
  const SplitSet lhs_splits = unrooted_splits(lhs, common);
  const SplitSet rhs_splits = unrooted_splits(rhs, common);

  int missing = 0;
  for (const Split& split : rhs_splits) {
    if (lhs_splits.find(split) == lhs_splits.end()) {
      ++missing;
    }
  }

  int extra = 0;
  for (const Split& split : lhs_splits) {
    if (rhs_splits.find(split) == rhs_splits.end()) {
      ++extra;
    }
  }

  const int rf = missing + extra;
  const int denom =
      static_cast<int>(lhs_splits.size() + rhs_splits.size());
  return {rf, denom == 0 ? 0.0 : static_cast<double>(rf) / denom};
}

}  // namespace qraft
