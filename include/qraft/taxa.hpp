#pragma once

#include <string>
#include <vector>

namespace qraft {

/// Taxon labels and their ids 0..size()-1. Lookups by label take O(log n).
class TaxonRegistry {
public:
  TaxonRegistry() = default;
  /// Keeps the distinct labels, numbered in lexicographic label order.
  explicit TaxonRegistry(std::vector<std::string> labels);
  /// Numbers `labels` in the given order. Throws std::invalid_argument when a
  /// label repeats.
  static TaxonRegistry from_ordered_unique(std::vector<std::string> labels);

  /// Throws std::invalid_argument for an unknown label.
  int id_for(const std::string &label) const;
  /// Throws std::out_of_range for an id outside 0..size()-1.
  const std::string &label_for(int id) const;
  int size() const { return static_cast<int>(labels_.size()); }
  /// Labels indexed by id.
  const std::vector<std::string> &labels() const { return labels_; }

private:
  /// Sorts id_by_label_order_ and rejects duplicate labels.
  void rebuild_label_order();

  std::vector<std::string> labels_;
  std::vector<int> id_by_label_order_;  ///< ids sorted by label, for id_for
};

} // namespace qraft
