#include <algorithm>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <qraft/taxa.hpp>

namespace qraft {

/// Deduplicates labels. Taxon ids follow lexicographic label order, so they
/// do not depend on the order in which labels were supplied.
TaxonRegistry::TaxonRegistry(std::vector<std::string> labels) {
  // An ordered map over the distinct labels: each input label, duplicate or
  // not, costs O(log n) comparisons, with no hashing.
  std::map<std::string_view, std::size_t, std::less<>> first_by_label;
  std::vector<std::size_t> first_occurrence;
  for (std::size_t i = 0; i < labels.size(); ++i) {
    first_by_label.try_emplace(labels[i], i);
  }
  first_occurrence.reserve(first_by_label.size());
  for (const auto& [label, index] : first_by_label) {
    (void)label;
    first_occurrence.push_back(index);
  }
  first_by_label.clear();
  if (first_occurrence.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("taxon registry exceeds its integer id domain");
  }
  labels_.reserve(first_occurrence.size());
  for (std::size_t index : first_occurrence) {
    labels_.push_back(std::move(labels[index]));
  }
  rebuild_label_order();
}

/// Taxon ids follow the given order; a repeated label throws.
TaxonRegistry TaxonRegistry::from_ordered_unique(std::vector<std::string> labels) {
  if (labels.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("taxon registry exceeds its integer id domain");
  }
  TaxonRegistry registry;
  registry.labels_ = std::move(labels);
  registry.rebuild_label_order();
  return registry;
}

/// Sorts ids by label for the binary search in id_for and rejects duplicate
/// labels.
void TaxonRegistry::rebuild_label_order() {
  id_by_label_order_.resize(labels_.size());
  std::iota(id_by_label_order_.begin(), id_by_label_order_.end(), 0);
  std::sort(id_by_label_order_.begin(), id_by_label_order_.end(),
            [&](int left, int right) {
              return labels_[static_cast<std::size_t>(left)] <
                     labels_[static_cast<std::size_t>(right)];
            });
  for (std::size_t i = 1; i < id_by_label_order_.size(); ++i) {
    const int previous = id_by_label_order_[i - 1];
    const int current = id_by_label_order_[i];
    if (labels_[static_cast<std::size_t>(previous)] ==
        labels_[static_cast<std::size_t>(current)]) {
      throw std::invalid_argument(
          "duplicate taxon label: " +
          labels_[static_cast<std::size_t>(current)]);
    }
  }
}

/// O(log n); throws std::invalid_argument for an unknown label.
int TaxonRegistry::id_for(const std::string &label) const {
  const auto it = std::lower_bound(
      id_by_label_order_.begin(), id_by_label_order_.end(), label,
      [&](int id, const std::string& sought) {
        return labels_[static_cast<std::size_t>(id)] < sought;
      });
  if (it == id_by_label_order_.end() ||
      labels_[static_cast<std::size_t>(*it)] != label) {
    throw std::invalid_argument("unknown taxon label: " + label);
  }
  return *it;
}

/// Throws std::out_of_range for an invalid id.
const std::string &TaxonRegistry::label_for(int id) const {
  if (id < 0 || id >= static_cast<int>(labels_.size())) {
    throw std::out_of_range("taxon id out of range");
  }
  return labels_[id];
}

} // namespace qraft
