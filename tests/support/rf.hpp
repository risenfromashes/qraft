#pragma once

#include <utility>

#include "qraft/tree.hpp"

namespace qraft {

/// Returns the unrooted Robinson-Foulds distance and normalized RF distance.
///
/// The comparison uses non-trivial splits induced on the taxa shared by both
/// trees. The normalized value is `rf / (lhs_splits + rhs_splits)`, as in
/// `tools/bench/rf_compare.py`, and is zero when neither tree has a split.
std::pair<int, double> rf_distance(const Tree& lhs, const Tree& rhs);

}  // namespace qraft
