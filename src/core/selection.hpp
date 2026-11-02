#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>

namespace qraft::internal {

// Score classes in maximization order: +Inf > finite > -Inf > NaN. The
// bipartition, regraft and NNI selections all rank scores this way.
enum class ScoreClass : unsigned char {
  NotANumber = 0,
  NegativeInfinity = 1,
  Finite = 2,
  PositiveInfinity = 3,
};

inline ScoreClass score_class(double score) noexcept {
  if (std::isnan(score)) return ScoreClass::NotANumber;
  if (std::isinf(score)) {
    return std::signbit(score) ? ScoreClass::NegativeInfinity
                               : ScoreClass::PositiveInfinity;
  }
  return ScoreClass::Finite;
}

// Strict weak order on scores with no epsilon: by class, then by value within
// the finite class. Two non-finite scores of one class are equivalent.
inline bool exact_score_less(double lhs, double rhs) noexcept {
  const ScoreClass lhs_class = score_class(lhs);
  const ScoreClass rhs_class = score_class(rhs);
  if (lhs_class != rhs_class) return lhs_class < rhs_class;
  if (lhs_class == ScoreClass::Finite) return lhs < rhs;
  return false;
}

// True when score has the class of exact_maximum and, for a finite maximum,
// lies within epsilon of it.
inline bool score_in_maximum_band(double score,
                                      double exact_maximum,
                                      double epsilon) noexcept {
  const ScoreClass value_class = score_class(score);
  const ScoreClass maximum_class = score_class(exact_maximum);
  if (value_class != maximum_class) return false;
  if (maximum_class != ScoreClass::Finite) return true;
  return std::fabs(score - exact_maximum) <= epsilon;
}

// exact_maximum is the first eligible candidate with the exact maximum score,
// selected the candidate with the smallest stable key in its band, and
// band_size the number of eligible candidates in the band. Both iterators are
// `last` when no candidate is eligible.
template <typename ForwardIt>
struct TieBandSelection {
  ForwardIt exact_maximum;
  ForwardIt selected;
  std::size_t band_size = 0;
};

// Selects in two independent passes. The first pass finds the exact numerical
// maximum. The second applies the stable key only within that maximum's epsilon
// band. stable_key_of must include a deterministic final fallback when distinct
// candidates can otherwise have the same key. An epsilon comparison is not a
// strict weak order, so it must never drive a sort.
template <typename ForwardIt,
          typename ScoreOf,
          typename StableKeyOf,
          typename Eligible>
TieBandSelection<ForwardIt> select_max_with_tie_band(
    ForwardIt first,
    ForwardIt last,
    double epsilon,
    ScoreOf score_of,
    StableKeyOf stable_key_of,
    Eligible eligible) {
  if (!std::isfinite(epsilon) || epsilon < 0.0) {
    throw std::invalid_argument(
        "score-selection epsilon must be finite and non-negative");
  }

  TieBandSelection<ForwardIt> result{last, last, 0};
  double exact_maximum = 0.0;
  for (ForwardIt candidate = first; candidate != last; ++candidate) {
    if (!std::invoke(eligible, *candidate)) continue;
    const double score = std::invoke(score_of, *candidate);
    if (result.exact_maximum == last ||
        exact_score_less(exact_maximum, score)) {
      result.exact_maximum = candidate;
      exact_maximum = score;
    }
  }

  if (result.exact_maximum == last) return result;

  for (ForwardIt candidate = first; candidate != last; ++candidate) {
    if (!std::invoke(eligible, *candidate)) continue;
    const double score = std::invoke(score_of, *candidate);
    if (!score_in_maximum_band(score, exact_maximum, epsilon)) continue;
    ++result.band_size;
    if (result.selected == last ||
        std::invoke(stable_key_of, *candidate) <
            std::invoke(stable_key_of, *result.selected)) {
      result.selected = candidate;
    }
  }
  return result;
}

// Every candidate is eligible.
template <typename ForwardIt, typename ScoreOf, typename StableKeyOf>
TieBandSelection<ForwardIt> select_max_with_tie_band(
    ForwardIt first,
    ForwardIt last,
    double epsilon,
    ScoreOf score_of,
    StableKeyOf stable_key_of) {
  return select_max_with_tie_band(
      first, last, epsilon, score_of, stable_key_of,
      [](const auto&) noexcept { return true; });
}

struct NniMove;
struct SolverTree;

// Test hook: applies the NNI gain threshold and tie-band selection to
// caller-supplied moves.
NniMove select_nni_move_for_tests(
    const SolverTree& tree,
    const std::vector<NniMove>& candidates,
    double min_gain);

}  // namespace qraft::internal
