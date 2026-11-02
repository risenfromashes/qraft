#pragma once

#include <cstddef>
#include <limits>
#include <stdexcept>

namespace qraft::internal {

// Overflow-checked size and index arithmetic. Each helper computes the exact
// result and throws std::length_error(message) when it does not fit the
// return type, so size computations fail instead of wrapping.

inline std::size_t checked_hdt_size_add(std::size_t left,
                                       std::size_t right,
                                       const char* message) {
  if (right > std::numeric_limits<std::size_t>::max() - left) {
    throw std::length_error(message);
  }
  return left + right;
}

inline std::size_t checked_hdt_size_mul(std::size_t left,
                                       std::size_t right,
                                       const char* message) {
  if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
    throw std::length_error(message);
  }
  return left * right;
}

// total += count * element_size.
inline void checked_hdt_accumulate_bytes(std::size_t& total,
                                         std::size_t count,
                                         std::size_t element_size,
                                         const char* message) {
  total = checked_hdt_size_add(
      total, checked_hdt_size_mul(count, element_size, message), message);
}

inline int checked_hdt_int_size(std::size_t value, const char* message) {
  if (value > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error(message);
  }
  return static_cast<int>(value);
}

inline long long checked_hdt_long_long_size(std::size_t value,
                                            const char* message) {
  if (value >
      static_cast<std::size_t>(std::numeric_limits<long long>::max())) {
    throw std::length_error(message);
  }
  return static_cast<long long>(value);
}

inline int checked_hdt_int_add(int left, int right, const char* message) {
  const long long value =
      static_cast<long long>(left) + static_cast<long long>(right);
  if (value < std::numeric_limits<int>::min() ||
      value > std::numeric_limits<int>::max()) {
    throw std::length_error(message);
  }
  return static_cast<int>(value);
}

inline int checked_hdt_int_mul(int left, int right, const char* message) {
  const long long value =
      static_cast<long long>(left) * static_cast<long long>(right);
  if (value < std::numeric_limits<int>::min() ||
      value > std::numeric_limits<int>::max()) {
    throw std::length_error(message);
  }
  return static_cast<int>(value);
}

// value * multiplier + addend; only the final result must fit in int.
inline int checked_hdt_int_mul_add(int value,
                                  int multiplier,
                                  int addend,
                                  const char* message) {
  const long long result = static_cast<long long>(value) *
                               static_cast<long long>(multiplier) +
                           static_cast<long long>(addend);
  if (result < std::numeric_limits<int>::min() ||
      result > std::numeric_limits<int>::max()) {
    throw std::length_error(message);
  }
  return static_cast<int>(result);
}

}  // namespace qraft::internal
