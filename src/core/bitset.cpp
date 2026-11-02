#include "internal.hpp"

#include <bit>

namespace qraft::internal {

Bitset::Bitset(int size) : size_(size), words_((size + 63) / 64, 0) {}

bool Bitset::empty() const {
  for (uint64_t word : words_) {
    if (word != 0) {
      return false;
    }
  }
  return true;
}

void Bitset::set(int pos) {
  if (pos < 0 || pos >= size_) {
    throw std::out_of_range("bitset set out of range");
  }
  words_[pos / 64] |= (uint64_t{1} << (pos % 64));
}

void Bitset::reset(int pos) {
  if (pos < 0 || pos >= size_) {
    throw std::out_of_range("bitset reset out of range");
  }
  words_[pos / 64] &= ~(uint64_t{1} << (pos % 64));
}

bool Bitset::test(int pos) const {
  if (pos < 0 || pos >= size_) {
    return false;
  }
  return (words_[pos / 64] & (uint64_t{1} << (pos % 64))) != 0;
}

int Bitset::count() const {
  int total = 0;
  for (uint64_t word : words_) {
    total += static_cast<int>(std::popcount(word));
  }
  return total;
}

void Bitset::clear() {
  std::fill(words_.begin(), words_.end(), uint64_t{0});
}

void Bitset::fill() {
  std::fill(words_.begin(), words_.end(), ~uint64_t{0});
  mask_unused_bits();
}

/// Padding bits are always zero, so comparing words compares the sets.
bool Bitset::operator==(const Bitset& other) const {
  return size_ == other.size_ && words_ == other.words_;
}

/// Clears the padding bits of the last word, which count() and operator==
/// rely on being zero.
void Bitset::mask_unused_bits() {
  if (size_ == 0 || words_.empty()) {
    return;
  }
  const int used = size_ % 64;
  if (used != 0) {
    words_.back() &= ((uint64_t{1} << used) - 1);
  }
}

}  // namespace qraft::internal
