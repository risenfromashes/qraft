#include "internal.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <mutex>
#include <optional>
#include <sstream>

namespace qraft::internal {
namespace {

// Whether stats_add, stats_max and stats_add_time record; see
// ScopedInferenceRun.
std::atomic<bool>& stats_gate() {
  static std::atomic<bool> enabled{true};
  return enabled;
}

}  // namespace

std::vector<Side> assignment_from_split(const Bitset& split, int m) {
  std::vector<Side> assignment(m, Side::Right);
  for (int i = 0; i < m; ++i) {
    if (split.test(i)) assignment[i] = Side::Left;
  }
  return assignment;
}

const std::vector<GeneIndex>& scoring_genes(const SolverContext& ctx) {
  if (ctx.owned_active_genes != nullptr) {
    return *ctx.owned_active_genes;
  }
  return ctx.active_genes == nullptr ? base_genes(ctx) : *ctx.active_genes;
}

const std::vector<GeneIndex>& base_genes(const SolverContext& ctx) {
  return ctx.genes_ref == nullptr ? ctx.genes : *ctx.genes_ref;
}

const TaxonRegistry& context_taxa(const SolverContext& ctx) {
  return ctx.taxa_ref == nullptr ? ctx.taxa : *ctx.taxa_ref;
}

const StartingTree& context_starting_tree(const SolverContext& ctx) {
  return ctx.starting_tree_ref == nullptr ? ctx.starting_tree : *ctx.starting_tree_ref;
}

RunStats& run_stats() {
  static RunStats stats;
  return stats;
}

bool run_stats_enabled() noexcept {
  return stats_gate().load(std::memory_order_acquire);
}

void set_run_stats_enabled(bool enabled) noexcept {
  stats_gate().store(enabled, std::memory_order_release);
}

void reset_run_stats() {
  RunStats& stats = run_stats();
#define QRAFT_RESET_STAT(name) stats.name.store(0, std::memory_order_relaxed);
  QRAFT_RUN_STATS(QRAFT_RESET_STAT)
#undef QRAFT_RESET_STAT
  stats.complexity_contract_satisfied.store(1, std::memory_order_relaxed);
  stats.nni_gain = WideCount{};
  stats.spr_gain = 0;
  stats.spr_gain_exact = true;
}

void wide_add(WideCount& sum, double value) {
  constexpr double kTwo64 = 18446744073709551616.0;
  const double rounded = std::nearbyint(value);
  if (!(rounded >= 0.0) || rounded >= kTwo64 * kTwo64) {
    throw std::overflow_error("wide count addend is negative, NaN or too large");
  }
  // Both parts are exact: `rounded` is an integer and kTwo64 a power of two.
  const double high_part = std::floor(rounded / kTwo64);
  const std::uint64_t add_high = static_cast<std::uint64_t>(high_part);
  const std::uint64_t add_low =
      static_cast<std::uint64_t>(rounded - high_part * kTwo64);
  const std::uint64_t low = sum.low + add_low;
  const std::uint64_t carry = low < add_low ? 1 : 0;
  constexpr std::uint64_t kMax = std::numeric_limits<std::uint64_t>::max();
  if (add_high > kMax - carry || sum.high > kMax - add_high - carry) {
    throw std::overflow_error("wide count sum reaches 2^128");
  }
  sum.low = low;
  sum.high += add_high + carry;
}

std::string wide_decimal(const WideCount& value) {
  // Divides four 32-bit limbs, most significant first, by 10 per digit.
  std::array<std::uint32_t, 4> limbs{
      static_cast<std::uint32_t>(value.high >> 32),
      static_cast<std::uint32_t>(value.high),
      static_cast<std::uint32_t>(value.low >> 32),
      static_cast<std::uint32_t>(value.low)};
  std::string digits;
  bool nonzero = true;
  while (nonzero) {
    std::uint64_t remainder = 0;
    nonzero = false;
    for (std::uint32_t& limb : limbs) {
      const std::uint64_t current = (remainder << 32) | limb;
      limb = static_cast<std::uint32_t>(current / 10);
      remainder = current % 10;
      nonzero = nonzero || limb != 0;
    }
    digits.push_back(static_cast<char>('0' + remainder));
  }
  std::reverse(digits.begin(), digits.end());
  return digits;
}

void stats_add(std::atomic<long long>& counter, long long value) {
  if (!run_stats_enabled()) return;
  counter.fetch_add(value, std::memory_order_relaxed);
}

void stats_max(std::atomic<long long>& counter, long long value) {
  if (!run_stats_enabled()) return;
  long long current = counter.load(std::memory_order_relaxed);
  while (current < value &&
         !counter.compare_exchange_weak(current, value,
                                        std::memory_order_relaxed,
                                        std::memory_order_relaxed)) {
  }
}

namespace {

// Counts a release of more than the category holds.
void note_release_underflow(long long previous, long long value) {
  if (value < 0 || previous < value) {
    stats_add(
        run_stats().accounting_violations);
  }
}

}  // namespace

void stats_acquire_payload(std::atomic<long long>& current,
                                    std::atomic<long long>& peak,
                                    long long value) {
  if (value <= 0) return;
  const long long category =
      current.fetch_add(value, std::memory_order_relaxed) + value;
  stats_max(peak, category);
  RunStats& stats = run_stats();
  const long long total = stats.live_payload.fetch_add(
                              value, std::memory_order_relaxed) +
                          value;
  stats_max(stats.peak_payload, total);
}

void stats_release_payload(std::atomic<long long>& current,
                                    long long value) {
  if (value <= 0) return;
  const long long previous =
      current.fetch_sub(value, std::memory_order_relaxed);
  note_release_underflow(previous, value);
  RunStats& stats = run_stats();
  const long long total_previous = stats.live_payload.fetch_sub(
      value, std::memory_order_relaxed);
  note_release_underflow(total_previous, value);
}

void stats_acquire_gauge(std::atomic<long long>& current,
                                  std::atomic<long long>& peak,
                                  long long value) {
  if (value <= 0) return;
  const long long observed =
      current.fetch_add(value, std::memory_order_relaxed) + value;
  stats_max(peak, observed);
}

void stats_release_gauge(std::atomic<long long>& current,
                                  long long value) {
  if (value <= 0) return;
  const long long previous =
      current.fetch_sub(value, std::memory_order_relaxed);
  note_release_underflow(previous, value);
}

void stats_add_time(std::atomic<long long>& counter,
                         std::chrono::nanoseconds elapsed) {
  if (!run_stats_enabled()) return;
  counter.fetch_add(elapsed.count(), std::memory_order_relaxed);
}

void stats_set_nni_gain(const WideCount& gain) {
  if (!run_stats_enabled()) return;
  run_stats().nni_gain = gain;
}

const char* progress_mode_name(ProgressMode mode) {
  switch (mode) {
    case ProgressMode::Off:
      return "off";
    case ProgressMode::Auto:
      return "auto";
    case ProgressMode::On:
      return "on";
  }
  return "unknown";
}

}  // namespace qraft::internal
