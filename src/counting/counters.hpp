#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <span>
#include <vector>

// Counters of HDT components for a fixed palette K = 2 or 4: the counter
// basis (unit, leaf, pair, triple and quartet counters), derived by
// tools/counting/generate_counters.py into generated/counters_k{2,4}.hpp.
// Every transition is a fixed linear or bilinear map, so it costs O(1) for a
// fixed K:
//   LEAF_TO_C  leaf -> C1 (U)
//   C_TO_G     C1 (U) -> G
//   GG_TO_G    G x G -> G
//   IG_FF      G -> Closed    joins a G to the vertex above it, giving a
//   IG_TF      G -> C1 (U)    component with no, one or two boundary
//   IG_TT      G -> C2        vertices
//   CLOSE      C1 (U) -> Closed
//   REVERSE    C1 U <-> D, and C2 -> C2
//   CC_TO_C    (C1 D or C2) x (C1 U or C2) -> Closed, C1 or C2
// plus the C2 identity element and the root readout of the channels:
//   K = 2: W(00|11), W(01|01)             (bipartition scoring)
//   K = 4: W(01|23), W(02|13), W(03|12)   (regrafting and NNI)
// Counters per interface (K=2 / K=4): Closed 2 / 3, C1 24 / 52 (both
// orientations), C2 23 / 56, G 15 / 32.
namespace qraft {

// Closed is the draft's R. C1 with orientation U is the draft's C1 and with
// orientation D its C3; C2 and G match. The largest interface has 56
// counters, so values are stored inline and every transition is
// allocation-free.
enum class CounterInterface : std::uint8_t {
  Closed,
  C1,
  C2,
  G,
};

enum class CounterOrientation : std::uint8_t {
  None,
  U,
  D,
};

inline constexpr std::size_t kMaxCounters = 56;

// Constant packed layout of one (palette, packed sort, colour mask): the
// stored full-coordinate bits, their count, and each coordinate's packed rank
// (kCounterNotPacked when the coordinate is not stored). The packed sorts are
// Closed, C1 with a U port, C1 with a D port only, C2 and G. The table is
// derived from the generated bit tables; hot paths read it through the
// unchecked accessor counters_packed_layout, and the validating functions
// below are the checked interface.
inline constexpr std::uint8_t kCounterNotPacked = 0xff;
struct PackedCounterLayout {
  std::uint64_t bits = 0;
  std::uint8_t width = 0;
  std::array<std::uint8_t, kMaxCounters> rank{};
};
using PackedCounterLayoutTable =
    std::array<std::array<std::array<PackedCounterLayout, 16>, 5>, 2>;
extern const PackedCounterLayoutTable kPackedCounterLayouts;

inline int counters_packed_sort(CounterInterface interface,
                                CounterOrientation orientation) noexcept {
  switch (interface) {
    case CounterInterface::Closed: return 0;
    case CounterInterface::C1:
      return orientation == CounterOrientation::D ? 2 : 1;
    case CounterInterface::C2: return 3;
    case CounterInterface::G: return 4;
  }
  return 0;
}

// Unchecked: palette must be 2 or 4, the mask inside the palette, and a C1
// orientation U or D. Callers pass values produced by the HDT itself.
inline const PackedCounterLayout& counters_packed_layout(
    int palette,
    CounterInterface interface,
    CounterOrientation orientation,
    std::uint8_t color_mask) noexcept {
  return kPackedCounterLayouts[palette == 4 ? 1 : 0]
                              [static_cast<std::size_t>(
                                  counters_packed_sort(interface, orientation))]
                              [color_mask & 15u];
}

// Full (unpacked) counters of one interface. Only the first `dimension`
// values are meaningful; the rest are uninitialized and never copied.
struct CounterVec {
  CounterInterface interface = CounterInterface::Closed;
  CounterOrientation orientation = CounterOrientation::None;
  std::uint8_t dimension = 0;
  std::array<double, kMaxCounters> values;

  CounterVec() = default;

  CounterVec(const CounterVec& other)
      : interface(other.interface),
        orientation(other.orientation),
        dimension(other.dimension) {
    std::copy_n(other.values.data(), dimension, values.data());
  }

  CounterVec& operator=(const CounterVec& other) {
    if (this == &other) return *this;
    interface = other.interface;
    orientation = other.orientation;
    dimension = other.dimension;
    std::copy_n(other.values.data(), dimension, values.data());
    return *this;
  }

  CounterVec(CounterVec&& other) noexcept
      : interface(other.interface),
        orientation(other.orientation),
        dimension(other.dimension) {
    std::copy_n(other.values.data(), dimension, values.data());
  }

  CounterVec& operator=(CounterVec&& other) noexcept {
    if (this == &other) return *this;
    interface = other.interface;
    orientation = other.orientation;
    dimension = other.dimension;
    std::copy_n(other.values.data(), dimension, values.data());
    return *this;
  }

  std::span<double> active() {
    return std::span<double>(values.data(), dimension);
  }
  std::span<const double> active() const {
    return std::span<const double>(values.data(), dimension);
  }
};

// Counters of one component under the current and the background colouring.
// When current_equals_background is set the two are equal, and the pair
// transitions below evaluate once.
struct CounterPair {
  CounterVec current;
  CounterVec background;
  bool current_equals_background = false;
};

// Non-owning view of stored counters: `dimension` is the interface's full
// dimension and `values` holds the counters packed under the owner's colour
// mask.
struct CounterView {
  CounterInterface interface = CounterInterface::Closed;
  CounterOrientation orientation = CounterOrientation::None;
  std::uint8_t dimension = 0;
  const double* values = nullptr;
};

/// Counters of `interface` for `palette`, or 0 for an unsupported palette.
int counters_dimension(int palette, CounterInterface interface);
int counters_max_dimension(int palette);
/// Generated scalar terms evaluated so far on the calling thread.
long long counters_thread_scalar_fmas();

CounterVec counters_leaf(int palette, int color, double weight);
CounterPair counters_leaf_dual(int palette,
                               int current_color,
                               int background_color,
                               double weight);
CounterVec counters_identity_c2(int palette);

CounterVec counters_c_to_g(int palette, const CounterVec& c1_u);
CounterVec counters_ig_ff(int palette, const CounterVec& g);
CounterVec counters_ig_tf(int palette, const CounterVec& g);
CounterVec counters_ig_tt(int palette, const CounterVec& g);
CounterVec counters_close(int palette, const CounterVec& c1_u);
CounterVec counters_reverse(int palette, const CounterVec& c);
CounterVec counters_gg(int palette,
                       const CounterVec& left,
                       const CounterVec& right);
CounterVec counters_cc(int palette,
                       const CounterVec& higher,
                       const CounterVec& lower);

/// The channels in the order listed above; the third is 0 for K = 2.
std::array<double, 3> counters_readout(int palette,
                                       const CounterVec& closed);

// Exact colour-graded kernels used by the HDT. A component's colour mask is
// the set of colours among its leaves; a counter involving a colour outside
// the mask is identically zero. The active counters are packed in increasing
// full-coordinate order. Packing changes storage only; every generated
// scalar reduction keeps its term order.
bool counters_coordinate_active(int palette,
                                CounterInterface interface,
                                int coordinate,
                                std::uint8_t color_mask);
std::uint64_t counters_active_coordinate_bits(
    int palette,
    CounterInterface interface,
    std::uint8_t color_mask);
// Packed HDT storage. Interface C1 is the direct sum of the D-only block
// (the draft's C3) and the U-only block (the draft's C1); a packed C1 payload
// stores only the block of its orientation, the other block being identically
// zero. Other interfaces are packed exactly as counters_active_coordinate_bits.
std::uint64_t counters_packed_coordinate_bits(int palette,
                                              CounterInterface interface,
                                              CounterOrientation orientation,
                                              std::uint8_t color_mask);
int counters_packed_coordinate_count(int palette,
                                     CounterInterface interface,
                                     CounterOrientation orientation,
                                     std::uint8_t color_mask);
std::uint8_t counters_value_color_mask(int palette,
                                       const CounterVec& value);
void counters_leaf_color_mask_into(int palette,
                                   int color,
                                   double weight,
                                   std::uint8_t color_mask,
                                   double* output);
void counters_c_to_g_color_mask_into(int palette,
                                     std::uint8_t color_mask,
                                     const double* input,
                                     double* output);
void counters_ig_ff_color_mask_into(int palette,
                                    std::uint8_t color_mask,
                                    const double* input,
                                    double* output);
void counters_ig_tf_color_mask_into(int palette,
                                    std::uint8_t color_mask,
                                    const double* input,
                                    double* output);
void counters_ig_tt_color_mask_into(int palette,
                                    std::uint8_t color_mask,
                                    const double* input,
                                    double* output);
void counters_close_color_mask_into(int palette,
                                    std::uint8_t color_mask,
                                    const double* input,
                                    double* output);
void counters_reverse_c1_color_mask_into(int palette,
                                         CounterOrientation input_orientation,
                                         std::uint8_t color_mask,
                                         const double* input,
                                         double* output);
void counters_reverse_c2_color_mask_into(int palette,
                                         std::uint8_t color_mask,
                                         const double* input,
                                         double* output);
void counters_gg_color_mask_into(int palette,
                                 std::uint8_t left_mask,
                                 std::uint8_t right_mask,
                                 const double* left,
                                 const double* right,
                                 double* output);
void counters_cc_color_mask_into(int palette,
                                 CounterInterface higher_interface,
                                 CounterInterface lower_interface,
                                 std::uint8_t higher_mask,
                                 std::uint8_t lower_mask,
                                 const double* higher,
                                 const double* lower,
                                 double* output);
std::array<double, 3> counters_color_mask_readout(
    int palette,
    std::uint8_t color_mask,
    const double* closed);

CounterPair counters_c_to_g(int palette, const CounterPair& c1_u);
CounterPair counters_ig_ff(int palette, const CounterPair& g);
CounterPair counters_ig_tf(int palette, const CounterPair& g);
CounterPair counters_ig_tt(int palette, const CounterPair& g);
CounterPair counters_close(int palette, const CounterPair& c1_u);
CounterPair counters_reverse(int palette, const CounterPair& c);
CounterPair counters_gg(int palette,
                        const CounterPair& left,
                        const CounterPair& right);
CounterPair counters_cc(int palette,
                        const CounterPair& higher,
                        const CounterPair& lower);

}  // namespace qraft
