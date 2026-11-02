#include "counting/counters.hpp"

#include "counting/generated/counters_k2.hpp"
#include "counting/generated/counters_k4.hpp"

#include <bit>
#include <stdexcept>
#include <string>
#include <utility>

namespace qraft {

namespace {

// Generated scalar terms evaluated on this thread, for run statistics.
thread_local long long g_counter_scalar_fmas = 0;

CounterVec make_full_descriptor(int palette,
                                CounterInterface interface,
                                CounterOrientation orientation) {
  CounterVec out;
  out.interface = interface;
  out.orientation = orientation;
  out.dimension = static_cast<std::uint8_t>(
      counters_dimension(palette, interface));
  if (out.dimension == 0) {
    throw std::invalid_argument("unsupported counter palette");
  }
  return out;
}

CounterVec make_full_value(int palette,
                           CounterInterface interface,
                           CounterOrientation orientation) {
  CounterVec out = make_full_descriptor(palette, interface, orientation);
  std::fill_n(out.values.data(), out.dimension, 0.0);
  return out;
}

void require_full_value(int palette,
                        const CounterVec& value,
                        CounterInterface interface,
                        CounterOrientation orientation,
                        const char* operation) {
  const int expected = counters_dimension(palette, interface);
  if (expected <= 0 || value.interface != interface ||
      value.orientation != orientation || value.dimension != expected) {
    throw std::invalid_argument(std::string(operation) +
                                ": invalid counter interface");
  }
}

template <auto Eval, class Terms>
CounterVec generated_full_unary(int palette,
                                const CounterVec& input,
                                CounterInterface output_interface,
                                CounterOrientation output_orientation,
                                const Terms& terms) {
  CounterVec output =
      make_full_descriptor(palette, output_interface, output_orientation);
  g_counter_scalar_fmas += static_cast<long long>(terms.size());
  Eval(input.values.data(), output.values.data());
  return output;
}

template <auto Eval, class Terms>
CounterVec generated_full_binary(int palette,
                                 const CounterVec& left,
                                 const CounterVec& right,
                                 CounterInterface output_interface,
                                 CounterOrientation output_orientation,
                                 const Terms& terms) {
  CounterVec output =
      make_full_descriptor(palette, output_interface, output_orientation);
  g_counter_scalar_fmas += static_cast<long long>(terms.size());
  Eval(left.values.data(), right.values.data(), output.values.data());
  return output;
}

template <auto Eval, class Terms>
CounterPair generated_full_unary_dual(
    int palette,
    const CounterPair& input,
    CounterInterface output_interface,
    CounterOrientation output_orientation,
    const Terms& terms) {
  CounterPair output;
  output.current =
      make_full_descriptor(palette, output_interface, output_orientation);
  output.background =
      make_full_descriptor(palette, output_interface, output_orientation);
  output.current_equals_background = false;
  g_counter_scalar_fmas +=
      2 * static_cast<long long>(terms.size());
  Eval(input.current.values.data(), input.background.values.data(),
       output.current.values.data(), output.background.values.data());
  return output;
}

template <auto Eval, class Terms>
CounterPair generated_full_binary_dual(
    int palette,
    const CounterPair& left,
    const CounterPair& right,
    CounterInterface output_interface,
    CounterOrientation output_orientation,
    const Terms& terms) {
  CounterPair output;
  output.current =
      make_full_descriptor(palette, output_interface, output_orientation);
  output.background =
      make_full_descriptor(palette, output_interface, output_orientation);
  output.current_equals_background = false;
  g_counter_scalar_fmas +=
      2 * static_cast<long long>(terms.size());
  Eval(left.current.values.data(), left.background.values.data(),
       right.current.values.data(), right.background.values.data(),
       output.current.values.data(), output.background.values.data());
  return output;
}

template <class Terms>
CounterVec sparse_full_constant(int palette,
                                CounterInterface interface,
                                CounterOrientation orientation,
                                const Terms& terms) {
  CounterVec output = make_full_value(palette, interface, orientation);
  g_counter_scalar_fmas += static_cast<long long>(terms.size());
  for (const auto& term : terms) {
    output.values[term.index] = term.sign > 0 ? 1.0 : -1.0;
  }
  return output;
}

// Terms are grouped by output channel; each group is summed in generated
// order before it is added to its channel.
template <class Terms>
std::array<double, 3> sparse_full_readout(const CounterVec& closed,
                                          const Terms& terms) {
  std::array<double, 3> output{0.0, 0.0, 0.0};
  g_counter_scalar_fmas += static_cast<long long>(terms.size());
  if (terms.empty()) return output;
  std::uint8_t active_out = terms.front().out;
  double sum = 0.0;
  for (const auto& term : terms) {
    if (term.out != active_out) {
      output[active_out] += sum;
      active_out = term.out;
      sum = 0.0;
    }
    const double value = closed.values[term.in];
    sum += term.sign > 0 ? value : -value;
  }
  output[active_out] += sum;
  return output;
}

}  // namespace

int counters_dimension(int palette, CounterInterface interface) {
  if (palette == 2) {
    namespace n = generated_counters_k2;
    switch (interface) {
      case CounterInterface::Closed: return n::dim_Closed;
      case CounterInterface::C1: return n::dim_C1;
      case CounterInterface::C2: return n::dim_C2;
      case CounterInterface::G: return n::dim_G;
    }
  }
  if (palette == 4) {
    namespace n = generated_counters_k4;
    switch (interface) {
      case CounterInterface::Closed: return n::dim_Closed;
      case CounterInterface::C1: return n::dim_C1;
      case CounterInterface::C2: return n::dim_C2;
      case CounterInterface::G: return n::dim_G;
    }
  }
  return 0;
}

int counters_max_dimension(int palette) {
  return std::max(
      {counters_dimension(palette, CounterInterface::Closed),
       counters_dimension(palette, CounterInterface::C1),
       counters_dimension(palette, CounterInterface::C2),
       counters_dimension(palette, CounterInterface::G)});
}

long long counters_thread_scalar_fmas() {
  return g_counter_scalar_fmas;
}

CounterVec counters_leaf(int palette, int color, double weight) {
  if (color < 0 || color >= palette) {
    throw std::invalid_argument("LEAF_TO_C: color outside palette");
  }
  CounterVec output;
  if (palette == 2) {
    namespace n = generated_counters_k2;
    output = sparse_full_constant(palette, CounterInterface::C1,
                                  CounterOrientation::U,
                                  n::LEAF_ABSENT_C1);
    const auto add = [&](const auto& terms) {
      g_counter_scalar_fmas += static_cast<long long>(terms.size());
      for (const auto& term : terms) {
        output.values[term.index] += term.sign > 0 ? weight : -weight;
      }
    };
    if (color == 0) add(n::LEAF_PRESENT_0_C1);
    else add(n::LEAF_PRESENT_1_C1);
    return output;
  }
  if (palette == 4) {
    namespace n = generated_counters_k4;
    output = sparse_full_constant(palette, CounterInterface::C1,
                                  CounterOrientation::U,
                                  n::LEAF_ABSENT_C1);
    const auto add = [&](const auto& terms) {
      g_counter_scalar_fmas += static_cast<long long>(terms.size());
      for (const auto& term : terms) {
        output.values[term.index] += term.sign > 0 ? weight : -weight;
      }
    };
    switch (color) {
      case 0: add(n::LEAF_PRESENT_0_C1); break;
      case 1: add(n::LEAF_PRESENT_1_C1); break;
      case 2: add(n::LEAF_PRESENT_2_C1); break;
      case 3: add(n::LEAF_PRESENT_3_C1); break;
      default: break;
    }
    return output;
  }
  throw std::invalid_argument("unsupported counter palette");
}

// A leaf of weight zero has only its unit counter, so both colourings agree.
CounterPair counters_leaf_dual(int palette,
                               int current_color,
                               int background_color,
                               double weight) {
  return CounterPair{counters_leaf(palette, current_color, weight),
                     counters_leaf(palette, background_color, weight),
                     current_color == background_color || weight == 0.0};
}

CounterVec counters_identity_c2(int palette) {
  if (palette == 2) {
    return sparse_full_constant(palette, CounterInterface::C2,
                                CounterOrientation::None,
                                generated_counters_k2::IDENTITY_C2);
  }
  if (palette == 4) {
    return sparse_full_constant(palette, CounterInterface::C2,
                                CounterOrientation::None,
                                generated_counters_k4::IDENTITY_C2);
  }
  throw std::invalid_argument("unsupported counter palette");
}

CounterVec counters_c_to_g(int palette, const CounterVec& c1_u) {
  require_full_value(palette, c1_u, CounterInterface::C1,
                     CounterOrientation::U, "C_TO_G");
  if (palette == 2) {
    return generated_full_unary<
        generated_counters_k2::eval_C_TO_G_C1_G>(
        palette, c1_u, CounterInterface::G, CounterOrientation::None,
        generated_counters_k2::C_TO_G_C1_G);
  }
  return generated_full_unary<generated_counters_k4::eval_C_TO_G_C1_G>(
      palette, c1_u, CounterInterface::G, CounterOrientation::None,
      generated_counters_k4::C_TO_G_C1_G);
}

CounterVec counters_ig_ff(int palette, const CounterVec& g) {
  require_full_value(palette, g, CounterInterface::G,
                     CounterOrientation::None, "IG_FF");
  if (palette == 2) {
    return generated_full_unary<
        generated_counters_k2::eval_IG_FF_G_CLOSED>(
        palette, g, CounterInterface::Closed, CounterOrientation::None,
        generated_counters_k2::IG_FF_G_CLOSED);
  }
  return generated_full_unary<generated_counters_k4::eval_IG_FF_G_CLOSED>(
      palette, g, CounterInterface::Closed, CounterOrientation::None,
      generated_counters_k4::IG_FF_G_CLOSED);
}

CounterVec counters_ig_tf(int palette, const CounterVec& g) {
  require_full_value(palette, g, CounterInterface::G,
                     CounterOrientation::None, "IG_TF");
  if (palette == 2) {
    return generated_full_unary<generated_counters_k2::eval_IG_TF_G_C1>(
        palette, g, CounterInterface::C1, CounterOrientation::U,
        generated_counters_k2::IG_TF_G_C1);
  }
  return generated_full_unary<generated_counters_k4::eval_IG_TF_G_C1>(
      palette, g, CounterInterface::C1, CounterOrientation::U,
      generated_counters_k4::IG_TF_G_C1);
}

CounterVec counters_ig_tt(int palette, const CounterVec& g) {
  require_full_value(palette, g, CounterInterface::G,
                     CounterOrientation::None, "IG_TT");
  if (palette == 2) {
    return generated_full_unary<generated_counters_k2::eval_IG_TT_G_C2>(
        palette, g, CounterInterface::C2, CounterOrientation::None,
        generated_counters_k2::IG_TT_G_C2);
  }
  return generated_full_unary<generated_counters_k4::eval_IG_TT_G_C2>(
      palette, g, CounterInterface::C2, CounterOrientation::None,
      generated_counters_k4::IG_TT_G_C2);
}

CounterVec counters_close(int palette, const CounterVec& c1_u) {
  require_full_value(palette, c1_u, CounterInterface::C1,
                     CounterOrientation::U, "CLOSE");
  if (palette == 2) {
    return generated_full_unary<
        generated_counters_k2::eval_CLOSE_C1_CLOSED>(
        palette, c1_u, CounterInterface::Closed, CounterOrientation::None,
        generated_counters_k2::CLOSE_C1_CLOSED);
  }
  return generated_full_unary<generated_counters_k4::eval_CLOSE_C1_CLOSED>(
      palette, c1_u, CounterInterface::Closed, CounterOrientation::None,
      generated_counters_k4::CLOSE_C1_CLOSED);
}

CounterVec counters_reverse(int palette, const CounterVec& c) {
  if (c.interface == CounterInterface::C1) {
    if (c.orientation != CounterOrientation::U &&
        c.orientation != CounterOrientation::D) {
      throw std::invalid_argument("REVERSE C1: missing orientation");
    }
    require_full_value(palette, c, CounterInterface::C1, c.orientation,
                       "REVERSE C1");
    const CounterOrientation output_orientation =
        c.orientation == CounterOrientation::U ? CounterOrientation::D
                                                : CounterOrientation::U;
    if (palette == 2) {
      return generated_full_unary<generated_counters_k2::eval_REV_C1_C1>(
          palette, c, CounterInterface::C1, output_orientation,
          generated_counters_k2::REV_C1_C1);
    }
    return generated_full_unary<generated_counters_k4::eval_REV_C1_C1>(
        palette, c, CounterInterface::C1, output_orientation,
        generated_counters_k4::REV_C1_C1);
  }
  require_full_value(palette, c, CounterInterface::C2,
                     CounterOrientation::None, "REVERSE C2");
  if (palette == 2) {
    return generated_full_unary<generated_counters_k2::eval_REV_C2_C2>(
        palette, c, CounterInterface::C2, CounterOrientation::None,
        generated_counters_k2::REV_C2_C2);
  }
  return generated_full_unary<generated_counters_k4::eval_REV_C2_C2>(
      palette, c, CounterInterface::C2, CounterOrientation::None,
      generated_counters_k4::REV_C2_C2);
}

CounterVec counters_gg(int palette,
                       const CounterVec& left,
                       const CounterVec& right) {
  require_full_value(palette, left, CounterInterface::G,
                     CounterOrientation::None, "GG_TO_G left");
  require_full_value(palette, right, CounterInterface::G,
                     CounterOrientation::None, "GG_TO_G right");
  if (palette == 2) {
    return generated_full_binary<generated_counters_k2::eval_GG_TO_G_G_G_G>(
        palette, left, right, CounterInterface::G, CounterOrientation::None,
        generated_counters_k2::GG_TO_G_G_G_G);
  }
  return generated_full_binary<generated_counters_k4::eval_GG_TO_G_G_G_G>(
      palette, left, right, CounterInterface::G, CounterOrientation::None,
      generated_counters_k4::GG_TO_G_G_G_G);
}

CounterVec counters_cc(int palette,
                       const CounterVec& higher,
                       const CounterVec& lower) {
  const bool higher_c1 = higher.interface == CounterInterface::C1;
  const bool lower_c1 = lower.interface == CounterInterface::C1;
  if (higher_c1) {
    require_full_value(palette, higher, CounterInterface::C1,
                       CounterOrientation::D, "CC_TO_C higher");
  } else {
    require_full_value(palette, higher, CounterInterface::C2,
                       CounterOrientation::None, "CC_TO_C higher");
  }
  if (lower_c1) {
    require_full_value(palette, lower, CounterInterface::C1,
                       CounterOrientation::U, "CC_TO_C lower");
  } else {
    require_full_value(palette, lower, CounterInterface::C2,
                       CounterOrientation::None, "CC_TO_C lower");
  }

  const CounterInterface output_interface =
      higher_c1 && lower_c1 ? CounterInterface::Closed
                            : (higher_c1 || lower_c1
                                   ? CounterInterface::C1
                                   : CounterInterface::C2);
  const CounterOrientation output_orientation =
      output_interface != CounterInterface::C1
          ? CounterOrientation::None
          : (higher_c1 ? CounterOrientation::D : CounterOrientation::U);

  if (palette == 2) {
    namespace n = generated_counters_k2;
    if (higher_c1 && lower_c1) {
      return generated_full_binary<n::eval_CC_TO_C_C1_C1_CLOSED>(
          palette, higher, lower, output_interface, output_orientation,
          n::CC_TO_C_C1_C1_CLOSED);
    }
    if (higher_c1) {
      return generated_full_binary<n::eval_CC_TO_C_C1_C2_C1>(
          palette, higher, lower, output_interface, output_orientation,
          n::CC_TO_C_C1_C2_C1);
    }
    if (lower_c1) {
      return generated_full_binary<n::eval_CC_TO_C_C2_C1_C1>(
          palette, higher, lower, output_interface, output_orientation,
          n::CC_TO_C_C2_C1_C1);
    }
    return generated_full_binary<n::eval_CC_TO_C_C2_C2_C2>(
        palette, higher, lower, output_interface, output_orientation,
        n::CC_TO_C_C2_C2_C2);
  }
  namespace n = generated_counters_k4;
  if (higher_c1 && lower_c1) {
    return generated_full_binary<n::eval_CC_TO_C_C1_C1_CLOSED>(
        palette, higher, lower, output_interface, output_orientation,
        n::CC_TO_C_C1_C1_CLOSED);
  }
  if (higher_c1) {
    return generated_full_binary<n::eval_CC_TO_C_C1_C2_C1>(
        palette, higher, lower, output_interface, output_orientation,
        n::CC_TO_C_C1_C2_C1);
  }
  if (lower_c1) {
    return generated_full_binary<n::eval_CC_TO_C_C2_C1_C1>(
        palette, higher, lower, output_interface, output_orientation,
        n::CC_TO_C_C2_C1_C1);
  }
  return generated_full_binary<n::eval_CC_TO_C_C2_C2_C2>(
      palette, higher, lower, output_interface, output_orientation,
      n::CC_TO_C_C2_C2_C2);
}

std::array<double, 3> counters_readout(int palette,
                                       const CounterVec& closed) {
  require_full_value(palette, closed, CounterInterface::Closed,
                     CounterOrientation::None, "READOUT");
  if (palette == 2) {
    return sparse_full_readout(closed,
                               generated_counters_k2::READOUT_CLOSED);
  }
  return sparse_full_readout(closed,
                             generated_counters_k4::READOUT_CLOSED);
}

namespace {

template <class Masks>
std::uint8_t value_color_mask_impl(const CounterVec& value,
                                   const Masks& masks) {
  std::uint8_t mask = 0;
  for (int coordinate = 0; coordinate < value.dimension; ++coordinate) {
    if (value.values[static_cast<std::size_t>(coordinate)] != 0.0) {
      mask |= masks[static_cast<std::size_t>(coordinate)];
    }
  }
  return mask;
}

template <class Terms>
std::array<double, 3> color_mask_readout_impl(
    int palette,
    std::uint8_t color_mask,
    const double* closed,
    const Terms& terms) {
  std::array<double, 3> output{0.0, 0.0, 0.0};
  if (color_mask >= (1u << palette)) {
    throw std::invalid_argument("counter color mask outside palette");
  }
  const PackedCounterLayout& layout = counters_packed_layout(
      palette, CounterInterface::Closed, CounterOrientation::None,
      color_mask);
  for (const auto& term : terms) {
    // A counter outside the packed layout is identically zero.
    const std::uint8_t packed = layout.rank[term.in];
    if (packed == kCounterNotPacked) continue;
    ++g_counter_scalar_fmas;
    output[term.out] +=
        term.sign > 0 ? closed[packed] : -closed[packed];
  }
  return output;
}

}  // namespace

std::uint64_t counters_active_coordinate_bits(
    int palette,
    CounterInterface interface,
    std::uint8_t color_mask) {
  if (palette != 2 && palette != 4) {
    throw std::invalid_argument("unsupported counter palette");
  }
  if (color_mask >= (1u << palette)) {
    throw std::invalid_argument("counter color mask outside palette");
  }
#define QRAFT_ACTIVE_BITS(IFACE)                                             \
  (palette == 2                                                             \
       ? generated_counters_k2::active_coordinate_bits_##IFACE[color_mask]   \
       : generated_counters_k4::active_coordinate_bits_##IFACE[color_mask])
  switch (interface) {
    case CounterInterface::Closed: return QRAFT_ACTIVE_BITS(Closed);
    case CounterInterface::C1: return QRAFT_ACTIVE_BITS(C1);
    case CounterInterface::C2: return QRAFT_ACTIVE_BITS(C2);
    case CounterInterface::G: return QRAFT_ACTIVE_BITS(G);
  }
#undef QRAFT_ACTIVE_BITS
  return 0;
}

namespace {

template <std::size_t N>
constexpr std::array<PackedCounterLayout, 16> packed_layout_row(
    const std::array<std::uint64_t, N>& bits) {
  std::array<PackedCounterLayout, 16> row{};
  // Masks beyond the palette (4..15 for K = 2) get an empty layout.
  for (std::size_t mask = 0; mask < 16; ++mask) {
    PackedCounterLayout layout;
    layout.bits = mask < N ? bits[mask] : 0;
    int rank = 0;
    for (std::size_t coordinate = 0; coordinate < kMaxCounters;
         ++coordinate) {
      layout.rank[coordinate] =
          ((layout.bits >> coordinate) & 1u) != 0
              ? static_cast<std::uint8_t>(rank++)
              : kCounterNotPacked;
    }
    layout.width = static_cast<std::uint8_t>(rank);
    row[mask] = layout;
  }
  return row;
}

constexpr PackedCounterLayoutTable make_packed_layouts() {
  PackedCounterLayoutTable table{};
  table[0] = {packed_layout_row(generated_counters_k2::active_coordinate_bits_Closed),
              packed_layout_row(generated_counters_k2::packed_coordinate_bits_C1U),
              packed_layout_row(generated_counters_k2::packed_coordinate_bits_C1D),
              packed_layout_row(generated_counters_k2::active_coordinate_bits_C2),
              packed_layout_row(generated_counters_k2::active_coordinate_bits_G)};
  table[1] = {packed_layout_row(generated_counters_k4::active_coordinate_bits_Closed),
              packed_layout_row(generated_counters_k4::packed_coordinate_bits_C1U),
              packed_layout_row(generated_counters_k4::packed_coordinate_bits_C1D),
              packed_layout_row(generated_counters_k4::active_coordinate_bits_C2),
              packed_layout_row(generated_counters_k4::active_coordinate_bits_G)};
  return table;
}

}  // namespace

extern const PackedCounterLayoutTable kPackedCounterLayouts =
    make_packed_layouts();

std::uint64_t counters_packed_coordinate_bits(int palette,
                                              CounterInterface interface,
                                              CounterOrientation orientation,
                                              std::uint8_t color_mask) {
  if (interface != CounterInterface::C1) {
    return counters_active_coordinate_bits(palette, interface, color_mask);
  }
  if (palette != 2 && palette != 4) {
    throw std::invalid_argument("unsupported counter palette");
  }
  if (color_mask >= (1u << palette)) {
    throw std::invalid_argument("counter color mask outside palette");
  }
  if (orientation == CounterOrientation::U) {
    return palette == 2 ? generated_counters_k2::packed_coordinate_bits_C1U[color_mask]
                        : generated_counters_k4::packed_coordinate_bits_C1U[color_mask];
  }
  if (orientation == CounterOrientation::D) {
    return palette == 2 ? generated_counters_k2::packed_coordinate_bits_C1D[color_mask]
                        : generated_counters_k4::packed_coordinate_bits_C1D[color_mask];
  }
  throw std::logic_error("packed C1 payload lacks an orientation");
}

int counters_packed_coordinate_count(int palette,
                                     CounterInterface interface,
                                     CounterOrientation orientation,
                                     std::uint8_t color_mask) {
  return std::popcount(counters_packed_coordinate_bits(
      palette, interface, orientation, color_mask));
}

bool counters_coordinate_active(int palette,
                                CounterInterface interface,
                                int coordinate,
                                std::uint8_t color_mask) {
  if (coordinate < 0 || coordinate >= counters_dimension(palette, interface)) {
    throw std::out_of_range("counter coordinate outside interface");
  }
  return (counters_active_coordinate_bits(palette, interface, color_mask) &
          (std::uint64_t{1} << coordinate)) != 0;
}

std::uint8_t counters_value_color_mask(int palette,
                                       const CounterVec& value) {
#define QRAFT_VALUE_MASK(IFACE)                                              \
  (palette == 2                                                             \
       ? value_color_mask_impl(                                             \
             value, generated_counters_k2::coordinate_color_mask_##IFACE)     \
       : value_color_mask_impl(                                             \
             value, generated_counters_k4::coordinate_color_mask_##IFACE))
  switch (value.interface) {
    case CounterInterface::Closed: return QRAFT_VALUE_MASK(Closed);
    case CounterInterface::C1: return QRAFT_VALUE_MASK(C1);
    case CounterInterface::C2: return QRAFT_VALUE_MASK(C2);
    case CounterInterface::G: return QRAFT_VALUE_MASK(G);
  }
#undef QRAFT_VALUE_MASK
  return 0;
}

void counters_leaf_color_mask_into(int palette,
                                   int color,
                                   double weight,
                                   std::uint8_t color_mask,
                                   double* output) {
  if (color < 0 || color >= palette || output == nullptr ||
      color_mask >= (1u << palette)) {
    throw std::invalid_argument("invalid color-graded leaf");
  }
  const PackedCounterLayout& layout = counters_packed_layout(
      palette, CounterInterface::C1, CounterOrientation::U, color_mask);
  std::fill_n(output, layout.width, 0.0);
  const auto add = [&](const auto& terms, double scale) {
    for (const auto& term : terms) {
      const std::uint8_t packed = layout.rank[term.index];
      if (packed == kCounterNotPacked) continue;
      ++g_counter_scalar_fmas;
      output[packed] +=
          term.sign > 0 ? scale : -scale;
    }
  };
  if (palette == 2) {
    namespace n = generated_counters_k2;
    add(n::LEAF_ABSENT_C1, 1.0);
    add(color == 0 ? n::LEAF_PRESENT_0_C1 : n::LEAF_PRESENT_1_C1,
        weight);
    return;
  }
  namespace n = generated_counters_k4;
  add(n::LEAF_ABSENT_C1, 1.0);
  switch (color) {
    case 0: add(n::LEAF_PRESENT_0_C1, weight); break;
    case 1: add(n::LEAF_PRESENT_1_C1, weight); break;
    case 2: add(n::LEAF_PRESENT_2_C1, weight); break;
    case 3: add(n::LEAF_PRESENT_3_C1, weight); break;
    default: break;
  }
}

#define QRAFT_COLOR_MASK_UNARY(NAME, ID)                                     \
  void NAME(int palette, std::uint8_t color_mask, const double* input,      \
            double* output) {                                               \
    int terms = 0;                                                          \
    if (palette == 2) {                                                     \
      terms = generated_counters_k2::eval_color_mask_##ID(                    \
          color_mask, input, output);                                       \
    } else {                                                                \
      terms = generated_counters_k4::eval_color_mask_##ID(                    \
          color_mask, input, output);                                       \
    }                                                                       \
    g_counter_scalar_fmas += terms;                                        \
  }

QRAFT_COLOR_MASK_UNARY(counters_c_to_g_color_mask_into, C_TO_G_C1_G)
QRAFT_COLOR_MASK_UNARY(counters_ig_ff_color_mask_into, IG_FF_G_CLOSED)
QRAFT_COLOR_MASK_UNARY(counters_ig_tf_color_mask_into, IG_TF_G_C1)
QRAFT_COLOR_MASK_UNARY(counters_ig_tt_color_mask_into, IG_TT_G_C2)
QRAFT_COLOR_MASK_UNARY(counters_close_color_mask_into, CLOSE_C1_CLOSED)
QRAFT_COLOR_MASK_UNARY(counters_reverse_c1_ud_color_mask_into, REV_C1_C1_UD)
QRAFT_COLOR_MASK_UNARY(counters_reverse_c1_du_color_mask_into, REV_C1_C1_DU)
QRAFT_COLOR_MASK_UNARY(counters_reverse_c2_color_mask_into, REV_C2_C2)

#undef QRAFT_COLOR_MASK_UNARY

void counters_reverse_c1_color_mask_into(int palette,
                                         CounterOrientation input_orientation,
                                         std::uint8_t color_mask,
                                         const double* input,
                                         double* output) {
  if (input_orientation == CounterOrientation::U) {
    counters_reverse_c1_ud_color_mask_into(palette, color_mask, input, output);
  } else if (input_orientation == CounterOrientation::D) {
    counters_reverse_c1_du_color_mask_into(palette, color_mask, input, output);
  } else {
    throw std::logic_error("C1 reversal input lacks an orientation");
  }
}

void counters_gg_color_mask_into(int palette,
                                 std::uint8_t left_mask,
                                 std::uint8_t right_mask,
                                 const double* left,
                                 const double* right,
                                 double* output) {
  const int terms = palette == 2
                        ? generated_counters_k2::eval_color_mask_GG_TO_G_G_G_G(
                              left_mask, right_mask, left, right, output)
                        : generated_counters_k4::eval_color_mask_GG_TO_G_G_G_G(
                              left_mask, right_mask, left, right, output);
  g_counter_scalar_fmas += terms;
}

void counters_cc_color_mask_into(int palette,
                                 CounterInterface higher_interface,
                                 CounterInterface lower_interface,
                                 std::uint8_t higher_mask,
                                 std::uint8_t lower_mask,
                                 const double* higher,
                                 const double* lower,
                                 double* output) {
  const bool higher_c1 = higher_interface == CounterInterface::C1;
  const bool lower_c1 = lower_interface == CounterInterface::C1;
  if ((!higher_c1 && higher_interface != CounterInterface::C2) ||
      (!lower_c1 && lower_interface != CounterInterface::C2)) {
    throw std::invalid_argument("invalid color-graded CC interface");
  }
#define QRAFT_COLOR_MASK_CC_DISPATCH(NS)                                     \
  (higher_c1 && lower_c1                                                    \
       ? NS::eval_color_mask_CC_TO_C_C1_C1_CLOSED(                         \
             higher_mask, lower_mask, higher, lower, output)                \
       : (higher_c1                                                         \
              ? NS::eval_color_mask_CC_TO_C_C1_C2_C1(                      \
                    higher_mask, lower_mask, higher, lower, output)         \
              : (lower_c1                                                  \
                     ? NS::eval_color_mask_CC_TO_C_C2_C1_C1(               \
                           higher_mask, lower_mask, higher, lower, output)  \
                     : NS::eval_color_mask_CC_TO_C_C2_C2_C2(               \
                           higher_mask, lower_mask, higher, lower, output))))
  const int terms = palette == 2
                        ? QRAFT_COLOR_MASK_CC_DISPATCH(generated_counters_k2)
                        : QRAFT_COLOR_MASK_CC_DISPATCH(generated_counters_k4);
#undef QRAFT_COLOR_MASK_CC_DISPATCH
  g_counter_scalar_fmas += terms;
}

std::array<double, 3> counters_color_mask_readout(
    int palette,
    std::uint8_t color_mask,
    const double* closed) {
  if (palette == 2) {
    return color_mask_readout_impl(
        palette, color_mask, closed,
        generated_counters_k2::READOUT_CLOSED);
  }
  return color_mask_readout_impl(
      palette, color_mask, closed,
      generated_counters_k4::READOUT_CLOSED);
}

#define QRAFT_FULL_BDP_DUAL_UNARY(NAME, INPUT_INTERFACE, INPUT_ORIENTATION,   \
                                 OUTPUT_INTERFACE, OUTPUT_ORIENTATION, ID)  \
  CounterPair NAME(int palette, const CounterPair& input) {                 \
    if (input.current_equals_background) {                                  \
      CounterVec current = NAME(palette, input.current);                     \
      return CounterPair{current, current, true};                            \
    }                                                                        \
    require_full_value(palette, input.current, INPUT_INTERFACE,              \
                       INPUT_ORIENTATION, #NAME " current");                \
    require_full_value(palette, input.background, INPUT_INTERFACE,           \
                       INPUT_ORIENTATION, #NAME " background");             \
    if (palette == 2) {                                                       \
      namespace n = generated_counters_k2;                                    \
      return generated_full_unary_dual<n::eval_dual_##ID>(                  \
          palette, input, OUTPUT_INTERFACE, OUTPUT_ORIENTATION, n::ID);      \
    }                                                                        \
    namespace n = generated_counters_k4;                                      \
    return generated_full_unary_dual<n::eval_dual_##ID>(                    \
        palette, input, OUTPUT_INTERFACE, OUTPUT_ORIENTATION, n::ID);        \
  }

QRAFT_FULL_BDP_DUAL_UNARY(
    counters_c_to_g, CounterInterface::C1, CounterOrientation::U,
    CounterInterface::G, CounterOrientation::None, C_TO_G_C1_G)
QRAFT_FULL_BDP_DUAL_UNARY(
    counters_ig_ff, CounterInterface::G, CounterOrientation::None,
    CounterInterface::Closed, CounterOrientation::None, IG_FF_G_CLOSED)
QRAFT_FULL_BDP_DUAL_UNARY(
    counters_ig_tf, CounterInterface::G, CounterOrientation::None,
    CounterInterface::C1, CounterOrientation::U, IG_TF_G_C1)
QRAFT_FULL_BDP_DUAL_UNARY(
    counters_ig_tt, CounterInterface::G, CounterOrientation::None,
    CounterInterface::C2, CounterOrientation::None, IG_TT_G_C2)
QRAFT_FULL_BDP_DUAL_UNARY(
    counters_close, CounterInterface::C1, CounterOrientation::U,
    CounterInterface::Closed, CounterOrientation::None, CLOSE_C1_CLOSED)

#undef QRAFT_FULL_BDP_DUAL_UNARY

CounterPair counters_reverse(int palette, const CounterPair& input) {
  if (input.current_equals_background) {
    CounterVec current = counters_reverse(palette, input.current);
    return CounterPair{current, current, true};
  }
  if (input.current.interface == CounterInterface::C1) {
    if (input.current.orientation != input.background.orientation ||
        (input.current.orientation != CounterOrientation::U &&
         input.current.orientation != CounterOrientation::D)) {
      throw std::invalid_argument("REVERSE C1 dual orientation mismatch");
    }
    require_full_value(palette, input.current, CounterInterface::C1,
                       input.current.orientation, "REVERSE C1 current");
    require_full_value(palette, input.background, CounterInterface::C1,
                       input.current.orientation, "REVERSE C1 background");
    const CounterOrientation output_orientation =
        input.current.orientation == CounterOrientation::U
            ? CounterOrientation::D
            : CounterOrientation::U;
    if (palette == 2) {
      namespace n = generated_counters_k2;
      return generated_full_unary_dual<n::eval_dual_REV_C1_C1>(
          palette, input, CounterInterface::C1, output_orientation,
          n::REV_C1_C1);
    }
    namespace n = generated_counters_k4;
    return generated_full_unary_dual<n::eval_dual_REV_C1_C1>(
        palette, input, CounterInterface::C1, output_orientation,
        n::REV_C1_C1);
  }
  require_full_value(palette, input.current, CounterInterface::C2,
                     CounterOrientation::None, "REVERSE C2 current");
  require_full_value(palette, input.background, CounterInterface::C2,
                     CounterOrientation::None, "REVERSE C2 background");
  if (palette == 2) {
    namespace n = generated_counters_k2;
    return generated_full_unary_dual<n::eval_dual_REV_C2_C2>(
        palette, input, CounterInterface::C2, CounterOrientation::None,
        n::REV_C2_C2);
  }
  namespace n = generated_counters_k4;
  return generated_full_unary_dual<n::eval_dual_REV_C2_C2>(
      palette, input, CounterInterface::C2, CounterOrientation::None,
      n::REV_C2_C2);
}

CounterPair counters_gg(int palette,
                        const CounterPair& left,
                        const CounterPair& right) {
  const bool identical = left.current_equals_background &&
                         right.current_equals_background;
  if (identical) {
    CounterVec current = counters_gg(palette, left.current, right.current);
    return CounterPair{current, current, true};
  }
  for (const CounterVec* value : {&left.current, &left.background,
                                  &right.current, &right.background}) {
    require_full_value(palette, *value, CounterInterface::G,
                       CounterOrientation::None, "GG_TO_G dual");
  }
  if (palette == 2) {
    namespace n = generated_counters_k2;
    return generated_full_binary_dual<n::eval_dual_GG_TO_G_G_G_G>(
        palette, left, right, CounterInterface::G, CounterOrientation::None,
        n::GG_TO_G_G_G_G);
  }
  namespace n = generated_counters_k4;
  return generated_full_binary_dual<n::eval_dual_GG_TO_G_G_G_G>(
      palette, left, right, CounterInterface::G, CounterOrientation::None,
      n::GG_TO_G_G_G_G);
}

CounterPair counters_cc(int palette,
                        const CounterPair& higher,
                        const CounterPair& lower) {
  const bool identical = higher.current_equals_background &&
                         lower.current_equals_background;
  if (identical) {
    CounterVec current =
        counters_cc(palette, higher.current, lower.current);
    return CounterPair{current, current, true};
  }
  const bool higher_c1 = higher.current.interface == CounterInterface::C1;
  const bool lower_c1 = lower.current.interface == CounterInterface::C1;
  const CounterInterface higher_interface =
      higher_c1 ? CounterInterface::C1 : CounterInterface::C2;
  const CounterInterface lower_interface =
      lower_c1 ? CounterInterface::C1 : CounterInterface::C2;
  const CounterOrientation higher_orientation =
      higher_c1 ? CounterOrientation::D : CounterOrientation::None;
  const CounterOrientation lower_orientation =
      lower_c1 ? CounterOrientation::U : CounterOrientation::None;
  require_full_value(palette, higher.current, higher_interface,
                     higher_orientation, "CC_TO_C higher current");
  require_full_value(palette, higher.background, higher_interface,
                     higher_orientation, "CC_TO_C higher background");
  require_full_value(palette, lower.current, lower_interface,
                     lower_orientation, "CC_TO_C lower current");
  require_full_value(palette, lower.background, lower_interface,
                     lower_orientation, "CC_TO_C lower background");
  const CounterInterface output_interface =
      higher_c1 && lower_c1
          ? CounterInterface::Closed
          : (higher_c1 || lower_c1 ? CounterInterface::C1
                                   : CounterInterface::C2);
  const CounterOrientation output_orientation =
      output_interface != CounterInterface::C1
          ? CounterOrientation::None
          : (higher_c1 ? CounterOrientation::D : CounterOrientation::U);
#define QRAFT_DUAL_CC_DISPATCH(NS)                                            \
  if (higher_c1 && lower_c1) {                                               \
    return generated_full_binary_dual<                                      \
        NS::eval_dual_CC_TO_C_C1_C1_CLOSED>(                                \
        palette, higher, lower, output_interface, output_orientation,        \
        NS::CC_TO_C_C1_C1_CLOSED);                                          \
  }                                                                          \
  if (higher_c1) {                                                           \
    return generated_full_binary_dual<NS::eval_dual_CC_TO_C_C1_C2_C1>(      \
        palette, higher, lower, output_interface, output_orientation,        \
        NS::CC_TO_C_C1_C2_C1);                                              \
  }                                                                          \
  if (lower_c1) {                                                            \
    return generated_full_binary_dual<NS::eval_dual_CC_TO_C_C2_C1_C1>(      \
        palette, higher, lower, output_interface, output_orientation,        \
        NS::CC_TO_C_C2_C1_C1);                                              \
  }                                                                          \
  return generated_full_binary_dual<NS::eval_dual_CC_TO_C_C2_C2_C2>(       \
      palette, higher, lower, output_interface, output_orientation,          \
      NS::CC_TO_C_C2_C2_C2)
  if (palette == 2) {
    namespace n = generated_counters_k2;
    QRAFT_DUAL_CC_DISPATCH(n);
  }
  namespace n = generated_counters_k4;
  QRAFT_DUAL_CC_DISPATCH(n);
#undef QRAFT_DUAL_CC_DISPATCH
}

}  // namespace qraft
