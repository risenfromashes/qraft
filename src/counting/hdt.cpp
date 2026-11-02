#include "counting/hdt.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <deque>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace qraft::internal {
namespace {

// Sort of a component in the tree being built: I, C (C1, C2 or Closed) or G.
enum class ActiveKind : std::uint8_t {
  I,
  C,
  G,
};

// At most two ports of a component. MutableActiveTree::external_ports holds
// those at the graph boundary.
struct ActivePortList {
  std::array<HdtPortKey, 2> values{};
  std::uint8_t count = 0;
};

void append_external_port(ActivePortList& ports,
                          const HdtPortKey& port,
                          bool unique = false) {
  if (unique) {
    for (int i = 0; i < ports.count; ++i) {
      if (ports.values[static_cast<std::size_t>(i)] == port) return;
    }
  }
  if (ports.count >= ports.values.size()) {
    throw std::logic_error("active component exceeds two external ports");
  }
  ports.values[ports.count++] = port;
}

struct MutableActiveLink {
  int child = -1;
  int previous = -1;
  int next = -1;
  HdtPortKey parent_port;
  HdtPortKey child_port;
};

struct MutableActiveNode {
  int hdt_root = -1;
  int source_center = -1;
  int parent = -1;
  int parent_link = -1;
  int first_child = -1;
  int last_child = -1;
  ActiveKind kind = ActiveKind::I;
};

// The tree of components during construction, indexed by active node. A
// node is active while hdt_root >= 0; its children are the links from
// first_child to last_child, and each link records the ports of both ends.
struct MutableActiveTree {
  std::vector<MutableActiveNode> nodes;
  std::vector<ActivePortList> external_ports;
  std::vector<MutableActiveLink> links;
  int root = -1;
  int active_nodes = 0;
};

enum class PlannedOperation : std::uint8_t {
  Survivor,
  GG,
  IG,
  CC,
};

// A composition planned for the current round; a is the I, the higher C or
// the first G.
struct PlannedGroup {
  PlannedOperation operation = PlannedOperation::Survivor;
  int a = -1;
  int b = -1;
};

struct RoundWorkspace {
  std::vector<int> parent;
  std::vector<int> parent_edge;
  std::vector<int> preorder;
  std::vector<int> traversal_stack;
  std::vector<std::uint8_t> used;
  std::vector<int> mutable_used_epoch;
  std::vector<PlannedGroup> compositions;
  std::vector<PlannedGroup> ig_compositions;
  std::vector<PlannedGroup> groups;
  std::vector<int> group_scratch;
  std::vector<int> group_of;
  std::vector<int> center_edge_seen;
};

template <typename T>
void add_vector_capacity_bytes(std::size_t& bytes,
                               const std::vector<T>& values) {
  checked_hdt_accumulate_bytes(
      bytes, values.capacity(), sizeof(T),
      "HDT builder scratch byte count exceeds size_t");
}

std::size_t live_builder_scratch_bytes(const MutableActiveTree& tree,
                                       const RoundWorkspace& workspace) {
  std::size_t bytes = 0;
  add_vector_capacity_bytes(bytes, tree.nodes);
  add_vector_capacity_bytes(bytes, tree.external_ports);
  add_vector_capacity_bytes(bytes, tree.links);
  add_vector_capacity_bytes(bytes, workspace.parent);
  add_vector_capacity_bytes(bytes, workspace.parent_edge);
  add_vector_capacity_bytes(bytes, workspace.preorder);
  add_vector_capacity_bytes(bytes, workspace.traversal_stack);
  add_vector_capacity_bytes(bytes, workspace.used);
  add_vector_capacity_bytes(bytes, workspace.mutable_used_epoch);
  add_vector_capacity_bytes(bytes, workspace.compositions);
  add_vector_capacity_bytes(bytes, workspace.ig_compositions);
  add_vector_capacity_bytes(bytes, workspace.groups);
  add_vector_capacity_bytes(bytes, workspace.group_scratch);
  add_vector_capacity_bytes(bytes, workspace.group_of);
  add_vector_capacity_bytes(bytes, workspace.center_edge_seen);
  return bytes;
}

HdtOutputKind output_kind_of(const CounterVec& value) {
  switch (value.interface) {
    case CounterInterface::Closed: return HdtOutputKind::Closed;
    case CounterInterface::C1:
      return value.orientation == CounterOrientation::U
                 ? HdtOutputKind::C1U
                 : HdtOutputKind::C1D;
    case CounterInterface::C2: return HdtOutputKind::C2;
    case CounterInterface::G: return HdtOutputKind::G;
  }
  throw std::logic_error("unknown counter output interface");
}

HdtOutputKind output_kind_of(const HdtComponentUnit& unit) {
  switch (unit.kind) {
    case HdtComponentKind::IAnchor:
      return HdtOutputKind::I;
    case HdtComponentKind::ClosedOpaque:
      return HdtOutputKind::Closed;
    case HdtComponentKind::LiveLeaf:
      return HdtOutputKind::C1U;
    case HdtComponentKind::C1Opaque:
      return unit.dp_orientation == CounterOrientation::D
                 ? HdtOutputKind::C1D
                 : HdtOutputKind::C1U;
    case HdtComponentKind::C2Opaque:
      return HdtOutputKind::C2;
    case HdtComponentKind::GOpaque:
      return HdtOutputKind::G;
  }
  throw std::logic_error("unknown component output kind");
}

[[noreturn]] [[gnu::noinline]] void throw_payloadless_kind() {
  throw std::logic_error("I node has no counters");
}

// Indexed by HdtOutputKind; the I entry is never returned.
constexpr std::array<std::pair<CounterInterface, CounterOrientation>, 6>
    kOutputKindDescriptors{{
        {CounterInterface::Closed, CounterOrientation::None},
        {CounterInterface::Closed, CounterOrientation::None},
        {CounterInterface::C1, CounterOrientation::U},
        {CounterInterface::C1, CounterOrientation::D},
        {CounterInterface::C2, CounterOrientation::None},
        {CounterInterface::G, CounterOrientation::None},
    }};

inline std::pair<CounterInterface, CounterOrientation> descriptor_for_kind(
    HdtOutputKind kind) {
  const auto index = static_cast<std::size_t>(kind);
  if (kind == HdtOutputKind::I || index >= kOutputKindDescriptors.size())
      [[unlikely]] {
    throw_payloadless_kind();
  }
  return kOutputKindDescriptors[index];
}

// Slot width of a packed value. A slot holds at least one value, which
// carries the free-list link once the slot is released.
inline int packed_payload_width(int palette,
                                CounterInterface interface,
                                CounterOrientation orientation,
                                std::uint8_t color_mask) noexcept {
  return std::max<int>(1, counters_packed_layout(palette, interface,
                                                 orientation, color_mask)
                              .width);
}

inline int packed_coordinate_count(int palette,
                                   CounterInterface interface,
                                   CounterOrientation orientation,
                                   std::uint8_t color_mask) noexcept {
  return counters_packed_layout(palette, interface, orientation, color_mask)
      .width;
}

// The 15 width classes. They must cover every packed width of an
// (interface, orientation, colour mask) after the round-up to one.
constexpr std::array<std::uint8_t, kMaxCounters + 1>
    kPackedPayloadClassByWidth = [] {
      std::array<std::uint8_t, kMaxCounters + 1> table{};
      table.fill(0xff);
      constexpr std::array<int, 15> widths{
          1, 2, 3, 4, 5, 7, 10, 12, 13, 15, 22, 23, 26, 32, 56};
      for (std::size_t i = 0; i < widths.size(); ++i) {
        table[static_cast<std::size_t>(widths[i])] =
            static_cast<std::uint8_t>(i);
      }
      return table;
    }();

[[noreturn]] [[gnu::noinline]] void throw_unreachable_width() {
  throw std::logic_error("unreachable packed-payload width");
}

inline std::uint8_t packed_payload_class_index(int width) {
  const std::uint8_t index =
      width >= 0 && width <= static_cast<int>(kMaxCounters)
          ? kPackedPayloadClassByWidth[static_cast<std::size_t>(width)]
          : std::uint8_t{0xff};
  if (index == 0xff) [[unlikely]] throw_unreachable_width();
  return index;
}

// Inverse of kPackedPayloadClassByWidth.
constexpr std::array<int, 15> kPackedPayloadWidths{
    1, 2, 3, 4, 5, 7, 10, 12, 13, 15, 22, 23, 26, 32, 56};

int packed_payload_width(std::uint8_t payload_class) {
  if (payload_class >= kPackedPayloadWidths.size()) {
    throw std::logic_error("invalid packed-payload class");
  }
  return kPackedPayloadWidths[payload_class];
}

int allocate_packed_payload(HdtPackedPayloadArena& arena,
                            std::uint8_t class_index) {
  const int width = packed_payload_width(class_index);
  int& free_head = arena.free_head[static_cast<std::size_t>(class_index)];
  if (free_head >= 0) {
    const int offset = free_head;
    free_head =
        static_cast<int>(arena.values[static_cast<std::size_t>(offset)]);
    return offset;
  }
  const int offset = checked_hdt_int_size(
      arena.values.size(),
      "packed payload exceeds its integer offset domain");
  const std::size_t next_size = checked_hdt_size_add(
      arena.values.size(), static_cast<std::size_t>(width),
      "packed payload capacity overflow");
  (void)checked_hdt_int_size(
      next_size, "packed payload exceeds its integer offset domain");
  arena.values.resize(next_size);
  return offset;
}

void release_packed_payload(HdtPackedPayloadArena& arena,
                            std::uint8_t class_index,
                            int offset) {
  if (offset < 0) return;
  const int width = packed_payload_width(class_index);
  if (static_cast<std::size_t>(offset) + width > arena.values.size()) {
    throw std::logic_error("released payload outside arena");
  }
  int& free_head = arena.free_head[static_cast<std::size_t>(class_index)];
  arena.values[static_cast<std::size_t>(offset)] =
      static_cast<double>(free_head);
  free_head = offset;
}

void pack_full_value(int palette,
                     const CounterVec& value,
                     std::uint8_t color_mask,
                     double* output) {
  std::uint64_t active = counters_packed_coordinate_bits(
      palette, value.interface, value.orientation, color_mask);
  int packed = 0;
  while (active != 0) {
    const int coordinate = std::countr_zero(active);
    output[packed++] = value.values[static_cast<std::size_t>(coordinate)];
    active &= active - 1;
  }
}

const double* node_current_data(const Hdt& hdt,
                                const HdtNode& node) noexcept {
  return hdt.current_payload.values.data() + node.dp_offset;
}

double* node_current_data(Hdt& hdt,
                          const HdtNode& node) noexcept {
  return hdt.current_payload.values.data() + node.dp_offset;
}

const double* node_background_data(const Hdt& hdt,
                                   const HdtNode& node) noexcept {
  const bool aliases_current =
      hdt.background_alias_current || node.background_dp_offset < 0;
  const HdtPackedPayloadArena& arena =
      aliases_current ? hdt.current_payload : hdt.background_payload;
  const int offset =
      aliases_current ? node.dp_offset : node.background_dp_offset;
  return arena.values.data() + offset;
}

double* node_background_data(Hdt& hdt,
                             const HdtNode& node) noexcept {
  const bool aliases_current =
      hdt.background_alias_current || node.background_dp_offset < 0;
  HdtPackedPayloadArena& arena =
      aliases_current ? hdt.current_payload : hdt.background_payload;
  const int offset =
      aliases_current ? node.dp_offset : node.background_dp_offset;
  return arena.values.data() + offset;
}

CounterView node_current_view(const Hdt& hdt,
                                   const HdtNode& node) {
  const auto [interface, orientation] = descriptor_for_kind(node.output_kind);
  return CounterView{
      interface, orientation, node.dp_dimension,
      node_current_data(hdt, node)};
}

CounterView node_background_view(const Hdt& hdt,
                                      const HdtNode& node) {
  const auto [interface, orientation] = descriptor_for_kind(node.output_kind);
  return CounterView{
      interface, orientation, node.dp_dimension,
      node_background_data(hdt, node)};
}

void allocate_node_background(Hdt& hdt, HdtNode& node) {
  if (hdt.background_alias_current) {
    throw std::logic_error("cannot allocate a frozen HDT background");
  }
  if (node.background_dp_offset >= 0 || node.dp_dimension == 0) return;
  const auto [interface, orientation] = descriptor_for_kind(node.output_kind);
  node.background_payload_class = packed_payload_class_index(
      packed_payload_width(hdt.palette, interface, orientation,
                           node.background_color_mask));
  node.background_dp_offset =
      allocate_packed_payload(hdt.background_payload,
                              node.background_payload_class);
}

// Gives a node whose background aliases its current value a background slot
// of its own. Callers promote before overwriting the current value, while it
// and its colour mask still equal the background.
void promote_node_background(Hdt& hdt, HdtNode& node) {
  if (node.background_dp_offset >= 0 || node.dp_dimension == 0) return;
  allocate_node_background(hdt, node);
  const auto [interface, orientation] = descriptor_for_kind(node.output_kind);
  const int active = packed_coordinate_count(
      hdt.palette, interface, orientation, node.current_color_mask);
  std::copy_n(node_current_data(hdt, node), active,
              node_background_data(hdt, node));
}

void allocate_node_payload(Hdt& hdt,
                           HdtNode& node,
                           HdtOutputKind output_kind) {
  node.output_kind = output_kind;
  if (output_kind == HdtOutputKind::I) {
    node.dp_offset = -1;
    node.background_dp_offset = -1;
    node.current_payload_class = 0xff;
    node.background_payload_class = 0xff;
    node.dp_dimension = 0;
    node.current_equals_background = true;
    return;
  }
  const auto [interface, orientation] = descriptor_for_kind(output_kind);
  const int dimension = counters_dimension(hdt.palette, interface);
  node.current_payload_class = packed_payload_class_index(
      packed_payload_width(hdt.palette, interface, orientation,
                           node.current_color_mask));
  node.dp_offset =
      allocate_packed_payload(hdt.current_payload,
                              node.current_payload_class);
  node.background_dp_offset = -1;
  node.background_payload_class = 0xff;
  node.dp_dimension = static_cast<std::uint8_t>(dimension);
}

void store_node_payload(Hdt& hdt,
                        HdtNode& node,
                        const CounterPair& value,
                        HdtOutputKind output_kind) {
  node.current_color_mask =
      counters_value_color_mask(hdt.palette, value.current);
  node.background_color_mask = value.current_equals_background
                                   ? node.current_color_mask
                                   : counters_value_color_mask(
                                         hdt.palette, value.background);
  node.current_equals_background = value.current_equals_background;
  allocate_node_payload(hdt, node, output_kind);
  if (output_kind == HdtOutputKind::I) return;
  pack_full_value(hdt.palette, value.current, node.current_color_mask,
                  node_current_data(hdt, node));
  if (hdt.background_alias_current) {
    if (!value.current_equals_background) {
      throw std::logic_error("non-frozen atom in frozen HDT fold");
    }
  } else if (!value.current_equals_background) {
    allocate_node_background(hdt, node);
    pack_full_value(hdt.palette, value.background, node.background_color_mask,
                    node_background_data(hdt, node));
  }
}

// The two evaluate_* helpers allocate the output node's slots and apply
// `eval` to the packed current values, and to the backgrounds when they
// differ from the current values.
template <class Eval>
void evaluate_unary_payload_color_mask(Hdt& hdt,
                                       HdtNode& output,
                                       HdtOutputKind output_kind,
                                       int input_index,
                                       Eval&& eval) {
  const HdtNode& input =
      hdt.nodes[static_cast<std::size_t>(input_index)];
  output.current_color_mask = input.current_color_mask;
  output.background_color_mask = input.background_color_mask;
  allocate_node_payload(hdt, output, output_kind);
  eval(output.current_color_mask,
       node_current_data(hdt, input), node_current_data(hdt, output));
  output.current_equals_background = input.current_equals_background;
  if (!input.current_equals_background) {
    if (hdt.background_alias_current) {
      throw std::logic_error("non-frozen unary node in frozen HDT fold");
    }
    allocate_node_background(hdt, output);
    eval(output.background_color_mask,
         node_background_data(hdt, input),
         node_background_data(hdt, output));
  }
}

template <class Eval>
void evaluate_binary_payload_color_mask(Hdt& hdt,
                                        HdtNode& output,
                                        HdtOutputKind output_kind,
                                        int left_index,
                                        int right_index,
                                        Eval&& eval) {
  const HdtNode& left =
      hdt.nodes[static_cast<std::size_t>(left_index)];
  const HdtNode& right =
      hdt.nodes[static_cast<std::size_t>(right_index)];
  output.current_color_mask =
      left.current_color_mask | right.current_color_mask;
  output.background_color_mask =
      left.background_color_mask | right.background_color_mask;
  allocate_node_payload(hdt, output, output_kind);
  eval(left.current_color_mask, right.current_color_mask,
       node_current_data(hdt, left), node_current_data(hdt, right),
       node_current_data(hdt, output));
  output.current_equals_background =
      left.current_equals_background && right.current_equals_background;
  if (!output.current_equals_background) {
    if (hdt.background_alias_current) {
      throw std::logic_error("non-frozen binary node in frozen HDT fold");
    }
    allocate_node_background(hdt, output);
    const double* left_background = node_background_data(hdt, left);
    const double* right_background = node_background_data(hdt, right);
    eval(left.background_color_mask, right.background_color_mask,
         left_background, right_background,
         node_background_data(hdt, output));
  }
}

// Applies a unary transform to a node in place. The output dimension may not
// exceed the input's; results go to scratch first, and a value moves to a new
// slot only when its width class changes.
template <class Eval>
void replace_unary_payload_color_mask(Hdt& hdt,
                                      HdtNode& node,
                                      HdtOutputKind output_kind,
                                      Eval&& eval) {
  const HdtNode input = node;
  const auto [output_interface, output_orientation] =
      descriptor_for_kind(output_kind);
  const int output_dimension =
      counters_dimension(hdt.palette, output_interface);
  if (output_dimension > input.dp_dimension) {
    throw std::logic_error("in-place masked unary transform grew its type");
  }
  std::array<double, kMaxCounters> current;
  eval(input.current_color_mask,
       node_current_data(hdt, input), current.data());
  std::array<double, kMaxCounters> background;
  if (!input.current_equals_background) {
    eval(input.background_color_mask,
         node_background_data(hdt, input), background.data());
    if (node.background_dp_offset < 0) {
      throw std::logic_error("masked unary node lacks background payload");
    }
  }

  const std::uint8_t next_current_class = packed_payload_class_index(
      packed_payload_width(hdt.palette, output_interface, output_orientation,
                           input.current_color_mask));
  const bool relocate_current =
      input.current_payload_class != next_current_class;
  const std::uint8_t next_background_class =
      input.current_equals_background
          ? std::uint8_t{0xff}
          : packed_payload_class_index(packed_payload_width(
                hdt.palette, output_interface, output_orientation,
                input.background_color_mask));
  const bool relocate_background =
      !input.current_equals_background &&
      input.background_payload_class != next_background_class;
  if (relocate_current) {
    node.dp_offset =
        allocate_packed_payload(hdt.current_payload, next_current_class);
    node.current_payload_class = next_current_class;
  }
  if (relocate_background) {
    node.background_dp_offset = allocate_packed_payload(
        hdt.background_payload, next_background_class);
    node.background_payload_class = next_background_class;
  }

  node.output_kind = output_kind;
  node.dp_dimension = static_cast<std::uint8_t>(output_dimension);
  node.current_equals_background = input.current_equals_background;
  node.current_color_mask = input.current_color_mask;
  node.background_color_mask = input.background_color_mask;
  const int current_count = packed_coordinate_count(
      hdt.palette, output_interface, output_orientation,
      input.current_color_mask);
  std::copy_n(current.data(), current_count,
              node_current_data(hdt, node));
  if (!input.current_equals_background) {
    const int background_count = packed_coordinate_count(
        hdt.palette, output_interface, output_orientation,
        input.background_color_mask);
    std::copy_n(background.data(), background_count,
                node_background_data(hdt, node));
  }
  if (relocate_current) {
    release_packed_payload(hdt.current_payload,
                           input.current_payload_class, input.dp_offset);
  }
  if (relocate_background) {
    release_packed_payload(hdt.background_payload,
                           input.background_payload_class,
                           input.background_dp_offset);
  }
}

// Expands packed values to the full coordinate vector; the omitted
// coordinates are identically zero.
CounterVec materialize_view(CounterView view,
                            int palette,
                            std::uint8_t color_mask) {
  CounterVec value;
  value.interface = view.interface;
  value.orientation = view.orientation;
  value.dimension = view.dimension;
  std::fill_n(value.values.data(), view.dimension, 0.0);
  std::uint64_t active = counters_packed_coordinate_bits(
      palette, view.interface, view.orientation, color_mask);
  int packed = 0;
  while (active != 0) {
    const int coordinate = std::countr_zero(active);
    value.values[static_cast<std::size_t>(coordinate)] =
        view.values[packed++];
    active &= active - 1;
  }
  return value;
}

// IG counters from the I's open ports: none gives Closed, an Up port C1U, a
// Down port C1D and an Up and a Down port C2.
CounterPair ig_output(int palette,
                      const CounterPair& g,
                      const ActivePortList& ports) {
  if (ports.count == 0) return counters_ig_ff(palette, g);
  if (ports.count == 1) {
    CounterPair value = counters_ig_tf(palette, g);
    if (ports.values[0].orientation == HdtPortOrientation::Down) {
      value = counters_reverse(palette, value);
    } else if (ports.values[0].orientation !=
               HdtPortOrientation::Up) {
      throw std::logic_error("IG one-port output is not ordinary");
    }
    return value;
  }
  if (ports.count == 2) {
    bool have_up = false;
    bool have_down = false;
    for (const HdtPortKey& port : ports.values) {
      have_up = have_up || port.orientation == HdtPortOrientation::Up;
      have_down = have_down ||
                  port.orientation == HdtPortOrientation::Down;
    }
    if (!have_up || !have_down) {
      throw std::logic_error("IG C2 output lacks U/D ports");
    }
    return counters_ig_tt(palette, g);
  }
  throw std::logic_error("IG output exceeds two ordinary ports");
}

void note_node_shape(Hdt& hdt, const HdtNode& node) {
  const int degree = (node.left >= 0 ? 1 : 0) + (node.right >= 0 ? 1 : 0);
  hdt.stats.max_recompute_degree =
      std::max(hdt.stats.max_recompute_degree, degree);
  hdt.stats.max_update_height =
      std::max(hdt.stats.max_update_height, node.height);
  if (node.atom_count <= 1) return;
  // A node over a atoms may have height at most 7*ceil(log2 a) + 1. Since
  // 2*9^7 < 10^7, this integer bound exceeds log_(10/9)(a) from the SODA'13
  // round constant, only slightly, and keeps the local update height
  // O(log a).
  const int ceil_log2 = static_cast<int>(
      std::bit_width(static_cast<unsigned int>(node.atom_count - 1)));
  const int allowed = 7 * ceil_log2 + 1;
  if (node.height > allowed) ++hdt.stats.local_balance_violations;
}

// Append an atom and its leaf node. add_atom stores the given counters;
// add_component_atom takes the unit's stored counters, or computes a live
// leaf's from its colour and weight.
int add_atom(Hdt& hdt,
             HdtExpandedAtom atom,
             const CounterPair& dp,
             HdtOutputKind output_kind,
             int source_center) {
  const int atom_index = checked_hdt_int_size(
      hdt.atoms.size(), "HDT exceeds its atom index domain");
  hdt.atoms.push_back(std::move(atom));
  HdtNode node;
  node.atom_descriptor = atom_index;
  node.operation = HdtOperation::Atom;
  node.output_kind = output_kind;
  node.source_center_key = source_center;
  store_node_payload(hdt, node, dp, output_kind);
  node.live_owner = hdt.atoms[static_cast<std::size_t>(atom_index)].owner;
  const int node_index = checked_hdt_int_size(
      hdt.nodes.size(), "HDT exceeds its node index domain");
  note_node_shape(hdt, node);
  hdt.nodes.push_back(std::move(node));
  return node_index;
}

int add_component_atom(Hdt& hdt,
                       HdtExpandedAtom atom,
                       const HdtComponentGraph& graph,
                       int component_index,
                       int source_center) {
  const HdtComponentUnit& component =
      graph.units[static_cast<std::size_t>(component_index)];
  const HdtOutputKind output_kind = output_kind_of(component);
  const int atom_index = checked_hdt_int_size(
      hdt.atoms.size(), "HDT exceeds its atom index domain");
  hdt.atoms.push_back(std::move(atom));
  HdtNode node;
  node.atom_descriptor = atom_index;
  node.operation = HdtOperation::Atom;
  node.source_center_key = source_center;
  node.current_color_mask = component.current_color_mask;
  node.background_color_mask = component.background_color_mask;
  allocate_node_payload(hdt, node, output_kind);
  if (output_kind != HdtOutputKind::I) {
    if (component.dp_dimension != node.dp_dimension) {
      throw std::logic_error("component/HDT atom payload mismatch");
    }
    if (component.kind == HdtComponentKind::LiveLeaf) {
      counters_leaf_color_mask_into(
          hdt.palette, component.current_color, component.weight,
          component.current_color_mask,
          node_current_data(hdt, node));
      node.current_equals_background = component.current_equals_background;
      if (hdt.background_alias_current) {
        if (!component.current_equals_background) {
          throw std::logic_error("non-frozen component in frozen HDT fold");
        }
      } else if (!component.current_equals_background) {
        allocate_node_background(hdt, node);
        counters_leaf_color_mask_into(
            hdt.palette, component.background_color, component.weight,
            component.background_color_mask,
            node_background_data(hdt, node));
      }
    } else {
      if (component.dp_offset < 0) {
        throw std::logic_error("component/HDT atom payload missing");
      }
      const auto [interface, orientation] = descriptor_for_kind(output_kind);
      CounterVec current;
      current.interface = interface;
      current.orientation = orientation;
      current.dimension = node.dp_dimension;
      std::copy_n(graph.current_payload.data() + component.dp_offset,
                  node.dp_dimension, current.values.data());
      pack_full_value(hdt.palette, current, component.current_color_mask,
                      node_current_data(hdt, node));
      if (hdt.background_alias_current) {
        if (!component.current_equals_background) {
          throw std::logic_error("non-frozen component in frozen HDT fold");
        }
      } else if (!component.current_equals_background) {
        node.current_equals_background = false;
        allocate_node_background(hdt, node);
        const double* background =
            graph.background_payload.data() +
            component.background_dp_offset;
        CounterVec background_value;
        background_value.interface = interface;
        background_value.orientation = orientation;
        background_value.dimension = node.dp_dimension;
        std::copy_n(background, node.dp_dimension,
                    background_value.values.data());
        pack_full_value(hdt.palette, background_value,
                        component.background_color_mask,
                        node_background_data(hdt, node));
      }
      node.current_equals_background = component.current_equals_background;
    }
  }
  node.live_owner = hdt.atoms[static_cast<std::size_t>(atom_index)].owner;
  const int node_index = checked_hdt_int_size(
      hdt.nodes.size(), "HDT exceeds its node index domain");
  note_node_shape(hdt, node);
  hdt.nodes.push_back(std::move(node));
  return node_index;
}

// Value of a transformed node before its post-transform, from the current
// values below it. Only live-leaf atoms, IG and CC nodes can reach here.
void evaluate_transformed_base_color_mask_into(const Hdt& hdt,
                                               const HdtNode& node,
                                               double* output) {
  if (node.operation == HdtOperation::Atom) {
    const HdtExpandedAtom& atom =
        hdt.atoms[static_cast<std::size_t>(node.atom_descriptor)];
    if (atom.component_kind == HdtComponentKind::LiveLeaf) {
      counters_leaf_color_mask_into(
          hdt.palette, atom.current_color, atom.weight,
          atom.current_color_mask, output);
      return;
    }
    throw std::logic_error("immutable transformed atom was marked dirty");
  }
  const HdtNode& left =
      hdt.nodes[static_cast<std::size_t>(node.left)];
  const HdtNode& right =
      hdt.nodes[static_cast<std::size_t>(node.right)];
  switch (node.operation) {
    case HdtOperation::GG:
      throw std::logic_error("G node has an invalid C post-transform");
    case HdtOperation::IG:
      counters_ig_tf_color_mask_into(
          hdt.palette, right.current_color_mask,
          node_current_data(hdt, right), output);
      return;
    case HdtOperation::CC:
      counters_cc_color_mask_into(
          hdt.palette, descriptor_for_kind(left.output_kind).first,
          descriptor_for_kind(right.output_kind).first,
          left.current_color_mask, right.current_color_mask,
          node_current_data(hdt, left), node_current_data(hdt, right), output);
      return;
    case HdtOperation::Atom: break;
  }
  throw std::logic_error("invalid transformed HDT operation");
}

// Recomputes a dirty node's current value from its children, or from the
// colour of a live-leaf atom. A new colour mask can change the width class,
// which moves the value to a new slot.
void recompute_node(Hdt& hdt, int node_index) {
  HdtNode& node = hdt.nodes[static_cast<std::size_t>(node_index)];
  bool next_equals_background = false;
  std::uint8_t next_color_mask = 0;
  if (node.operation == HdtOperation::Atom) {
    const HdtExpandedAtom& atom =
        hdt.atoms[static_cast<std::size_t>(node.atom_descriptor)];
    next_equals_background =
        atom.current_color == atom.background_color || atom.weight == 0.0;
    next_color_mask = atom.current_color_mask;
  } else {
    next_equals_background =
        hdt.nodes[static_cast<std::size_t>(node.left)]
            .current_equals_background &&
        hdt.nodes[static_cast<std::size_t>(node.right)]
            .current_equals_background;
    next_color_mask =
        hdt.nodes[static_cast<std::size_t>(node.left)].current_color_mask |
        hdt.nodes[static_cast<std::size_t>(node.right)].current_color_mask;
  }
  // Backgrounds never change; copy an aliased one out before the current
  // value moves away from it.
  if (!next_equals_background && node.background_dp_offset < 0) {
    promote_node_background(hdt, node);
  }
  if (node.output_kind == HdtOutputKind::I) {
    node.current_equals_background = next_equals_background;
    node.current_color_mask = next_color_mask;
    return;
  }

  const auto [output_interface, output_orientation] =
      descriptor_for_kind(node.output_kind);
  const int old_dp_offset = node.dp_offset;
  const std::uint8_t old_payload_class = node.current_payload_class;
  const std::uint8_t next_payload_class = packed_payload_class_index(
      packed_payload_width(hdt.palette, output_interface, output_orientation,
                           next_color_mask));
  const bool relocated = old_payload_class != next_payload_class;
  if (relocated) {
    node.dp_offset =
        allocate_packed_payload(hdt.current_payload, next_payload_class);
    node.current_payload_class = next_payload_class;
  }
  node.current_color_mask = next_color_mask;
  double* const output = node_current_data(hdt, node);

  if (node.post_transform != HdtTransform::None) {
    std::array<double, kMaxCounters> base;
    evaluate_transformed_base_color_mask_into(hdt, node, base.data());
    if (node.post_transform == HdtTransform::CToG) {
      counters_c_to_g_color_mask_into(
          hdt.palette, next_color_mask, base.data(), output);
    } else if (node.post_transform == HdtTransform::CloseRoot) {
      counters_close_color_mask_into(
          hdt.palette, next_color_mask, base.data(), output);
    } else {
      throw std::logic_error("unknown HDT transform");
    }
  } else if (node.operation == HdtOperation::Atom) {
    const HdtExpandedAtom& atom =
        hdt.atoms[static_cast<std::size_t>(node.atom_descriptor)];
    if (atom.component_kind == HdtComponentKind::LiveLeaf) {
      counters_leaf_color_mask_into(
          hdt.palette, atom.current_color, atom.weight, next_color_mask,
          output);
    }
  } else {
    const HdtNode& left =
        hdt.nodes[static_cast<std::size_t>(node.left)];
    const HdtNode& right =
        hdt.nodes[static_cast<std::size_t>(node.right)];
    switch (node.operation) {
      case HdtOperation::GG:
        counters_gg_color_mask_into(
            hdt.palette, left.current_color_mask, right.current_color_mask,
            node_current_data(hdt, left), node_current_data(hdt, right),
            output);
        break;
      case HdtOperation::IG:
        if (node.output_kind == HdtOutputKind::Closed) {
          counters_ig_ff_color_mask_into(
              hdt.palette, right.current_color_mask,
              node_current_data(hdt, right), output);
        } else if (node.output_kind == HdtOutputKind::C2) {
          counters_ig_tt_color_mask_into(
              hdt.palette, right.current_color_mask,
              node_current_data(hdt, right), output);
        } else if (node.output_kind == HdtOutputKind::C1U) {
          counters_ig_tf_color_mask_into(
              hdt.palette, right.current_color_mask,
              node_current_data(hdt, right), output);
        } else {
          std::array<double, kMaxCounters> c1_u;
          counters_ig_tf_color_mask_into(
              hdt.palette, right.current_color_mask,
              node_current_data(hdt, right), c1_u.data());
          counters_reverse_c1_color_mask_into(
              hdt.palette, CounterOrientation::U, right.current_color_mask,
              c1_u.data(), output);
        }
        break;
      case HdtOperation::CC:
        counters_cc_color_mask_into(
            hdt.palette,
            descriptor_for_kind(left.output_kind).first,
            descriptor_for_kind(right.output_kind).first,
            left.current_color_mask, right.current_color_mask,
            node_current_data(hdt, left), node_current_data(hdt, right),
            output);
        break;
      case HdtOperation::Atom: break;
    }
  }
  if (relocated) {
    release_packed_payload(hdt.current_payload, old_payload_class,
                           old_dp_offset);
  }
  node.current_equals_background = next_equals_background;
}

void set_mutable_hdt_output(Hdt& hdt,
                            const MutableActiveTree& tree,
                            int active_index);

// Builds the atoms and the rooted mutable tree from the component graph.
// Active node u < U is component unit u; node U + j is the local background
// G of the j-th anchor (in unit order) that has one, joined by edge E + j
// after the graph's E edges. Each node's incidence is its unit's ascending
// edge list followed by its local edge. The DFS pushes it in reverse, so
// preorder and every child list follow incidence order.
MutableActiveTree build_mutable_active_tree(const HdtComponentGraph& graph,
                                            std::size_t expanded_atom_count,
                                            Hdt& hdt,
                                            RoundWorkspace& workspace,
                                            std::size_t& scratch_peak_bytes) {
  const int unit_count = graph.component_units();
  if (expanded_atom_count < graph.units.size()) {
    throw std::logic_error("active atom capacity underflow");
  }
  const int node_count = checked_hdt_int_size(
      expanded_atom_count, "active graph exceeds node domain");
  const int edge_count = checked_hdt_int_size(
      graph.edges.size(), "active graph exceeds edge domain");
  (void)checked_hdt_int_size(
      checked_hdt_size_add(graph.edges.size(),
                           expanded_atom_count - graph.units.size(),
                           "active edge capacity overflow"),
      "active graph exceeds edge domain");

  MutableActiveTree tree;
  tree.nodes.resize(static_cast<std::size_t>(node_count));
  tree.external_ports.resize(static_cast<std::size_t>(node_count));
  for (const HdtExternalPort& external : graph.external_ports) {
    append_external_port(
        tree.external_ports[static_cast<std::size_t>(external.unit)],
        external.port);
  }

  for (int unit_index = 0; unit_index < unit_count; ++unit_index) {
    const HdtComponentUnit& unit =
        graph.units[static_cast<std::size_t>(unit_index)];
    HdtExpandedAtom atom;
    atom.component_kind = unit.kind;
    atom.component_unit = unit_index;
    atom.owner = unit.owner;
    atom.source_vertex_id = unit.source_vertex_id;
    atom.weight = unit.weight;
    atom.current_color = unit.current_color;
    atom.background_color = unit.background_color;
    atom.current_color_mask = unit.current_color_mask;
    atom.background_color_mask = unit.background_color_mask;
    atom.port_begin = checked_hdt_int_size(
        hdt.atom_ports.size(), "HDT exceeds its port index domain");
    for (int i = 0; i < unit.incidence_count; ++i) {
      const int edge_index = graph.incidence_edges[static_cast<std::size_t>(
          unit.incidence_begin + i)];
      const HdtIncidenceEdge& edge =
          graph.edges[static_cast<std::size_t>(edge_index)];
      hdt.atom_ports.push_back(edge.a == unit_index ? edge.a_port : edge.b_port);
    }
    const ActivePortList& unit_external =
        tree.external_ports[static_cast<std::size_t>(unit_index)];
    for (int i = 0; i < unit_external.count; ++i) {
      hdt.atom_ports.push_back(
          unit_external.values[static_cast<std::size_t>(i)]);
    }
    atom.port_count = checked_hdt_int_size(
                          hdt.atom_ports.size(),
                          "HDT exceeds its port index domain") -
                      atom.port_begin;

    MutableActiveNode& node = tree.nodes[static_cast<std::size_t>(unit_index)];
    node.source_center = unit.source_vertex_id;
    if (unit.kind == HdtComponentKind::IAnchor) {
      node.kind = ActiveKind::I;
    } else if (unit.kind == HdtComponentKind::GOpaque) {
      node.kind = ActiveKind::G;
    } else {
      node.kind = ActiveKind::C;
    }
    node.hdt_root = add_component_atom(
        hdt, std::move(atom), graph, unit_index, node.source_center);
    if (unit_index == graph.root_unit &&
        unit.kind == HdtComponentKind::IAnchor) {
      hdt.root_anchor_atom_node = node.hdt_root;
    }
  }

  // local_of_unit[u] is the active node of u's local background G, or -1;
  // local_anchor[j] is the unit that owns local node U + j.
  std::vector<int> local_of_unit(static_cast<std::size_t>(unit_count), -1);
  std::vector<int> local_anchor;
  local_anchor.reserve(expanded_atom_count - graph.units.size());
  for (int unit_index = 0; unit_index < unit_count; ++unit_index) {
    const HdtComponentUnit& unit =
        graph.units[static_cast<std::size_t>(unit_index)];
    if (unit.kind != HdtComponentKind::IAnchor ||
        unit.local_background_g == nullptr) {
      continue;
    }
    HdtExpandedAtom atom;
    atom.component_kind = HdtComponentKind::GOpaque;
    atom.component_unit = unit_index;
    atom.source_vertex_id = unit.source_vertex_id;
    atom.port_begin = checked_hdt_int_size(
        hdt.atom_ports.size(), "HDT exceeds its port index domain");
    hdt.atom_ports.push_back(HdtPortKey{
        HdtPortKind::Center, unit.source_vertex_id,
        HdtPortOrientation::Center});
    atom.port_count = 1;
    const int local_index = unit_count +
                            static_cast<int>(local_anchor.size());
    if (local_index >= node_count) {
      throw std::logic_error("HDT atom count changed during expansion");
    }
    MutableActiveNode& local =
        tree.nodes[static_cast<std::size_t>(local_index)];
    local.kind = ActiveKind::G;
    local.source_center = unit.source_vertex_id;
    local.hdt_root = add_atom(hdt, std::move(atom), *unit.local_background_g,
                              HdtOutputKind::G, unit.source_vertex_id);
    local_of_unit[static_cast<std::size_t>(unit_index)] = local_index;
    local_anchor.push_back(unit_index);
  }
  if (unit_count + static_cast<int>(local_anchor.size()) != node_count) {
    throw std::logic_error("HDT atom count changed during expansion");
  }

  // Endpoints and ports of active edge e: graph edge e < E, else the
  // center edge of local node U + (e - E).
  const auto edge_endpoint = [&](int edge_index, int component, int& other,
                                 HdtPortKey& own_port,
                                 HdtPortKey& other_port) {
    if (edge_index < edge_count) {
      const HdtIncidenceEdge& edge =
          graph.edges[static_cast<std::size_t>(edge_index)];
      const bool is_a = edge.a == component;
      other = is_a ? edge.b : edge.a;
      own_port = is_a ? edge.a_port : edge.b_port;
      other_port = is_a ? edge.b_port : edge.a_port;
      return;
    }
    const int local = unit_count + (edge_index - edge_count);
    const int anchor = local_anchor[static_cast<std::size_t>(
        edge_index - edge_count)];
    const HdtPortKey center{HdtPortKind::Center,
                            graph.units[static_cast<std::size_t>(anchor)]
                                .source_vertex_id,
                            HdtPortOrientation::Center};
    other = component == anchor ? local : anchor;
    own_port = center;
    other_port = center;
  };

  std::vector<int>& parent = workspace.parent;
  std::vector<int>& parent_edge = workspace.parent_edge;
  std::vector<int>& preorder = workspace.preorder;
  std::vector<int>& stack = workspace.traversal_stack;
  parent.assign(static_cast<std::size_t>(node_count), -2);
  parent_edge.assign(static_cast<std::size_t>(node_count), -1);
  preorder.clear();
  preorder.reserve(static_cast<std::size_t>(node_count));
  tree.root = graph.root_unit;
  tree.active_nodes = node_count;
  if (tree.root >= 0 && tree.root < node_count) {
    parent[static_cast<std::size_t>(tree.root)] = -1;
    stack.clear();
    stack.reserve(static_cast<std::size_t>(node_count));
    stack.push_back(tree.root);
    const auto visit = [&](int component, int edge_index) {
      int next = -1;
      HdtPortKey own_port;
      HdtPortKey other_port;
      edge_endpoint(edge_index, component, next, own_port, other_port);
      if (parent[static_cast<std::size_t>(next)] != -2) return;
      parent[static_cast<std::size_t>(next)] = component;
      parent_edge[static_cast<std::size_t>(next)] = edge_index;
      stack.push_back(next);
    };
    while (!stack.empty()) {
      const int component = stack.back();
      stack.pop_back();
      preorder.push_back(component);
      if (component >= unit_count) {
        visit(component, edge_count + (component - unit_count));
        continue;
      }
      const int local = local_of_unit[static_cast<std::size_t>(component)];
      if (local >= 0) visit(component, edge_count + (local - unit_count));
      const HdtComponentUnit& unit =
          graph.units[static_cast<std::size_t>(component)];
      for (int i = unit.incidence_count - 1; i >= 0; --i) {
        visit(component, graph.incidence_edges[static_cast<std::size_t>(
                             unit.incidence_begin + i)]);
      }
    }
  }
  if (static_cast<int>(preorder.size()) != node_count) {
    throw std::logic_error("HDT active graph disconnected");
  }

  tree.links.reserve(static_cast<std::size_t>(node_count));
  for (int node = 0; node < node_count; ++node) {
    tree.nodes[static_cast<std::size_t>(node)].parent =
        parent[static_cast<std::size_t>(node)];
  }
  for (int child : preorder) {
    const int up = parent[static_cast<std::size_t>(child)];
    if (up < 0) continue;
    int other = -1;
    HdtPortKey child_port;
    HdtPortKey parent_port;
    edge_endpoint(parent_edge[static_cast<std::size_t>(child)], child, other,
                  child_port, parent_port);
    const int link_index = checked_hdt_int_size(
        tree.links.size(), "mutable tree exceeds link domain");
    tree.links.emplace_back();
    MutableActiveLink& link = tree.links.back();
    link.child = child;
    link.parent_port = parent_port;
    link.child_port = child_port;
    MutableActiveNode& parent_node = tree.nodes[static_cast<std::size_t>(up)];
    if (parent_node.first_child < 0) {
      parent_node.first_child = link_index;
    } else {
      const int previous = parent_node.last_child;
      tree.links[static_cast<std::size_t>(previous)].next = link_index;
      tree.links[static_cast<std::size_t>(link_index)].previous = previous;
    }
    parent_node.last_child = link_index;
    tree.nodes[static_cast<std::size_t>(child)].parent_link = link_index;
  }

  for (int node = 0; node < node_count; ++node) {
    set_mutable_hdt_output(hdt, tree, node);
  }
  std::size_t initial_scratch = live_builder_scratch_bytes(tree, workspace);
  add_vector_capacity_bytes(initial_scratch, local_of_unit);
  add_vector_capacity_bytes(initial_scratch, local_anchor);
  scratch_peak_bytes = std::max(scratch_peak_bytes, initial_scratch);
  return tree;
}

void mutable_unlink_child(MutableActiveTree& tree, int parent, int link) {
  MutableActiveLink& removed = tree.links[static_cast<std::size_t>(link)];
  const int previous = removed.previous;
  const int next = removed.next;
  if (previous < 0) {
    tree.nodes[static_cast<std::size_t>(parent)].first_child = next;
  } else {
    tree.links[static_cast<std::size_t>(previous)].next = next;
  }
  if (next < 0) {
    tree.nodes[static_cast<std::size_t>(parent)].last_child = previous;
  } else {
    tree.links[static_cast<std::size_t>(next)].previous = previous;
  }
  removed.previous = -1;
  removed.next = -1;
}

void set_mutable_hdt_output(Hdt& hdt,
                            const MutableActiveTree& tree,
                            int active_index) {
  const MutableActiveNode& active_node =
      tree.nodes[static_cast<std::size_t>(active_index)];
  HdtNode& node =
      hdt.nodes[static_cast<std::size_t>(active_node.hdt_root)];
  node.output_port_count = 0;
  node.source_center_key = active_node.source_center;
  if (active_node.kind == ActiveKind::I) {
    node.output_kind = HdtOutputKind::I;
    return;
  }
  auto append = [&](const HdtPortKey& port) {
    if (node.output_port_count >= node.output_ports.size()) {
      throw std::logic_error("mutable HDT component exceeds two ports");
    }
    node.output_ports[node.output_port_count++] = port;
  };
  const ActivePortList& external_ports =
      tree.external_ports[static_cast<std::size_t>(active_index)];
  for (int i = 0; i < external_ports.count; ++i) {
    append(external_ports.values[static_cast<std::size_t>(i)]);
  }
  if (active_node.parent_link >= 0) {
    append(tree.links[static_cast<std::size_t>(active_node.parent_link)]
               .child_port);
  }
  for (int link = active_node.first_child; link >= 0;
       link = tree.links[static_cast<std::size_t>(link)].next) {
    append(tree.links[static_cast<std::size_t>(link)].parent_port);
  }
}

// One round of SODA'13, Sec. 4, on the components present at its start:
// every C without a Down port whose parent is an I becomes a G; each I pairs
// its G children in child order (GG), and an I with one G child and fewer
// than two other Down arms absorbs it (IG); each maximal downward path of C
// components pairs them from the top (CC). A round must remove at least a
// ninth of the components.
void execute_mutable_round(MutableActiveTree& tree,
                           Hdt& hdt,
                           int round,
                           RoundWorkspace& workspace) {
  std::vector<int>& preorder = workspace.preorder;
  std::size_t active_end = 0;
  const std::size_t previous_size = preorder.size();
  for (std::size_t position = 0; position < previous_size; ++position) {
    const int component = preorder[position];
    MutableActiveNode& current =
        tree.nodes[static_cast<std::size_t>(component)];
    if (current.hdt_root < 0) continue;
    if (component != tree.root &&
        (current.parent < 0 ||
         tree.nodes[static_cast<std::size_t>(current.parent)].hdt_root < 0)) {
      throw std::logic_error("mutable HDT tree disconnected");
    }
    preorder[active_end++] = component;
    if (current.kind != ActiveKind::C || current.parent < 0 ||
        tree.nodes[static_cast<std::size_t>(current.parent)].kind !=
            ActiveKind::I) {
      continue;
    }
    bool down_open = current.first_child >= 0;
    const ActivePortList& current_external_ports =
        tree.external_ports[static_cast<std::size_t>(component)];
    for (int i = 0; i < current_external_ports.count && !down_open; ++i) {
      const HdtPortKey& port = current_external_ports.values[
          static_cast<std::size_t>(i)];
      if (port.kind == HdtPortKind::OrdinaryEdge &&
          port.orientation == HdtPortOrientation::Down) {
        down_open = true;
      }
    }
    if (down_open) continue;

    current.kind = ActiveKind::G;
    current.source_center =
        tree.nodes[static_cast<std::size_t>(current.parent)].source_center;
    HdtNode& root = hdt.nodes[static_cast<std::size_t>(
        current.hdt_root)];
    if (root.post_transform != HdtTransform::None) {
      throw std::logic_error("multiple fused transforms on HDT node");
    }
    root.post_transform = HdtTransform::CToG;
    replace_unary_payload_color_mask(
        hdt, root, HdtOutputKind::G,
        [&](std::uint8_t mask, const double* input, double* output) {
          counters_c_to_g_color_mask_into(
              hdt.palette, mask, input, output);
        });
    root.output_kind = HdtOutputKind::G;
    root.source_center_key = current.source_center;
    const HdtPortKey center{HdtPortKind::Center,
                                  current.source_center,
                                  HdtPortOrientation::Center};
    root.output_port_count = 1;
    root.output_ports[0] = center;
    MutableActiveLink& parent_link =
        tree.links[static_cast<std::size_t>(current.parent_link)];
    parent_link.parent_port = center;
    parent_link.child_port = center;
    ++hdt.stats.c_to_g_ops;
  }
  preorder.resize(active_end);
  if (static_cast<int>(active_end) != tree.active_nodes) {
    throw std::logic_error("mutable HDT preorder lost an active node");
  }

  const int node_slots = checked_hdt_int_size(
      tree.nodes.size(), "mutable tree exceeds node domain");
  std::vector<int>& used_epoch = workspace.mutable_used_epoch;
  std::vector<PlannedGroup>& compositions = workspace.compositions;
  std::vector<PlannedGroup>& ig_compositions = workspace.ig_compositions;
  if (used_epoch.size() != static_cast<std::size_t>(node_slots)) {
    used_epoch.assign(static_cast<std::size_t>(node_slots), 0);
  }
  const auto is_used = [&](int node) {
    return used_epoch[static_cast<std::size_t>(node)] == round;
  };
  auto mark_used = [&](int node) {
    used_epoch[static_cast<std::size_t>(node)] = round;
  };
  compositions.clear();
  compositions.reserve(static_cast<std::size_t>(tree.active_nodes));
  ig_compositions.clear();
  ig_compositions.reserve(static_cast<std::size_t>(tree.active_nodes));

  for (int i_component : preorder) {
    const MutableActiveNode& i_node =
        tree.nodes[static_cast<std::size_t>(i_component)];
    if (i_node.kind != ActiveKind::I) continue;
    int pending_g = -1;
    int g_count = 0;
    int down_open = 0;
    const ActivePortList& i_external_ports =
        tree.external_ports[static_cast<std::size_t>(i_component)];
    for (int i = 0; i < i_external_ports.count; ++i) {
      const HdtPortKey& port = i_external_ports.values[
          static_cast<std::size_t>(i)];
      if (port.kind == HdtPortKind::OrdinaryEdge &&
          port.orientation == HdtPortOrientation::Down) {
        ++down_open;
      }
    }
    for (int link = i_node.first_child; link >= 0;
         link = tree.links[static_cast<std::size_t>(link)].next) {
      const int child = tree.links[static_cast<std::size_t>(link)].child;
      if (tree.nodes[static_cast<std::size_t>(child)].kind != ActiveKind::G) {
        ++down_open;
        continue;
      }
      ++g_count;
      if (pending_g < 0) {
        pending_g = child;
      } else {
        mark_used(pending_g);
        mark_used(child);
        compositions.push_back(
            PlannedGroup{PlannedOperation::GG, pending_g, child});
        pending_g = -1;
      }
    }
    if (g_count == 1 && down_open < 2 && !is_used(pending_g)) {
      mark_used(i_component);
      mark_used(pending_g);
      ig_compositions.push_back(
          PlannedGroup{PlannedOperation::IG, i_component, pending_g});
    }
  }
  compositions.insert(compositions.end(), ig_compositions.begin(),
                      ig_compositions.end());

  for (int top : preorder) {
    if (tree.nodes[static_cast<std::size_t>(top)].kind != ActiveKind::C ||
        is_used(top)) {
      continue;
    }
    const int parent = tree.nodes[static_cast<std::size_t>(top)].parent;
    if (parent >= 0 &&
        tree.nodes[static_cast<std::size_t>(parent)].kind == ActiveKind::C) {
      continue;
    }
    workspace.group_scratch.clear();
    int cursor = top;
    while (cursor >= 0 &&
           tree.nodes[static_cast<std::size_t>(cursor)].kind == ActiveKind::C) {
      workspace.group_scratch.push_back(cursor);
      int next_c = -1;
      for (int link = tree.nodes[static_cast<std::size_t>(cursor)].first_child;
           link >= 0; link = tree.links[static_cast<std::size_t>(link)].next) {
        const int child = tree.links[static_cast<std::size_t>(link)].child;
        if (tree.nodes[static_cast<std::size_t>(child)].kind == ActiveKind::C) {
          next_c = child;
          break;
        }
      }
      cursor = next_c;
    }
    for (int i = 0; i + 1 < static_cast<int>(workspace.group_scratch.size());
         i += 2) {
      const int higher = workspace.group_scratch[static_cast<std::size_t>(i)];
      const int lower =
          workspace.group_scratch[static_cast<std::size_t>(i + 1)];
      if (is_used(higher) || is_used(lower)) {
        continue;
      }
      mark_used(higher);
      mark_used(lower);
      compositions.push_back(
          PlannedGroup{PlannedOperation::CC, higher, lower});
    }
  }

  if (compositions.empty()) {
    throw std::logic_error("mutable SODA round made no progress");
  }

  for (const PlannedGroup& group : compositions) {
    MutableActiveNode& a_node =
        tree.nodes[static_cast<std::size_t>(group.a)];
    MutableActiveNode& b_node =
        tree.nodes[static_cast<std::size_t>(group.b)];
    const int a_hdt_root = a_node.hdt_root;
    const int b_hdt_root = b_node.hdt_root;
    const int a_source_center = a_node.source_center;
    const ActivePortList a_external_ports =
        tree.external_ports[static_cast<std::size_t>(group.a)];
    const ActivePortList b_external_ports =
        tree.external_ports[static_cast<std::size_t>(group.b)];
    ActiveKind result_kind = ActiveKind::I;
    int result_source_center = -1;
    ActivePortList result_external_ports;
    const int result_hdt_root = checked_hdt_int_size(
        hdt.nodes.size(), "HDT exceeds its node index domain");
    hdt.nodes.emplace_back();
    HdtNode& node = hdt.nodes.back();
    node.atom_count = hdt.nodes[static_cast<std::size_t>(a_hdt_root)].atom_count +
                      hdt.nodes[static_cast<std::size_t>(b_hdt_root)].atom_count;
    node.height = 1 + std::max(
        hdt.nodes[static_cast<std::size_t>(a_hdt_root)].height,
        hdt.nodes[static_cast<std::size_t>(b_hdt_root)].height);

    if (group.operation == PlannedOperation::GG) {
      if (a_node.parent < 0 || a_node.parent != b_node.parent ||
          a_source_center != b_node.source_center) {
        throw std::logic_error("mutable GG operands do not share a center");
      }
      result_kind = ActiveKind::G;
      result_source_center = a_source_center;
      evaluate_binary_payload_color_mask(
          hdt, node, HdtOutputKind::G, a_hdt_root, b_hdt_root,
          [&](std::uint8_t left_mask, std::uint8_t right_mask,
              const double* left, const double* right, double* output) {
            counters_gg_color_mask_into(
                hdt.palette, left_mask, right_mask, left, right, output);
          });
      result_external_ports = a_external_ports;
      for (int i = 0; i < b_external_ports.count; ++i) {
        append_external_port(
            result_external_ports,
            b_external_ports.values[static_cast<std::size_t>(i)], true);
      }
      node.operation = HdtOperation::GG;
      node.left = a_hdt_root;
      node.right = b_hdt_root;
      ++hdt.stats.gg_ops;
      mutable_unlink_child(tree, a_node.parent, b_node.parent_link);
    } else if (group.operation == PlannedOperation::IG) {
      if (b_node.parent != group.a) {
        throw std::logic_error("mutable IG operands are not parent/child");
      }
      result_kind = ActiveKind::C;
      result_source_center = a_source_center;
      ActivePortList open_ports = a_external_ports;
      if (a_node.parent_link >= 0) {
        append_external_port(
            open_ports,
            tree.links[static_cast<std::size_t>(a_node.parent_link)]
                .child_port);
      }
      for (int link = a_node.first_child; link >= 0;
           link = tree.links[static_cast<std::size_t>(link)].next) {
        if (tree.links[static_cast<std::size_t>(link)].child == group.b) {
          continue;
        }
        append_external_port(
            open_ports,
            tree.links[static_cast<std::size_t>(link)].parent_port);
      }
      if (open_ports.count == 0) {
        evaluate_unary_payload_color_mask(
            hdt, node, HdtOutputKind::Closed, b_hdt_root,
            [&](std::uint8_t mask, const double* input, double* output) {
              counters_ig_ff_color_mask_into(
                  hdt.palette, mask, input, output);
            });
      } else if (open_ports.count == 1 &&
                 open_ports.values[0].orientation ==
                     HdtPortOrientation::Up) {
        evaluate_unary_payload_color_mask(
            hdt, node, HdtOutputKind::C1U, b_hdt_root,
            [&](std::uint8_t mask, const double* input, double* output) {
              counters_ig_tf_color_mask_into(
                  hdt.palette, mask, input, output);
            });
      } else if (open_ports.count == 2) {
        evaluate_unary_payload_color_mask(
            hdt, node, HdtOutputKind::C2, b_hdt_root,
            [&](std::uint8_t mask, const double* input, double* output) {
              counters_ig_tt_color_mask_into(
                  hdt.palette, mask, input, output);
            });
      } else {
        const CounterPair value = ig_output(
            hdt.palette, hdt_node_dp(hdt, b_hdt_root), open_ports);
        store_node_payload(hdt, node, value,
                           output_kind_of(value.current));
      }
      result_external_ports = a_external_ports;
      node.operation = HdtOperation::IG;
      node.left = a_hdt_root;
      node.right = b_hdt_root;
      ++hdt.stats.ig_ops;
      mutable_unlink_child(tree, group.a, b_node.parent_link);
    } else {
      if (b_node.parent != group.a) {
        throw std::logic_error("mutable CC operands are not parent/child");
      }
      const MutableActiveLink& shared =
          tree.links[static_cast<std::size_t>(b_node.parent_link)];
      const bool a_higher =
          shared.parent_port.orientation == HdtPortOrientation::Down &&
          shared.child_port.orientation == HdtPortOrientation::Up;
      if (!a_higher) {
        throw std::logic_error("mutable CC ordinary-port orientation mismatch");
      }
      // CC joins a value open at its bottom (C3 or C2) above one open at its
      // top (C1 or C2); the counter kernels read the operands in those
      // orientations.
      const HdtOutputKind higher_kind =
          hdt.nodes[static_cast<std::size_t>(a_hdt_root)].output_kind;
      const HdtOutputKind lower_kind =
          hdt.nodes[static_cast<std::size_t>(b_hdt_root)].output_kind;
      if ((higher_kind != HdtOutputKind::C1D &&
           higher_kind != HdtOutputKind::C2) ||
          (lower_kind != HdtOutputKind::C1U &&
           lower_kind != HdtOutputKind::C2)) {
        throw std::logic_error("mutable CC operands have the wrong orientation");
      }
      result_kind = ActiveKind::C;
      const HdtOutputKind cc_output =
          hdt.nodes[static_cast<std::size_t>(a_hdt_root)].output_kind ==
                      HdtOutputKind::C1D &&
                  hdt.nodes[static_cast<std::size_t>(b_hdt_root)].output_kind ==
                      HdtOutputKind::C1U
              ? HdtOutputKind::Closed
              : (hdt.nodes[static_cast<std::size_t>(a_hdt_root)].output_kind ==
                             HdtOutputKind::C1D
                     ? HdtOutputKind::C1D
                     : (hdt.nodes[static_cast<std::size_t>(b_hdt_root)]
                                    .output_kind ==
                                HdtOutputKind::C1U
                            ? HdtOutputKind::C1U
                            : HdtOutputKind::C2));
      const CounterInterface higher_interface = descriptor_for_kind(
          hdt.nodes[static_cast<std::size_t>(a_hdt_root)].output_kind).first;
      const CounterInterface lower_interface = descriptor_for_kind(
          hdt.nodes[static_cast<std::size_t>(b_hdt_root)].output_kind).first;
      evaluate_binary_payload_color_mask(
          hdt, node, cc_output, a_hdt_root, b_hdt_root,
          [&](std::uint8_t higher_mask, std::uint8_t lower_mask,
              const double* higher, const double* lower, double* output) {
            counters_cc_color_mask_into(
                hdt.palette, higher_interface, lower_interface,
                higher_mask, lower_mask, higher, lower, output);
          });
      result_external_ports = a_external_ports;
      for (int i = 0; i < b_external_ports.count; ++i) {
        append_external_port(
            result_external_ports,
            b_external_ports.values[static_cast<std::size_t>(i)]);
      }
      result_source_center = a_source_center;
      node.operation = HdtOperation::CC;
      node.left = a_hdt_root;
      node.right = b_hdt_root;
      ++hdt.stats.cc_ops;

      const int removed_link = b_node.parent_link;
      MutableActiveLink& removed =
          tree.links[static_cast<std::size_t>(removed_link)];
      const int previous = removed.previous;
      const int after = removed.next;
      const int adopted = b_node.first_child;
      if (adopted < 0) {
        if (previous < 0) a_node.first_child = after;
        else tree.links[static_cast<std::size_t>(previous)].next = after;
        if (after < 0) a_node.last_child = previous;
        else tree.links[static_cast<std::size_t>(after)].previous = previous;
      } else {
        if (previous < 0) a_node.first_child = adopted;
        else tree.links[static_cast<std::size_t>(previous)].next = adopted;
        tree.links[static_cast<std::size_t>(adopted)].previous = previous;
        const int tail = b_node.last_child;
        for (int link = adopted; link >= 0;
             link = tree.links[static_cast<std::size_t>(link)].next) {
          const int child = tree.links[static_cast<std::size_t>(link)].child;
          tree.nodes[static_cast<std::size_t>(child)].parent = group.a;
        }
        tree.links[static_cast<std::size_t>(tail)].next = after;
        if (after < 0) a_node.last_child = tail;
        else tree.links[static_cast<std::size_t>(after)].previous = tail;
      }
      removed.previous = -1;
      removed.next = -1;
      b_node.first_child = -1;
      b_node.last_child = -1;
    }

    hdt.nodes[static_cast<std::size_t>(node.left)].parent = result_hdt_root;
    hdt.nodes[static_cast<std::size_t>(node.right)].parent = result_hdt_root;
    note_node_shape(hdt, node);
    a_node.hdt_root = result_hdt_root;
    a_node.source_center = result_source_center;
    a_node.kind = result_kind;
    tree.external_ports[static_cast<std::size_t>(group.a)] =
        result_external_ports;
    b_node.hdt_root = -1;
    b_node.parent = -1;
    b_node.parent_link = -1;
    --tree.active_nodes;
    set_mutable_hdt_output(hdt, tree, group.a);
  }

  const long long removed = static_cast<long long>(compositions.size());
  const long long before = tree.active_nodes + removed;
  if (9LL * removed < before) {
    throw std::logic_error(
        "mutable SODA round violated the 1-in-9 shrink bound");
  }
}

// Builds the recolouring metadata: dirty flags and the owner index, a
// counting sort of the live-leaf atoms by owner slot.
void finalize_hdt_metadata(Hdt& hdt,
                           int owner_count,
                           HdtOwnerIndexView owner_index) {
  const bool by_position = owner_index.owner_position != nullptr;
  if (by_position &&
      (owner_index.position_begin < 0 ||
       owner_index.position_end < owner_index.position_begin)) {
    throw std::invalid_argument("invalid owner-position index");
  }
  const int slot_count = by_position
                             ? owner_index.position_end -
                                   owner_index.position_begin
                             : owner_count;
  if (slot_count < 0) {
    throw std::invalid_argument("negative owner index domain");
  }
  hdt.owner_index_by_position = by_position;
  hdt.owner_position_begin =
      by_position ? owner_index.position_begin : 0;
  auto owner_slot = [&](int owner) {
    if (!by_position) return owner;
    if (owner < 0 ||
        static_cast<std::size_t>(owner) >=
            owner_index.owner_position->size()) {
      return -1;
    }
    const int position = (*owner_index.owner_position)[
        static_cast<std::size_t>(owner)];
    return position >= owner_index.position_begin &&
                   position < owner_index.position_end
               ? position - owner_index.position_begin
               : -1;
  };
  hdt.dirty.assign(hdt.nodes.size(), 0);
  hdt.recompute_stack.clear();
  hdt.recompute_stack.reserve(static_cast<std::size_t>(checked_hdt_int_mul_add(
      hdt.stats.max_update_height, 2, 4,
      "recompute stack capacity overflow")));
  hdt.owner_offsets.assign(
      checked_hdt_size_add(static_cast<std::size_t>(slot_count), 1u,
                           "owner-offset capacity overflow"),
      0);
  const int atom_count = checked_hdt_int_size(
      hdt.atoms.size(), "HDT exceeds its atom index domain");
  for (int node = 0; node < atom_count; ++node) {
    if (hdt.nodes[static_cast<std::size_t>(node)].operation !=
            HdtOperation::Atom ||
        hdt.nodes[static_cast<std::size_t>(node)].atom_descriptor != node) {
      throw std::logic_error("atom/node prefix invariant failed");
    }
    const int owner = hdt.nodes[static_cast<std::size_t>(node)].live_owner;
    const int slot = owner_slot(owner);
    if (slot >= 0 && slot < slot_count) {
      ++hdt.owner_offsets[static_cast<std::size_t>(slot) + 1u];
    }
  }
  for (int slot = 0; slot < slot_count; ++slot) {
    hdt.owner_offsets[static_cast<std::size_t>(slot) + 1u] +=
        hdt.owner_offsets[static_cast<std::size_t>(slot)];
  }
  hdt.owner_atom_nodes.assign(
      static_cast<std::size_t>(hdt.owner_offsets.back()), -1);
  for (int node = atom_count - 1; node >= 0; --node) {
    const int owner = hdt.nodes[static_cast<std::size_t>(node)].live_owner;
    const int slot = owner_slot(owner);
    if (slot < 0 || slot >= slot_count) continue;
    const int position = --hdt.owner_offsets[
        static_cast<std::size_t>(slot) + 1u];
    hdt.owner_atom_nodes[static_cast<std::size_t>(position)] = node;
  }
  // owner_offsets[s + 1] now holds the begin of slot s; shift into place.
  for (int slot = 1; slot < slot_count; ++slot) {
    hdt.owner_offsets[static_cast<std::size_t>(slot)] =
        hdt.owner_offsets[static_cast<std::size_t>(slot) + 1u];
  }
  hdt.owner_offsets[static_cast<std::size_t>(slot_count)] =
      checked_hdt_int_size(hdt.owner_atom_nodes.size(),
                           "owner index exceeds integer domain");
}

void set_error(std::string* error, const std::string& message) {
  if (error != nullptr) *error = message;
}

}  // namespace

static void build_hdt_impl(Hdt& hdt,
                              const HdtComponentGraph& graph,
                              bool require_closed_root,
                              bool materialize_update_metadata,
                              bool background_alias_current,
                              HdtOwnerIndexView owner_index) {
  const long long scalar_fmas_before = counters_thread_scalar_fmas();
#ifndef NDEBUG
  std::string graph_error;
  if (!validate_component_graph(graph, &graph_error)) {
    throw std::invalid_argument("invalid component graph: " + graph_error);
  }
#endif
  hdt.root = -1;
  hdt.root_anchor_atom_node = -1;
  hdt.nodes.clear();
  hdt.current_payload.clear();
  hdt.background_payload.clear();
  hdt.atoms.clear();
  hdt.atom_ports.clear();
  hdt.owner_offsets.clear();
  hdt.owner_atom_nodes.clear();
  hdt.recompute_stack.clear();
  hdt.owner_position_begin = 0;
  hdt.owner_index_by_position = false;
  hdt.dirty.clear();
  hdt.stats = HdtBuildStats{};
  hdt.palette = graph.palette;
  hdt.closed_root_required = require_closed_root;
  hdt.background_alias_current = background_alias_current;
  hdt.background_restricted = graph.background_restricted;
  hdt.stats.component_units = graph.component_units();
  if (graph.units.empty()) {
    // A gene without a live leaf gives the empty graph, A = H = R = 0. Both
    // structural predicates hold vacuously; H = 2A - 1 applies only to a
    // nonempty graph.
    hdt.stats.node_count_ok = true;
    hdt.stats.atom_bound_ok = true;
    hdt.stats.scalar_fmas =
        counters_thread_scalar_fmas() - scalar_fmas_before;
    return;
  }

  const std::size_t component_units = graph.units.size();
  std::size_t initial_atoms = component_units;
  for (const HdtComponentUnit& unit : graph.units) {
    if (unit.kind == HdtComponentKind::IAnchor &&
        unit.local_background_g != nullptr) {
      initial_atoms = checked_hdt_size_add(
          initial_atoms, 1u, "HDT atom count overflow");
    }
  }
  const int initial_atom_count = checked_hdt_int_size(
      initial_atoms, "HDT exceeds its atom index domain");
  const int node_capacity = checked_hdt_int_mul_add(
      initial_atom_count, 2, -1,
      "HDT exceeds its node index domain");
  const std::size_t local_background_count = initial_atoms - component_units;
  const std::size_t atom_port_capacity = checked_hdt_size_add(
      checked_hdt_size_add(
          checked_hdt_size_mul(graph.edges.size(), 2u,
                               "HDT port capacity overflow"),
          graph.external_ports.size(), "HDT port capacity overflow"),
      local_background_count, "HDT port capacity overflow");
  (void)checked_hdt_int_size(
      atom_port_capacity, "HDT exceeds its port index domain");
  // The builder emits exactly A atoms and A-1 binary nodes. Payload arenas
  // grow independently by exact active-coordinate width.
  hdt.atoms.reserve(initial_atoms);
  hdt.nodes.reserve(static_cast<std::size_t>(node_capacity));
  hdt.atom_ports.reserve(atom_port_capacity);

  std::size_t build_scratch_peak_bytes = 0;
  RoundWorkspace round_workspace;
  MutableActiveTree mutable_tree = build_mutable_active_tree(
      graph, initial_atoms, hdt, round_workspace, build_scratch_peak_bytes);
  hdt.stats.initial_atoms = mutable_tree.active_nodes;
  if (hdt.stats.initial_atoms != initial_atom_count) {
    throw std::logic_error("HDT atom count changed during expansion");
  }

  // The round bound is a guard; the one-ninth shrink gives O(log A) rounds.
  int round = 0;
  while (mutable_tree.active_nodes > 1) {
    ++round;
    if (static_cast<long long>(round) >
        4LL * static_cast<long long>(hdt.stats.initial_atoms) + 16LL) {
      throw std::logic_error("SODA builder exceeded round bound");
    }
    execute_mutable_round(mutable_tree, hdt, round, round_workspace);
  }
  build_scratch_peak_bytes = std::max(
      build_scratch_peak_bytes,
      live_builder_scratch_bytes(mutable_tree, round_workspace));
  hdt.stats.build_scratch_peak_bytes = build_scratch_peak_bytes;
  hdt.stats.rounds = round;
  hdt.root = mutable_tree.nodes[static_cast<std::size_t>(mutable_tree.root)]
                 .hdt_root;
  HdtNode& root = hdt.nodes[static_cast<std::size_t>(hdt.root)];
  if (require_closed_root && root.output_kind == HdtOutputKind::C1U) {
    root.post_transform = HdtTransform::CloseRoot;
    replace_unary_payload_color_mask(
        hdt, root, HdtOutputKind::Closed,
        [&](std::uint8_t mask, const double* input, double* output) {
          counters_close_color_mask_into(
              hdt.palette, mask, input, output);
        });
    root.output_kind = HdtOutputKind::Closed;
    root.output_port_count = 0;
  }
  if (require_closed_root &&
      root.output_kind != HdtOutputKind::Closed) {
    throw std::logic_error("HDT root is not Closed");
  }
  if (materialize_update_metadata) {
    finalize_hdt_metadata(hdt, graph.owner_count, owner_index);
  }

  hdt.stats.reachable_nodes = checked_hdt_int_size(
      hdt.nodes.size(), "HDT exceeds its node index domain");
  hdt.stats.node_count_ok =
      static_cast<long long>(hdt.stats.reachable_nodes) ==
      2LL * static_cast<long long>(hdt.stats.initial_atoms) - 1LL;
  hdt.stats.atom_bound_ok =
      static_cast<long long>(hdt.stats.initial_atoms) <=
      2LL * static_cast<long long>(hdt.stats.component_units);
  if (!hdt.stats.node_count_ok || !hdt.stats.atom_bound_ok ||
      hdt.stats.max_recompute_degree > 2 ||
      hdt.stats.local_balance_violations != 0) {
    throw std::logic_error(
        "HDT construction violated a structural invariant");
  }
  if (materialize_update_metadata) {
#ifndef NDEBUG
    std::string error;
    if (!validate_hdt(graph, hdt, &error)) {
      throw std::logic_error("invalid HDT: " + error);
    }
#endif
  }
  hdt.stats.scalar_fmas =
      counters_thread_scalar_fmas() - scalar_fmas_before;
}

void build_hdt_into(Hdt& hdt,
                              const HdtComponentGraph& graph,
                              bool require_closed_root,
                              HdtOwnerIndexView owner_index) {
  build_hdt_impl(hdt, graph, require_closed_root, true, false,
                           owner_index);
}

void build_frozen_fold_hdt_into(
    Hdt& hdt,
    const HdtComponentGraph& graph,
    bool require_closed_root) {
  build_hdt_impl(hdt, graph, require_closed_root, false, true, {});
}

Hdt build_hdt(const HdtComponentGraph& graph,
                                 bool require_closed_root,
                                 HdtOwnerIndexView owner_index) {
  Hdt hdt;
  build_hdt_into(hdt, graph, require_closed_root, owner_index);
  return hdt;
}

HdtState build_hdt_state(
    HdtComponentGraph graph,
    HdtOwnerIndexView owner_index) {
  HdtState state;
  state.components = std::move(graph);
  state.hdt = build_hdt(state.components, true, owner_index);
  return state;
}

namespace {
int flush_pending_recompute(Hdt& hdt);
}  // namespace

CounterVec hdt_node_background(const Hdt& hdt, int node) {
  flush_pending_recompute(const_cast<Hdt&>(hdt));
  return hdt_node_background_stored(hdt, node);
}

CounterVec hdt_node_background_stored(const Hdt& hdt, int node) {
  if (node < 0 || node >= static_cast<int>(hdt.nodes.size())) {
    throw std::out_of_range(
        "HDT background payload node out of range");
  }
  const HdtNode& value =
      hdt.nodes[static_cast<std::size_t>(node)];
  return materialize_view(node_background_view(hdt, value), hdt.palette,
                          value.background_color_mask);
}

CounterPair hdt_node_dp(const Hdt& hdt, int node) {
  flush_pending_recompute(const_cast<Hdt&>(hdt));
  return hdt_node_dp_stored(hdt, node);
}

CounterPair hdt_node_dp_stored(const Hdt& hdt, int node) {
  if (node < 0 || node >= static_cast<int>(hdt.nodes.size())) {
    throw std::out_of_range("HDT dual payload node out of range");
  }
  const HdtNode& value = hdt.nodes[static_cast<std::size_t>(node)];
  if (value.output_kind == HdtOutputKind::I) return CounterPair{};
  return CounterPair{materialize_view(node_current_view(hdt, value),
                                     hdt.palette,
                                     value.current_color_mask),
                     materialize_view(node_background_view(hdt, value),
                                     hdt.palette,
                                     value.background_color_mask),
                     value.current_equals_background};
}

bool validate_hdt(const HdtComponentGraph& graph,
                            const Hdt& hdt,
                            std::string* error) {
  flush_pending_recompute(const_cast<Hdt&>(hdt));
  if (graph.units.empty()) return hdt.root == -1 && hdt.nodes.empty();
  if (hdt.root < 0 || hdt.root >= static_cast<int>(hdt.nodes.size())) {
    set_error(error, "invalid root");
    return false;
  }
  if (!hdt.stats.node_count_ok) {
    set_error(error, "H != 2A-1");
    return false;
  }
  if (!hdt.stats.atom_bound_ok) {
    set_error(error, "A > 2R");
    return false;
  }
  if (hdt.stats.max_recompute_degree > 2) {
    set_error(error, "nonbinary recompute node");
    return false;
  }
  if (hdt.stats.local_balance_violations != 0) {
    set_error(error, "local balance violation");
    return false;
  }
  if (hdt.root_anchor_atom_node >= 0) {
    if (hdt.root_anchor_atom_node >= static_cast<int>(hdt.nodes.size()) ||
        hdt.nodes[static_cast<std::size_t>(hdt.root_anchor_atom_node)]
                .operation != HdtOperation::Atom) {
      set_error(error, "invalid stored root-anchor atom");
      return false;
    }
    const HdtExpandedAtom& atom = hdt.atoms[static_cast<std::size_t>(
        hdt.nodes[static_cast<std::size_t>(hdt.root_anchor_atom_node)]
            .atom_descriptor)];
    if (atom.component_kind != HdtComponentKind::IAnchor ||
        atom.component_unit != graph.root_unit) {
      set_error(error, "stored root-anchor atom does not name graph root");
      return false;
    }
  }
  const HdtNode& root =
      hdt.nodes[static_cast<std::size_t>(hdt.root)];
  if (hdt.closed_root_required) {
    if (root.output_kind != HdtOutputKind::Closed) {
      set_error(error, "root output is not Closed");
      return false;
    }
  } else if (root.output_kind == HdtOutputKind::I ||
             root.output_port_count != graph.external_ports.size()) {
    set_error(error, "open root output does not match external ports");
    return false;
  }
  if (root.parent != -1) {
    set_error(error, "root has a parent");
    return false;
  }
  std::vector<std::uint8_t> seen(hdt.nodes.size(), 0);
  std::vector<int> stack{hdt.root};
  int reached = 0;
  while (!stack.empty()) {
    const int node = stack.back();
    stack.pop_back();
    if (node < 0 || node >= static_cast<int>(hdt.nodes.size()) ||
        seen[static_cast<std::size_t>(node)] != 0) {
      set_error(error, "node reachability is not one-to-one");
      return false;
    }
    seen[static_cast<std::size_t>(node)] = 1;
    ++reached;
    const HdtNode& value = hdt.nodes[static_cast<std::size_t>(node)];
    if (value.output_kind == HdtOutputKind::I) {
      if (value.dp_offset != -1 || value.background_dp_offset != -1 ||
          value.dp_dimension != 0 || value.current_payload_class != 0xff ||
          value.background_payload_class != 0xff) {
        set_error(error, "I node owns a packed payload");
        return false;
      }
    } else {
      const auto [interface, orientation] =
          descriptor_for_kind(value.output_kind);
      if (value.dp_dimension != counters_dimension(hdt.palette, interface)) {
        set_error(error, "packed payload dimension mismatch");
        return false;
      }
      const std::uint8_t expected_current_class = packed_payload_class_index(
          packed_payload_width(hdt.palette, interface, orientation,
                               value.current_color_mask));
      if (value.current_payload_class != expected_current_class ||
          value.dp_offset < 0) {
        set_error(error, "packed current payload class mismatch");
        return false;
      }
      const int current_width = packed_payload_width(expected_current_class);
      const auto& current_values = hdt.current_payload.values;
      if (static_cast<std::size_t>(value.dp_offset) + current_width >
          current_values.size()) {
        set_error(error, "packed current payload offset mismatch");
        return false;
      }
      if (value.background_dp_offset < 0) {
        if (value.background_payload_class != 0xff) {
          set_error(error, "aliased background has a payload class");
          return false;
        }
      } else {
        const std::uint8_t expected_background_class =
            packed_payload_class_index(packed_payload_width(
                hdt.palette, interface, orientation,
                value.background_color_mask));
        if (hdt.background_alias_current ||
            value.background_payload_class != expected_background_class) {
          set_error(error, "packed background payload class mismatch");
          return false;
        }
        const int background_width =
            packed_payload_width(expected_background_class);
        const auto& background_values = hdt.background_payload.values;
        if (static_cast<std::size_t>(value.background_dp_offset) +
                background_width >
            background_values.size()) {
          set_error(error, "packed background payload offset mismatch");
          return false;
        }
      }
    }
    if (value.operation == HdtOperation::Atom) {
      if (value.left >= 0 || value.right >= 0 || value.atom_descriptor < 0 ||
          value.atom_descriptor >= static_cast<int>(hdt.atoms.size())) {
        set_error(error, "invalid atom node");
        return false;
      }
    } else {
      if (value.left < 0 || value.right < 0 ||
          value.left >= static_cast<int>(hdt.nodes.size()) ||
          value.right >= static_cast<int>(hdt.nodes.size()) ||
          hdt.nodes[static_cast<std::size_t>(value.left)].parent != node ||
          hdt.nodes[static_cast<std::size_t>(value.right)].parent != node) {
        set_error(error, "invalid binary parent link");
        return false;
      }
      stack.push_back(value.right);
      stack.push_back(value.left);
    }
  }
  if (reached != static_cast<int>(hdt.nodes.size())) {
    set_error(error, "HDT contains unreachable nodes");
    return false;
  }
  return true;
}

int hdt_recolor_owner_ranges(
    HdtState& state,
    const std::vector<int>& owner_order,
    std::span<const std::pair<int, int>> owner_ranges,
    std::uint8_t color) {
  const long long scalar_fmas_before = counters_thread_scalar_fmas();
  Hdt& hdt = state.hdt;
  if (color >= hdt.palette) {
    throw std::invalid_argument("recolor outside palette");
  }
  const int owner_order_count = checked_hdt_int_size(
      owner_order.size(), "recolor exceeds integer owner domain");
  const int component_count = state.components.component_units();
  const int owner_offset_count = checked_hdt_int_size(
      hdt.owner_offsets.size(),
      "owner index exceeds its integer domain");
  for (const auto& [begin_owner, end_owner] : owner_ranges) {
    if (begin_owner < 0 || end_owner < begin_owner ||
        end_owner > owner_order_count) {
      throw std::invalid_argument("invalid recolor owner range");
    }
  }
  auto recolor_component_unit = [&](int unit_index) {
    if (unit_index < 0 || unit_index >= component_count) {
      throw std::logic_error("live HDT atom lacks its component unit");
    }
    HdtComponentUnit& unit =
        state.components.units[static_cast<std::size_t>(unit_index)];
    if (unit.kind != HdtComponentKind::LiveLeaf) {
      throw std::logic_error("mutable HDT atom is not a live component");
    }
    unit.current_color = color;
    unit.current_color_mask = unit.weight == 0.0
                                  ? 0
                                  : static_cast<std::uint8_t>(1u << color);
    unit.current_equals_background =
        color == unit.background_color || unit.weight == 0.0;
  };
  int changed = 0;
  auto recolor_slot = [&](int slot) {
    if (slot < 0 || slot >= owner_offset_count - 1) {
      return;
    }
    const int begin = hdt.owner_offsets[static_cast<std::size_t>(slot)];
    const int end = hdt.owner_offsets[static_cast<std::size_t>(slot) + 1u];
    for (int p = begin; p < end; ++p) {
      const int node_index =
          hdt.owner_atom_nodes[static_cast<std::size_t>(p)];
      HdtNode& node =
          hdt.nodes[static_cast<std::size_t>(node_index)];
      HdtExpandedAtom& atom = hdt.atoms[static_cast<std::size_t>(
          node.atom_descriptor)];
      if (atom.current_color == color) continue;
      atom.current_color = color;
      atom.current_color_mask = atom.weight == 0.0
                                    ? 0
                                    : static_cast<std::uint8_t>(1u << color);
      recolor_component_unit(atom.component_unit);
      ++changed;
      // Mark the path to the root; above a dirty node it is already marked.
      for (int mark = node_index; mark >= 0 &&
                                      hdt.dirty[static_cast<std::size_t>(mark)] == 0;
           mark = hdt.nodes[static_cast<std::size_t>(mark)].parent) {
        hdt.dirty[static_cast<std::size_t>(mark)] = 1;
      }
    }
  };

  if (hdt.owner_index_by_position) {
    if (owner_offset_count == 0) {
      throw std::logic_error("owner-position index has no offsets");
    }
    const int indexed_end = checked_hdt_int_add(
        hdt.owner_position_begin, owner_offset_count - 1,
        "owner-position interval overflow");
    for (const auto& [begin_owner, end_owner] : owner_ranges) {
      const int begin_position =
          std::max(begin_owner, hdt.owner_position_begin);
      const int end_position = std::min(end_owner, indexed_end);
      for (int position = begin_position; position < end_position;
           ++position) {
        recolor_slot(position - hdt.owner_position_begin);
      }
    }
  } else {
    for (const auto& [begin_owner, end_owner] : owner_ranges) {
      for (int owner_index = begin_owner; owner_index < end_owner;
           ++owner_index) {
        recolor_slot(owner_order[static_cast<std::size_t>(owner_index)]);
      }
    }
  }
  if (changed != 0) hdt.recompute_pending = true;
  state.recolor_changed_leaves += changed;
  state.scalar_fmas +=
      counters_thread_scalar_fmas() - scalar_fmas_before;
  return changed;
}

namespace {

// Recomputes the dirty nodes once each, in post-order. They form a union of
// paths to the root, so the walk descends only into dirty children; a token
// ~node schedules a node after its children. Returns the recomputed count.
int flush_pending_recompute(Hdt& hdt) {
  if (!hdt.recompute_pending) return 0;
  hdt.recompute_pending = false;
  int recomputed = 0;
  std::vector<int>& stack = hdt.recompute_stack;
  stack.clear();
  stack.push_back(hdt.root);
  while (!stack.empty()) {
    const int token = stack.back();
    stack.pop_back();
    if (token < 0) {
      const int node = ~token;
      recompute_node(hdt, node);
      hdt.dirty[static_cast<std::size_t>(node)] = 0;
      ++recomputed;
      continue;
    }
    const int node = token;
    if (node < 0 || hdt.dirty[static_cast<std::size_t>(node)] == 0) continue;
    stack.push_back(~node);
    const HdtNode& value =
        hdt.nodes[static_cast<std::size_t>(node)];
    if (value.right >= 0 &&
        hdt.dirty[static_cast<std::size_t>(value.right)] != 0) {
      stack.push_back(value.right);
    }
    if (value.left >= 0 &&
        hdt.dirty[static_cast<std::size_t>(value.left)] != 0) {
      stack.push_back(value.left);
    }
  }
  return recomputed;
}

}  // namespace

void hdt_flush_recompute(HdtState& state) {
  if (!state.hdt.recompute_pending) return;
  const long long scalar_fmas_before = counters_thread_scalar_fmas();
  const int recomputed = flush_pending_recompute(state.hdt);
  state.recolor_marked_nodes += recomputed;
  state.transition_evals += recomputed;
  state.scalar_fmas +=
      counters_thread_scalar_fmas() - scalar_fmas_before;
}

int hdt_recolor_owner_range(HdtState& state,
                                  const std::vector<int>& owner_order,
                                  int begin_owner,
                                  int end_owner,
                                  std::uint8_t color) {
  const std::pair<int, int> range{begin_owner, end_owner};
  return hdt_recolor_owner_ranges(
      state, owner_order, std::span<const std::pair<int, int>>(&range, 1),
      color);
}

int hdt_recolor_owners(HdtState& state,
                             const std::vector<int>& owners,
                             std::uint8_t color) {
  return hdt_recolor_owner_range(
      state, owners, 0,
      checked_hdt_int_size(owners.size(),
                           "recolor exceeds integer owner domain"),
      color);
}

int hdt_recolor_owner(HdtState& state,
                            int owner,
                            std::uint8_t color) {
  return hdt_recolor_owners(state, std::vector<int>{owner}, color);
}

std::array<double, 3> hdt_root_counts(
    const HdtState& state) {
  if (state.hdt.root < 0) return {0.0, 0.0, 0.0};
  if (state.hdt.recompute_pending) {
    // A pending recompute exists only on a state that was recolored, which
    // is never a const object.
    hdt_flush_recompute(const_cast<HdtState&>(state));
  }
  const HdtNode& root =
      state.hdt.nodes[static_cast<std::size_t>(state.hdt.root)];
  return counters_color_mask_readout(
      state.hdt.palette, root.current_color_mask,
      node_current_data(state.hdt, root));
}

int hdt_live_atom_count(const HdtState& state) {
  return state.hdt.owner_offsets.empty() ? 0
                                         : state.hdt.owner_offsets.back();
}

}  // namespace qraft::internal
