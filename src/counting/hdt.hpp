#pragma once

#include "counting/units.hpp"

#include <array>
#include <cstdint>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace qraft::internal {

// How an HDT node combines its children. GG merges two G components at one
// centre; IG absorbs the single G child of an I (left) into a C or Closed
// component; CC composes a C (left) with the C below its Down port (right).
enum class HdtOperation : std::uint8_t {
  Atom,
  GG,
  IG,
  CC,
};

// Unary transform fused onto a node's output. CToG turns a C without a Down
// port whose parent is an I into a G at the I's centre; CloseRoot closes a
// root C1 with an Up port.
enum class HdtTransform : std::uint8_t {
  None,
  CToG,
  CloseRoot,
};

// Counter interface of a node's output. I has no counters; C1U and C1D are
// C1 with an Up or a Down port (the draft's C1 and C3).
enum class HdtOutputKind : std::uint8_t {
  I,
  Closed,
  C1U,
  C1D,
  C2,
  G,
};

// An HDT leaf: a component unit, or the local background G of an anchor
// (component_kind GOpaque, component_unit the anchor). Its ports are
// Hdt::atom_ports[port_begin, +port_count).
struct HdtExpandedAtom {
  HdtComponentKind component_kind = HdtComponentKind::IAnchor;
  int component_unit = -1;
  int owner = -1;
  int source_vertex_id = -1;
  int port_begin = 0;
  int port_count = 0;
  double weight = 0.0;
  std::uint8_t current_color = 0;
  std::uint8_t background_color = 0;
  std::uint8_t current_color_mask = 0;
  std::uint8_t background_color_mask = 0;
};

struct HdtNode {
  int left = -1;
  int right = -1;
  int parent = -1;
  int atom_descriptor = -1;  // index into Hdt::atoms for an atom node
  // Open ports of the component the node represents.
  std::array<HdtPortKey, 2> output_ports{};
  // Gene vertex keying the component; a G component's centre.
  int source_center_key = -1;
  // Slots of the packed counters in the current and background arenas. With
  // background_dp_offset -1 the background equals the current value and is
  // read from the current slot.
  int dp_offset = -1;
  int background_dp_offset = -1;
  int atom_count = 1;  // atoms below the node
  int height = 0;
  int live_owner = -1;  // owner of a live-leaf atom, else -1
  HdtOperation operation = HdtOperation::Atom;
  HdtTransform post_transform = HdtTransform::None;
  HdtOutputKind output_kind = HdtOutputKind::I;
  std::uint8_t output_port_count = 0;
  std::uint8_t dp_dimension = 0;  // full counter dimension; 0 for I
  bool current_equals_background = false;
  // Colours with nonzero mass below the node; they select the packed
  // coordinates and so the width class of each slot.
  std::uint8_t current_color_mask = 0;
  std::uint8_t background_color_mask = 0;
  // Width classes of the two slots; 0xff when the slot is absent.
  std::uint8_t current_payload_class = 0xff;
  std::uint8_t background_payload_class = 0xff;
};

// Structural record of one build over R component units, A atoms (units plus
// local background G atoms) and H nodes. The build requires H = 2A - 1,
// A <= 2R, recompute degree at most two and no local balance violation.
struct HdtBuildStats {
  int component_units = 0;
  int initial_atoms = 0;
  int reachable_nodes = 0;
  int max_recompute_degree = 0;
  int local_balance_violations = 0;
  int max_update_height = 0;
  int rounds = 0;
  long long gg_ops = 0;
  long long ig_ops = 0;
  long long cc_ops = 0;
  long long c_to_g_ops = 0;
  long long scalar_fmas = 0;
  std::size_t build_scratch_peak_bytes = 0;
  bool node_count_ok = false;
  bool atom_bound_ok = false;
};

// Optional restriction of the owner index. Without it each owner id is a
// slot. With it only owners whose owner_position lies in
// [position_begin, position_end) are indexed, at slot
// position - position_begin.
struct HdtOwnerIndexView {
  const std::vector<int>* owner_position = nullptr;
  int position_begin = 0;
  int position_end = 0;
};

struct HdtPackedPayloadArena {
  // One slab holds every active-coordinate width. A released slot keeps its
  // width class and is reused only by that class, through a per-class free
  // list threaded through the slot's first value.
  static constexpr std::size_t kClassCount = 15;
  std::vector<double> values;
  std::array<int, kClassCount> free_head = make_empty_free_heads();

  static constexpr std::array<int, kClassCount> make_empty_free_heads() {
    std::array<int, kClassCount> heads{};
    for (int& head : heads) head = -1;
    return heads;
  }

  void clear() {
    values.clear();
    free_head = make_empty_free_heads();
  }

  std::size_t capacity_bytes() const {
    std::size_t bytes = 0;
    checked_hdt_accumulate_bytes(
        bytes, values.capacity(), sizeof(double),
        "packed payload byte count exceeds size_t");
    return bytes;
  }
};

struct Hdt {
  int palette = 0;
  int root = -1;
  // Atom node of the graph's root unit when that unit is an I anchor.
  int root_anchor_atom_node = -1;
  // Atom nodes form the prefix [0, atoms.size()), node i holding atom i.
  std::vector<HdtNode> nodes;
  HdtPackedPayloadArena current_payload;
  HdtPackedPayloadArena background_payload;
  std::vector<HdtExpandedAtom> atoms;
  std::vector<HdtPortKey> atom_ports;
  // Owner index: slot s lists its live-leaf atom nodes in increasing order in
  // owner_atom_nodes[owner_offsets[s], owner_offsets[s + 1]). Every duplicate
  // leaf of an owner has its own entry. Slots are owner ids, or positions
  // minus owner_position_begin when owner_index_by_position.
  std::vector<int> owner_offsets;
  std::vector<int> owner_atom_nodes;
  std::vector<int> recompute_stack;
  int owner_position_begin = 0;
  bool owner_index_by_position = false;
  std::vector<std::uint8_t> dirty;
  // Recolors mark dirty paths only; the dirty nodes are recomputed once, in
  // post-order, at the next read of a current value (hdt_flush_recompute).
  bool recompute_pending = false;
  HdtBuildStats stats;
  bool closed_root_required = true;
  // A frozen fold: every background equals its current value and is read
  // from the current arena. It has no owner index or dirty flags.
  bool background_alias_current = false;
  bool background_restricted = false;
};

// A component graph with its HDT. Recolouring updates both, so the graph's
// units always carry the current colours.
struct HdtState {
  HdtState() = default;
  HdtState(const HdtState&) = delete;
  HdtState& operator=(const HdtState&) = delete;
  HdtState(HdtState&&) noexcept = default;
  HdtState& operator=(HdtState&&) noexcept = default;

  HdtComponentGraph components;
  Hdt hdt;
  // Work charged to this state: recoloured leaf atoms, recomputed nodes,
  // transition evaluations and scalar FMAs.
  long long recolor_changed_leaves = 0;
  long long recolor_marked_nodes = 0;
  long long transition_evals = 0;
  long long scalar_fmas = 0;
  // The initial HDT's cost is already in the run statistics.
  bool initial_hdt_accounted = false;
  // Components are valid but the HDT has not been built. Only a traversal
  // entry that immediately contracts or builds may receive such a state.
  bool hdt_deferred = false;
};

// Builds the HDT of a component graph by the rounds of Brodal et al.
// (SODA'13, Sec. 4). Each round removes at least a ninth of the components,
// and a node over a atoms has height O(log a). With require_closed_root a
// root C1 with an Up port is closed, and any other non-Closed root throws.
Hdt build_hdt(const HdtComponentGraph& graph,
              bool require_closed_root = true,
              HdtOwnerIndexView owner_index = {});

void build_hdt_into(Hdt& output,
                    const HdtComponentGraph& graph,
                    bool require_closed_root = true,
                    HdtOwnerIndexView owner_index = {});

// Builds a frozen fold (see Hdt::background_alias_current) for reading its
// root once. Every unit's current value must equal its background.
void build_frozen_fold_hdt_into(
    Hdt& output,
    const HdtComponentGraph& graph,
    bool require_closed_root = true);

HdtState build_hdt_state(HdtComponentGraph graph,
                         HdtOwnerIndexView owner_index = {});

// Checks the structural invariants of an HDT built from `graph`: the build
// predicates, the root's interface, parent links, reachability of every node
// and each slot's width class. Completes pending recomputation first.
bool validate_hdt(const HdtComponentGraph& graph,
                  const Hdt& hdt,
                  std::string* error = nullptr);

// Recolour every live-leaf atom of the selected owners and return how many
// atoms changed colour. The owners are ranges of positions in owner_order;
// the wrappers select all of `owners`, or the single `owner`. The changed
// atoms' paths to the root are marked dirty and recomputed at the next read.
// In a position-indexed HDT position p is slot p - owner_position_begin, so
// owner_order must be the order that owner_position indexes; positions
// outside the index are skipped.
int hdt_recolor_owners(HdtState& state,
                       const std::vector<int>& owners,
                       std::uint8_t color);

int hdt_recolor_owner_range(HdtState& state,
                            const std::vector<int>& owner_order,
                            int begin,
                            int end,
                            std::uint8_t color);

int hdt_recolor_owner_ranges(
    HdtState& state,
    const std::vector<int>& owner_order,
    std::span<const std::pair<int, int>> owner_ranges,
    std::uint8_t color);

int hdt_recolor_owner(HdtState& state,
                      int owner,
                      std::uint8_t color);

// Channels of the root's current counters; zeros for an empty HDT. Completes
// pending recomputation first.
std::array<double, 3> hdt_root_counts(const HdtState& state);

// These accessors first complete any pending recomputation.
CounterPair hdt_node_dp(const Hdt& hdt, int node);
CounterVec hdt_node_background(const Hdt& hdt, int node);

// Stored values without completing pending recomputation. Valid for values
// pending recolors cannot change: every atom that is not a live leaf, and
// every background (an unpromoted node's stored current equals it).
CounterPair hdt_node_dp_stored(const Hdt& hdt, int node);
CounterVec hdt_node_background_stored(const Hdt& hdt, int node);

// Recomputes every node marked dirty since the last flush, each once in
// post-order, and charges the work to the state's counters.
void hdt_flush_recompute(HdtState& state);

// Number of live-leaf atoms in the owner index.
int hdt_live_atom_count(const HdtState& state);

}  // namespace qraft::internal
