#pragma once

#include "counting/hdt.hpp"

#include <cstdint>
#include <vector>

namespace qraft::internal {

// LightCopy may return an unchanged rebuild when the input already meets the
// Contract bound. HeavyReplace (a state that replaces its input for a heavy
// child) and InitialCompact (removal of absent leaves before the first HDT)
// always contract. The mode is also recorded for the run statistics.
enum class HdtContractMode : std::uint8_t {
  LightCopy,
  HeavyReplace,
  InitialCompact,
};

// Work, size and timing of one Contract. size_bound is 4*max(1,s)+3 for s
// selected live atoms. persistent_scratch_bytes is held for the whole call,
// ephemeral_scratch_peak_bytes is the largest transient allocation on top of
// it, and scratch_peak_bytes is their sum.
struct HdtContractStats {
  HdtContractMode mode = HdtContractMode::LightCopy;
  int input_units = 0;
  int output_units = 0;
  int selected_live_atoms = 0;
  int anchors = 0;
  int corridors_emitted = 0;
  int background_units_absorbed = 0;
  int polytomy_arms_folded = 0;
  int size_bound = 0;
  bool size_bound_ok = false;
  long long fragment_hdt_nodes = 0;
  long long gg_ops = 0;
  long long ig_ops = 0;
  long long cc_ops = 0;
  long long c_to_g_ops = 0;
  long long mark_ns = 0;
  long long background_fold_ns = 0;
  long long corridor_fold_ns = 0;
  long long emit_ns = 0;
  long long hdt_rebuild_ns = 0;
  long long scalar_fmas = 0;
  std::size_t persistent_scratch_bytes = 0;
  std::size_t ephemeral_scratch_peak_bytes = 0;
  std::size_t scratch_peak_bytes = 0;
};

struct HdtContractResult {
  HdtState state;
  HdtContractStats stats;
};

// `contract` is meaningful only when `contracted` is true.
struct HdtInitialStateResult {
  HdtState state;
  HdtContractStats contract;
  bool contracted = false;
};

// intermediate_units counts the emitted graph before its Contract.
struct HdtExtractStats {
  int source_hdt_nodes = 0;
  int marked_hdt_nodes = 0;
  int hanging_components = 0;
  int intermediate_units = 0;
  int selected_live_atoms = 0;
};

struct HdtExtractContractResult {
  HdtState state;
  HdtExtractStats extract;
  HdtContractStats contract;
};

// Scratch reused by the Extract calls of one worker. An entry is valid only
// when its stamp equals `epoch`, so no call clears the arrays. The HDT marks
// are indexed by HDT node, the ordinary entries by source edge key and the
// center entries by source vertex. peak_bytes is the largest total capacity.
struct HdtExtractWorkspace {
  std::uint64_t epoch = 0;
  std::vector<std::uint64_t> hdt_mark_epoch;
  std::vector<std::uint64_t> ordinary_epoch;
  std::vector<int> ordinary_unit;
  std::vector<HdtPortKey> ordinary_port;
  std::vector<std::uint64_t> center_epoch;
  std::vector<int> center_anchor;
  std::size_t peak_bytes = 0;
};

// The owners whose live leaves a Contract keeps. The interval form, used by
// the traversals, borrows two arrays of an HdtQueryOwnerIndex: owners are
// iterated as owner_order[owner_begin, owner_end), and membership is O(1)
// through owner_position, the query leaf position, tested against
// [membership_begin, membership_end). Building it is O(1). The owning dense
// form serves standalone callers and tests. A non-negative
// exact_selected_atom_count is the number of selected live atoms, which
// spares a scan of the component graph.
struct HdtLiveSelector {
  const std::vector<int>* owner_order = nullptr;
  const std::vector<int>* owner_position = nullptr;
  int owner_begin = 0;
  int owner_end = 0;
  int membership_begin = 0;
  int membership_end = 0;
  int exact_selected_atom_count = -1;
  std::vector<int> owned_owners;
  std::vector<std::uint8_t> dense_membership;

  int selected_owner_count() const;
  int selected_owner_at(int index) const;
  bool contains_owner(int owner) const;
  std::uint8_t operator[](std::size_t owner) const {
    return contains_owner(static_cast<int>(owner)) ? 1u : 0u;
  }
};

// Contract: a new state whose root counters equal the input's for every
// colouring of the selected owners, with every other unit frozen at its
// background value. It keeps the root, the selected live leaves and the
// branch points of their spanning tree, folds each corridor between two kept
// anchors into one C2 unit and every other hanging branch into its I
// anchor's local background G, then builds the HDT. It reads only the
// component graph, so a deferred state is accepted. The output has at most
// 4*max(1,s)+3 units for s selected live atoms; a larger one throws. The
// input is not modified.
HdtContractResult hdt_contract(
    const HdtState& input,
    const HdtLiveSelector& selector,
    HdtContractMode mode = HdtContractMode::LightCopy);

// Builds the first state of a gene. A graph with absent leaves is first
// compacted by an InitialCompact Contract that keeps every live owner.
HdtInitialStateResult hdt_build_initial_state(HdtComponentGraph graph);
// True when the graph contains a weight-zero absent-leaf unit, which requires
// the InitialCompact Contract before an HDT is built.
bool hdt_component_graph_has_absent_leaf(const HdtComponentGraph& graph);

// Owning selector; owners out of range or repeated are ignored.
HdtLiveSelector hdt_owner_selector(int owner_count,
                                   const std::vector<int>& owners);

// Borrowing selector; both vectors must outlive it. Pass -1 as
// exact_selected_atom_count when the count is unknown.
HdtLiveSelector hdt_interval_owner_selector(
    const std::vector<int>& represented_owner_order,
    const std::vector<int>& owner_position,
    int owner_begin,
    int owner_end,
    int membership_begin,
    int membership_end,
    int exact_selected_atom_count);

// Live leaves of selected owners, each duplicate leaf counted. O(1) when the
// selector carries the exact count, otherwise a scan of the units.
int hdt_selected_atom_count(const HdtState& state,
                            const HdtLiveSelector& selector);

// 4*max(1,s)+3, the Contract output bound for s selected live atoms.
int hdt_contract_size_bound(int selected_live_atoms);

// Extract/Contract: the child state for the selected owners, built from the
// parent's HDT without a rescan of its component graph. It marks the HDT
// paths from every selected atom, and from the root I anchor, to the root;
// walks the marked nodes and emits each atom it reaches as a unit and each
// unmarked internal node as one opaque unit, all but the selected live
// leaves frozen at their background values; reconnects the emitted units by
// their source ports; and contracts the result unless it already meets the
// Contract bound with an anchor root. With an interval selector the work
// is linear in the selected owners and the marked nodes, once the workspace
// has grown to the source's key domains. The parent is not modified.
HdtExtractContractResult hdt_extract_contract(
    const HdtState& parent,
    const HdtLiveSelector& selector,
    HdtExtractWorkspace& workspace);

// As above with a temporary workspace.
HdtExtractContractResult hdt_extract_contract(
    const HdtState& parent,
    const HdtLiveSelector& selector);

}  // namespace qraft::internal
