#pragma once

#include "counting/counters.hpp"
#include "core/checked.hpp"
#include "counting/query_taxa.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace qraft::internal {

// An ordinary port is one end of a gene edge, keyed by the edge's child
// vertex. A centre port joins a G component to the I anchor at its centre,
// keyed by the anchor's gene vertex.
enum class HdtPortKind : std::uint8_t {
  OrdinaryEdge,
  Center,
};

// Direction of the edge seen from the port's own component: Up leads to the
// parent side, Down to the child side. Centre ports use Center.
enum class HdtPortOrientation : std::uint8_t {
  Up,
  Down,
  Center,
};

struct HdtPortKey {
  int source_id = -1;
  HdtPortKind kind = HdtPortKind::OrdinaryEdge;
  HdtPortOrientation orientation = HdtPortOrientation::Up;

  constexpr HdtPortKey() = default;
  constexpr HdtPortKey(HdtPortKind kind_value,
                             int source_id_value,
                             HdtPortOrientation orientation_value)
      : source_id(source_id_value),
        kind(kind_value),
        orientation(orientation_value) {}
};

inline bool operator==(const HdtPortKey& left,
                       const HdtPortKey& right) {
  return left.kind == right.kind && left.source_id == right.source_id &&
         left.orientation == right.orientation;
}

// LiveLeaf is a gene leaf with an owner; its counters follow from its colour
// and weight. IAnchor is a gene vertex kept as an I component, which has no
// counters. The opaque kinds store counters of the named interface; they are
// contracted components, and C1Opaque also holds weight-zero absent leaves.
enum class HdtComponentKind : std::uint8_t {
  LiveLeaf,
  IAnchor,
  ClosedOpaque,
  C1Opaque,
  C2Opaque,
  GOpaque,
};

// Domains of the gene vertex and edge keys carried by a component graph's
// ports. Ordinary-edge keys are child vertex ids, so both domains are the
// gene's node count.
struct HdtSourceIndex {
  int root_vertex = -1;
  int vertex_slots = 0;
  int edge_slots = 0;
};

struct HdtComponentUnit {
  HdtComponentKind kind = HdtComponentKind::IAnchor;
  std::uint64_t stable_id = 0;
  int owner = -1;  // -1 unless kind is LiveLeaf
  int source_vertex_id = -1;
  // The unit's edges are incidence_edges[incidence_begin, +incidence_count).
  int incidence_begin = 0;
  int incidence_count = 0;
  double weight = 0.0;  // leaf mass; zero for an absent leaf
  std::uint8_t current_color = 0;
  std::uint8_t background_color = 0;
  // Colours with nonzero mass in the unit; they select packed coordinates.
  std::uint8_t current_color_mask = 0;
  std::uint8_t background_color_mask = 0;
  // Offsets of stored counters in the graph's payload arrays; -1 for a live
  // leaf or an anchor, and for a background equal to the current value.
  int dp_offset = -1;
  int background_dp_offset = -1;
  std::uint8_t dp_dimension = 0;
  CounterOrientation dp_orientation = CounterOrientation::None;
  bool current_equals_background = false;
  // For an anchor: the contracted branches hanging at its centre, folded
  // into one G value whose current equals its background.
  std::shared_ptr<const CounterPair> local_background_g;
};

struct HdtIncidenceEdge {
  int a = -1;
  int b = -1;
  HdtPortKey a_port;
  HdtPortKey b_port;
};

struct HdtExternalPort {
  int unit = -1;
  HdtPortKey port;
};

// Current and background colour of each owner. The defaults apply to owners
// outside the vectors and to leaves without an owner.
struct HdtColorPolicy {
  std::uint8_t default_current = 1;
  std::uint8_t default_background = 0;
  std::vector<std::uint8_t> current_by_owner;
  std::vector<std::uint8_t> background_by_owner;
};

// A tree of component units joined by edges with matching port keys. Stored
// counters of opaque units live in the two payload arrays.
struct HdtComponentGraph {
  int palette = 0;
  int root_unit = -1;
  std::shared_ptr<const HdtSourceIndex> source;
  std::vector<HdtComponentUnit> units;
  std::vector<double> current_payload;
  std::vector<double> background_payload;
  std::vector<HdtIncidenceEdge> edges;
  // Ports left open at the graph's boundary; empty for a whole gene.
  std::vector<HdtExternalPort> external_ports;
  // Edge indices grouped by unit, in increasing edge order within a unit.
  std::vector<int> incidence_edges;
  int owner_count = 0;
  long long construction_scalar_fmas = 0;
  // True when every leaf's background colour is 0, or 0 or 3 for K=4.
  bool background_restricted = false;

  int component_units() const {
    return checked_hdt_int_size(
        units.size(), "component graph exceeds its integer domain");
  }
};

// One unit per gene vertex (unit id = vertex id) and one edge per gene edge,
// with the parent's port Down and the child's Up. Leaves with an owner are
// live; every other leaf is a weight-zero C1 atom that keeps the gene's
// topology. Every duplicate leaf of an owner is kept.
HdtComponentGraph build_base_component_graph(
    const GeneIndex& gene,
    const OwnerProjection& owners,
    int palette,
    const HdtColorPolicy& colors);

// Builds the base graph of the gene restricted to leaves with an owner. This
// avoids building weight-zero atoms for absent leaves only to contract them.
HdtComponentGraph build_induced_base_component_graph(
    const GeneIndex& gene,
    const OwnerProjection& owners,
    int palette,
    const HdtColorPolicy& colors);

// Restricts `gene` to its leaves with an owner: branches without such a leaf
// are dropped and the remaining degree-two vertices suppressed. Retained
// vertices keep their gene order; with no live leaf the result is empty.
// `changed` is false, and the result empty, when every leaf is live; callers
// then keep using `gene`.
GeneIndex induce_gene_for_owner_projection(
    const GeneIndex& gene,
    const OwnerProjection& owners,
    bool& changed);

// Rebuilds every unit's incidence range from graph.edges.
void rebuild_component_graph_indices(HdtComponentGraph& graph);

// Checks that the graph is a tree of R units with R - 1 edges whose port keys
// match, and that each unit's degree and payload metadata fit its kind.
bool validate_component_graph(const HdtComponentGraph& graph,
                              std::string* error = nullptr);

// Store or read a unit's counters. A live leaf stores only metadata; its
// counters are computed from its colours and weight. An anchor has none.
void hdt_component_store_dp(HdtComponentGraph& graph,
                            HdtComponentUnit& unit,
                            const CounterPair& value);
CounterPair hdt_component_dp(const HdtComponentGraph& graph,
                             int unit);


}  // namespace qraft::internal
