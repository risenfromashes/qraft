#include "counting/units.hpp"

#include <algorithm>
#include <deque>
#include <limits>
#include <stdexcept>

namespace qraft::internal {
namespace {

std::uint8_t policy_color(const std::vector<std::uint8_t>& colors,
                          int owner,
                          std::uint8_t fallback) {
  return owner >= 0 && static_cast<std::size_t>(owner) < colors.size()
             ? colors[static_cast<std::size_t>(owner)]
             : fallback;
}

bool is_leaf_node(const GeneIndex::Node& node) {
  return node.key != -1;
}

void set_error(std::string* error, const std::string& message) {
  if (error != nullptr) *error = message;
}

bool ordinary_pair(const HdtPortKey& a, const HdtPortKey& b) {
  return a.kind == HdtPortKind::OrdinaryEdge &&
         b.kind == HdtPortKind::OrdinaryEdge &&
         a.source_id == b.source_id && a.source_id >= 0 &&
         ((a.orientation == HdtPortOrientation::Up &&
           b.orientation == HdtPortOrientation::Down) ||
          (a.orientation == HdtPortOrientation::Down &&
           b.orientation == HdtPortOrientation::Up));
}

bool center_pair(const HdtPortKey& a, const HdtPortKey& b) {
  return a.kind == HdtPortKind::Center &&
         b.kind == HdtPortKind::Center && a.source_id == b.source_id &&
         a.source_id >= 0 &&
         a.orientation == HdtPortOrientation::Center &&
         b.orientation == HdtPortOrientation::Center;
}

CounterInterface component_interface(HdtComponentKind kind) {
  switch (kind) {
    case HdtComponentKind::LiveLeaf:
    case HdtComponentKind::C1Opaque:
      return CounterInterface::C1;
    case HdtComponentKind::C2Opaque:
      return CounterInterface::C2;
    case HdtComponentKind::GOpaque:
      return CounterInterface::G;
    case HdtComponentKind::ClosedOpaque:
      return CounterInterface::Closed;
    case HdtComponentKind::IAnchor:
      break;
  }
  throw std::logic_error("I anchor has no component payload interface");
}

bool valid_component_orientation(HdtComponentKind kind,
                                 CounterOrientation orientation) {
  if (kind == HdtComponentKind::LiveLeaf) {
    return orientation == CounterOrientation::U;
  }
  if (kind == HdtComponentKind::C1Opaque) {
    return orientation == CounterOrientation::U ||
           orientation == CounterOrientation::D;
  }
  return kind != HdtComponentKind::IAnchor &&
         orientation == CounterOrientation::None;
}

// Preorder and parents of `gene` from `root`; the root is its own parent.
GeneTraversal induced_gene_traversal(const GeneIndex& gene, int root) {
  GeneTraversal traversal;
  const int node_count = checked_hdt_int_size(
      gene.nodes.size(), "induced source exceeds its integer node domain");
  traversal.parent.assign(gene.nodes.size(), -1);
  if (root < 0 || root >= node_count) return traversal;
  traversal.preorder.reserve(gene.nodes.size());

  std::vector<std::pair<int, int>> stack{{root, 0}};
  traversal.parent[static_cast<std::size_t>(root)] = root;
  while (!stack.empty()) {
    const int node = stack.back().first;
    int& next = stack.back().second;
    if (next == 0) traversal.preorder.push_back(node);
    if (next < gene_degree(gene, node)) {
      const int neighbor =
          gene_neighbors(gene, node)[static_cast<std::size_t>(next++)];
      if (neighbor < 0 || neighbor >= node_count ||
          neighbor == traversal.parent[static_cast<std::size_t>(node)]) {
        continue;
      }
      if (traversal.parent[static_cast<std::size_t>(neighbor)] >= 0) {
        throw std::logic_error("induced source is not a tree");
      }
      traversal.parent[static_cast<std::size_t>(neighbor)] = node;
      stack.push_back({neighbor, 0});
      continue;
    }
    stack.pop_back();
  }
  if (traversal.preorder.size() != gene.nodes.size()) {
    throw std::logic_error("induced source is disconnected");
  }
  return traversal;
}

// Keeps the live leaves and the vertices where three or more kept edges meet.
// An edge is kept when both of its sides hold a live leaf; each path of
// suppressed degree-two vertices becomes one edge.
GeneIndex induce_gene_for_owners_impl(const GeneIndex& gene,
                                      const OwnerProjection& owners,
                                      bool& changed) {
  changed = false;
  if (gene.nodes.empty()) return GeneIndex{};
  const int node_count = checked_hdt_int_size(
      gene.nodes.size(), "source exceeds its integer node domain");
  if (gene.root < 0 || gene.root >= node_count ||
      gene.root_traversal.parent.size() != gene.nodes.size() ||
      gene.root_traversal.preorder.size() != gene.nodes.size()) {
    throw std::logic_error("source traversal is unavailable for induction");
  }

  std::vector<std::uint8_t> live_leaf(gene.nodes.size(), 0);
  std::vector<int> subtree_live(gene.nodes.size(), 0);
  int live_count = 0;
  int physical_leaf_count = 0;
  for (int node = 0; node < node_count; ++node) {
    const GeneIndex::Node& source =
        gene.nodes[static_cast<std::size_t>(node)];
    if (!is_leaf_node(source)) continue;
    ++physical_leaf_count;
    if (owner_for_gene_leaf(source, owners) >= 0) {
      live_leaf[static_cast<std::size_t>(node)] = 1;
      ++live_count;
    }
  }
  if (live_count == physical_leaf_count) return GeneIndex{};
  changed = true;
  if (live_count == 0) return GeneIndex{};

  const std::vector<int>& preorder = gene.root_traversal.preorder;
  for (auto it = preorder.rbegin(); it != preorder.rend(); ++it) {
    const int node = *it;
    int count = live_leaf[static_cast<std::size_t>(node)] != 0 ? 1 : 0;
    for (int neighbor : gene_neighbors(gene, node)) {
      if (neighbor >= 0 && neighbor < node_count &&
          gene.root_traversal.parent[static_cast<std::size_t>(neighbor)] ==
              node) {
        count += subtree_live[static_cast<std::size_t>(neighbor)];
      }
    }
    subtree_live[static_cast<std::size_t>(node)] = count;
  }
  if (subtree_live[static_cast<std::size_t>(gene.root)] != live_count) {
    throw std::logic_error("source traversal lost a live leaf");
  }

  auto kept_edge = [&](int a, int b) {
    int child = -1;
    if (gene.root_traversal.parent[static_cast<std::size_t>(b)] == a) {
      child = b;
    } else if (gene.root_traversal.parent[static_cast<std::size_t>(a)] == b) {
      child = a;
    } else {
      throw std::logic_error("source adjacency disagrees with its traversal");
    }
    const int below = subtree_live[static_cast<std::size_t>(child)];
    return below > 0 && below < live_count;
  };

  std::vector<int> kept_degree(gene.nodes.size(), 0);
  if (live_count > 1) {
    for (int node = 0; node < node_count; ++node) {
      for (int neighbor : gene_neighbors(gene, node)) {
        if (neighbor >= 0 && neighbor < node_count &&
            kept_edge(node, neighbor)) {
          ++kept_degree[static_cast<std::size_t>(node)];
        }
      }
    }
  }

  // A non-leaf vertex with kept degree two lies inside a path and is
  // suppressed; one with kept degree one cannot occur.
  std::vector<int> source_to_induced(gene.nodes.size(), -1);
  GeneIndex induced;
  induced.nodes.reserve(checked_hdt_size_mul(
      static_cast<std::size_t>(live_count), 2u,
      "induced source reserve exceeds its size domain"));
  for (int node = 0; node < node_count; ++node) {
    const bool retain =
        live_leaf[static_cast<std::size_t>(node)] != 0 ||
        (kept_degree[static_cast<std::size_t>(node)] > 0 &&
         kept_degree[static_cast<std::size_t>(node)] != 2);
    if (!retain) continue;
    source_to_induced[static_cast<std::size_t>(node)] =
        checked_hdt_int_size(
            induced.nodes.size(),
            "induced source exceeds its integer node domain");
    GeneIndex::Node compact;
    double compact_mass = 1.0;
    if (live_leaf[static_cast<std::size_t>(node)] != 0) {
      compact = gene.nodes[static_cast<std::size_t>(node)];
      compact_mass = gene_leaf_mass(gene, node);
    }
    induced.nodes.push_back(compact);
    if (!gene.leaf_mass.empty()) induced.leaf_mass.push_back(compact_mass);
  }
  if (induced.nodes.empty()) {
    throw std::logic_error("source induction lost every live leaf");
  }

  // Retained sources are numbered in source order, so rows are emitted in
  // induced node order.
  induced.adjacency_offsets.reserve(induced.nodes.size() + 1u);
  induced.adjacency.reserve(2u * (induced.nodes.size() - 1u));
  induced.adjacency_offsets.push_back(0);
  for (int source = 0; source < node_count; ++source) {
    const int compact_source =
        source_to_induced[static_cast<std::size_t>(source)];
    if (compact_source < 0) continue;
    for (int first : gene_neighbors(gene, source)) {
      if (first < 0 || first >= node_count ||
          (live_count > 1 && !kept_edge(source, first))) {
        continue;
      }
      if (live_count == 1) continue;
      // Follow the path of suppressed vertices to the next retained one.
      int previous = source;
      int current = first;
      int steps = 0;
      while (source_to_induced[static_cast<std::size_t>(current)] < 0) {
        int next = -1;
        for (int neighbor : gene_neighbors(gene, current)) {
          if (neighbor == previous || neighbor < 0 ||
              neighbor >= node_count || !kept_edge(current, neighbor)) {
            continue;
          }
          if (next >= 0) {
            throw std::logic_error(
                "induced source corridor has more than two arms");
          }
          next = neighbor;
        }
        if (next < 0 || ++steps > node_count) {
          throw std::logic_error("induced source corridor is incomplete");
        }
        previous = current;
        current = next;
      }
      induced.adjacency.push_back(
          source_to_induced[static_cast<std::size_t>(current)]);
    }
    induced.adjacency_offsets.push_back(checked_hdt_int_size(
        induced.adjacency.size(),
        "induced source exceeds its half-edge domain"));
  }

  // Root at the first retained internal vertex, else at node 0 (one or two
  // live leaves).
  induced.root = 0;
  for (int node = 0; node < static_cast<int>(induced.nodes.size()); ++node) {
    if (!is_leaf_node(induced.nodes[static_cast<std::size_t>(node)]) &&
        gene_degree(induced, node) >= 2) {
      induced.root = node;
      break;
    }
  }
  induced.reference_root = induced.root;
  induced.root_traversal = induced_gene_traversal(induced, induced.root);
  return induced;
}

}  // namespace

void hdt_component_store_dp(HdtComponentGraph& graph,
                            HdtComponentUnit& unit,
                            const CounterPair& value) {
  const int dimension = value.current.dimension;
  if (unit.kind == HdtComponentKind::IAnchor || dimension <= 0 ||
      value.background.dimension != dimension ||
      value.current.interface != value.background.interface ||
      value.current.orientation != value.background.orientation ||
      value.current.interface != component_interface(unit.kind) ||
      !valid_component_orientation(unit.kind, value.current.orientation) ||
      dimension != counters_dimension(graph.palette, value.current.interface)) {
    throw std::invalid_argument("invalid component dual payload");
  }
  unit.current_color_mask =
      counters_value_color_mask(graph.palette, value.current);
  unit.background_color_mask = value.current_equals_background
                                   ? unit.current_color_mask
                                   : counters_value_color_mask(
                                         graph.palette, value.background);
  if (unit.kind == HdtComponentKind::LiveLeaf) {
    unit.dp_offset = -1;
    unit.background_dp_offset = -1;
    unit.dp_dimension = static_cast<std::uint8_t>(dimension);
    unit.dp_orientation = value.current.orientation;
    unit.current_equals_background = value.current_equals_background;
    return;
  }
  const bool dimension_changed = unit.dp_dimension != dimension;
  const std::size_t old_end = unit.dp_offset < 0
                                  ? 0
                                  : checked_hdt_size_add(
                                        static_cast<std::size_t>(unit.dp_offset),
                                        static_cast<std::size_t>(dimension),
                                        "component payload offset overflow");
  if (unit.dp_offset < 0 || unit.dp_dimension != dimension ||
      old_end > graph.current_payload.size()) {
    unit.dp_offset = checked_hdt_int_size(
        graph.current_payload.size(),
        "component payload exceeds its integer offset domain");
    unit.dp_dimension = static_cast<std::uint8_t>(dimension);
    const std::size_t next_size = checked_hdt_size_add(
        graph.current_payload.size(), static_cast<std::size_t>(dimension),
        "component payload capacity overflow");
    (void)checked_hdt_int_size(
        next_size, "component payload exceeds its integer offset domain");
    graph.current_payload.resize(next_size);
  }
  std::copy_n(value.current.values.data(), dimension,
              graph.current_payload.data() + unit.dp_offset);
  if (value.current_equals_background) {
    unit.background_dp_offset = -1;
  } else {
    const std::size_t background_end =
        unit.background_dp_offset < 0
            ? 0
            : checked_hdt_size_add(
                  static_cast<std::size_t>(unit.background_dp_offset),
                  static_cast<std::size_t>(dimension),
                  "component background payload offset overflow");
    if (dimension_changed || unit.background_dp_offset < 0 ||
        background_end > graph.background_payload.size()) {
      unit.background_dp_offset = checked_hdt_int_size(
          graph.background_payload.size(),
          "component background payload exceeds its integer offset domain");
      const std::size_t next_size = checked_hdt_size_add(
          graph.background_payload.size(), static_cast<std::size_t>(dimension),
          "component background payload capacity overflow");
      (void)checked_hdt_int_size(
          next_size,
          "component background payload exceeds its integer offset domain");
      graph.background_payload.resize(next_size);
    }
    std::copy_n(value.background.values.data(), dimension,
                graph.background_payload.data() +
                    unit.background_dp_offset);
  }
  unit.dp_orientation = value.current.orientation;
  unit.current_equals_background = value.current_equals_background;
}

CounterPair hdt_component_dp(const HdtComponentGraph& graph,
                             int unit_index) {
  if (unit_index < 0 || unit_index >= graph.component_units()) {
    throw std::out_of_range("component payload unit out of range");
  }
  const HdtComponentUnit& unit =
      graph.units[static_cast<std::size_t>(unit_index)];
  if (unit.kind == HdtComponentKind::LiveLeaf) {
    return counters_leaf_dual(graph.palette, unit.current_color,
                              unit.background_color, unit.weight);
  }
  if (unit.dp_offset < 0 || unit.dp_dimension == 0) return CounterPair{};
  const std::size_t end = static_cast<std::size_t>(unit.dp_offset) +
                          static_cast<std::size_t>(unit.dp_dimension);
  const std::size_t background_end =
      unit.background_dp_offset < 0
          ? 0
          : static_cast<std::size_t>(unit.background_dp_offset) +
                static_cast<std::size_t>(unit.dp_dimension);
  if (end > graph.current_payload.size() ||
      (!unit.current_equals_background &&
       (unit.background_dp_offset < 0 ||
        background_end > graph.background_payload.size()))) {
    throw std::logic_error("component payload offset out of range");
  }
  CounterPair value;
  value.current.dimension = unit.dp_dimension;
  value.background.dimension = unit.dp_dimension;
  value.current.interface = value.background.interface =
      component_interface(unit.kind);
  value.current.orientation = value.background.orientation =
      unit.dp_orientation;
  std::copy_n(graph.current_payload.data() + unit.dp_offset,
              unit.dp_dimension, value.current.values.data());
  const double* background = unit.current_equals_background
                                 ? graph.current_payload.data() + unit.dp_offset
                                 : graph.background_payload.data() +
                                       unit.background_dp_offset;
  std::copy_n(background, unit.dp_dimension,
              value.background.values.data());
  value.current_equals_background = unit.current_equals_background;
  return value;
}

HdtComponentGraph build_base_component_graph(
    const GeneIndex& gene,
    const OwnerProjection& owners,
    int palette,
    const HdtColorPolicy& colors) {
  const long long scalar_fmas_before = counters_thread_scalar_fmas();
  HdtComponentGraph graph;
  graph.palette = palette;
  graph.owner_count = owners.owner_count;
  graph.background_restricted = true;
  if (gene.nodes.empty()) return graph;

  const int gene_node_count = checked_hdt_int_size(
      gene.nodes.size(), "source exceeds its integer node domain");

  auto source = std::make_shared<HdtSourceIndex>();
  source->root_vertex = gene.root;
  source->vertex_slots = gene_node_count;
  source->edge_slots = gene_node_count;
  const std::vector<int>& source_parent = gene.root_traversal.parent;

  // Only absent leaves store a payload, so the payload arrays grow on demand.
  graph.units.reserve(gene.nodes.size());
  for (int vertex = 0; vertex < gene_node_count; ++vertex) {
    const GeneIndex::Node& node = gene.nodes[static_cast<std::size_t>(vertex)];
    HdtComponentUnit unit;
    unit.stable_id = static_cast<std::uint64_t>(vertex);
    unit.source_vertex_id = vertex;
    CounterPair payload;
    const bool have_payload = is_leaf_node(node);
    bool store_payload = false;
    if (have_payload) {
      unit.owner = owner_for_gene_leaf(node, owners);
      unit.kind = unit.owner >= 0 ? HdtComponentKind::LiveLeaf
                                  : HdtComponentKind::C1Opaque;
      // A leaf without an owner is outside the subproblem. As a weight-zero
      // C1 atom it keeps the gene's topology and contributes no mass.
      unit.weight = unit.owner >= 0 ? gene_leaf_mass(gene, vertex) : 0.0;
      unit.current_color = policy_color(colors.current_by_owner, unit.owner,
                                        colors.default_current);
      unit.background_color = policy_color(
          colors.background_by_owner, unit.owner, colors.default_background);
      unit.current_color_mask = unit.weight == 0.0
                                    ? 0
                                    : static_cast<std::uint8_t>(
                                          1u << unit.current_color);
      unit.background_color_mask = unit.weight == 0.0
                                       ? 0
                                       : static_cast<std::uint8_t>(
                                             1u << unit.background_color);
      const bool allowed_background =
          palette == 2 ? unit.background_color == 0
                       : (unit.background_color == 0 ||
                          unit.background_color == 3);
      graph.background_restricted =
          graph.background_restricted && allowed_background;
      if (unit.kind == HdtComponentKind::LiveLeaf) {
        unit.dp_dimension = static_cast<std::uint8_t>(
            counters_dimension(palette, CounterInterface::C1));
        unit.dp_orientation = CounterOrientation::U;
        unit.current_equals_background =
            unit.current_color == unit.background_color || unit.weight == 0.0;
      } else {
        payload = counters_leaf_dual(palette, unit.current_color,
                                     unit.background_color, unit.weight);
        store_payload = true;
      }
    } else {
      unit.kind = HdtComponentKind::IAnchor;
    }
    graph.units.push_back(std::move(unit));
    if (store_payload) {
      hdt_component_store_dp(graph, graph.units.back(), payload);
    }
  }
  graph.root_unit = gene.root;

  graph.edges.reserve(gene.nodes.size() - 1u);
  for (int parent : gene.root_traversal.preorder) {
    for (int child : gene_neighbors(gene, parent)) {
      if (child < 0 || child >= gene_node_count ||
          source_parent[static_cast<std::size_t>(child)] != parent) {
        continue;
      }
      const int edge_id = child;
      HdtIncidenceEdge edge;
      edge.a = parent;
      edge.b = child;
      edge.a_port = HdtPortKey{HdtPortKind::OrdinaryEdge, edge_id,
                                     HdtPortOrientation::Down};
      edge.b_port = HdtPortKey{HdtPortKind::OrdinaryEdge, edge_id,
                                     HdtPortOrientation::Up};
      graph.edges.push_back(edge);
    }
  }
  graph.source = std::move(source);
  rebuild_component_graph_indices(graph);
  graph.construction_scalar_fmas =
      counters_thread_scalar_fmas() - scalar_fmas_before;
  return graph;
}

HdtComponentGraph build_induced_base_component_graph(
    const GeneIndex& gene,
    const OwnerProjection& owners,
    int palette,
    const HdtColorPolicy& colors) {
  bool changed = false;
  GeneIndex induced = induce_gene_for_owners_impl(gene, owners, changed);
  return changed ? build_base_component_graph(induced, owners, palette, colors)
                 : build_base_component_graph(gene, owners, palette, colors);
}

GeneIndex induce_gene_for_owner_projection(
    const GeneIndex& gene,
    const OwnerProjection& owners,
    bool& changed) {
  return induce_gene_for_owners_impl(gene, owners, changed);
}

void rebuild_component_graph_indices(HdtComponentGraph& graph) {
  const int unit_count = graph.component_units();
  const int edge_count = checked_hdt_int_size(
      graph.edges.size(), "component graph exceeds its edge domain");
  for (HdtComponentUnit& unit : graph.units) {
    unit.incidence_count = 0;
  }
  for (const HdtIncidenceEdge& edge : graph.edges) {
    if (edge.a >= 0 && edge.a < unit_count) {
      HdtComponentUnit& unit =
          graph.units[static_cast<std::size_t>(edge.a)];
      unit.incidence_count = checked_hdt_int_add(
          unit.incidence_count, 1,
          "component incidence count overflow");
    }
    if (edge.b >= 0 && edge.b < unit_count) {
      HdtComponentUnit& unit =
          graph.units[static_cast<std::size_t>(edge.b)];
      unit.incidence_count = checked_hdt_int_add(
          unit.incidence_count, 1,
          "component incidence count overflow");
    }
  }
  int total = 0;
  for (int unit = 0; unit < unit_count; ++unit) {
    HdtComponentUnit& component =
        graph.units[static_cast<std::size_t>(unit)];
    component.incidence_begin = total;
    total = checked_hdt_int_add(
        total, component.incidence_count,
        "component incidence storage overflow");
    component.incidence_count = 0;
  }
  graph.incidence_edges.assign(static_cast<std::size_t>(total), -1);
  for (int edge_index = 0; edge_index < edge_count; ++edge_index) {
    const HdtIncidenceEdge& edge =
        graph.edges[static_cast<std::size_t>(edge_index)];
    for (int unit : {edge.a, edge.b}) {
      if (unit < 0 || unit >= unit_count) continue;
      HdtComponentUnit& component =
          graph.units[static_cast<std::size_t>(unit)];
      const int position =
          component.incidence_begin + component.incidence_count++;
      graph.incidence_edges[static_cast<std::size_t>(position)] = edge_index;
    }
  }

}

bool validate_component_graph(const HdtComponentGraph& graph,
                              std::string* error) {
  if (graph.units.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      graph.edges.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      graph.external_ports.size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    set_error(error, "component graph exceeds its integer domain");
    return false;
  }
  const int unit_count = static_cast<int>(graph.units.size());
  if (graph.palette != 2 && graph.palette != 4) {
    set_error(error, "invalid palette");
    return false;
  }
  if (unit_count == 0) {
    if (graph.root_unit != -1 || !graph.edges.empty() ||
        !graph.current_payload.empty() || !graph.background_payload.empty()) {
      set_error(error, "nonempty metadata on empty graph");
      return false;
    }
    return true;
  }
  if (graph.root_unit < 0 || graph.root_unit >= unit_count) {
    set_error(error, "invalid root unit");
    return false;
  }
  if (static_cast<int>(graph.edges.size()) != unit_count - 1) {
    set_error(error, "component graph edge count is not R-1");
    return false;
  }
  for (const HdtIncidenceEdge& edge : graph.edges) {
    if (edge.a < 0 || edge.a >= unit_count || edge.b < 0 ||
        edge.b >= unit_count || edge.a == edge.b) {
      set_error(error, "invalid incidence endpoint");
      return false;
    }
    if (!ordinary_pair(edge.a_port, edge.b_port) &&
        !center_pair(edge.a_port, edge.b_port)) {
      set_error(error, "incompatible incidence port keys");
      return false;
    }
  }
  std::vector<int> external_count(static_cast<std::size_t>(unit_count), 0);
  for (const HdtExternalPort& external : graph.external_ports) {
    if (external.unit < 0 || external.unit >= unit_count ||
        external.port.source_id < 0 ||
        (external.port.kind == HdtPortKind::Center &&
         external.port.orientation != HdtPortOrientation::Center) ||
        (external.port.kind == HdtPortKind::OrdinaryEdge &&
         external.port.orientation == HdtPortOrientation::Center)) {
      set_error(error, "invalid external port");
      return false;
    }
    ++external_count[static_cast<std::size_t>(external.unit)];
  }

  std::vector<std::uint8_t> seen(static_cast<std::size_t>(unit_count), 0);
  std::vector<int> stack{graph.root_unit};
  int reached = 0;
  while (!stack.empty()) {
    const int unit = stack.back();
    stack.pop_back();
    if (seen[static_cast<std::size_t>(unit)] != 0) continue;
    seen[static_cast<std::size_t>(unit)] = 1;
    ++reached;
    const HdtComponentUnit& component =
        graph.units[static_cast<std::size_t>(unit)];
    for (int i = component.incidence_count - 1; i >= 0; --i) {
      const int edge_index = graph.incidence_edges[static_cast<std::size_t>(
          component.incidence_begin + i)];
      const HdtIncidenceEdge& edge =
          graph.edges[static_cast<std::size_t>(edge_index)];
      stack.push_back(edge.a == unit ? edge.b : edge.a);
    }
  }
  if (reached != unit_count) {
    set_error(error, "component graph is disconnected");
    return false;
  }

  for (int unit_index = 0; unit_index < unit_count; ++unit_index) {
    const HdtComponentUnit& unit =
        graph.units[static_cast<std::size_t>(unit_index)];
    const int degree = unit.incidence_count +
                       external_count[static_cast<std::size_t>(unit_index)];
    const bool should_have_payload = unit.kind != HdtComponentKind::IAnchor &&
                                     unit.kind != HdtComponentKind::LiveLeaf;
    const std::size_t payload_end =
        unit.dp_offset < 0
            ? 0
            : static_cast<std::size_t>(unit.dp_offset) + unit.dp_dimension;
    const std::size_t background_end =
        unit.background_dp_offset < 0
            ? 0
            : static_cast<std::size_t>(unit.background_dp_offset) +
                  unit.dp_dimension;
    if ((should_have_payload &&
         (unit.dp_offset < 0 || unit.dp_dimension != counters_dimension(
                                                    graph.palette,
                                                    component_interface(
                                                        unit.kind)) ||
          payload_end > graph.current_payload.size() ||
          (unit.current_equals_background
               ? (unit.background_dp_offset >= 0 &&
                  background_end > graph.background_payload.size())
               : (unit.background_dp_offset < 0 ||
                  background_end > graph.background_payload.size())) ||
          !valid_component_orientation(unit.kind, unit.dp_orientation))) ||
        (unit.kind == HdtComponentKind::LiveLeaf &&
         (unit.dp_offset != -1 || unit.background_dp_offset != -1 ||
          unit.dp_dimension != counters_dimension(
                                   graph.palette, CounterInterface::C1) ||
          unit.dp_orientation != CounterOrientation::U)) ||
        (unit.kind == HdtComponentKind::IAnchor &&
         (unit.dp_offset != -1 || unit.background_dp_offset != -1 ||
          unit.dp_dimension != 0 ||
          unit.dp_orientation != CounterOrientation::None))) {
      set_error(error, "invalid component payload metadata");
      return false;
    }
    if (unit.local_background_g != nullptr &&
        (unit.kind != HdtComponentKind::IAnchor ||
         unit.local_background_g->current.interface != CounterInterface::G ||
         unit.local_background_g->background.interface !=
             CounterInterface::G)) {
      set_error(error, "invalid local background G payload");
      return false;
    }
    switch (unit.kind) {
      case HdtComponentKind::LiveLeaf:
        if (unit.owner < 0 || degree > 1) {
          set_error(error, "invalid live-leaf unit");
          return false;
        }
        break;
      case HdtComponentKind::IAnchor:
        if (unit.owner != -1) {
          set_error(error, "I anchor has a live owner");
          return false;
        }
        break;
      case HdtComponentKind::ClosedOpaque:
        if (degree != 0 || unit.owner != -1) {
          set_error(error, "invalid Closed opaque unit");
          return false;
        }
        break;
      case HdtComponentKind::C1Opaque:
        if (degree > 1 || unit.owner != -1) {
          set_error(error, "invalid C1 opaque unit");
          return false;
        }
        break;
      case HdtComponentKind::C2Opaque:
        if (degree != 2 || unit.owner != -1) {
          set_error(error, "invalid C2 opaque unit");
          return false;
        }
        break;
      case HdtComponentKind::GOpaque:
        if (degree != 1 || unit.owner != -1) {
          set_error(error, "invalid G opaque unit");
          return false;
        }
        break;
    }
  }
  return true;
}

}  // namespace qraft::internal
