#include "internal.hpp"

#include <limits>

namespace qraft::internal {
namespace {

/// The empty gene: no nodes, both roots -1.
GeneIndex canonical_empty_gene_index() {
  GeneIndex result;
  result.root = -1;
  result.reference_root = -1;
  result.root_traversal = traversal(result, result.root);
  return result;
}

/// Checks the CSR framing before any neighbor span is formed.
void validate_gene_adjacency_shape(const GeneIndex& gene) {
  constexpr std::size_t kMaxIndexedEdges =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  if (gene.adjacency.size() > kMaxIndexedEdges) {
    throw std::length_error("gene adjacency exceeds its half-edge domain");
  }
  if (gene.nodes.empty()) {
    if (!gene.adjacency.empty() || gene.adjacency_offsets.size() > 1u) {
      throw std::invalid_argument("empty gene has adjacency");
    }
    return;
  }
  if (gene.adjacency_offsets.size() != gene.nodes.size() + 1u ||
      gene.adjacency_offsets.front() != 0 ||
      static_cast<std::size_t>(gene.adjacency_offsets.back()) !=
          gene.adjacency.size()) {
    throw std::invalid_argument("gene adjacency offsets are inconsistent");
  }
  for (std::size_t i = 1; i < gene.adjacency_offsets.size(); ++i) {
    if (gene.adjacency_offsets[i] < gene.adjacency_offsets[i - 1u]) {
      throw std::invalid_argument("gene adjacency offsets are inconsistent");
    }
  }
}

/// Appends one node's CSR row; rows are emitted in node order.
class GeneAdjacencyWriter {
 public:
  GeneAdjacencyWriter(GeneIndex& gene, std::size_t node_count,
                      std::size_t directed_edges)
      : gene_(gene) {
    gene_.adjacency_offsets.clear();
    gene_.adjacency.clear();
    gene_.adjacency_offsets.reserve(node_count + 1u);
    gene_.adjacency.reserve(directed_edges);
    gene_.adjacency_offsets.push_back(0);
  }

  void push(int neighbor) { gene_.adjacency.push_back(neighbor); }

  void end_node() {
    if (gene_.adjacency.size() >
        static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::length_error("gene adjacency exceeds its half-edge domain");
    }
    gene_.adjacency_offsets.push_back(
        static_cast<int>(gene_.adjacency.size()));
  }

 private:
  GeneIndex& gene_;
};

/// A tree on n nodes has 2(n - 1) half-edges.
std::size_t tree_directed_edges(std::size_t node_count) {
  return node_count == 0 ? 0 : 2u * (node_count - 1u);
}

/// Mutable tree adjacency with stable incidence order and O(1) edge unlink.
/// Each node's half-edges form a doubly linked list; a new edge is appended
/// at the end of both endpoints' lists. Construction rejects input that is
/// not a symmetric, connected tree.
class MutableTreeAdjacency {
 public:
  explicit MutableTreeAdjacency(const GeneIndex& gene,
                                std::size_t extra_directed_edges = 0)
      : nodes_(gene.nodes.size()) {
    constexpr std::size_t kMaxIndexedEdges =
        static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (gene.nodes.size() > kMaxIndexedEdges) {
      throw std::length_error("gene adjacency exceeds its integer node domain");
    }
    validate_gene_adjacency_shape(gene);
    const std::size_t directed_edges = gene.adjacency.size();
    if (extra_directed_edges > kMaxIndexedEdges - directed_edges) {
      throw std::length_error("gene adjacency exceeds its half-edge domain");
    }
    edges_.reserve(directed_edges + extra_directed_edges);
    for (int node = 0; node < static_cast<int>(gene.nodes.size()); ++node) {
      for (int neighbor : gene_neighbors(gene, node)) {
        append_directed(node, neighbor);
      }
    }
    pair_original_edges();
  }

  void add_node() { nodes_.push_back(NodeLinks{}); }

  int degree(int node) const {
    return nodes_[static_cast<std::size_t>(node)].degree;
  }

  int first_edge(int node) const {
    return nodes_[static_cast<std::size_t>(node)].first;
  }

  int next_edge(int edge) const {
    return edges_[static_cast<std::size_t>(edge)].next;
  }

  int from(int edge) const {
    return edges_[static_cast<std::size_t>(edge)].from;
  }

  int neighbor(int edge) const {
    return edges_[static_cast<std::size_t>(edge)].to;
  }

  template <typename Fn>
  void for_each_edge(int node, Fn&& fn) const {
    for (int edge = first_edge(node); edge >= 0; edge = next_edge(edge)) {
      fn(edge);
    }
  }

  int add_edge(int a, int b) {
    const int ab = append_directed(a, b);
    const int ba = append_directed(b, a);
    edges_[static_cast<std::size_t>(ab)].reverse = ba;
    edges_[static_cast<std::size_t>(ba)].reverse = ab;
    return ab;
  }

  void remove_edge(int edge) {
    HalfEdge& forward = edges_[static_cast<std::size_t>(edge)];
    if (!forward.active) return;
    const int reverse = forward.reverse;
    if (reverse < 0 ||
        reverse >= static_cast<int>(edges_.size()) ||
        !edges_[static_cast<std::size_t>(reverse)].active) {
      throw std::invalid_argument("gene adjacency must be a symmetric tree");
    }
    unlink_directed(edge);
    unlink_directed(reverse);
  }

 private:
  struct NodeLinks {
    int first = -1;
    int last = -1;
    int degree = 0;
  };

  struct HalfEdge {
    int from = -1;
    int to = -1;
    int reverse = -1;
    int previous = -1;
    int next = -1;
    bool active = true;
  };

  int append_directed(int from, int to) {
    if (edges_.size() >=
        static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::length_error("gene adjacency exceeds its half-edge domain");
    }
    const int edge = static_cast<int>(edges_.size());
    NodeLinks& links = nodes_[static_cast<std::size_t>(from)];
    edges_.push_back(HalfEdge{from, to, -1, links.last, -1, true});
    if (links.last >= 0) {
      edges_[static_cast<std::size_t>(links.last)].next = edge;
    } else {
      links.first = edge;
    }
    links.last = edge;
    ++links.degree;
    return edge;
  }

  void unlink_directed(int edge) {
    HalfEdge& current = edges_[static_cast<std::size_t>(edge)];
    NodeLinks& links = nodes_[static_cast<std::size_t>(current.from)];
    if (current.previous >= 0) {
      edges_[static_cast<std::size_t>(current.previous)].next = current.next;
    } else {
      links.first = current.next;
    }
    if (current.next >= 0) {
      edges_[static_cast<std::size_t>(current.next)].previous = current.previous;
    } else {
      links.last = current.previous;
    }
    current.active = false;
    --links.degree;
  }

  /// Pairs every half-edge with its reverse during a DFS. Each discovered
  /// node's list is scanned once for the reverse, so pairing is linear. An
  /// unpaired half-edge means a cycle or asymmetry; a second DFS root means a
  /// forest.
  void pair_original_edges() {
    const int n = static_cast<int>(nodes_.size());
    std::vector<unsigned char> discovered(static_cast<std::size_t>(n), 0);
    std::vector<int> stack;
    stack.reserve(static_cast<std::size_t>(n));
    int component_count = 0;
    for (int start = 0; start < n; ++start) {
      if (discovered[static_cast<std::size_t>(start)] != 0) continue;
      ++component_count;
      discovered[static_cast<std::size_t>(start)] = 1;
      stack.push_back(start);
      while (!stack.empty()) {
        const int node = stack.back();
        stack.pop_back();
        for (int edge = first_edge(node); edge >= 0; edge = next_edge(edge)) {
          HalfEdge& forward = edges_[static_cast<std::size_t>(edge)];
          if (forward.reverse >= 0) continue;
          const int child = forward.to;
          if (child < 0 || child >= n) {
            throw std::invalid_argument("gene adjacency contains an invalid node");
          }
          if (discovered[static_cast<std::size_t>(child)] != 0) continue;

          int reverse = -1;
          for (int candidate = first_edge(child); candidate >= 0;
               candidate = next_edge(candidate)) {
            const HalfEdge& back = edges_[static_cast<std::size_t>(candidate)];
            if (back.to == node && back.reverse < 0) {
              reverse = candidate;
              break;
            }
          }
          if (reverse < 0) {
            throw std::invalid_argument("gene adjacency must be symmetric");
          }
          forward.reverse = reverse;
          edges_[static_cast<std::size_t>(reverse)].reverse = edge;
          discovered[static_cast<std::size_t>(child)] = 1;
          stack.push_back(child);
        }
      }
    }
    for (const HalfEdge& edge : edges_) {
      if (edge.reverse < 0) {
        throw std::invalid_argument("gene adjacency must be acyclic and symmetric");
      }
    }
    if (component_count > 1) {
      throw std::invalid_argument("gene adjacency must be connected");
    }
  }

  std::vector<NodeLinks> nodes_;
  std::vector<HalfEdge> edges_;
};

}  // namespace

/// Removes a keyless node of degree one (the unary Newick parse root, as in
/// `((A,C));`), and its neighbor in turn if that becomes one, so that no
/// dangling vertex lengthens a path; then splices out every degree-two node
/// that carries no leaf key. The edge that replaces a spliced node goes to the
/// end of both neighbors' incidence lists.
/// When any node is removed, the kept nodes are renumbered in order,
/// taxon_to_node is rebuilt, and the gene is rooted at its first keyless node
/// of degree at least two.
void suppress_gene_degree_two(GeneIndex& gene, int taxon_count) {
  if (gene.nodes.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("gene index exceeds its integer node domain");
  }
  const int n = static_cast<int>(gene.nodes.size());
  if (n == 0) return;
  MutableTreeAdjacency adjacency(gene, static_cast<std::size_t>(n) * 2);
  std::vector<unsigned char> active(static_cast<std::size_t>(n), 1);

  for (int i = 0; i < n; ++i) {
    int stub = i;
    while (active[static_cast<std::size_t>(stub)] != 0 &&
           gene.nodes[static_cast<std::size_t>(stub)].key == -1 &&
           adjacency.degree(stub) == 1) {
      const int edge = adjacency.first_edge(stub);
      const int next = adjacency.neighbor(edge);
      adjacency.remove_edge(edge);
      active[static_cast<std::size_t>(stub)] = 0;
      stub = next;
    }
  }
  for (int i = 0; i < n; ++i) {
    if (active[static_cast<std::size_t>(i)] == 0 ||
        gene.nodes[static_cast<std::size_t>(i)].key != -1 ||
        adjacency.degree(i) != 2) {
      continue;
    }
    const int first = adjacency.first_edge(i);
    const int second = adjacency.next_edge(first);
    const int a = adjacency.neighbor(first);
    const int b = adjacency.neighbor(second);
    adjacency.remove_edge(first);
    adjacency.remove_edge(second);
    active[static_cast<std::size_t>(i)] = 0;
    adjacency.add_edge(a, b);
  }

  std::vector<int> remap(n, -1);
  int next = 0;
  for (int i = 0; i < n; ++i) {
    if (active[static_cast<std::size_t>(i)] != 0) remap[i] = next++;
  }
  if (next == n) return;

  std::vector<GeneIndex::Node> compact(next);
  std::vector<double> compact_mass;
  if (!gene.leaf_mass.empty()) compact_mass.resize(next);
  GeneIndex compact_adjacency;
  GeneAdjacencyWriter writer(compact_adjacency, static_cast<std::size_t>(next),
                             tree_directed_edges(static_cast<std::size_t>(next)));
  for (int i = 0; i < n; ++i) {
    const int ni = remap[i];
    if (ni < 0) continue;
    compact[ni] = gene.nodes[i];
    if (!compact_mass.empty()) compact_mass[ni] = gene.leaf_mass[i];
    adjacency.for_each_edge(i, [&](int edge) {
      const int nb = adjacency.neighbor(edge);
      if (nb >= 0 && nb < n && remap[nb] >= 0) writer.push(remap[nb]);
    });
    writer.end_node();
  }
  gene.nodes = std::move(compact);
  gene.leaf_mass = std::move(compact_mass);
  gene.adjacency_offsets = std::move(compact_adjacency.adjacency_offsets);
  gene.adjacency = std::move(compact_adjacency.adjacency);
  gene.taxon_to_node.assign(taxon_count, -1);
  for (int i = 0; i < static_cast<int>(gene.nodes.size()); ++i) {
    const int taxon = gene.nodes[i].taxon();
    if (taxon >= 0 && taxon < taxon_count) gene.taxon_to_node[taxon] = i;
  }
  gene.root = -1;
  for (int i = 0; i < static_cast<int>(gene.nodes.size()); ++i) {
    if (gene.nodes[i].key == -1 && gene_degree(gene, i) >= 2) {
      gene.root = i;
      break;
    }
  }
  if (gene.root < 0 && !gene.nodes.empty()) gene.root = 0;
}

/// Indexes a parsed gene tree: checks leaf labels, builds the CSR adjacency in
/// the tree's incidence order and suppresses degree-two nodes. Every duplicate
/// leaf of a taxon is kept; taxon_to_node names its last occurrence.
GeneIndex index_gene(const GeneTree& gene, const TaxonRegistry& taxa) {
  if (gene.nodes.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("gene tree exceeds its integer node domain");
  }
  if (gene.nodes.empty()) return canonical_empty_gene_index();
  if (gene.root < 0 ||
      gene.root >= static_cast<int>(gene.nodes.size())) {
    throw std::invalid_argument("gene tree has an invalid root");
  }
  for (int i = 0; i < static_cast<int>(gene.nodes.size()); ++i) {
    const Tree::Node& node = gene.nodes[static_cast<std::size_t>(i)];
    const bool labeled = !node.label.empty();
    const std::size_t degree = node.neighbors.size();
    if (labeled && degree > 1u) {
      throw std::invalid_argument("gene tree has a labeled internal node");
    }
    // Newick permits a unary, unlabeled parse root (for example `(A);`);
    // suppress_gene_degree_two removes it. Every other degree-zero/one vertex
    // must be a physical labeled leaf.
    if (!labeled && degree <= 1u && !(i == gene.root && degree == 1u)) {
      throw std::invalid_argument("gene tree has an unlabeled leaf");
    }
  }
  GeneIndex out;
  out.nodes.resize(gene.nodes.size());
  out.root = gene.root;
  out.taxon_to_node.assign(taxa.size(), -1);
  std::size_t directed_edges = 0;
  for (const Tree::Node& node : gene.nodes) {
    directed_edges += node.neighbors.size();
  }
  GeneAdjacencyWriter writer(out, gene.nodes.size(), directed_edges);
  for (int i = 0; i < static_cast<int>(gene.nodes.size()); ++i) {
    for (int neighbor : gene.nodes[i].neighbors) writer.push(neighbor);
    writer.end_node();
    if (!gene.nodes[i].label.empty()) {
      const int taxon = taxa.id_for(gene.nodes[i].label);
      out.nodes[i].set_taxon(taxon);
      out.taxon_to_node[taxon] = i;
    }
  }
  suppress_gene_degree_two(out, taxa.size());
  out.root_traversal = traversal(out, out.root);
  // The reference root is the leaf of the smallest taxon present (its last
  // occurrence). A split key describes the side away from this leaf, so keys
  // do not depend on the Newick parse root, even when taxon 0 is absent.
  out.reference_root = out.root;
  for (int taxon = 0; taxon < static_cast<int>(out.taxon_to_node.size());
       ++taxon) {
    if (out.taxon_to_node[static_cast<std::size_t>(taxon)] >= 0) {
      out.reference_root =
          out.taxon_to_node[static_cast<std::size_t>(taxon)];
      break;
    }
  }
  return out;
}

void set_gene_adjacency(GeneIndex& gene,
                        const std::vector<std::vector<int>>& lists) {
  std::size_t directed_edges = 0;
  for (const std::vector<int>& list : lists) directed_edges += list.size();
  GeneAdjacencyWriter writer(gene, lists.size(), directed_edges);
  for (const std::vector<int>& list : lists) {
    for (int neighbor : list) writer.push(neighbor);
    writer.end_node();
  }
}

/// Iterative DFS from root (gene.root when root is negative) in incidence
/// order. parent[root] is root itself.
GeneTraversal traversal(const GeneIndex& gene, int root) {
  GeneTraversal tr;
  const int n = static_cast<int>(gene.nodes.size());
  tr.parent.assign(n, -1);
  if (n == 0) return tr;
  tr.preorder.reserve(static_cast<std::size_t>(n));
  root = root >= 0 ? root : gene.root;
  std::vector<std::pair<int, int>> stack{{root, 0}};
  tr.parent[root] = root;
  while (!stack.empty()) {
    const int node = stack.back().first;
    int& next = stack.back().second;
    if (next == 0) tr.preorder.push_back(node);
    const std::span<const int> neighbors = gene_neighbors(gene, node);
    if (next < static_cast<int>(neighbors.size())) {
      const int nb = neighbors[static_cast<std::size_t>(next++)];
      if (nb == tr.parent[node]) continue;
      tr.parent[nb] = node;
      stack.push_back({nb, 0});
      continue;
    }
    stack.pop_back();
  }
  return tr;
}

/// The tree induced on the leaves that `local_of(node)` maps to a local
/// index: other leaves are pruned and unlabelled nodes of degree at most two
/// are suppressed. Kept leaves carry their local index and mass, and every
/// duplicate leaf is kept. Pruning a node reschedules only its neighbor, so
/// the restriction is linear in the input gene. The result is rooted at its
/// first internal node (node 0 when it has none), which is also its
/// reference root.
template <typename LocalOf>
GeneIndex restrict_gene(const GeneIndex& original, LocalOf&& local_of) {
  struct MNode {
    int local = -1;
    double mass = 0.0;
    bool is_leaf = false;
    bool deleted = false;
    bool queued = false;
  };

  if (original.nodes.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("gene restriction exceeds its integer node domain");
  }
  const int n = static_cast<int>(original.nodes.size());
  MutableTreeAdjacency adjacency(
      original, static_cast<std::size_t>(n) * 2 + 8);
  if (n == 0) return canonical_empty_gene_index();
  std::vector<MNode> nodes(n);
  int mapped_leaf_count = 0;
  for (int i = 0; i < n; ++i) {
    const int local = local_of(original.nodes[static_cast<std::size_t>(i)]);
    if (local >= 0) {
      nodes[i].is_leaf = true;
      nodes[i].local = local;
      nodes[i].mass = gene_leaf_mass(original, i);
      ++mapped_leaf_count;
    }
  }
  if (mapped_leaf_count == 0) return canonical_empty_gene_index();

  std::vector<int> queue;
  queue.reserve(nodes.size());
  auto schedule = [&](int node) {
    if (node < 0 || node >= static_cast<int>(nodes.size())) return;
    MNode& current = nodes[static_cast<std::size_t>(node)];
    if (current.deleted || current.is_leaf || current.queued) return;
    current.queued = true;
    queue.push_back(node);
  };

  for (int i = 0; i < static_cast<int>(nodes.size()); ++i) schedule(i);
  for (std::size_t qpos = 0; qpos < queue.size(); ++qpos) {
    const int node = queue[qpos];
    MNode& current = nodes[static_cast<std::size_t>(node)];
    current.queued = false;
    if (current.deleted || current.is_leaf) continue;

    const int degree = adjacency.degree(node);
    if (degree == 2) {
      const int first = adjacency.first_edge(node);
      const int second = adjacency.next_edge(first);
      const int a = adjacency.neighbor(first);
      const int b = adjacency.neighbor(second);
      adjacency.remove_edge(first);
      adjacency.remove_edge(second);
      current.deleted = true;
      adjacency.add_edge(a, b);
    } else if (degree <= 1) {
      int neighbor = -1;
      if (degree == 1) {
        const int edge = adjacency.first_edge(node);
        neighbor = adjacency.neighbor(edge);
        adjacency.remove_edge(edge);
      }
      current.deleted = true;
      schedule(neighbor);
    }
  }

  std::vector<int> old_to_new(nodes.size(), -1);
  int active = 0;
  for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
    if (!nodes[i].deleted) old_to_new[i] = active++;
  }
  if (active == 0) return canonical_empty_gene_index();

  GeneIndex result;
  result.nodes.resize(active);
  result.root = 0;
  GeneAdjacencyWriter writer(result, static_cast<std::size_t>(active),
                             tree_directed_edges(static_cast<std::size_t>(active)));
  for (int i = 0; i < static_cast<int>(nodes.size()); ++i) {
    const int ni = old_to_new[i];
    if (ni < 0) continue;
    adjacency.for_each_edge(i, [&](int edge) {
      const int nb = adjacency.neighbor(edge);
      const int nn = nb >= 0 && nb < static_cast<int>(old_to_new.size())
                         ? old_to_new[nb]
                         : -1;
      if (nn >= 0) writer.push(nn);
    });
    writer.end_node();
    if (nodes[i].is_leaf) {
      result.nodes[ni].set_local(nodes[i].local);
      set_gene_leaf_mass(result, ni, nodes[i].mass);
    }
  }
  for (int i = 0; i < static_cast<int>(result.nodes.size()); ++i) {
    if (result.nodes[i].local() < 0 && gene_degree(result, i) >= 2) {
      result.root = i;
      break;
    }
  }
  result.reference_root = result.root;
  result.root_traversal = traversal(result, result.root);
  return result;
}

GeneIndex restrict_gene_by_taxon(const GeneIndex& gene,
                                 const std::vector<int>& taxon_to_local) {
  return restrict_gene(gene, [&](const GeneIndex::Node& node) {
    const int taxon = node.taxon();
    return taxon >= 0 && taxon < static_cast<int>(taxon_to_local.size())
               ? taxon_to_local[static_cast<std::size_t>(taxon)]
               : -1;
  });
}

GeneIndex restrict_gene_by_local(const GeneIndex& gene,
                                 const std::vector<int>& local_to_local) {
  return restrict_gene(gene, [&](const GeneIndex::Node& node) {
    const int local = node.local();
    return local >= 0 && local < static_cast<int>(local_to_local.size())
               ? local_to_local[static_cast<std::size_t>(local)]
               : -1;
  });
}

}  // namespace qraft::internal
