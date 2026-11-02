#include "solver/solver.hpp"

#include <algorithm>
#include <stdexcept>

namespace qraft::internal {

TaxonIndex make_taxon_index(const std::vector<int>& taxa,
                                   int taxon_count) {
  TaxonIndex index;
  std::vector<std::pair<int, int>> first_positions;
  first_positions.reserve(taxa.size());
  for (int position = 0; position < static_cast<int>(taxa.size()); ++position) {
    const int taxon = taxa[static_cast<std::size_t>(position)];
    if (taxon >= 0 && taxon < taxon_count) {
      first_positions.push_back({taxon, position});
    }
  }
  std::sort(first_positions.begin(), first_positions.end());
  first_positions.erase(
      std::unique(first_positions.begin(), first_positions.end(),
                  [](const auto& a, const auto& b) {
                    return a.first == b.first;
                  }),
      first_positions.end());
  std::sort(first_positions.begin(), first_positions.end(),
            [](const auto& a, const auto& b) {
              return a.second < b.second;
            });

  index.taxa.reserve(first_positions.size());
  index.by_taxon.reserve(first_positions.size());
  for (const auto& [taxon, unused_position] : first_positions) {
    (void)unused_position;
    const int local = static_cast<int>(index.taxa.size());
    index.taxa.push_back(taxon);
    index.by_taxon.push_back({taxon, local});
  }
  std::sort(index.by_taxon.begin(), index.by_taxon.end());
  stats_add(run_stats().active_taxa_entries,
                        static_cast<long long>(taxa.size()));
  stats_add(run_stats().active_taxa_index_entries,
                        static_cast<long long>(index.by_taxon.size()));
  return index;
}

int active_index(const ActiveTaxa& active, int taxon) {
  return active.active_by_taxon.find(taxon);
}

ActiveTaxa make_active_taxa(const std::vector<int>& taxa, int taxon_count) {
  ActiveTaxa active;
  active.active_by_taxon = make_taxon_index(taxa, taxon_count);
  active.taxa = active.active_by_taxon.taxa;
  return active;
}

void note_gene_set_use(const std::vector<GeneIndex>& genes) {
  long long source_nodes = 0;
  for (const GeneIndex& gene : genes) {
    source_nodes += static_cast<long long>(gene.nodes.size());
  }
  stats_add(run_stats().gene_set_scoring_calls);
  stats_add(run_stats().gene_set_scoring_nodes,
                        source_nodes);
}

SolverTree base_tree(const std::vector<int>& taxa) {
  SolverTree tree;
  if (taxa.empty()) return tree;
  if (taxa.size() == 1) {
    tree.nodes.push_back(SolverTree::Node{{}, taxa[0], true});
    return tree;
  }
  if (taxa.size() == 2) {
    tree.nodes.push_back(SolverTree::Node{{1}, taxa[0], true});
    tree.nodes.push_back(SolverTree::Node{{0}, taxa[1], true});
    return tree;
  }
  tree.nodes.push_back(SolverTree::Node{{}, -1, true});
  for (int taxon : taxa) {
    const int leaf = static_cast<int>(tree.nodes.size());
    tree.nodes.push_back(SolverTree::Node{{0}, taxon, true});
    tree.nodes[0].neighbors.push_back(leaf);
  }
  return tree;
}

std::vector<int> active_neighbors(const SolverTree& tree, int node) {
  std::vector<int> out;
  if (node < 0 || node >= static_cast<int>(tree.nodes.size()) ||
      !tree.nodes[static_cast<std::size_t>(node)].active) {
    return out;
  }
  for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
    if (nb >= 0 && nb < static_cast<int>(tree.nodes.size()) &&
        tree.nodes[static_cast<std::size_t>(nb)].active) {
      out.push_back(nb);
    }
  }
  return out;
}

bool has_solver_edge(const SolverTree& tree, int a, int b) {
  if (a < 0 || b < 0 || a >= static_cast<int>(tree.nodes.size()) ||
      b >= static_cast<int>(tree.nodes.size())) {
    return false;
  }
  const auto& nbs = tree.nodes[static_cast<std::size_t>(a)].neighbors;
  return tree.nodes[static_cast<std::size_t>(a)].active &&
         tree.nodes[static_cast<std::size_t>(b)].active &&
         std::find(nbs.begin(), nbs.end(), b) != nbs.end();
}

void add_solver_edge(SolverTree& tree, int a, int b) {
  if (a < 0 || b < 0 || a >= static_cast<int>(tree.nodes.size()) ||
      b >= static_cast<int>(tree.nodes.size())) {
    return;
  }
  auto& an = tree.nodes[static_cast<std::size_t>(a)].neighbors;
  if (std::find(an.begin(), an.end(), b) == an.end()) an.push_back(b);
  auto& bn = tree.nodes[static_cast<std::size_t>(b)].neighbors;
  if (std::find(bn.begin(), bn.end(), a) == bn.end()) bn.push_back(a);
}

void remove_solver_edge(SolverTree& tree, int a, int b) {
  if (a < 0 || b < 0 || a >= static_cast<int>(tree.nodes.size()) ||
      b >= static_cast<int>(tree.nodes.size())) {
    return;
  }
  auto& an = tree.nodes[static_cast<std::size_t>(a)].neighbors;
  an.erase(std::remove(an.begin(), an.end(), b), an.end());
  auto& bn = tree.nodes[static_cast<std::size_t>(b)].neighbors;
  bn.erase(std::remove(bn.begin(), bn.end(), a), bn.end());
}

void suppress_solver_degree_two(SolverTree& tree) {
  bool changed = true;
  while (changed) {
    changed = false;
    for (int i = 0; i < static_cast<int>(tree.nodes.size()); ++i) {
      if (!tree.nodes[static_cast<std::size_t>(i)].active ||
          tree.nodes[static_cast<std::size_t>(i)].taxon >= 0) {
        continue;
      }
      const std::vector<int> nbs = active_neighbors(tree, i);
      if (nbs.size() != 2) continue;
      remove_solver_edge(tree, i, nbs[0]);
      remove_solver_edge(tree, i, nbs[1]);
      tree.nodes[static_cast<std::size_t>(i)].active = false;
      tree.nodes[static_cast<std::size_t>(i)].neighbors.clear();
      add_solver_edge(tree, nbs[0], nbs[1]);
      changed = true;
    }
  }
}

namespace {

/// Checks that `edge` is an edge of `tree`, or the sentinel of a one-vertex
/// or empty tree.
void validate_regraft_edge(const SolverTree& tree,
                              const RegraftEdge& edge) {
  if (edge.a >= 0 && edge.b >= 0) {
    if (has_solver_edge(tree, edge.a, edge.b) &&
        has_solver_edge(tree, edge.b, edge.a)) {
      return;
    }
    throw std::logic_error("invalid regraft edge");
  }

  int active_count = 0;
  int singleton = -1;
  for (int a = 0; a < static_cast<int>(tree.nodes.size()); ++a) {
    if (tree.nodes[static_cast<std::size_t>(a)].active) {
      ++active_count;
      singleton = a;
    }
  }

  if (active_count == 0) {
    if (edge.a == -1 && edge.b == -1) return;
    throw std::logic_error("nonempty regraft endpoint for an empty tree");
  }
  if (active_count == 1) {
    if (edge.a == singleton && edge.b == -1) return;
    throw std::logic_error("invalid regraft endpoint for a one-vertex tree");
  }
  throw std::logic_error("invalid regraft edge sentinel");
}

/// The vertex of `tree` that connects to the other side at `edge`.
int subdivide_regraft_edge(SolverTree& tree, const RegraftEdge& edge) {
  if (edge.a < 0) {
    if (edge.b != -1) {
      throw std::logic_error("invalid regraft endpoint for an empty tree");
    }
    const int node = static_cast<int>(tree.nodes.size());
    tree.nodes.push_back(SolverTree::Node{{}, -1, true});
    return node;
  }
  if (edge.b == -1) {
    // A one-vertex tree connects at its vertex. A new connector there would
    // be an active unlabelled vertex of degree two: the tree would not be
    // binary, and regrafting, which scores the three arms of each internal
    // vertex, would give that vertex's two edges incomplete scores, which
    // could let a wrong edge win.
    return edge.a;
  }
  if (!has_solver_edge(tree, edge.a, edge.b)) {
    throw std::logic_error("regraft edge disappeared before connecting");
  }
  remove_solver_edge(tree, edge.a, edge.b);
  const int node = static_cast<int>(tree.nodes.size());
  tree.nodes.push_back(SolverTree::Node{{}, -1, true});
  add_solver_edge(tree, node, edge.a);
  add_solver_edge(tree, node, edge.b);
  return node;
}

}  // namespace

SolverTree connect_at_edges(SolverTree left,
                    SolverTree right,
                    RegraftEdge left_edge,
                    RegraftEdge right_edge) {
  validate_regraft_edge(left, left_edge);
  validate_regraft_edge(right, right_edge);
  const int right_offset = static_cast<int>(left.nodes.size());
  for (const SolverTree::Node& node : right.nodes) {
    SolverTree::Node copy;
    copy.taxon = node.taxon;
    copy.active = node.active;
    copy.neighbors.reserve(node.neighbors.size());
    for (int nb : node.neighbors) copy.neighbors.push_back(nb + right_offset);
    left.nodes.push_back(std::move(copy));
  }
  if (right_edge.a >= 0) right_edge.a += right_offset;
  if (right_edge.b >= 0) right_edge.b += right_offset;
  const int left_connector = subdivide_regraft_edge(left, left_edge);
  const int right_connector = subdivide_regraft_edge(left, right_edge);
  add_solver_edge(left, left_connector, right_connector);
  return left;
}

}  // namespace qraft::internal
