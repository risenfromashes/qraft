#pragma once

#include "solver/solver.hpp"

#include <array>

namespace qraft::internal {

/// An internal vertex of a solved side with the number of active taxa in each
/// of its three arms. Arms are numbered by the vertex's active neighbors in
/// adjacency order, as are VertexScores.
struct VertexColorPlan {
  int vertex_node = -1;
  std::array<int, 3> branch_size{0, 0, 0};
};

/// A vertex's regraft score for each arm, for one gene or summed over genes.
/// Of the quartets with one taxon in each arm and one on the other side, it
/// counts those whose gene topology pairs the other side with that arm, minus
/// twice those that pair it with another arm.
struct VertexScores {
  std::array<double, 3> raw{0.0, 0.0, 0.0};
};

/// The scoring input for regrafting the other side into one solved side.
struct RegraftJob {
  const SolverTree* child_tree = nullptr;
  const std::vector<int>* active_taxa_ids = nullptr;   ///< the side's taxa
  const std::vector<int>* context_taxa_ids = nullptr;  ///< the other side
  QueryTree query;  ///< the side over its active taxa
  std::vector<VertexColorPlan> plans;
  std::vector<int> source_to_plan;  ///< child-tree node -> plan index or -1
  /// The arm number of each of a query vertex's three branches, in
  /// vertex_arm_coloring order; -1 for a vertex that is not scored.
  std::vector<std::array<int, 3>> direction_by_query_vertex;
  /// The owner of each active taxon (by position in active_taxa_ids), the
  /// owners in query leaf order, and the other side's owners. All three hold
  /// taxon ids until the scoring pass assigns owners.
  std::vector<int> query_active_to_owner;
  std::vector<int> query_leaf_owner_order;
  std::vector<int> context_owner_ids;
  /// The most scored query vertices that share one plan.
  int max_vertices_per_plan = 0;
  /// Exact integer totals, when requested: valid only if every gene's
  /// counting stayed in the exactly representable range.
  bool exact_requested = false;
  bool exact_valid = false;
  std::vector<std::array<RegraftCount, 3>> exact_totals;
  std::vector<VertexScores> totals;  ///< per plan, over all scoring genes
};

}  // namespace qraft::internal
