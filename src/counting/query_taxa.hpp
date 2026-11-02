#pragma once

#include "solver/solver.hpp"

#include <cstdint>
#include <utility>
#include <vector>

namespace qraft::internal {

// A query's owners are its active taxa, numbered in query order. A gene leaf
// maps to an owner through its taxon (input genes) or through its local index
// (restricted genes); every leaf that maps nowhere has owner -1. Several
// leaves may map to one owner (a duplicated taxon), and each stays a separate
// weighted leaf.

// Missing: the map is not built and every lookup misses.
enum class OwnerKeyMode : std::uint8_t {
  Missing = 0,
  Dense = 1,
  Sorted = 2,
};

/// A key -> owner map: a dense table when the key domain is linear in the
/// owner count, and otherwise a sorted list.
struct OwnerKeyMap {
  OwnerKeyMode mode = OwnerKeyMode::Missing;
  std::vector<int> dense;
  std::vector<std::pair<int, int>> sorted;

  /// The owner of `key`, or -1.
  int find(int key) const;
  /// One more than the largest owner stored, or 0.
  int inferred_owner_count() const;
};

struct OwnerProjection {
  OwnerKeyMap taxon_keys;
  OwnerKeyMap local_keys;
  // Direct taxon -> owner table, read when taxon_keys is Missing.
  std::vector<int> taxon_to_owner;
  // Owners that currently count; empty means all of them.
  std::vector<std::uint8_t> live_owner;
  std::vector<int> live_taxa;
  int owner_count = 0;
};

/// The projection of the scoring genes of `ctx` onto the taxa of
/// `local_by_taxon`, whose order numbers the owners. Local keys are built
/// only when the scoring genes are a restricted gene set.
OwnerProjection make_owner_projection(const SolverContext& ctx,
                                      const TaxonIndex& local_by_taxon);

int owner_for_taxon(const OwnerProjection& owners, int taxon);
/// The live owner of a gene leaf, or -1.
int owner_for_gene_leaf(const GeneIndex::Node& node,
                        const OwnerProjection& owners);
/// Whether `gene` has at least four leaves with an owner, the fewest that
/// form a quartet. A gene with fewer adds nothing to any channel, and the
/// scorers skip it: a gene restricted to two of a subproblem's taxa is a
/// single edge, which no HDT composition represents.
bool gene_has_owner_quartet(const GeneIndex& gene,
                            const OwnerProjection& owners);

}  // namespace qraft::internal
