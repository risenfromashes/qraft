#pragma once

#include <vector>

#include "qraft/options.hpp"
#include "qraft/tree.hpp"

namespace qraft {

/// Infers a species tree from `genes` by the QRAFT pipeline: the MRE starting
/// tree, the quartet-guided recursive divide with contextual regrafting, and
/// NNI refinement, each configurable through `options`.
///
/// The taxa are the distinct leaf labels. Gene trees are treated as unrooted;
/// a gene may lack taxa, in which case the bipartition test ranks taxa along
/// the anchor path instead of measuring path lengths, or repeat one, in which
/// case the test is skipped. Scores are summed in gene order, so reordering the genes may
/// change near-tied choices. For the same genes and options the result does
/// not depend on the thread count. The returned tree is unrooted (`root` is
/// -1) and has one leaf per taxon.
///
/// Calls with progress Off and no `stats_json` may run concurrently when they
/// use the same thread count; any other call waits to run alone.
///
/// Throws std::invalid_argument for an empty gene list, a malformed gene or
/// starting tree, or bipartition test options out of range, and
/// std::runtime_error when the starting tree or statistics file cannot be
/// opened.
SpeciesTree infer_species_tree(std::vector<GeneTree> genes, const Options& options = {});

/// Throws std::invalid_argument when `options` hold a value that
/// infer_species_tree rejects (bipartition-test and SPR settings), so a
/// caller can check them before reading any input. infer_species_tree runs
/// the same checks.
void validate_options(const Options& options);

}  // namespace qraft
