#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace qraft::internal {

// The key of a gene split: the size of its side away from the gene's
// reference root and that side's exact identifier in h[0], below 2^32, with
// h[1..3] zero (exact_split_ids.hpp). Equal sides have equal keys and
// distinct sides distinct keys.
struct SplitKey256 {
  std::uint32_t size = 0;
  std::array<std::uint64_t, 4> h{};

  bool operator==(const SplitKey256&) const = default;
};

// Lexicographic order on (size, h[0], h[1], h[2], h[3]).
bool split_key_less(const SplitKey256& left, const SplitKey256& right);

// One eligible gene edge: the key of the side away from the gene's reference
// root, the gene, the edge's preorder ordinal in that gene, the gene node
// below the edge, and one taxon on that side. Sorting by (key, gene, edge)
// makes the first record of every (key, gene) run the deterministic
// representative.
struct SplitOccurrence256 {
  SplitKey256 key;
  std::uint32_t gene = 0;
  std::uint32_t edge = 0;
  std::uint32_t node = 0;
  std::uint32_t one_taxon = 0;
};

// One record per distinct key after per-gene deduplication. split_id is the
// zero-based rank of key in SplitKey256 lexicographic order. count is the
// number of genes that contain the split. first_tree and first_edge are the
// smallest (gene, edge) occurrence, and rep_tree, rep_node and one_taxon come
// from that same occurrence.
struct SplitAggregate256 {
  SplitKey256 key;
  std::uint32_t split_id = 0;
  std::uint32_t count = 0;
  std::uint32_t first_tree = 0;
  std::uint32_t first_edge = 0;
  std::uint32_t rep_tree = 0;
  std::uint32_t rep_node = 0;
  std::uint32_t one_taxon = 0;
};

// Pass, item and record counts of one aggregation, reported to the
// complexity tests.
struct SplitRadixStats {
  long long occurrence_passes = 0;
  long long occurrence_items = 0;
  long long priority_passes = 0;
  long long priority_items = 0;
  long long emitted_occurrences = 0;
  long long deduplicated_gene_splits = 0;
  long long distinct_splits = 0;

  long long total_passes() const {
    return occurrence_passes + priority_passes;
  }
  long long total_items() const {
    return occurrence_items + priority_items;
  }
};

// With 16-bit digits, the priority key (count, first_tree, first_edge,
// split_id) takes 4 * 2 passes.
inline constexpr int kSplitPriorityRadixPasses = 8;

// Comparator form of the MRE priority order; the tests use it as the oracle
// for radix_sort_split_priority256.
bool split_priority_less_by_key(const SplitAggregate256& left,
                                const SplitAggregate256& right);

// Sorts the records of `chunks`, taken in order, by (key, gene, edge),
// removes repeated observations of a key in one gene, aggregates support, and
// assigns dense key-rank split ids. The returned vector is in key/split-id
// order. The chunks are freed. The passes and the aggregation run on up to
// one part per worker, and the result does not depend on the workers.
// Every key must have h = {id, 0, 0, 0} with id < 2^32 (SplitKey256);
// another key throws std::logic_error.
std::vector<SplitAggregate256> aggregate_split_occurrences256(
    std::vector<std::vector<SplitOccurrence256>>& chunks,
    SplitRadixStats& stats);

// The same for one vector of records, which is left sorted by
// (key, gene, edge).
std::vector<SplitAggregate256> aggregate_split_occurrences256(
    std::vector<SplitOccurrence256>& occurrences,
    SplitRadixStats& stats);

// Produces the MRE priority order
// (count descending, first gene, first edge, key) using split_id for the final
// key tie-break. split_id is order-isomorphic to the complete key. The passes
// run in parallel as in aggregate_split_occurrences256.
void radix_sort_split_priority256(
    std::vector<SplitAggregate256>& aggregates,
    SplitRadixStats& stats);

}  // namespace qraft::internal
