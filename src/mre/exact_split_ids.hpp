#pragma once

#include "mre/split_records.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace qraft::internal {

struct GeneIndex;
struct GeneTraversal;

// The copies of one gene rooted at its reference root, for the exact copy
// test. scan_taxa records each taxon's first and last copy in preorder;
// scan_nodes then gives every node's number of distinct taxa and whether it
// holds all or none of every taxon's copies (its subtree is the preorder
// interval [tin, tin + size), so it is closed when the first and last copy
// of every taxon below it fall inside). The taxon-indexed arrays are reset on
// the taxa the previous gene touched, so a scan costs O(gene) after an O(n)
// construction.
class GeneCopyScan {
 public:
  explicit GeneCopyScan(int taxon_count);

  // Scans the taxa of `gene` in the preorder of `tr` (its reference
  // traversal) and returns the number of distinct taxa.
  int scan_taxa(const GeneIndex& gene, const GeneTraversal& tr);
  // Computes the node data; call after scan_taxa on the same gene.
  void scan_nodes(const GeneIndex& gene, const GeneTraversal& tr);

  int32_t tin(int node) const { return tin_[static_cast<std::size_t>(node)]; }
  // Preorder positions of a taxon's first and last copy, -1 when absent.
  int32_t first(int taxon) const {
    return first_[static_cast<std::size_t>(taxon)];
  }
  int32_t last(int taxon) const { return last_[static_cast<std::size_t>(taxon)]; }
  // The distinct taxa of the gene, in order of first copy.
  const std::vector<int32_t>& taxa() const { return touched_; }
  int32_t distinct(int node) const {
    return distinct_[static_cast<std::size_t>(node)];
  }
  // False at the reference root.
  bool closed(int node) const;

 private:
  int taxon_count_ = 0;
  int root_ = -1;
  std::vector<int32_t> first_;
  std::vector<int32_t> last_;
  std::vector<int32_t> touched_;
  std::vector<int32_t> tin_;
  std::vector<int32_t> size_;
  std::vector<int32_t> min_first_;
  std::vector<int32_t> max_last_;
  std::vector<int32_t> distinct_;
};

// A forest listed children first: parent[i] > i, or -1 for a root. A vertex
// with taxon[i] >= 0 is a leaf; the leaves of one tree carry distinct taxa,
// and every vertex has a leaf below it.
struct ClusterForest {
  std::vector<int32_t> parent;
  std::vector<int32_t> taxon;
};

// Exact identifiers of the forest's clusters: ids[i] in [1, size], and
// ids[i] == ids[j] exactly when vertices i and j have the same leaf taxa.
// Halves the taxon range [0, taxon_count) recursively and restricts the forest
// to each half, so the time is linear in the forest plus
// O(leaves * log(taxon_count)) and the recursion depth is logarithmic. The
// forest is consumed.
std::vector<uint32_t> exact_cluster_ids(ClusterForest&& forest,
                                        int taxon_count);

// For every node of every gene, rooted at its reference root: the exact
// identifier of the node's distinct taxa when the node holds all or none of
// the copies of every taxon and at least two distinct taxa, and 0 otherwise.
// Genes with fewer than four distinct taxa get 0 throughout. Identifiers are
// comparable across genes: equal exactly when the taxon sets are equal.
struct GeneClusterIds {
  std::vector<std::size_t> offset;
  std::vector<uint32_t> id;

  uint32_t at(int gene, int node) const {
    return id[offset[static_cast<std::size_t>(gene)] +
              static_cast<std::size_t>(node)];
  }
};

GeneClusterIds exact_gene_cluster_ids(const std::vector<GeneIndex>& genes,
                                      int taxon_count);

// The MRE's distinct splits in priority order; for the tests.
std::vector<SplitAggregate256> mre_split_entries_for_tests(
    const std::vector<GeneIndex>& genes,
    int taxon_count);

}  // namespace qraft::internal
