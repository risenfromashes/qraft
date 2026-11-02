#pragma once

#include "internal.hpp"

// The bipartition test. For a proposed bipartition A|B of the active taxa and
// anchors a in A, b in B, every active taxon z gets
//   J(z) = (1/k) sum_G [d_G(a,z) - d_G(b,z)],
// where d_G counts the edges between two leaves of input gene G. The anchor
// separation is (min_{y in B} J(y) - max_{x in A} J(x)) / 2, and the
// bipartition is accepted when the separation of some tested anchor pair
// exceeds t_k (bipartition_test_threshold).
//
// Genes that miss taxa (but hold none twice) use J on shared genes: with h_G(z)
// the position, counted in edges from a, at which z meets the a-b path of
// gene G, two taxa x in A and y in B are compared only on the k_xy genes that
// hold a, b, x and y,
//   Delta(x, y) = (1/k_xy) sum_{G holding a, b, x, y} [h_G(y) - h_G(x)],
// and the bipartition is accepted when Delta(x, y) > t(k_xy, R_xy) for the
// pairs tested (test_bipartition), with R_xy = 2(mean taxa of those genes - 1).
// On complete genes Delta(x, y) = (J(y) - J(x)) / 2 and this is the J test.

namespace qraft::internal {

// Leaf-to-leaf distances of one input gene in O(1) time per query and O(L)
// words for L leaves. Distances count edges of the input gene after
// degree-two suppression; a restricted gene set does not preserve them, so
// the index is built once from the input genes.
//
// The leaves are ranked in the preorder of gene.root_traversal, and
// lca_depth[r] is the depth, from gene.root, of the LCA of the leaves of
// ranks r and r + 1. Every subtree is a contiguous range of ranks, so the LCA
// of the leaves of ranks i < j has depth min lca_depth[i..j-1]. That minimum
// takes O(1) time from three levels: within a block of 32 entries, a mask per
// entry marks the suffix minima of its block up to it; within a superblock of
// 32 blocks, the same masks over the block minima; and a sparse table over
// the superblock minima, O(L / 1024 log L) words. A leaf's depth is one more
// than the larger LCA depth with its neighbours in rank order: below an
// internal root every internal node has at least two children, so one of the
// neighbours meets the leaf at its parent. A root that is a leaf (a gene of
// at most two leaves, or a tree given with a leaf root) has depth 0.
struct GeneDistanceIndex {
  // The rank of each leaf, keyed by taxon or by node, whichever domain is
  // smaller, and -1 elsewhere.
  std::vector<int> rank;
  bool rank_by_taxon = false;
  // gene.root is a leaf; it has rank 0.
  bool leaf_root = false;
  std::vector<int> lca_depth;
  std::vector<std::uint32_t> lca_masks;
  std::vector<int> block_min;
  std::vector<std::uint32_t> block_masks;
  // Level j holds the minima of 2^j superblocks from each superblock, at
  // offset j * superblocks.
  std::vector<int> sparse;
  int superblocks = 0;

  // The rank of the leaf of `taxon` at `node`, or -1 when the gene lacks it.
  int leaf_rank(int taxon, int node) const {
    const int key = rank_by_taxon ? taxon : node;
    return key >= 0 && static_cast<std::size_t>(key) < rank.size()
               ? rank[static_cast<std::size_t>(key)]
               : -1;
  }
  // The depth of the leaf of rank r.
  int leaf_depth(int r) const;
  // The depth of the LCA of the leaves of ranks r and s (r == s allowed).
  int common_depth(int r, int s) const;
  // The distance between the leaves of ranks r and s.
  int distance(int r, int s) const;
};

// Distance indices of every input gene, shared read-only by all subproblems.
// The test applies only when no gene holds a taxon twice (single_copy);
// otherwise `genes` is empty. complete_single_copy also requires every gene
// to hold every taxon.
struct BipartitionTestIndex {
  bool single_copy = false;
  bool complete_single_copy = false;
  int taxon_count = 0;
  std::vector<GeneDistanceIndex> genes;
  // The number of genes holding each taxon, when single_copy.
  std::vector<int> taxon_genes;
  // The number of taxa each gene holds, when single_copy.
  std::vector<int> gene_taxa;
};

struct BipartitionTestResult {
  // False when there is no usable index; the bipartition is then not tested.
  bool available = false;
  bool accepted = false;
  // The separation, threshold, anchors and k of the threshold of one tested
  // anchor pair: the accepting pair, or else the pair whose separation
  // exceeds its threshold by the most (on complete genes every pair has the
  // same t_k). k is the gene count; on genes that miss taxa the separation,
  // threshold and k are those of the pair of taxa that decided, k being its
  // shared genes k_xy.
  double gap = 0.0;
  double threshold = 0.0;
  int anchor_a = -1;
  int anchor_b = -1;
  int threshold_genes = 0;
};

// O(N + kn) time and O(N + n + k) space for N gene nodes, k genes and n taxa,
// parallel over genes.
std::shared_ptr<const BipartitionTestIndex> make_bipartition_test_index(
    const std::vector<GeneIndex>& genes, int taxon_count);

// t_k = C R / sqrt(k) * (sqrt(log(k + e)) - sqrt(log(K0 + e))), with
// R = 2(n - 1), k = gene_count, n = taxon_count, C = coefficient and
// K0 = crossover. It is negative below K0, zero at K0 and positive above.
// Throws when an argument is invalid or rounding would change the sign.
double bipartition_test_threshold(int gene_count, int taxon_count,
                                       double coefficient, double crossover);

// The same threshold with the range R given directly (R = 2(n - 1) above).
// J on shared genes uses R_xy = 2(mean taxa of the shared genes - 1).
double bipartition_test_threshold_range(int gene_count, double range,
                                        double coefficient, double crossover);

// Tests `split` (bit i set: active_taxa[i] is in A) with the first
// options.bipartition_test.anchor_pairs of the anchor pairs (first A, first
// B), (last A, last B), (first A, last B), (last A, first B), by position in
// active_taxa; a repeated pair is skipped. t_k uses the k input genes and the
// n taxa of the whole run. Each pair costs O(k m) time for m active taxa,
// and on genes that miss taxa O(k m) space.
// On genes that miss taxa, a1, a2 and b1, b2 are the two taxa of each side
// held by the most genes (ties by position), the pairs are (a1, b1),
// (a2, b2), (a1, b2), (a2, b1), and each tests (x*, y) for y in B \ {b} and
// (x, y*) for x in A \ {a}, where x* and y* are the running maximum and
// minimum of the attachment order on shared genes; the pair of taxa of
// smallest margin decides, and pairs that share no gene are skipped. A side
// holding only its anchor accepts. An anchor pair is skipped when a side has
// no taxon in a gene with both anchors or x* and y* share no gene; with no
// usable pair the result is not available.
BipartitionTestResult test_bipartition(
    const SolverContext& ctx, const std::vector<int>& active_taxa,
    const Bitset& split);

// For tests: the test of genes that miss taxa, run whatever the genes hold.
// Needs a single-copy index.
BipartitionTestResult test_bipartition_on_shared_genes_for_tests(
    const SolverContext& ctx, const std::vector<int>& active_taxa,
    const Bitset& split);

}  // namespace qraft::internal
