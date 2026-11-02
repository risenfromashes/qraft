#include "divide/bipartition_test.hpp"
#include "core/checked.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace qraft::internal {
namespace {

// A tile of the bipartition test's (gene, taxon) grid spans at least this many
// taxa, so that a large subproblem is cut along its taxa and needs partial sums
// of only a few gene ranges.
constexpr int kTestTaxaPerTile = 64;

}  // namespace

namespace {

// Entries per block, and blocks per superblock, of the range-minimum index.
constexpr int kRmqGroupBits = 5;
constexpr int kRmqGroup = 1 << kRmqGroupBits;

// The minimum of values[l..r] for l <= r in one group of kRmqGroup entries.
// masks[r] has bit t set when entry (group start + t) is smaller than every
// later entry up to r, so the lowest such bit at or after l is the minimum.
int group_min(const int* values, const std::uint32_t* masks, int l, int r) {
  const int start = r & ~(kRmqGroup - 1);
  const int offset = l - start;
  return values[start + offset + std::countr_zero(masks[r] >> offset)];
}

// The suffix-minimum masks of `values`, in groups of kRmqGroup entries: an
// entry clears the entries of its group at least as large, then adds itself.
std::vector<std::uint32_t> group_masks(const std::vector<int>& values) {
  std::vector<std::uint32_t> masks(values.size());
  for (std::size_t p = 0; p < values.size(); ++p) {
    const std::size_t start = p & ~static_cast<std::size_t>(kRmqGroup - 1);
    std::uint32_t mask = p == start ? 0u : masks[p - 1];
    while (mask != 0u) {
      const int top = 31 - std::countl_zero(mask);
      if (values[start + static_cast<std::size_t>(top)] < values[p]) break;
      mask &= ~(1u << top);
    }
    masks[p] = mask | (1u << (p - start));
  }
  return masks;
}

// The minima of consecutive groups of kRmqGroup entries.
std::vector<int> group_minima(const std::vector<int>& values) {
  std::vector<int> minima((values.size() + kRmqGroup - 1) / kRmqGroup);
  for (std::size_t p = 0; p < values.size(); ++p) {
    const std::size_t group = p / kRmqGroup;
    minima[group] = p % kRmqGroup == 0 ? values[p]
                                       : std::min(minima[group], values[p]);
  }
  return minima;
}

}  // namespace

int GeneDistanceIndex::leaf_depth(int r) const {
  if (leaf_root && r == 0) return 0;
  const int count = static_cast<int>(lca_depth.size());
  const int left = r > 0 ? lca_depth[static_cast<std::size_t>(r - 1)] : -1;
  const int right = r < count ? lca_depth[static_cast<std::size_t>(r)] : -1;
  return std::max(left, right) + 1;
}

int GeneDistanceIndex::common_depth(int r, int s) const {
  if (r == s) return leaf_depth(r);
  if (r > s) std::swap(r, s);
  // The minimum of lca_depth[r..s-1].
  const int l = r;
  const int h = s - 1;
  const int lb = l >> kRmqGroupBits;
  const int hb = h >> kRmqGroupBits;
  const int* values = lca_depth.data();
  const std::uint32_t* masks = lca_masks.data();
  if (lb == hb) return group_min(values, masks, l, h);
  int best = std::min(group_min(values, masks, l, (lb << kRmqGroupBits) | (kRmqGroup - 1)),
                      group_min(values, masks, hb << kRmqGroupBits, h));
  if (lb + 1 == hb) return best;
  // Blocks lb + 1 .. hb - 1.
  const int bl = lb + 1;
  const int bh = hb - 1;
  const int ls = bl >> kRmqGroupBits;
  const int hs = bh >> kRmqGroupBits;
  const int* minima = block_min.data();
  const std::uint32_t* block = block_masks.data();
  if (ls == hs) return std::min(best, group_min(minima, block, bl, bh));
  best = std::min({best,
                   group_min(minima, block, bl, (ls << kRmqGroupBits) | (kRmqGroup - 1)),
                   group_min(minima, block, hs << kRmqGroupBits, bh)});
  if (ls + 1 == hs) return best;
  // Superblocks ls + 1 .. hs - 1.
  const int sl = ls + 1;
  const int sh = hs - 1;
  const int level = std::bit_width(static_cast<unsigned>(sh - sl + 1)) - 1;
  const int* row = sparse.data() + static_cast<std::size_t>(level) * superblocks;
  return std::min({best, row[sl], row[sh - (1 << level) + 1]});
}

int GeneDistanceIndex::distance(int r, int s) const {
  if (r == s) return 0;
  // At most the gene's node count minus one, so the result fits in int.
  return static_cast<int>(static_cast<long long>(leaf_depth(r)) + leaf_depth(s) -
                          2LL * common_depth(r, s));
}

namespace {

// The node of `taxon` in `gene`, or -1 when the gene lacks it. An empty gene
// has no taxon map.
int taxon_node(const GeneIndex& gene, int taxon) {
  return static_cast<std::size_t>(taxon) < gene.taxon_to_node.size()
             ? gene.taxon_to_node[static_cast<std::size_t>(taxon)]
             : -1;
}

// The distance index of one input gene (GeneDistanceIndex), O(|G|) time.
// Checks the leaf-depth rule against the depths of the traversal.
GeneDistanceIndex build_distance_index(const GeneIndex& gene) {
  GeneDistanceIndex out;
  const int size = checked_hdt_int_size(gene.nodes.size(), "gene too large");
  if (size == 0) return out;
  const std::vector<int>& preorder = gene.root_traversal.preorder;
  const std::vector<int>& parent = gene.root_traversal.parent;
  if (preorder.size() != gene.nodes.size() || parent.size() != gene.nodes.size()) {
    throw std::logic_error("bipartition test needs a full traversal");
  }
  out.rank_by_taxon = gene.taxon_to_node.size() <= gene.nodes.size();
  out.rank.assign(out.rank_by_taxon ? gene.taxon_to_node.size() : gene.nodes.size(), -1);
  out.leaf_root = gene.nodes[static_cast<std::size_t>(gene.root)].taxon() >= 0;
  std::size_t leaves = 0;
  for (const GeneIndex::Node& node : gene.nodes) leaves += node.taxon() >= 0 ? 1u : 0u;
  out.lca_depth.reserve(leaves == 0 ? 0 : leaves - 1);
  std::vector<int> depth(gene.nodes.size(), 0);
  std::vector<int> leaf_depths;
  leaf_depths.reserve(leaves);
  // The smallest depth since the last leaf: the depth of the next leaf's LCA
  // with it, plus one.
  int since_leaf = std::numeric_limits<int>::max();
  for (const int u : preorder) {
    if (u != gene.root) {
      const int p = parent[static_cast<std::size_t>(u)];
      if (p < 0 || p >= size) {
        throw std::logic_error("bipartition test needs a full traversal");
      }
      depth[static_cast<std::size_t>(u)] = depth[static_cast<std::size_t>(p)] + 1;
    }
    const int d = depth[static_cast<std::size_t>(u)];
    since_leaf = std::min(since_leaf, d);
    const int taxon = gene.nodes[static_cast<std::size_t>(u)].taxon();
    if (taxon < 0) continue;
    const int r = static_cast<int>(leaf_depths.size());
    if (r > 0) out.lca_depth.push_back(since_leaf - 1);
    since_leaf = std::numeric_limits<int>::max();
    const int key = out.rank_by_taxon ? taxon : u;
    if (static_cast<std::size_t>(key) >= out.rank.size()) {
      throw std::logic_error("bipartition test leaf outside the taxon map");
    }
    out.rank[static_cast<std::size_t>(key)] = r;
    leaf_depths.push_back(d);
  }
  for (int r = 0; r < static_cast<int>(leaf_depths.size()); ++r) {
    if (out.leaf_depth(r) != leaf_depths[static_cast<std::size_t>(r)]) {
      throw std::logic_error("bipartition test gene has a unary internal node");
    }
  }
  out.lca_masks = group_masks(out.lca_depth);
  out.block_min = group_minima(out.lca_depth);
  out.block_masks = group_masks(out.block_min);
  std::vector<int> level = group_minima(out.block_min);
  out.superblocks = static_cast<int>(level.size());
  const int levels = out.superblocks == 0 ? 0 : std::bit_width(static_cast<unsigned>(out.superblocks));
  out.sparse.resize(checked_hdt_size_mul(static_cast<std::size_t>(levels),
                                         level.size(), "bipartition test sparse table overflow"));
  std::copy(level.begin(), level.end(), out.sparse.begin());
  for (int j = 1; j < levels; ++j) {
    const int* below = out.sparse.data() + static_cast<std::size_t>(j - 1) * out.superblocks;
    int* row = out.sparse.data() + static_cast<std::size_t>(j) * out.superblocks;
    const int half = 1 << (j - 1);
    for (int s = 0; s + 2 * half <= out.superblocks; ++s) {
      row[s] = std::min(below[s], below[s + half]);
    }
  }
  return out;
}

}  // namespace

std::shared_ptr<const BipartitionTestIndex> make_bipartition_test_index(
    const std::vector<GeneIndex>& genes, int taxon_count) {
  const auto started = std::chrono::steady_clock::now();
  auto index = std::make_shared<BipartitionTestIndex>();
  index->taxon_count = taxon_count;
  if (taxon_count < 4 || genes.empty()) return index;
  const int gene_count = checked_hdt_int_size(genes.size(), "too many genes");
  std::vector<unsigned char> single(genes.size(), 1);
  std::vector<unsigned char> complete(genes.size(), 1);
  std::vector<int> held(genes.size(), 0);
  parallel_for_inner(0, gene_count, [&](int g) {
    const GeneIndex& gene = genes[g];
    std::vector<unsigned char> seen(static_cast<std::size_t>(taxon_count), 0);
    int count = 0;
    for (int v = 0; v < static_cast<int>(gene.nodes.size()); ++v) {
      const int taxon = gene.nodes[static_cast<std::size_t>(v)].taxon();
      if (taxon < 0) continue;
      // The test needs at most one leaf for each taxon in every gene (a gene
      // of one taxon is a single node of degree 0). A duplicated taxon
      // disables it; picking one of its leaves through taxon_to_node would
      // not be sound.
      if (taxon >= taxon_count || seen[taxon] ||
          gene_degree(gene, v) > 1) {
        single[g] = 0;
        break;
      }
      seen[taxon] = 1;
      ++count;
    }
    held[g] = count;
    if (count != taxon_count) complete[g] = 0;
  }, 1);
  const auto all_set = [](const std::vector<unsigned char>& flags) {
    return std::all_of(flags.begin(), flags.end(), [](auto x) { return x != 0; });
  };
  index->single_copy = all_set(single);
  index->complete_single_copy = index->single_copy && all_set(complete);
  if (index->single_copy) {
    index->gene_taxa = std::move(held);
    index->taxon_genes.assign(static_cast<std::size_t>(taxon_count), 0);
    for (const GeneIndex& gene : genes) {
      for (int taxon = 0; taxon < taxon_count; ++taxon) {
        if (taxon_node(gene, taxon) >= 0) {
          ++index->taxon_genes[static_cast<std::size_t>(taxon)];
        }
      }
    }
    index->genes.resize(genes.size());
    parallel_for_inner(0, gene_count, [&](int g) {
      index->genes[g] = build_distance_index(genes[g]);
    }, 1);
    std::int64_t bytes = 0;
    for (const GeneDistanceIndex& gene : index->genes) {
      bytes += static_cast<std::int64_t>(
          gene.rank.capacity() * sizeof(int) + gene.lca_depth.capacity() * sizeof(int) +
          gene.lca_masks.capacity() * sizeof(std::uint32_t) +
          gene.block_min.capacity() * sizeof(int) +
          gene.block_masks.capacity() * sizeof(std::uint32_t) +
          gene.sparse.capacity() * sizeof(int));
    }
    stats_add(run_stats().bipartition_test_index_bytes, bytes);
  }
  stats_add_time(run_stats().bipartition_test_index_ns,
                      std::chrono::steady_clock::now() - started);
  return index;
}

double bipartition_test_threshold(int gene_count, int taxon_count,
                                  double coefficient, double crossover) {
  if (taxon_count < 2) {
    throw std::invalid_argument("bipartition test needs at least two taxa");
  }
  return bipartition_test_threshold_range(
      gene_count, 2.0 * (taxon_count - 1), coefficient, crossover);
}

double bipartition_test_threshold_range(int gene_count, double range,
                                        double coefficient, double crossover) {
  if (gene_count <= 0 || !std::isfinite(range) || range <= 0.0 ||
      !std::isfinite(coefficient) || coefficient <= 0.0 ||
      !std::isfinite(crossover) || crossover < 0.0) {
    throw std::invalid_argument(
        "bipartition test needs genes, a finite positive range, a finite "
        "positive coefficient, and a finite nonnegative crossover");
  }
  if (gene_count == crossover) return 0.0;
  const long double k = gene_count;
  const long double k0 = crossover;
  const long double e = std::numbers::e_v<long double>;
  const long double log_k = std::log(k + e);
  const long double log_k0 = std::log(k0 + e);
  // sqrt(a) - sqrt(b) = (a - b) / (sqrt(a) + sqrt(b)). Near the crossover,
  // log1p computes a - b without cancellation and keeps its sign on either
  // side of K0.
  const long double log_ratio = std::abs(k - k0) <= (k0 + e) / 2.0L
      ? std::log1p((k - k0) / (k0 + e)) : log_k - log_k0;
  const double threshold = static_cast<double>(
      static_cast<long double>(coefficient) * range / std::sqrt(k) *
      log_ratio / (std::sqrt(log_k) + std::sqrt(log_k0)));
  // Zero is valid only at the crossover. A threshold that rounds to zero or
  // to the wrong sign would change the configured decision rule.
  if (!std::isfinite(threshold) || threshold == 0.0 ||
      (threshold > 0.0) != (gene_count > crossover)) {
    throw std::invalid_argument("bipartition test threshold unrepresentable");
  }
  return threshold;
}

namespace {

// The outcome of one anchor pair on genes that miss taxa: whether it could be
// tested, and the separation, threshold and shared-gene count of the pair of
// taxa that decided it.
struct PairOutcome {
  bool usable = false;
  double gap = 0.0;
  double threshold = 0.0;
  int genes = 0;
};

// Shared inputs of the per-pair statistics.
struct IncompleteInput {
  const BipartitionTestIndex& index;
  const std::vector<GeneIndex>& genes;
  const std::vector<int>& active_taxa;
  const Bitset& split;
  const BipartitionTestOptions& options;
  int k;
  int m;
};

// Gene ranges of about four per worker, so per-gene work is split evenly and
// partial rows are summed in gene order.
int incomplete_gene_blocks(int k) {
  return std::clamp(4 * current_inner_parallel_worker_limit(), 1, k);
}

int block_first(int count, int block, int blocks) {
  return static_cast<int>(static_cast<std::int64_t>(count) * block / blocks);
}

// The position, counted in edges from a, at which x's path to the a-b path
// meets it in gene g: (d(a,x) + d(a,b) - d(b,x)) / 2. With q(u,v) the depth
// of LCA(u,v) and A the depth of a, that is A - q(a,x) - q(a,b) + q(b,x), for
// the leaf ranks ra, rb and rx and a_minus_ab = A - q(a,b).
int attachment(const GeneDistanceIndex& distances, int ra, int rb,
               int a_minus_ab, int rx) {
  return a_minus_ab - distances.common_depth(ra, rx) +
         distances.common_depth(rb, rx);
}

// J on shared genes. In a gene holding the anchors a and b, h_g(z) is the
// position, counted in edges from a, at which z meets the a-b path. Two taxa
// x and y are compared only on the k_xy genes that hold a, b, x and y:
//   Delta(x, y) = (1/k_xy) sum over those genes of [h_g(y) - h_g(x)].
// Given the taxa a gene holds, its tree is an MSC gene tree on the species
// tree restricted to them, so the expected sign of h_g(y) - h_g(x) is the
// species quartet of a, b, x and y, whatever else the gene misses.
//
// A pair tested passes when Delta(x, y) > t(k_xy, R_xy), with R_xy = 2(mean
// taxa of the shared genes - 1) and t = bipartition_test_threshold_range. A
// running maximum over A and a running minimum over B, in active order and
// keeping the champion on ties, find x* and y*: when the comparisons are
// accurate, the taxon of A attaching last and the taxon of B attaching first.
// The pairs (x*, y) and (x, y*) that
// share a gene are tested, and the pair of smallest margin decides. On
// complete genes Delta(x, y) is J(y) - J(x) halved, every pair has the
// threshold t_k, and the decision is that of the J test and of testing every
// pair. On genes that miss taxa the comparisons use different genes and need
// not form one order, so a pair outside the star can fail while the star
// passes; the consistency argument is for the star. With only the anchor on
// one side there is no pair and the bipartition is accepted. The anchor pair
// is unusable when a side has no taxon in a gene with both anchors, or x* and
// y* share no gene.
//
// O(k m) time for the positions, O(k) per comparison and O(k m) space, for m
// active taxa.
PairOutcome shared_gene_pair(const IncompleteInput& in, int a, int b) {
  const int m = in.m;
  const int k = in.k;
  PairOutcome out;
  bool pairs_a = false;
  bool pairs_b = false;
  for (int i = 0; i < m; ++i) {
    const int taxon = in.active_taxa[i];
    if (in.split.test(i)) pairs_a = pairs_a || taxon != a;
    else pairs_b = pairs_b || taxon != b;
  }
  if (!pairs_a || !pairs_b) {
    // A side of one taxon: the bipartition is in every tree.
    out.usable = true;
    out.gap = std::numeric_limits<double>::infinity();
    return out;
  }
  // position[i * k + g] = h_g(active_taxa[i]), or -1 when gene g lacks a, b or
  // the taxon, or the taxon is an anchor.
  std::vector<int> position(checked_hdt_size_mul(
      static_cast<std::size_t>(m), static_cast<std::size_t>(k),
      "bipartition test positions overflow"), -1);
  const int blocks = incomplete_gene_blocks(k);
  parallel_for_dynamic(0, blocks, [&](int block) {
    for (int g = block_first(k, block, blocks);
         g < block_first(k, block + 1, blocks); ++g) {
      const GeneIndex& gene = in.genes[g];
      const GeneDistanceIndex& distances = in.index.genes[g];
      const int an = taxon_node(gene, a);
      const int bn = taxon_node(gene, b);
      if (an < 0 || bn < 0) continue;
      const int ra = distances.leaf_rank(a, an);
      const int rb = distances.leaf_rank(b, bn);
      const int a_minus_ab = distances.leaf_depth(ra) - distances.common_depth(ra, rb);
      for (int i = 0; i < m; ++i) {
        const int taxon = in.active_taxa[i];
        const int x = taxon_node(gene, taxon);
        if (x < 0 || taxon == a || taxon == b) continue;
        position[static_cast<std::size_t>(i) * k + g] =
            attachment(distances, ra, rb, a_minus_ab, distances.leaf_rank(taxon, x));
      }
    }
  }, 1);
  stats_add(run_stats().bipartition_test_distance_queries, 2LL * k * m);
  // Integer sums over the genes holding a, b, x and y, in gene order.
  struct Shared {
    std::int64_t difference = 0;  // sum of h_g(y) - h_g(x)
    std::int64_t genes = 0;       // k_xy
    std::int64_t taxa = 0;        // sum of the taxa those genes hold
  };
  const auto compare = [&](int x, int y) {
    Shared s;
    const int* px = position.data() + static_cast<std::size_t>(x) * k;
    const int* py = position.data() + static_cast<std::size_t>(y) * k;
    for (int g = 0; g < k; ++g) {
      if (px[g] < 0 || py[g] < 0) continue;
      s.difference += py[g] - px[g];
      ++s.genes;
      s.taxa += in.index.gene_taxa[static_cast<std::size_t>(g)];
    }
    return s;
  };
  // The running maximum (side A) or minimum (side B) of the attachment
  // order. The first champion is the first taxon held by a gene with both
  // anchors; a later taxon replaces it when it attaches strictly later (A) or
  // earlier (B) on their shared genes.
  const auto champion = [&](bool side_a) {
    int c = -1;
    for (int i = 0; i < m; ++i) {
      const int taxon = in.active_taxa[i];
      if (in.split.test(i) != side_a || taxon == a || taxon == b) continue;
      if (c < 0) {
        const int* p = position.data() + static_cast<std::size_t>(i) * k;
        if (std::any_of(p, p + k, [](int h) { return h >= 0; })) c = i;
        continue;
      }
      const std::int64_t difference = compare(c, i).difference;
      if (side_a ? difference > 0 : difference < 0) c = i;
    }
    return c;
  };
  const int xs = champion(true);
  const int ys = champion(false);
  if (xs < 0 || ys < 0) return out;
  // (x*, y*), then (x*, y) for the rest of B, then (x, y*) for the rest of A.
  std::vector<std::pair<int, int>> pairs{{xs, ys}};
  for (int i = 0; i < m; ++i) {
    if (i != ys && !in.split.test(i) && in.active_taxa[i] != b) {
      pairs.emplace_back(xs, i);
    }
  }
  for (int i = 0; i < m; ++i) {
    if (i != xs && in.split.test(i) && in.active_taxa[i] != a) {
      pairs.emplace_back(i, ys);
    }
  }
  double best_margin = std::numeric_limits<double>::infinity();
  for (const auto& [x, y] : pairs) {
    const Shared s = compare(x, y);
    if (s.genes == 0) {
      if (x == xs && y == ys) return out;
      continue;
    }
    const double genes = static_cast<double>(s.genes);
    const double gap = static_cast<double>(s.difference) / genes;
    // Every shared gene holds a, b, x and y, so the range is at least 6.
    const double range = 2.0 * (static_cast<double>(s.taxa) / genes - 1.0);
    const double threshold = bipartition_test_threshold_range(
        checked_hdt_int_size(static_cast<std::size_t>(s.genes), "too many genes"),
        range, in.options.coefficient, in.options.crossover);
    if (gap - threshold < best_margin) {
      best_margin = gap - threshold;
      out.gap = gap;
      out.threshold = threshold;
      out.genes = static_cast<int>(s.genes);
    }
  }
  out.usable = true;
  return out;
}

// The test on genes that miss taxa (test_bipartition), with J on shared genes
// (shared_gene_pair). The anchors are the two taxa of each side held by the
// most genes, ties by position; an unusable anchor pair is skipped. All sums
// are integers in gene order, so the result does not depend on the workers.
BipartitionTestResult test_bipartition_incomplete(
    const SolverContext& ctx, const std::vector<int>& active_taxa,
    const Bitset& split) {
  BipartitionTestResult result;
  const auto& index = *ctx.bipartition_test_index;
  const auto& genes = base_genes(ctx);
  const int k = checked_hdt_int_size(genes.size(), "too many genes");
  const int m = checked_hdt_int_size(active_taxa.size(), "too many active taxa");
  if (index.genes.size() != genes.size() || split.size() != m) {
    throw std::logic_error("bipartition test index or split size mismatch");
  }
  std::array<int, 2> side_a{-1, -1};
  std::array<int, 2> side_b{-1, -1};
  for (int i = 0; i < m; ++i) {
    const int taxon = active_taxa[i];
    if (taxon < 0 || taxon >= index.taxon_count) {
      throw std::logic_error("bipartition test has an invalid active taxon");
    }
    std::array<int, 2>& best = split.test(i) ? side_a : side_b;
    const int held = index.taxon_genes[static_cast<std::size_t>(taxon)];
    const auto held_by = [&](int t) {
      return t < 0 ? -1 : index.taxon_genes[static_cast<std::size_t>(t)];
    };
    if (held > held_by(best[0])) {
      best[1] = best[0];
      best[0] = taxon;
    } else if (held > held_by(best[1])) {
      best[1] = taxon;
    }
  }
  if (side_a[0] < 0 || side_b[0] < 0) {
    throw std::logic_error("bipartition test needs two nonempty sides");
  }
  const std::array<std::pair<int, int>, 4> anchor_pairs{{
      {side_a[0], side_b[0]}, {side_a[1], side_b[1]},
      {side_a[0], side_b[1]}, {side_a[1], side_b[0]}}};
  const IncompleteInput in{index, genes, active_taxa, split,
                           ctx.options.bipartition_test, k, m};
  result.gap = -std::numeric_limits<double>::infinity();
  double best_margin = -std::numeric_limits<double>::infinity();
  int tested = 0;
  bool negative_accept = false;
  for (int pair = 0; pair < static_cast<int>(anchor_pairs.size()) &&
                     tested < ctx.options.bipartition_test.anchor_pairs;
       ++pair) {
    const auto anchors = anchor_pairs[pair];
    if (anchors.first < 0 || anchors.second < 0 ||
        std::find(anchor_pairs.begin(), anchor_pairs.begin() + pair, anchors) !=
            anchor_pairs.begin() + pair) {
      continue;
    }
    const PairOutcome outcome = shared_gene_pair(in, anchors.first, anchors.second);
    if (!outcome.usable) {
      stats_add(run_stats().bipartition_test_uncovered_pairs);
      continue;
    }
    ++tested;
    result.available = true;
    // Report the accepting pair, or else the pair of largest margin.
    const double margin = outcome.gap - outcome.threshold;
    if (tested == 1 || margin > best_margin) {
      best_margin = margin;
      result.anchor_a = anchors.first;
      result.anchor_b = anchors.second;
      result.threshold_genes = outcome.genes;
      result.threshold = outcome.threshold;
      result.gap = outcome.gap;
    }
    stats_add(run_stats().bipartition_test_shared_genes, outcome.genes);
    if (outcome.gap > outcome.threshold) {
      result.accepted = true;
      negative_accept = outcome.gap < 0.0;
      break;
    }
  }
  if (!result.available) {
    stats_add(run_stats().bipartition_test_untestable);
    return result;
  }
  stats_add(run_stats().bipartition_tests);
  stats_add(run_stats().bipartition_test_incomplete);
  stats_add(result.accepted ? run_stats().bipartition_test_accepts
                            : run_stats().bipartition_test_rejects);
  if (negative_accept) {
    stats_add(run_stats().bipartition_test_negative_accepts);
  }
  return result;
}

}  // namespace

BipartitionTestResult test_bipartition_on_shared_genes_for_tests(
    const SolverContext& ctx, const std::vector<int>& active_taxa,
    const Bitset& split) {
  if (!ctx.bipartition_test_index || !ctx.bipartition_test_index->single_copy) {
    throw std::logic_error("the shared-gene test needs a single-copy index");
  }
  return test_bipartition_incomplete(ctx, active_taxa, split);
}

BipartitionTestResult test_bipartition(
    const SolverContext& ctx, const std::vector<int>& active_taxa,
    const Bitset& split) {
  BipartitionTestResult result;
  if (!ctx.bipartition_test_index ||
      !ctx.bipartition_test_index->single_copy) return result;
  if (ctx.options.bipartition_test.anchor_pairs < 1 ||
      ctx.options.bipartition_test.anchor_pairs > 4) {
    throw std::invalid_argument("bipartition test requires 1..4 anchor pairs");
  }
  if (!ctx.bipartition_test_index->complete_single_copy) {
    const auto started = std::chrono::steady_clock::now();
    result = test_bipartition_incomplete(ctx, active_taxa, split);
    stats_add_time(run_stats().bipartition_test_ns,
                   std::chrono::steady_clock::now() - started);
    return result;
  }
  const auto started = std::chrono::steady_clock::now();
  const auto& index = *ctx.bipartition_test_index;
  const auto& genes = base_genes(ctx);
  const int k = checked_hdt_int_size(genes.size(), "too many genes");
  const int m = checked_hdt_int_size(active_taxa.size(), "too many active taxa");
  if (index.genes.size() != genes.size() || split.size() != m) {
    throw std::logic_error("bipartition test index or split size mismatch");
  }
  // The first and last taxon of each side, in active order.
  int a = -1;
  int b = -1;
  int last_a = -1;
  int last_b = -1;
  for (int i = 0; i < m; ++i) {
    const int taxon = active_taxa[i];
    if (taxon < 0 || taxon >= index.taxon_count) {
      throw std::logic_error("bipartition test has an invalid active taxon");
    }
    int& anchor = split.test(i) ? a : b;
    if (anchor < 0) anchor = taxon;
    (split.test(i) ? last_a : last_b) = taxon;
  }
  if (a < 0 || b < 0) {
    throw std::logic_error("bipartition test needs two nonempty sides");
  }
  // Integer sums over the genes: no rounding before the final division and no
  // dependence on how the work is split. The (gene, taxon) grid is cut into
  // tiles of whole gene ranges and taxon ranges, about four per worker; each
  // tile adds into its gene range's partial row, and the partial rows are
  // summed per taxon.
  const int workers = current_inner_parallel_worker_limit();
  const int target_tiles = 4 * workers;
  const int taxon_blocks = std::clamp(m / kTestTaxaPerTile, 1, target_tiles);
  const int gene_blocks =
      std::clamp((target_tiles + taxon_blocks - 1) / taxon_blocks, 1, k);
  std::vector<std::int64_t> partial(checked_hdt_size_mul(
      static_cast<std::size_t>(gene_blocks), active_taxa.size(),
      "bipartition test partial sums overflow"));
  const std::array<std::pair<int, int>, 4> anchor_pairs{{
      {a, b}, {last_a, last_b}, {a, last_b}, {last_a, b}}};
  result.available = true;
  result.gap = -std::numeric_limits<double>::infinity();
  result.threshold_genes = k;
  result.threshold = bipartition_test_threshold(
      k, index.taxon_count, ctx.options.bipartition_test.coefficient,
      ctx.options.bipartition_test.crossover);
  std::int64_t best_gap_numerator = std::numeric_limits<std::int64_t>::min();
  for (int pair = 0; pair < ctx.options.bipartition_test.anchor_pairs; ++pair) {
    const auto anchors = anchor_pairs[pair];
    if (std::find(anchor_pairs.begin(), anchor_pairs.begin() + pair, anchors) !=
        anchor_pairs.begin() + pair) {
      continue;
    }
    std::fill(partial.begin(), partial.end(), 0);
    parallel_for_dynamic(0, gene_blocks * taxon_blocks, [&](int tile) {
      const int gene_block = tile / taxon_blocks;
      const int taxon_block = tile % taxon_blocks;
      const int g_first = static_cast<int>(
          static_cast<std::int64_t>(k) * gene_block / gene_blocks);
      const int g_last = static_cast<int>(
          static_cast<std::int64_t>(k) * (gene_block + 1) / gene_blocks);
      const int i_first = static_cast<int>(
          static_cast<std::int64_t>(m) * taxon_block / taxon_blocks);
      const int i_last = static_cast<int>(
          static_cast<std::int64_t>(m) * (taxon_block + 1) / taxon_blocks);
      std::int64_t* row =
          partial.data() + static_cast<std::size_t>(gene_block) * m;
      for (int g = g_first; g < g_last; ++g) {
        const GeneIndex& gene = genes[g];
        const GeneDistanceIndex& distances = index.genes[g];
        const auto rank_of = [&](int taxon) {
          return distances.leaf_rank(
              taxon, distances.rank_by_taxon ? -1 : gene.taxon_to_node[taxon]);
        };
        const int ra = rank_of(anchors.first);
        const int rb = rank_of(anchors.second);
        // d(a,x) - d(b,x) = A - B + 2(q(b,x) - q(a,x)), with A and B the
        // anchors' depths and q(u,v) the depth of LCA(u,v).
        const std::int64_t anchor_difference =
            static_cast<std::int64_t>(distances.leaf_depth(ra)) - distances.leaf_depth(rb);
        for (int i = i_first; i < i_last; ++i) {
          const int rx = rank_of(active_taxa[i]);
          row[i] += anchor_difference +
                    2 * (static_cast<std::int64_t>(distances.common_depth(rb, rx)) -
                         distances.common_depth(ra, rx));
        }
      }
    }, 1);
    std::vector<std::int64_t> sums(active_taxa.size(), 0);
    parallel_for_inner(0, taxon_blocks, [&](int taxon_block) {
      const int i_first = static_cast<int>(
          static_cast<std::int64_t>(m) * taxon_block / taxon_blocks);
      const int i_last = static_cast<int>(
          static_cast<std::int64_t>(m) * (taxon_block + 1) / taxon_blocks);
      for (int gene_block = 0; gene_block < gene_blocks; ++gene_block) {
        const std::int64_t* row =
            partial.data() + static_cast<std::size_t>(gene_block) * m;
        for (int i = i_first; i < i_last; ++i) sums[i] += row[i];
      }
    }, 1);
    // Each distance difference is at most |G| - 1 <= INT_MAX in magnitude
    // and k <= INT_MAX, so the sums and their difference fit in 64 bits.
    std::int64_t max_a = std::numeric_limits<std::int64_t>::min();
    std::int64_t min_b = std::numeric_limits<std::int64_t>::max();
    for (int i = 0; i < m; ++i) {
      if (split.test(i)) {
        max_a = std::max(max_a, sums[i]);
      } else {
        min_b = std::min(min_b, sums[i]);
      }
    }
    // The anchor separation (min_B J - max_A J) / 2, where J is the sum of
    // a taxon's row entries divided by k.
    const std::int64_t gap_numerator = min_b - max_a;
    const double gap = static_cast<double>(gap_numerator) / (2.0 * k);
    if (!std::isfinite(gap)) break;
    best_gap_numerator = std::max(best_gap_numerator, gap_numerator);
    if (gap > result.gap) {
      result.gap = gap;
      result.anchor_a = anchors.first;
      result.anchor_b = anchors.second;
    }
    stats_add(run_stats().bipartition_test_distance_queries,
                          2LL * k * m);
    // Accepting on any of a fixed number of anchor pairs is sound: in the
    // population limit the separation is positive for every pair a in A,
    // b in B exactly when A|B is a species-tree edge, however the pairs are
    // chosen.
    if (std::isfinite(result.threshold) && result.gap > result.threshold) {
      result.accepted = true;
      break;
    }
  }
  stats_add(run_stats().bipartition_tests);
  stats_add(result.accepted
                            ? run_stats().bipartition_test_accepts
                            : run_stats().bipartition_test_rejects);
  if (best_gap_numerator != std::numeric_limits<std::int64_t>::min() &&
      best_gap_numerator < 0) {
    stats_max(
        run_stats().bipartition_test_max_negative_gap_twice_k,
        -best_gap_numerator);
    if (result.accepted) {
      stats_add(run_stats().bipartition_test_negative_accepts);
    }
  }
  stats_add_time(run_stats().bipartition_test_ns,
                      std::chrono::steady_clock::now() - started);
  return result;
}

}  // namespace qraft::internal
