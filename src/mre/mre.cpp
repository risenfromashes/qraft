#include "internal.hpp"
#include "genes/gene_lca.hpp"
#include "mre/exact_split_ids.hpp"
#include "mre/split_records.hpp"

#include <array>
#include <atomic>
#include <climits>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <queue>
#include <sstream>

// MRE starting tree: the majority-rule extended consensus of the gene trees.
//
// Counting emits one record per eligible gene edge, keyed by the side away
// from the gene's reference root, and radix-sorts the distinct splits into
// priority order: support (the number of genes that contain the split)
// descending, then first gene, first edge and key. The key is the side's size
// and the exact identifier of its taxon set (exact_split_ids.hpp).
//
// Assembly starts from a star and accepts, in priority order, each split
// compatible with those already accepted, up to n - 3 splits. For every
// (starting-tree node, gene) pair an AVL index holds the node's children
// keyed by the preorder time of their LCA in that gene, so a candidate is
// tested against one node's children by a range query. An accepted split
// re-parents only the smaller of that node's two child groups, at most half
// of its children, at O(log n) index work per moved child and gene.

namespace qraft::internal {
namespace {

std::size_t checked_mre_size_product(std::size_t first,
                                     std::size_t second,
                                     const char* what) {
  if (first != 0 && second > std::numeric_limits<std::size_t>::max() / first) {
    throw std::length_error(what);
  }
  return first * second;
}

// first * second as a flat index; it must also fit in int.
std::size_t checked_mre_flat_slots(int first,
                                   int second,
                                   const char* what) {
  if (first < 0 || second < 0) throw std::length_error(what);
  const std::size_t slots = checked_mre_size_product(
      static_cast<std::size_t>(first), static_cast<std::size_t>(second), what);
  if (slots > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error(what);
  }
  return slots;
}

// The starting tree has at most max(2, 2n + 1) nodes, and the per-(node,
// gene) arrays over that many nodes must be int-indexable.
void validate_mre_dimensions(const std::vector<GeneIndex>& genes,
                             int taxon_count) {
  if (taxon_count < 0 ||
      taxon_count > (std::numeric_limits<int>::max() - 1) / 2) {
    throw std::length_error("MRE taxon count exceeds its integer domain");
  }
  if (genes.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("MRE gene count exceeds its integer domain");
  }
  const int max_nodes = std::max(2, 2 * taxon_count + 1);
  (void)checked_mre_flat_slots(
      max_nodes, static_cast<int>(genes.size()),
      "MRE flat state exceeds its integer index domain");
  constexpr std::size_t kMaxEulerNodes =
      (static_cast<std::size_t>(std::numeric_limits<int>::max()) + 1u) / 2u;
  for (const GeneIndex& gene : genes) {
    if (gene.nodes.size() > kMaxEulerNodes ||
        gene.taxon_to_node.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::length_error("MRE gene index exceeds its integer domain");
    }
  }
}

using MreEntry = SplitAggregate256;

bool entry_less(const MreEntry& a, const MreEntry& b) {
  return split_priority_less_by_key(a, b);
}

// Counts the splits of genes [begin_gene, end_gene) and returns them in
// priority order. A gene edge is eligible when the gene has at least four
// distinct taxa, the side away from the reference root holds at least two of
// them and leaves at least two on the other side, and that side holds all or
// none of the copies of every duplicated taxon; `ids` gives the side's exact
// identifier, 0 when it splits some taxon's copies. Records are sorted by
// (key, gene, edge), a total order, so the entries do not depend on how the
// genes are chunked across workers.
std::vector<MreEntry> count_mre_entries_cpu_range(
    const std::vector<GeneIndex>& genes,
    int taxon_count,
    const GeneClusterIds& ids,
    int begin_gene,
    int end_gene) {
  begin_gene = std::max(0, std::min(begin_gene,
                                    static_cast<int>(genes.size())));
  end_gene = std::max(begin_gene,
                      std::min(end_gene, static_cast<int>(genes.size())));
  const int gene_count = end_gene - begin_gene;
  if (gene_count <= 0) return {};
  const int chunks =
      std::max(1, std::min(gene_count, current_inner_parallel_worker_limit()));
  std::vector<std::vector<SplitOccurrence256>> local_records(
      static_cast<std::size_t>(chunks));

  parallel_for_inner(0, chunks, [&](int chunk) {
    const int begin = begin_gene + static_cast<int>(
        static_cast<int64_t>(gene_count) * chunk / chunks);
    const int end = begin_gene + static_cast<int>(
        static_cast<int64_t>(gene_count) * (chunk + 1) / chunks);
    auto& local = local_records[static_cast<std::size_t>(chunk)];
    std::size_t capacity = 0;
    for (int g = begin; g < end; ++g) {
      const std::size_t nodes = genes[static_cast<std::size_t>(g)].nodes.size();
      if (nodes > std::numeric_limits<std::size_t>::max() - capacity) {
        throw std::length_error(
            "MRE split records exceed the addressable size domain");
      }
      capacity += nodes;
    }
    local.reserve(capacity);
    std::vector<uint32_t> subtree_size;
    std::vector<uint32_t> one_taxon;
    std::vector<int> representative;
    GeneTraversal reference_scratch;

    for (int g = begin; g < end; ++g) {
      const GeneIndex& gene = genes[static_cast<std::size_t>(g)];
      if (gene.nodes.empty()) {
        progress_mre_count_gene();
        continue;
      }
      const GeneTraversal& tr =
          gene_reference_traversal(gene, reference_scratch);
      const int node_count = static_cast<int>(gene.nodes.size());
      // The ids certify all-or-none copies, so only the representative
      // counts and the first taxon in incidence order are needed here. The
      // representative of a taxon is its last node, as in the gene's
      // reference root. These passes take O(N + kn) time and words over the
      // chunks.
      representative.assign(static_cast<std::size_t>(taxon_count), -1);
      for (int node = 0; node < node_count; ++node) {
        const int taxon = gene.nodes[static_cast<std::size_t>(node)].taxon();
        if (taxon >= 0) representative[static_cast<std::size_t>(taxon)] = node;
      }
      uint32_t distinct_taxa = 0;
      int reference = -1;
      for (int node : representative) {
        if (node < 0) continue;
        if (reference < 0) reference = node;
        ++distinct_taxa;
      }
      if (reference >= 0 && gene.reference_root != reference) {
        throw std::logic_error("gene split traversal is not canonically rooted");
      }
      subtree_size.assign(static_cast<std::size_t>(node_count), 0);
      one_taxon.assign(static_cast<std::size_t>(node_count), UINT32_MAX);
      for (auto it = tr.preorder.rbegin(); it != tr.preorder.rend(); ++it) {
        const int node = *it;
        const std::size_t v = static_cast<std::size_t>(node);
        const int taxon = gene.nodes[v].taxon();
        if (taxon >= 0 &&
            representative[static_cast<std::size_t>(taxon)] == node) {
          subtree_size[v] = 1;
          one_taxon[v] = static_cast<uint32_t>(taxon);
        }
        for (int nb : gene_neighbors(gene, node)) {
          const std::size_t c = static_cast<std::size_t>(nb);
          if (tr.parent[c] != node) continue;
          subtree_size[v] += subtree_size[c];
          if (one_taxon[v] == UINT32_MAX) one_taxon[v] = one_taxon[c];
        }
      }
      uint32_t edge_id = 0;
      for (int child : tr.preorder) {
        if (child == gene.reference_root) continue;
        const uint32_t id = ids.at(g, child);
        const std::size_t v = static_cast<std::size_t>(child);
        const uint32_t size = subtree_size[v];
        if (id != 0 && distinct_taxa >= 4 && size >= 2 &&
            size <= distinct_taxa - 2) {
          SplitKey256 key;
          key.size = size;
          key.h = {id, 0, 0, 0};
          local.push_back(SplitOccurrence256{
              key, static_cast<uint32_t>(g), edge_id,
              static_cast<uint32_t>(child), one_taxon[v]});
        }
        ++edge_id;
      }
      progress_mre_count_gene();
    }
  }, 1);

  SplitRadixStats radix_stats;
  std::vector<MreEntry> entries =
      aggregate_split_occurrences256(local_records, radix_stats);
  radix_sort_split_priority256(entries, radix_stats);
  return entries;
}

std::vector<MreEntry> count_mre_entries(
    const std::vector<GeneIndex>& genes,
    int taxon_count) {
  const GeneClusterIds ids = exact_gene_cluster_ids(genes, taxon_count);
  return count_mre_entries_cpu_range(genes, taxon_count, ids, 0,
                                     static_cast<int>(genes.size()));
}

// Link-cut tree mirroring the starting tree, used for root-to-leaf path
// queries: path_size and kth_on_path (the k-th node from root, 0-based), each
// O(log n) amortized.
class LinkCutTree {
 public:
  explicit LinkCutTree(int n = 0) { reset(n); }

  void reset(int n) {
    if (n < 0) throw std::invalid_argument("negative Link-Cut tree size");
    ch_.assign(static_cast<std::size_t>(n), {-1, -1});
    parent_.assign(static_cast<std::size_t>(n), -1);
    rev_.assign(static_cast<std::size_t>(n), 0);
    size_.assign(static_cast<std::size_t>(n), 1);
    push_stack_.clear();
    push_stack_.reserve(static_cast<std::size_t>(n));
  }

  void reserve(int n) {
    if (n < 0) throw std::invalid_argument("negative Link-Cut tree capacity");
    const std::size_t capacity = static_cast<std::size_t>(n);
    ch_.reserve(capacity);
    parent_.reserve(capacity);
    rev_.reserve(capacity);
    size_.reserve(capacity);
    push_stack_.reserve(capacity);
  }

  void add_node() {
    ch_.push_back({-1, -1});
    parent_.push_back(-1);
    rev_.push_back(0);
    size_.push_back(1);
  }

  void link(int a, int b) {
    makeroot(a);
    if (findroot(b) != a) parent_[static_cast<std::size_t>(a)] = b;
  }

  void cut(int a, int b) {
    makeroot(a);
    access(b);
    splay(b);
    if (ch_[static_cast<std::size_t>(b)][0] == a) {
      ch_[static_cast<std::size_t>(b)][0] = -1;
      parent_[static_cast<std::size_t>(a)] = -1;
      pull(b);
    }
  }

  int path_size(int root, int node) {
    makeroot(root);
    access(node);
    splay(node);
    return size_[static_cast<std::size_t>(node)];
  }

  int kth_on_path(int root, int node, int k) {
    makeroot(root);
    access(node);
    splay(node);
    int cur = node;
    while (cur >= 0) {
      push(cur);
      const int left = ch_[static_cast<std::size_t>(cur)][0];
      const int left_size =
          left >= 0 ? size_[static_cast<std::size_t>(left)] : 0;
      if (k < left_size) {
        cur = left;
      } else if (k == left_size) {
        splay(cur);
        return cur;
      } else {
        k -= left_size + 1;
        cur = ch_[static_cast<std::size_t>(cur)][1];
      }
    }
    return node;
  }

 private:
  std::vector<std::array<int, 2>> ch_;
  std::vector<int> parent_;
  std::vector<uint8_t> rev_;
  std::vector<int> size_;
  std::vector<int> push_stack_;

  bool is_root(int x) const {
    const int p = parent_[static_cast<std::size_t>(x)];
    return p < 0 ||
           (ch_[static_cast<std::size_t>(p)][0] != x &&
            ch_[static_cast<std::size_t>(p)][1] != x);
  }

  void pull(int x) {
    int size = 1;
    for (int side = 0; side < 2; ++side) {
      const int c = ch_[static_cast<std::size_t>(x)][side];
      if (c >= 0) size += size_[static_cast<std::size_t>(c)];
    }
    size_[static_cast<std::size_t>(x)] = size;
  }

  void push(int x) {
    if (x < 0 || rev_[static_cast<std::size_t>(x)] == 0) return;
    std::swap(ch_[static_cast<std::size_t>(x)][0],
              ch_[static_cast<std::size_t>(x)][1]);
    for (int side = 0; side < 2; ++side) {
      const int c = ch_[static_cast<std::size_t>(x)][side];
      if (c >= 0) rev_[static_cast<std::size_t>(c)] ^= 1;
    }
    rev_[static_cast<std::size_t>(x)] = 0;
  }

  void push_path(int x) {
    push_stack_.clear();
    push_stack_.push_back(x);
    while (!is_root(x)) {
      x = parent_[static_cast<std::size_t>(x)];
      push_stack_.push_back(x);
    }
    for (auto it = push_stack_.rbegin(); it != push_stack_.rend(); ++it) {
      push(*it);
    }
  }

  void rotate(int x) {
    const int p = parent_[static_cast<std::size_t>(x)];
    const int g = parent_[static_cast<std::size_t>(p)];
    const int side = ch_[static_cast<std::size_t>(p)][1] == x ? 1 : 0;
    const int b = ch_[static_cast<std::size_t>(x)][side ^ 1];
    if (!is_root(p)) {
      ch_[static_cast<std::size_t>(g)]
         [ch_[static_cast<std::size_t>(g)][1] == p ? 1 : 0] = x;
    }
    parent_[static_cast<std::size_t>(x)] = g;
    ch_[static_cast<std::size_t>(x)][side ^ 1] = p;
    parent_[static_cast<std::size_t>(p)] = x;
    ch_[static_cast<std::size_t>(p)][side] = b;
    if (b >= 0) parent_[static_cast<std::size_t>(b)] = p;
    pull(p);
    pull(x);
  }

  void splay(int x) {
    push_path(x);
    while (!is_root(x)) {
      const int p = parent_[static_cast<std::size_t>(x)];
      const int g = parent_[static_cast<std::size_t>(p)];
      if (!is_root(p)) {
        const bool zigzig =
            (ch_[static_cast<std::size_t>(p)][0] == x) ==
            (ch_[static_cast<std::size_t>(g)][0] == p);
        rotate(zigzig ? p : x);
      }
      rotate(x);
    }
  }

  void access(int x) {
    int last = -1;
    for (int y = x; y >= 0; y = parent_[static_cast<std::size_t>(y)]) {
      splay(y);
      ch_[static_cast<std::size_t>(y)][1] = last;
      pull(y);
      last = y;
    }
    splay(x);
  }

  void makeroot(int x) {
    access(x);
    rev_[static_cast<std::size_t>(x)] ^= 1;
    push(x);
  }

  int findroot(int x) {
    access(x);
    while (true) {
      push(x);
      const int left = ch_[static_cast<std::size_t>(x)][0];
      if (left < 0) break;
      x = left;
    }
    splay(x);
    return x;
  }
};

// Aggregate over a set of children of one starting-tree node, for one gene.
struct ChildAgg {
  int count = 0;
  // The children's taxa present in the gene. Children of one starting-tree
  // node hold disjoint taxa, so the sum is at most the taxon count.
  uint32_t leaf_sum = 0;
  // The children's gene LCAs with the smallest and largest preorder time. The
  // LCA of these two is the LCA of the whole set.
  int min_node = -1;
  int max_node = -1;
};

// (preorder time of the child's LCA in the gene, child). A child with no taxa
// in the gene has time -1 and so lies outside every subtree range.
struct ChildKey {
  int tin = 0;
  int child = 0;
};

bool key_less(const ChildKey& a, const ChildKey& b) {
  return a.tin != b.tin ? a.tin < b.tin : a.child < b.child;
}

bool key_equal(const ChildKey& a, const ChildKey& b) {
  return a.tin == b.tin && a.child == b.child;
}

// For every (starting-tree node, gene) pair, an AVL tree of the node's
// children ordered by ChildKey with subtree aggregates. The children whose
// gene LCA lies in the subtree of a gene node form the key range
// [tin, tout) of that node: its aggregate costs O(log d) and its enumeration
// O(log d + output), for d children.
class ChildRangeForest {
 public:
  struct InitialChild {
    int child = -1;
    int key_tin = 0;
    uint32_t leaf_count = 0;
    int lca_node = -1;
  };

  ChildRangeForest(int max_nodes,
                   int genes,
                   const std::vector<StaticLcaIndex>& lca)
      : max_nodes_(max_nodes), genes_(genes), lca_(lca) {
    const std::size_t flat_slots = checked_mre_flat_slots(
        max_nodes_, genes_,
        "MRE child index exceeds its integer index domain");
    roots_.assign(flat_slots, -1);
    // Flat node pool partitioned into fixed per-gene segments
    // [gene*max_nodes_, (gene+1)*max_nodes_). Pre-sized and never reallocated,
    // so distinct genes own disjoint index ranges and can be updated in
    // parallel without data races. Every starting-tree node but the root is
    // the child of one parent, so a gene has fewer than max_nodes_ live AVL
    // nodes and its segment never overflows.
    nodes_.assign(flat_slots, Node{});
    free_.assign(static_cast<std::size_t>(genes_), {});
    next_.assign(static_cast<std::size_t>(genes_), 0);
    operations_.assign(static_cast<std::size_t>(genes_), 0);
    rotations_.assign(static_cast<std::size_t>(genes_), 0);
    max_height_.assign(static_cast<std::size_t>(genes_), 0);
    height_violations_.assign(static_cast<std::size_t>(genes_), 0);
    enumeration_visits_.assign(static_cast<std::size_t>(genes_), 0);
    for (int gene = 0; gene < genes_; ++gene) {
      next_[static_cast<std::size_t>(gene)] = gene * max_nodes_;
    }
  }

  void insert_child(int parent,
                    int gene,
                    int child,
                    int key_tin,
                    uint32_t leaf_count,
                    int lca_node) {
    const int node = allocate(gene);
    Node& n = nodes_[static_cast<std::size_t>(node)];
    n.left = n.right = -1;
    n.height = 1;
    n.key = ChildKey{key_tin, child};
    n.leaf_count = leaf_count;
    n.lca_node = lca_node;
    update(node, gene);
    int& root = roots_[root_index(parent, gene)];
    root = insert(root, node, gene);
    ++operations_[static_cast<std::size_t>(gene)];
    note_root(root, gene);
  }

  // Builds the index of a parent that has none yet by sorting and a balanced
  // build, O(d log d).
  void build_children(int parent,
                      int gene,
                      std::vector<InitialChild>& children) {
    std::sort(children.begin(), children.end(),
              [](const InitialChild& a, const InitialChild& b) {
                if (a.key_tin != b.key_tin) return a.key_tin < b.key_tin;
                return a.child < b.child;
              });
    int& root = roots_[root_index(parent, gene)];
    if (root >= 0) {
      throw std::logic_error("MRE child index rebuilt over a live root");
    }
    std::vector<int> ordered;
    ordered.reserve(children.size());
    for (const InitialChild& child : children) {
      const int node = allocate(gene);
      Node& n = nodes_[static_cast<std::size_t>(node)];
      n.left = n.right = -1;
      n.height = 1;
      n.key = ChildKey{child.key_tin, child.child};
      n.leaf_count = child.leaf_count;
      n.lca_node = child.lca_node;
      n.agg = ChildAgg{1, child.leaf_count, child.lca_node, child.lca_node};
      ordered.push_back(node);
    }
    root = build_balanced(ordered, 0, static_cast<int>(ordered.size()), gene);
    ++operations_[static_cast<std::size_t>(gene)];
    note_root(root, gene);
  }

  void erase_child(int parent, int gene, int child, int key_tin) {
    int& root = roots_[root_index(parent, gene)];
    bool erased = false;
    root = erase(root, ChildKey{key_tin, child}, gene, erased);
    if (!erased) {
      throw std::logic_error("MRE child index erased a missing key");
    }
    ++operations_[static_cast<std::size_t>(gene)];
    note_root(root, gene);
  }

  ChildAgg range_query(int parent, int gene, int tin_lo, int tin_hi) {
    ++operations_[static_cast<std::size_t>(gene)];
    return range_aggregate(root(parent, gene), ChildKey{tin_lo, INT_MIN},
                           ChildKey{tin_hi, INT_MIN}, gene);
  }

  std::vector<int> enumerate_range(int parent,
                                   int gene,
                                   int tin_lo,
                                   int tin_hi) {
    ++operations_[static_cast<std::size_t>(gene)];
    std::vector<int> out;
    enumerate_range_impl(root(parent, gene), ChildKey{tin_lo, INT_MIN},
                         ChildKey{tin_hi, INT_MIN}, gene, out);
    return out;
  }

  std::vector<int> enumerate_complement(int parent,
                                        int gene,
                                        int tin_lo,
                                        int tin_hi) {
    ++operations_[static_cast<std::size_t>(gene)];
    std::vector<int> out;
    const int parent_root = root(parent, gene);
    enumerate_prefix_impl(parent_root, ChildKey{tin_lo, INT_MIN}, gene, out);
    enumerate_suffix_impl(parent_root, ChildKey{tin_hi, INT_MIN}, gene, out);
    return out;
  }

  uint32_t aggregate_leaf_sum(int parent, int gene) const {
    const int r = root(parent, gene);
    return r < 0 ? 0 : nodes_[static_cast<std::size_t>(r)].agg.leaf_sum;
  }

  int aggregate_lca(int parent, int gene) const {
    const int r = root(parent, gene);
    if (r < 0) return -1;
    const ChildAgg& agg = nodes_[static_cast<std::size_t>(r)].agg;
    return lca_[static_cast<std::size_t>(gene)].lca(agg.min_node,
                                                    agg.max_node);
  }

  long long operation_count() const {
    return std::accumulate(operations_.begin(), operations_.end(), 0LL);
  }

  long long rotation_count() const {
    return std::accumulate(rotations_.begin(), rotations_.end(), 0LL);
  }

  int maximum_height() const {
    return max_height_.empty()
               ? 0
               : *std::max_element(max_height_.begin(), max_height_.end());
  }

  long long height_violation_count() const {
    return std::accumulate(height_violations_.begin(),
                           height_violations_.end(), 0LL);
  }

  long long enumeration_visit_count() const {
    return std::accumulate(enumeration_visits_.begin(),
                           enumeration_visits_.end(), 0LL);
  }

  // 2 * ceil(log2(count + 1)) + 1, above the AVL height bound
  // 1.44 * log2(count + 2); checked after every update.
  static int height_limit_for_count(int count) {
    if (count <= 0) return 0;
    uint64_t capacity = 1;
    int ceil_log = 0;
    const uint64_t required = static_cast<uint64_t>(count) + 1;
    while (capacity < required) {
      capacity <<= 1;
      ++ceil_log;
    }
    return 2 * ceil_log + 1;
  }

  bool validate_root_for_tests(int parent, int gene) const {
    return validate_tree(root(parent, gene), gene).valid;
  }

 private:
  // key.child is the starting-tree child this node stands for.
  struct Node {
    int left = -1;
    int right = -1;
    int height = 1;
    ChildKey key;
    uint32_t leaf_count = 0;
    int lca_node = -1;
    ChildAgg agg;
  };

  int max_nodes_ = 0;
  int genes_ = 0;
  const std::vector<StaticLcaIndex>& lca_;
  std::vector<int> roots_;
  std::vector<Node> nodes_;
  std::vector<std::vector<int>> free_;
  std::vector<int> next_;
  std::vector<long long> operations_;
  std::vector<long long> rotations_;
  std::vector<int> max_height_;
  std::vector<long long> height_violations_;
  std::vector<long long> enumeration_visits_;

  std::size_t root_index(int parent, int gene) const {
    // Parent-major: accepting a split visits every gene's root for a fixed
    // parent, so those roots are contiguous. The per-gene pools in nodes_ are
    // what is updated in parallel, and they are disjoint.
    return static_cast<std::size_t>(parent) * static_cast<std::size_t>(genes_) +
           static_cast<std::size_t>(gene);
  }

  int root(int parent, int gene) const {
    return roots_[root_index(parent, gene)];
  }

  int allocate(int gene) {
    std::vector<int>& fr = free_[static_cast<std::size_t>(gene)];
    if (!fr.empty()) {
      const int node = fr.back();
      fr.pop_back();
      return node;
    }
    const int node = next_[static_cast<std::size_t>(gene)]++;
    const int segment_end = (gene + 1) * max_nodes_;
    if (node >= segment_end) {
      throw std::logic_error("MRE child index exceeded its flat segment");
    }
    return node;
  }

  void release(int gene, int node) {
    if (node >= 0) free_[static_cast<std::size_t>(gene)].push_back(node);
  }

  int combine_min_tin_node(int gene, int a, int b) const {
    if (a < 0) return b;
    if (b < 0) return a;
    const StaticLcaIndex& lca = lca_[static_cast<std::size_t>(gene)];
    return lca.tin[static_cast<std::size_t>(a)] <=
                   lca.tin[static_cast<std::size_t>(b)]
               ? a
               : b;
  }

  int combine_max_tin_node(int gene, int a, int b) const {
    if (a < 0) return b;
    if (b < 0) return a;
    const StaticLcaIndex& lca = lca_[static_cast<std::size_t>(gene)];
    return lca.tin[static_cast<std::size_t>(a)] >=
                   lca.tin[static_cast<std::size_t>(b)]
               ? a
               : b;
  }

  void update(int node, int gene) {
    if (node < 0) return;
    Node& n = nodes_[static_cast<std::size_t>(node)];
    n.height = 1 + std::max(node_height(n.left), node_height(n.right));
    n.agg.count = 1;
    n.agg.leaf_sum = n.leaf_count;
    n.agg.min_node = n.lca_node;
    n.agg.max_node = n.lca_node;
    if (n.left >= 0) {
      const ChildAgg& left = nodes_[static_cast<std::size_t>(n.left)].agg;
      n.agg.count += left.count;
      n.agg.leaf_sum += left.leaf_sum;
      n.agg.min_node = combine_min_tin_node(gene, n.agg.min_node,
                                            left.min_node);
      n.agg.max_node = combine_max_tin_node(gene, n.agg.max_node,
                                            left.max_node);
    }
    if (n.right >= 0) {
      const ChildAgg& right = nodes_[static_cast<std::size_t>(n.right)].agg;
      n.agg.count += right.count;
      n.agg.leaf_sum += right.leaf_sum;
      n.agg.min_node = combine_min_tin_node(gene, n.agg.min_node,
                                            right.min_node);
      n.agg.max_node = combine_max_tin_node(gene, n.agg.max_node,
                                            right.max_node);
    }
  }

  ChildAgg aggregate(int node) const {
    return node >= 0 ? nodes_[static_cast<std::size_t>(node)].agg : ChildAgg{};
  }

  int node_height(int node) const {
    return node >= 0 ? nodes_[static_cast<std::size_t>(node)].height : 0;
  }

  int balance_factor(int node) const {
    if (node < 0) return 0;
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    return node_height(n.left) - node_height(n.right);
  }

  ChildAgg combine_aggregates(const ChildAgg& a,
                              const ChildAgg& b,
                              int gene) const {
    ChildAgg out;
    out.count = a.count + b.count;
    out.leaf_sum = a.leaf_sum + b.leaf_sum;
    out.min_node = combine_min_tin_node(gene, a.min_node, b.min_node);
    out.max_node = combine_max_tin_node(gene, a.max_node, b.max_node);
    return out;
  }

  ChildAgg self_aggregate(int node) const {
    if (node < 0) return ChildAgg{};
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    return ChildAgg{1, n.leaf_count, n.lca_node, n.lca_node};
  }

  int rotate_left(int root, int gene) {
    Node& r = nodes_[static_cast<std::size_t>(root)];
    const int next_root = r.right;
    if (next_root < 0) {
      throw std::logic_error("MRE invalid left AVL rotation");
    }
    Node& next = nodes_[static_cast<std::size_t>(next_root)];
    r.right = next.left;
    next.left = root;
    update(root, gene);
    update(next_root, gene);
    ++rotations_[static_cast<std::size_t>(gene)];
    return next_root;
  }

  int rotate_right(int root, int gene) {
    Node& r = nodes_[static_cast<std::size_t>(root)];
    const int next_root = r.left;
    if (next_root < 0) {
      throw std::logic_error("MRE invalid right AVL rotation");
    }
    Node& next = nodes_[static_cast<std::size_t>(next_root)];
    r.left = next.right;
    next.right = root;
    update(root, gene);
    update(next_root, gene);
    ++rotations_[static_cast<std::size_t>(gene)];
    return next_root;
  }

  int rebalance(int root, int gene) {
    if (root < 0) return -1;
    update(root, gene);
    const int balance = balance_factor(root);
    if (balance > 1) {
      Node& r = nodes_[static_cast<std::size_t>(root)];
      if (balance_factor(r.left) < 0) {
        r.left = rotate_left(r.left, gene);
      }
      root = rotate_right(root, gene);
    } else if (balance < -1) {
      Node& r = nodes_[static_cast<std::size_t>(root)];
      if (balance_factor(r.right) > 0) {
        r.right = rotate_right(r.right, gene);
      }
      root = rotate_left(root, gene);
    }
    if (std::abs(balance_factor(root)) > 1) {
      ++height_violations_[static_cast<std::size_t>(gene)];
      throw std::logic_error("MRE AVL rebalance invariant failed");
    }
    return root;
  }

  int insert(int root, int node, int gene) {
    if (root < 0) return node;
    Node& r = nodes_[static_cast<std::size_t>(root)];
    const ChildKey key = nodes_[static_cast<std::size_t>(node)].key;
    if (key_equal(key, r.key)) {
      throw std::logic_error("MRE child index inserted a duplicate key");
    }
    if (key_less(key, r.key)) {
      r.left = insert(r.left, node, gene);
    } else {
      r.right = insert(r.right, node, gene);
    }
    return rebalance(root, gene);
  }

  int build_balanced(const std::vector<int>& ordered,
                     int begin,
                     int end,
                     int gene) {
    if (begin >= end) return -1;
    const int middle = begin + (end - begin) / 2;
    const int root = ordered[static_cast<std::size_t>(middle)];
    Node& r = nodes_[static_cast<std::size_t>(root)];
    r.left = build_balanced(ordered, begin, middle, gene);
    r.right = build_balanced(ordered, middle + 1, end, gene);
    update(root, gene);
    return root;
  }

  int detach_min(int root, int& minimum, int gene) {
    Node& r = nodes_[static_cast<std::size_t>(root)];
    if (r.left < 0) {
      minimum = root;
      return r.right;
    }
    r.left = detach_min(r.left, minimum, gene);
    return rebalance(root, gene);
  }

  int erase(int root,
            const ChildKey& key,
            int gene,
            bool& erased) {
    if (root < 0) return -1;
    Node& r = nodes_[static_cast<std::size_t>(root)];
    if (key_equal(r.key, key)) {
      erased = true;
      const int left = r.left;
      const int right = r.right;
      if (left < 0 || right < 0) {
        const int replacement = left >= 0 ? left : right;
        release(gene, root);
        return replacement;
      }
      int successor = -1;
      const int remaining_right = detach_min(right, successor, gene);
      Node& replacement = nodes_[static_cast<std::size_t>(successor)];
      replacement.left = left;
      replacement.right = remaining_right;
      release(gene, root);
      return rebalance(successor, gene);
    }
    if (key_less(key, r.key)) {
      r.left = erase(r.left, key, gene, erased);
    } else {
      r.right = erase(r.right, key, gene, erased);
    }
    return erased ? rebalance(root, gene) : root;
  }

  ChildAgg suffix_aggregate(int node,
                            const ChildKey& lo,
                            int gene) const {
    if (node < 0) return ChildAgg{};
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    if (key_less(n.key, lo)) {
      return suffix_aggregate(n.right, lo, gene);
    }
    ChildAgg out = aggregate(n.right);
    out = combine_aggregates(out, self_aggregate(node), gene);
    return combine_aggregates(out, suffix_aggregate(n.left, lo, gene),
                              gene);
  }

  ChildAgg prefix_aggregate(int node,
                            const ChildKey& hi,
                            int gene) const {
    if (node < 0) return ChildAgg{};
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    if (!key_less(n.key, hi)) {
      return prefix_aggregate(n.left, hi, gene);
    }
    ChildAgg out = aggregate(n.left);
    out = combine_aggregates(out, self_aggregate(node), gene);
    return combine_aggregates(out, prefix_aggregate(n.right, hi, gene),
                              gene);
  }

  ChildAgg range_aggregate(int root,
                           const ChildKey& lo,
                           const ChildKey& hi,
                           int gene) const {
    int split_node = root;
    while (split_node >= 0) {
      const Node& n = nodes_[static_cast<std::size_t>(split_node)];
      if (key_less(n.key, lo)) {
        split_node = n.right;
      } else if (!key_less(n.key, hi)) {
        split_node = n.left;
      } else {
        ChildAgg out = suffix_aggregate(n.left, lo, gene);
        out = combine_aggregates(out, self_aggregate(split_node), gene);
        return combine_aggregates(out,
                                  prefix_aggregate(n.right, hi, gene), gene);
      }
    }
    return ChildAgg{};
  }

  struct ValidationResult {
    bool valid = true;
    bool empty = true;
    ChildKey min_key;
    ChildKey max_key;
    int height = 0;
    ChildAgg agg;
  };

  static bool aggregate_equal(const ChildAgg& a, const ChildAgg& b) {
    return a.count == b.count && a.leaf_sum == b.leaf_sum &&
           a.min_node == b.min_node && a.max_node == b.max_node;
  }

  ValidationResult validate_tree(int node, int gene) const {
    if (node < 0) return ValidationResult{};
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    const ValidationResult left = validate_tree(n.left, gene);
    const ValidationResult right = validate_tree(n.right, gene);
    ValidationResult out;
    out.empty = false;
    out.height = 1 + std::max(left.height, right.height);
    out.min_key = left.empty ? n.key : left.min_key;
    out.max_key = right.empty ? n.key : right.max_key;
    ChildAgg expected = combine_aggregates(left.agg, self_aggregate(node),
                                           gene);
    expected = combine_aggregates(expected, right.agg, gene);
    out.agg = expected;
    const bool ordered_left = left.empty || key_less(left.max_key, n.key);
    const bool ordered_right = right.empty || key_less(n.key, right.min_key);
    const bool balanced = std::abs(left.height - right.height) <= 1;
    out.valid = left.valid && right.valid && ordered_left && ordered_right &&
                balanced && n.height == out.height &&
                aggregate_equal(n.agg, expected) &&
                out.height <= height_limit_for_count(expected.count);
    return out;
  }

  void note_root(int root, int gene) {
    const int height = node_height(root);
    int& maximum = max_height_[static_cast<std::size_t>(gene)];
    maximum = std::max(maximum, height);
    const int count = root < 0
                          ? 0
                          : nodes_[static_cast<std::size_t>(root)].agg.count;
    if (height > height_limit_for_count(count)) {
      ++height_violations_[static_cast<std::size_t>(gene)];
      throw std::logic_error("MRE AVL height check failed");
    }
  }

  void enumerate_range_impl(int node,
                            const ChildKey& lo,
                            const ChildKey& hi,
                            int gene,
                            std::vector<int>& out) {
    if (node < 0) return;
    ++enumeration_visits_[static_cast<std::size_t>(gene)];
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    if (key_less(n.key, lo)) {
      enumerate_range_impl(n.right, lo, hi, gene, out);
      return;
    }
    if (!key_less(n.key, hi)) {
      enumerate_range_impl(n.left, lo, hi, gene, out);
      return;
    }
    enumerate_range_impl(n.left, lo, hi, gene, out);
    out.push_back(n.key.child);
    enumerate_range_impl(n.right, lo, hi, gene, out);
  }

  void enumerate_all_impl(int node,
                          int gene,
                          std::vector<int>& out) {
    if (node < 0) return;
    ++enumeration_visits_[static_cast<std::size_t>(gene)];
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    enumerate_all_impl(n.left, gene, out);
    out.push_back(n.key.child);
    enumerate_all_impl(n.right, gene, out);
  }

  void enumerate_prefix_impl(int node,
                             const ChildKey& bound,
                             int gene,
                             std::vector<int>& out) {
    if (node < 0) return;
    ++enumeration_visits_[static_cast<std::size_t>(gene)];
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    if (!key_less(n.key, bound)) {
      enumerate_prefix_impl(n.left, bound, gene, out);
      return;
    }
    enumerate_all_impl(n.left, gene, out);
    out.push_back(n.key.child);
    enumerate_prefix_impl(n.right, bound, gene, out);
  }

  void enumerate_suffix_impl(int node,
                             const ChildKey& bound,
                             int gene,
                             std::vector<int>& out) {
    if (node < 0) return;
    ++enumeration_visits_[static_cast<std::size_t>(gene)];
    const Node& n = nodes_[static_cast<std::size_t>(node)];
    if (key_less(n.key, bound)) {
      enumerate_suffix_impl(n.right, bound, gene, out);
      return;
    }
    enumerate_suffix_impl(n.left, bound, gene, out);
    out.push_back(n.key.child);
    enumerate_all_impl(n.right, gene, out);
  }
};

// Mutable starting-tree topology. Each node's children form a linked list of
// slots, one slot per live parent-child edge, drawn from a shared pool of
// max_nodes slots and recycled on detach, so space stays O(n) however
// accepted clusters nest. Detach moves the last child into the removed
// child's slot (swap-with-last); this fixes the child order that
// materialize_mre_tree emits.
class FlatChildTopology {
 public:
  explicit FlatChildTopology(int max_nodes = 0) { reset(max_nodes); }

  void reset(int max_nodes) {
    if (max_nodes < 0) {
      throw std::invalid_argument("negative MRE topology size");
    }
    max_nodes_ = max_nodes;
    head_.assign(static_cast<std::size_t>(max_nodes_), -1);
    tail_.assign(static_cast<std::size_t>(max_nodes_), -1);
    size_.assign(static_cast<std::size_t>(max_nodes_), 0);
    child_slot_.assign(static_cast<std::size_t>(max_nodes_), -1);
    slots_.assign(static_cast<std::size_t>(max_nodes_), Slot{});
    free_.clear();
    free_.reserve(static_cast<std::size_t>(max_nodes_));
    next_slot_ = 0;
    live_slots_ = 0;
    peak_live_slots_ = 0;
  }

  void reset_node(int node) {
    require_node(node);
    if (head_[static_cast<std::size_t>(node)] >= 0 ||
        tail_[static_cast<std::size_t>(node)] >= 0 ||
        size_[static_cast<std::size_t>(node)] != 0 ||
        child_slot_[static_cast<std::size_t>(node)] >= 0) {
      throw std::logic_error("MRE topology reused a live node");
    }
  }

  void attach(int parent, int child) {
    require_node(parent);
    require_node(child);
    if (parent == child || child_slot_[static_cast<std::size_t>(child)] >= 0) {
      throw std::logic_error("MRE topology attached a non-free child");
    }
    const int slot = allocate_slot();
    const int old_tail = tail_[static_cast<std::size_t>(parent)];
    Slot& record = slots_[static_cast<std::size_t>(slot)];
    record.parent = parent;
    record.child = child;
    record.prev = old_tail;
    record.next = -1;
    if (old_tail >= 0) {
      slots_[static_cast<std::size_t>(old_tail)].next = slot;
    } else {
      head_[static_cast<std::size_t>(parent)] = slot;
    }
    tail_[static_cast<std::size_t>(parent)] = slot;
    child_slot_[static_cast<std::size_t>(child)] = slot;
    ++size_[static_cast<std::size_t>(parent)];
    ++live_slots_;
    peak_live_slots_ = std::max(peak_live_slots_, live_slots_);
  }

  void detach(int parent, int child) {
    require_node(parent);
    require_node(child);
    const int slot = child_slot_[static_cast<std::size_t>(child)];
    const int last = tail_[static_cast<std::size_t>(parent)];
    if (slot < 0 || last < 0 ||
        slots_[static_cast<std::size_t>(slot)].parent != parent) {
      throw std::logic_error("MRE topology detached a missing child");
    }

    if (slot != last) {
      const int replacement_child =
          slots_[static_cast<std::size_t>(last)].child;
      slots_[static_cast<std::size_t>(slot)].child = replacement_child;
      child_slot_[static_cast<std::size_t>(replacement_child)] = slot;
    }
    const int new_tail = slots_[static_cast<std::size_t>(last)].prev;
    if (new_tail >= 0) {
      slots_[static_cast<std::size_t>(new_tail)].next = -1;
    } else {
      head_[static_cast<std::size_t>(parent)] = -1;
    }
    tail_[static_cast<std::size_t>(parent)] = new_tail;
    child_slot_[static_cast<std::size_t>(child)] = -1;
    --size_[static_cast<std::size_t>(parent)];
    --live_slots_;
    release_slot(last);
  }

  int size(int parent) const {
    require_node(parent);
    return size_[static_cast<std::size_t>(parent)];
  }

  template <class F>
  void for_each_child(int parent, F&& body) const {
    require_node(parent);
    for (int slot = head_[static_cast<std::size_t>(parent)]; slot >= 0;
         slot = slots_[static_cast<std::size_t>(slot)].next) {
      body(slots_[static_cast<std::size_t>(slot)].child);
    }
  }

  std::vector<int> children(int parent) const {
    std::vector<int> out;
    out.reserve(static_cast<std::size_t>(size(parent)));
    for_each_child(parent, [&](int child) { out.push_back(child); });
    return out;
  }

  int allocated_slots() const { return next_slot_; }
  int live_slots() const { return live_slots_; }
  int peak_live_slots() const { return peak_live_slots_; }
  int slot_limit() const { return max_nodes_; }

  bool validate_for_tests() const {
    std::vector<uint8_t> seen_slot(static_cast<std::size_t>(next_slot_), 0);
    std::vector<uint8_t> seen_child(static_cast<std::size_t>(max_nodes_), 0);
    int listed_slots = 0;
    for (int parent = 0; parent < max_nodes_; ++parent) {
      int previous = -1;
      int listed_for_parent = 0;
      for (int slot = head_[static_cast<std::size_t>(parent)]; slot >= 0;
           slot = slots_[static_cast<std::size_t>(slot)].next) {
        if (slot >= next_slot_ ||
            seen_slot[static_cast<std::size_t>(slot)] != 0) {
          return false;
        }
        seen_slot[static_cast<std::size_t>(slot)] = 1;
        const Slot& record = slots_[static_cast<std::size_t>(slot)];
        if (record.parent != parent || record.prev != previous ||
            record.child < 0 || record.child >= max_nodes_ ||
            seen_child[static_cast<std::size_t>(record.child)] != 0 ||
            child_slot_[static_cast<std::size_t>(record.child)] != slot) {
          return false;
        }
        seen_child[static_cast<std::size_t>(record.child)] = 1;
        previous = slot;
        ++listed_for_parent;
        ++listed_slots;
        if (listed_slots > live_slots_) return false;
      }
      if (previous != tail_[static_cast<std::size_t>(parent)] ||
          listed_for_parent != size_[static_cast<std::size_t>(parent)] ||
          ((listed_for_parent == 0) !=
           (head_[static_cast<std::size_t>(parent)] < 0))) {
        return false;
      }
    }
    if (listed_slots != live_slots_) return false;
    std::vector<uint8_t> seen_free(static_cast<std::size_t>(next_slot_), 0);
    for (int slot : free_) {
      if (slot < 0 || slot >= next_slot_ ||
          seen_slot[static_cast<std::size_t>(slot)] != 0 ||
          seen_free[static_cast<std::size_t>(slot)] != 0) {
        return false;
      }
      seen_free[static_cast<std::size_t>(slot)] = 1;
    }
    if (static_cast<int>(free_.size()) + live_slots_ != next_slot_) {
      return false;
    }
    for (int child = 0; child < max_nodes_; ++child) {
      const int slot = child_slot_[static_cast<std::size_t>(child)];
      if ((slot >= 0) !=
          (seen_child[static_cast<std::size_t>(child)] != 0)) {
        return false;
      }
    }
    return peak_live_slots_ <= max_nodes_ && next_slot_ <= max_nodes_;
  }

 private:
  struct Slot {
    int parent = -1;
    int child = -1;
    int prev = -1;
    int next = -1;
  };

  int max_nodes_ = 0;
  std::vector<int> head_;
  std::vector<int> tail_;
  std::vector<int> size_;
  std::vector<int> child_slot_;
  std::vector<Slot> slots_;
  std::vector<int> free_;
  int next_slot_ = 0;
  int live_slots_ = 0;
  int peak_live_slots_ = 0;

  void require_node(int node) const {
    if (node < 0 || node >= max_nodes_) {
      throw std::out_of_range("MRE topology node outside flat domain");
    }
  }

  int allocate_slot() {
    if (!free_.empty()) {
      const int slot = free_.back();
      free_.pop_back();
      return slot;
    }
    if (next_slot_ >= max_nodes_) {
      throw std::logic_error("MRE topology exceeded its linear slot bound");
    }
    return next_slot_++;
  }

  void release_slot(int slot) {
    slots_[static_cast<std::size_t>(slot)] = Slot{};
    free_.push_back(slot);
  }
};

// Starting tree under construction. Node 0 is the initial root and nodes
// 1..n the taxon leaves; each accepted split adds one node. leaf_count(v) is
// the number of v's taxa; lca_in_gene(v, g) is the LCA in gene g of v's taxa
// present in g, or -1 when there are none, and present_in_gene counts those
// taxa. dyn mirrors the topology for path queries.
struct MreState {
  int n = 0;
  int k = 0;
  int root = 0;
  int node_count = 0;
  int max_nodes = 0;
  std::vector<int> parent;
  std::vector<int> taxon;
  std::vector<uint32_t> leaf_count;
  std::vector<uint32_t> degree;
  FlatChildTopology topology;
  std::vector<int> taxon_to_leaf;
  std::vector<int> lca_flat;
  std::vector<StaticLcaIndex> static_lca;
  LinkCutTree dyn;
  std::unique_ptr<ChildRangeForest> child_index;

  int& lca_in_gene(int node, int gene) {
    return lca_flat[static_cast<std::size_t>(node) *
                        static_cast<std::size_t>(k) +
                    static_cast<std::size_t>(gene)];
  }

  int lca_in_gene(int node, int gene) const {
    return lca_flat[static_cast<std::size_t>(node) *
                        static_cast<std::size_t>(k) +
                    static_cast<std::size_t>(gene)];
  }

  int new_node() {
    if (node_count >= max_nodes) {
      throw std::runtime_error("MRE consensus exceeded laminar node bound");
    }
    const int node = node_count++;
    if (node >= static_cast<int>(dyn_size_initialized)) {
      dyn.add_node();
      ++dyn_size_initialized;
    }
    parent[static_cast<std::size_t>(node)] = -1;
    taxon[static_cast<std::size_t>(node)] = -1;
    leaf_count[static_cast<std::size_t>(node)] = 0;
    degree[static_cast<std::size_t>(node)] = 0;
    topology.reset_node(node);
    return node;
  }

 private:
  int dyn_size_initialized = 0;
};

MreState initialize_state(const std::vector<GeneIndex>& genes,
                            int taxon_count) {
  validate_mre_dimensions(genes, taxon_count);
  MreState st;
  st.n = taxon_count;
  st.k = static_cast<int>(genes.size());
  st.max_nodes = std::max(2, 2 * taxon_count + 1);
  st.parent.assign(static_cast<std::size_t>(st.max_nodes), -1);
  st.taxon.assign(static_cast<std::size_t>(st.max_nodes), -1);
  st.leaf_count.assign(static_cast<std::size_t>(st.max_nodes), 0);
  st.degree.assign(static_cast<std::size_t>(st.max_nodes), 0);
  st.topology.reset(st.max_nodes);
  st.taxon_to_leaf.assign(static_cast<std::size_t>(taxon_count), -1);
  st.lca_flat.assign(
      checked_mre_size_product(
          static_cast<std::size_t>(st.max_nodes),
          static_cast<std::size_t>(std::max(1, st.k)),
          "MRE LCA state exceeds the addressable size domain"),
      -1);
  st.static_lca.resize(static_cast<std::size_t>(st.k));
  parallel_for_inner(0, st.k, [&](int g) {
    st.static_lca[static_cast<std::size_t>(g)] =
        build_static_lca_index(genes[static_cast<std::size_t>(g)], taxon_count);
  }, 1);
  st.dyn.reset(0);
  st.dyn.reserve(st.max_nodes);
  st.child_index =
      std::make_unique<ChildRangeForest>(st.max_nodes, st.k, st.static_lca);

  st.root = st.new_node();
  st.leaf_count[static_cast<std::size_t>(st.root)] =
      static_cast<uint32_t>(taxon_count);
  for (int g = 0; g < st.k; ++g) {
    st.lca_in_gene(st.root, g) =
        genes[static_cast<std::size_t>(g)].reference_root;
  }

  for (int taxon = 0; taxon < taxon_count; ++taxon) {
    const int leaf = st.new_node();
    st.taxon[static_cast<std::size_t>(leaf)] = taxon;
    st.leaf_count[static_cast<std::size_t>(leaf)] = 1;
    st.taxon_to_leaf[static_cast<std::size_t>(taxon)] = leaf;
    st.parent[static_cast<std::size_t>(leaf)] = st.root;
    st.topology.attach(st.root, leaf);
    ++st.degree[static_cast<std::size_t>(st.root)];
    st.dyn.link(st.root, leaf);
  }
  // Builds each gene's child index of the initial star in parallel. A gene
  // writes only its own LCA entries and pool segment, so the build is
  // race-free and identical to a serial one.
  parallel_for_inner(0, st.k, [&](int g) {
    const StaticLcaIndex& lca = st.static_lca[static_cast<std::size_t>(g)];
    std::vector<ChildRangeForest::InitialChild> initial_children;
    initial_children.reserve(static_cast<std::size_t>(taxon_count));
    for (int taxon = 0; taxon < taxon_count; ++taxon) {
      const int leaf = st.taxon_to_leaf[static_cast<std::size_t>(taxon)];
      const int lca_node = lca.node_for_taxon(taxon);
      st.lca_in_gene(leaf, g) = lca_node;
      const int tin =
          lca_node >= 0 ? lca.tin[static_cast<std::size_t>(lca_node)] : -1;
      initial_children.push_back(
          ChildRangeForest::InitialChild{leaf, tin, lca_node >= 0 ? 1u : 0u,
                                         lca_node});
    }
    st.child_index->build_children(st.root, g, initial_children);
  }, 1);
  return st;
}

// True when node has taxa in gene and they all lie in rep_node's subtree,
// i.e. node's cluster, restricted to gene, is inside the candidate.
bool subset_node_in_candidate(const MreState& st,
                              int node,
                              int gene,
                              int rep_node) {
  const int lca = st.lca_in_gene(node, gene);
  return st.static_lca[static_cast<std::size_t>(gene)].in_subtree(lca,
                                                                  rep_node);
}

// The deepest node on the root-to-leaf path that is not inside the
// candidate. Containment is monotone along the path and the root is never
// inside, so a binary search takes O(log n) path queries.
int find_u0(MreState& st, int leaf, int gene, int rep_node) {
  const int depth = st.dyn.path_size(st.root, leaf) - 1;
  int lo = 0;
  int hi = depth;
  while (lo < hi) {
    const int mid = (lo + hi + 1) >> 1;
    const int u = st.dyn.kth_on_path(st.root, leaf, mid);
    if (subset_node_in_candidate(st, u, gene, rep_node)) {
      hi = mid - 1;
    } else {
      lo = mid;
    }
  }
  return st.dyn.kth_on_path(st.root, leaf, lo);
}

// The number of v's taxa present in gene g: 0 or 1 for a leaf, and for an
// internal node the sum its child index holds for g.
uint32_t present_in_gene(const MreState& st, int v, int g) {
  if (st.taxon[static_cast<std::size_t>(v)] >= 0) {
    return st.lca_in_gene(v, g) >= 0 ? 1u : 0u;
  }
  return st.child_index->aggregate_leaf_sum(v, g);
}

void detach_topology(MreState& st, int parent, int child) {
  st.topology.detach(parent, child);
  st.parent[static_cast<std::size_t>(child)] = -1;
  --st.degree[static_cast<std::size_t>(parent)];
  st.dyn.cut(parent, child);
}

void attach_topology(MreState& st, int parent, int child) {
  st.parent[static_cast<std::size_t>(child)] = parent;
  st.topology.attach(parent, child);
  ++st.degree[static_cast<std::size_t>(parent)];
  st.dyn.link(parent, child);
}

int child_key_tin(const MreState& st, int child, int gene) {
  const int lca = st.lca_in_gene(child, gene);
  if (lca < 0) return -1;
  return st.static_lca[static_cast<std::size_t>(gene)]
      .tin[static_cast<std::size_t>(lca)];
}

void erase_child_index(MreState& st, int parent, int gene, int child) {
  st.child_index->erase_child(parent, gene, child,
                              child_key_tin(st, child, gene));
}

void insert_child_index(MreState& st, int parent, int gene, int child) {
  st.child_index->insert_child(
      parent, gene, child, child_key_tin(st, child, gene),
      present_in_gene(st, child, gene), st.lca_in_gene(child, gene));
}

// The per-gene index work of one accepted split. Each gene owns a disjoint
// AVL pool (see ChildRangeForest's per-gene segments) and LCA slots, and an
// update reads only that gene's state and fields of the nodes it names that
// never change afterwards, so genes can be updated in any interleaving as long
// as each gene sees the updates in acceptance order. AVL updates and bulk
// builds are deterministic, so the result equals a serial run.
struct MreGeneUpdate {
  // A complement move; otherwise the contained children move to z.
  bool complement = false;
  // z's index is bulk-built from its children, from eight moved children on.
  bool batch = false;
  int u0 = -1;
  int z = -1;
  int old_parent = -1;
  std::vector<int> moved;
};

// Builds a new parent's child index for one gene from `children`, its
// topology children when the split was accepted.
void build_parent_gene_index(MreState& st,
                             int parent,
                             int g,
                             const std::vector<int>& children,
                             int extra_child) {
  std::vector<ChildRangeForest::InitialChild> initial;
  initial.reserve(children.size() + (extra_child >= 0 ? 1 : 0));
  const auto add = [&](int child) {
    initial.push_back(ChildRangeForest::InitialChild{
        child,
        child_key_tin(st, child, g),
        present_in_gene(st, child, g),
        st.lca_in_gene(child, g)});
  };
  for (int child : children) add(child);
  if (extra_child >= 0) add(extra_child);
  st.child_index->build_children(parent, g, initial);
}

void apply_gene_update(MreState& st, const MreGeneUpdate& update, int g) {
  if (!update.complement) {
    for (int child : update.moved) erase_child_index(st, update.u0, g, child);
    if (update.batch) {
      build_parent_gene_index(st, update.z, g, update.moved, -1);
    } else {
      for (int child : update.moved) insert_child_index(st, update.z, g, child);
    }
    st.lca_in_gene(update.z, g) = st.child_index->aggregate_lca(update.z, g);
    insert_child_index(st, update.u0, g, update.z);
    return;
  }
  // z inherits u0's pre-move LCA (same taxa). u0 leaves its parent's index
  // under that pre-move key, then its LCA is recomputed from the children it
  // keeps, and z's index is built before z enters the parent's, since z's
  // present-taxa count is read from it.
  const int z_lca = st.lca_in_gene(update.u0, g);
  if (update.old_parent >= 0) {
    erase_child_index(st, update.old_parent, g, update.u0);
  }
  st.lca_in_gene(update.z, g) = z_lca;
  for (int child : update.moved) erase_child_index(st, update.u0, g, child);
  st.lca_in_gene(update.u0, g) = st.child_index->aggregate_lca(update.u0, g);
  if (update.batch) {
    build_parent_gene_index(st, update.z, g, update.moved, update.u0);
  } else {
    for (int child : update.moved) insert_child_index(st, update.z, g, child);
    insert_child_index(st, update.z, g, update.u0);
  }
  if (update.old_parent >= 0) {
    insert_child_index(st, update.old_parent, g, update.z);
  }
}

// Applies the gene updates of accepted splits on helper workers while the
// assembly tests the next candidates. The genes are split into contiguous
// groups, each with a lock and the number of updates applied to it. A
// candidate reads only its representative gene, so the assembly holds that
// gene's group, applying its pending updates itself if no helper has, and
// never waits for a helper to start. Helpers take any group with pending
// updates, and are joined before the stream is destroyed.
//
// An update is freed once every group has applied it. When the updates not
// yet freed hold more than kBacklogPerTaxon * n moved children, the assembly
// applies every group's pending updates itself, so the stream holds O(n)
// words however many moves the run makes.
class MreGeneUpdateStream {
 public:
  MreGeneUpdateStream(MreState& st, int max_updates)
      : st_(st), control_(std::make_shared<Control>()) {
    updates_.resize(static_cast<std::size_t>(std::max(0, max_updates)));
    backlog_limit_ = kBacklogPerTaxon * static_cast<std::size_t>(
                                            std::max(1, st.n));
    const int workers = current_inner_parallel_worker_limit();
    const int group_count =
        std::max(1, std::min(st.k, kGroupsPerWorker * workers));
    groups_ = std::make_unique<Group[]>(static_cast<std::size_t>(group_count));
    group_count_ = group_count;
    group_of_.resize(static_cast<std::size_t>(st.k));
    for (int j = 0; j < group_count; ++j) {
      Group& group = groups_[static_cast<std::size_t>(j)];
      group.first = static_cast<int>(static_cast<int64_t>(st.k) * j / group_count);
      group.last =
          static_cast<int>(static_cast<int64_t>(st.k) * (j + 1) / group_count);
      for (int g = group.first; g < group.last; ++g) {
        group_of_[static_cast<std::size_t>(g)] = j;
      }
    }
    control_->stream = this;
    const int helpers = std::max(0, std::min(workers - 1, group_count));
    join_ = std::make_shared<ParallelLoopJoinState>(helpers);
    auto& exec = get_executor();
    for (int h = 0; h < helpers; ++h) {
      try {
        exec.silent_async(
            [control = control_, join = join_, start = h]() noexcept {
              ParallelLoopCompletion completion(*join);
              run_helper(*control, start);
            });
      } catch (...) {
        // Helpers that could not be submitted leave their work to the others
        // and to the assembly.
        join_->complete(helpers - h);
        break;
      }
    }
  }

  MreGeneUpdateStream(const MreGeneUpdateStream&) = delete;
  MreGeneUpdateStream& operator=(const MreGeneUpdateStream&) = delete;

  ~MreGeneUpdateStream() { stop_helpers(); }

  // Locks gene g's group with every update issued so far applied to it.
  std::unique_lock<std::mutex> hold(int gene) {
    Group& group = groups_[static_cast<std::size_t>(
        group_of_[static_cast<std::size_t>(gene)])];
    std::unique_lock<std::mutex> lock = lock_group(group);
    apply_pending(group, issued_);
    return lock;
  }

  void issue(MreGeneUpdate update) {
    if (issued_ >= static_cast<long long>(updates_.size())) {
      throw std::logic_error("MRE accepted more splits than its bound");
    }
    backlog_ += update.moved.size() + 1;
    updates_[static_cast<std::size_t>(issued_)] =
        std::make_unique<MreGeneUpdate>(std::move(update));
    ++issued_;
    control_->issued.store(issued_, std::memory_order_release);
    wake_helpers();
    release_applied();
    if (backlog_ > backlog_limit_) {
      for (int j = 0; j < group_count_; ++j) {
        Group& group = groups_[static_cast<std::size_t>(j)];
        const std::unique_lock<std::mutex> lock = lock_group(group);
        apply_pending(group, issued_);
      }
      release_applied();
    }
  }

  // Applies every issued update and stops the helpers.
  void finish() {
    stop_helpers();
    rethrow_helper_failure();
    for (int j = 0; j < group_count_; ++j) {
      Group& group = groups_[static_cast<std::size_t>(j)];
      std::lock_guard<std::mutex> lock(group.mutex);
      apply_pending(group, issued_);
    }
    release_applied();
  }

 private:
  static constexpr int kGroupsPerWorker = 4;
  // Moved children per taxon the unapplied updates may hold.
  static constexpr std::size_t kBacklogPerTaxon = 2;
  // Pause iterations before an idle helper sleeps.
  static constexpr int kHelperSpins = 256;

  struct alignas(64) Group {
    std::mutex mutex;
    std::atomic<long long> applied{0};
    int first = 0;
    int last = 0;
  };

  // Shared with the helpers.
  struct Control {
    MreGeneUpdateStream* stream = nullptr;
    std::atomic<long long> issued{0};
    std::atomic<std::uint32_t> wake{0};
    std::atomic<int> sleepers{0};
    std::atomic<bool> closing{false};
    std::mutex failure_mutex;
    std::exception_ptr failure;
    std::atomic<bool> failed{false};
  };

  // Applies updates [applied, upto) to the genes of `group`, whose lock the
  // caller holds.
  void apply_pending(Group& group, long long upto) {
    const long long from = group.applied.load(std::memory_order_relaxed);
    if (from >= upto) return;
    for (int g = group.first; g < group.last; ++g) {
      for (long long i = from; i < upto; ++i) {
        apply_gene_update(st_, *updates_[static_cast<std::size_t>(i)], g);
      }
    }
    // Release: the reads of these updates happen before release_applied
    // frees them.
    group.applied.store(upto, std::memory_order_release);
  }

  // Locks `group` for the assembly thread. A helper that failed while applying
  // a group's updates records the failure before it unlocks the group, so
  // checking after the lock never retries a partly applied batch.
  std::unique_lock<std::mutex> lock_group(Group& group) {
    std::unique_lock<std::mutex> lock(group.mutex);
    rethrow_helper_failure();
    return lock;
  }

  // Frees the updates every group has applied. Only the assembly thread
  // calls it; a group never reads an update below its applied count again.
  void release_applied() {
    long long applied = issued_;
    for (int j = 0; j < group_count_; ++j) {
      applied = std::min(applied, groups_[static_cast<std::size_t>(j)]
                                      .applied.load(std::memory_order_acquire));
    }
    for (; released_ < applied; ++released_) {
      std::unique_ptr<MreGeneUpdate>& update =
          updates_[static_cast<std::size_t>(released_)];
      backlog_ -= update->moved.size() + 1;
      update.reset();
    }
  }

  static void run_helper(Control& control, int start) {
    MreGeneUpdateStream& stream = *control.stream;
    while (true) {
      // The wake count is read before `closing`: stop_helpers sets closing
      // before it bumps the count, so a helper that still sees closing unset
      // sleeps on a count the bump changes and cannot miss the stop.
      const std::uint32_t wake = control.wake.load(std::memory_order_acquire);
      if (control.closing.load(std::memory_order_acquire) ||
          control.failed.load(std::memory_order_relaxed)) {
        return;
      }
      const long long issued = control.issued.load(std::memory_order_acquire);
      bool worked = false;
      for (int step = 0; step < stream.group_count_; ++step) {
        Group& group = stream.groups_[static_cast<std::size_t>(
            (start + step) % stream.group_count_)];
        if (group.applied.load(std::memory_order_relaxed) >= issued) continue;
        std::unique_lock<std::mutex> lock(group.mutex, std::try_to_lock);
        if (!lock.owns_lock()) continue;
        // A helper that failed in this group recorded the failure before it
        // unlocked the group; its batch may be partly applied, so stop.
        if (control.failed.load(std::memory_order_relaxed)) return;
        try {
          stream.apply_pending(group, issued);
        } catch (...) {
          std::lock_guard<std::mutex> failure_lock(control.failure_mutex);
          if (!control.failure) control.failure = std::current_exception();
          control.failed.store(true, std::memory_order_relaxed);
        }
        worked = true;
      }
      if (worked) continue;
      bool woken = false;
      for (int spin = 0; spin < kHelperSpins && !woken; ++spin) {
        cpu_relax();
        woken = control.wake.load(std::memory_order_acquire) != wake;
      }
      if (woken) continue;
      control.sleepers.fetch_add(1, std::memory_order_seq_cst);
      control.wake.wait(wake, std::memory_order_acquire);
      control.sleepers.fetch_sub(1, std::memory_order_relaxed);
    }
  }

  static void cpu_relax() {
#if defined(__x86_64__) || defined(__i386__)
    __builtin_ia32_pause();
#endif
  }

  void wake_helpers() {
    control_->wake.fetch_add(1, std::memory_order_seq_cst);
    if (control_->sleepers.load(std::memory_order_seq_cst) > 0) {
      control_->wake.notify_all();
    }
  }

  // Joins every submitted helper; one that starts after closing returns at
  // once. Idempotent.
  void stop_helpers() {
    if (join_ == nullptr) return;
    control_->closing.store(true, std::memory_order_seq_cst);
    wake_helpers();
    wait_for_parallel_loop(get_executor(), join_);
    join_.reset();
  }

  void rethrow_helper_failure() {
    if (!control_->failed.load(std::memory_order_relaxed)) return;
    std::lock_guard<std::mutex> lock(control_->failure_mutex);
    std::rethrow_exception(control_->failure);
  }

  MreState& st_;
  std::shared_ptr<Control> control_;
  std::shared_ptr<ParallelLoopJoinState> join_;
  std::vector<std::unique_ptr<MreGeneUpdate>> updates_;
  long long issued_ = 0;
  // Updates [0, released_) are freed; backlog_ counts the moved children
  // (plus one per update) of the others.
  long long released_ = 0;
  std::size_t backlog_ = 0;
  std::size_t backlog_limit_ = 0;
  std::unique_ptr<Group[]> groups_;
  int group_count_ = 1;
  std::vector<int> group_of_;
};

// Accepts the split as a new child z of u0 that adopts the contained
// children, updating the topology and returning the genes' update.
MreGeneUpdate accept_move_contained(MreState& st,
                                    int u0,
                                    std::vector<int> contained) {
  const int z = st.new_node();
  uint32_t taxa = 0;
  for (int child : contained) {
    taxa += st.leaf_count[static_cast<std::size_t>(child)];
  }
  st.leaf_count[static_cast<std::size_t>(z)] = taxa;

  // The link-cut tree and the children, parent and degree arrays are shared
  // by all genes.
  for (int child : contained) {
    detach_topology(st, u0, child);
    attach_topology(st, z, child);
  }
  attach_topology(st, u0, z);

  MreGeneUpdate update;
  update.batch = contained.size() >= 8;
  update.u0 = u0;
  update.z = z;
  update.moved = std::move(contained);
  return update;
}

// Accepts the split by moving the complement instead: a new node z takes
// u0's place under its parent and adopts u0's other children, and u0, left
// with the contained children, becomes z's child.
MreGeneUpdate accept_move_complement(MreState& st,
                                     int u0,
                                     std::vector<int> complement) {
  const int old_parent = st.parent[static_cast<std::size_t>(u0)];
  const uint32_t old_leaf_count = st.leaf_count[static_cast<std::size_t>(u0)];

  const int z = st.new_node();
  st.leaf_count[static_cast<std::size_t>(z)] = old_leaf_count;

  if (old_parent >= 0) {
    detach_topology(st, old_parent, u0);
    attach_topology(st, old_parent, z);
  } else {
    st.root = z;
  }
  uint32_t moved_taxa = 0;
  for (int child : complement) {
    moved_taxa += st.leaf_count[static_cast<std::size_t>(child)];
    detach_topology(st, u0, child);
    attach_topology(st, z, child);
  }
  st.leaf_count[static_cast<std::size_t>(u0)] = old_leaf_count - moved_taxa;
  attach_topology(st, z, u0);

  MreGeneUpdate update;
  update.complement = true;
  update.batch = complement.size() >= 8;
  update.u0 = u0;
  update.z = z;
  update.old_parent = old_parent;
  update.moved = std::move(complement);
  return update;
}

// Fills the preorder and, for each node, its preorder interval [tin, tout).
void fill_traversal(StartingTree& starting_tree) {
  starting_tree.preorder.clear();
  starting_tree.tin.assign(starting_tree.nodes.size(), 0);
  starting_tree.tout.assign(starting_tree.nodes.size(), 0);
  if (starting_tree.root < 0) return;
  std::vector<std::pair<int, int>> stack{{starting_tree.root, 0}};
  while (!stack.empty()) {
    const int node = stack.back().first;
    int& next = stack.back().second;
    if (next == 0) {
      starting_tree.tin[static_cast<std::size_t>(node)] =
          static_cast<int>(starting_tree.preorder.size());
      starting_tree.preorder.push_back(node);
    }
    if (next < static_cast<int>(
                   starting_tree.nodes[static_cast<std::size_t>(node)].children.size())) {
      stack.push_back({starting_tree.nodes[static_cast<std::size_t>(node)]
                           .children[static_cast<std::size_t>(next++)],
                       0});
      continue;
    }
    starting_tree.tout[static_cast<std::size_t>(node)] =
        static_cast<int>(starting_tree.preorder.size());
    stack.pop_back();
  }
}

// Children keep the topology's slot order.
StartingTree materialize_mre_tree(const MreState& st) {
  StartingTree out;
  out.root = st.root;
  out.nodes.resize(static_cast<std::size_t>(st.node_count));
  out.taxon_to_leaf.assign(static_cast<std::size_t>(st.n), -1);
  for (int node = 0; node < st.node_count; ++node) {
    out.nodes[static_cast<std::size_t>(node)].parent =
        st.parent[static_cast<std::size_t>(node)];
    out.nodes[static_cast<std::size_t>(node)].taxon =
        st.taxon[static_cast<std::size_t>(node)];
    out.nodes[static_cast<std::size_t>(node)].children =
        st.topology.children(node);
    if (st.taxon[static_cast<std::size_t>(node)] >= 0) {
      out.taxon_to_leaf[static_cast<std::size_t>(
          st.taxon[static_cast<std::size_t>(node)])] = node;
    }
  }
  fill_traversal(out);
  return out;
}

}  // namespace

// Test probe on synthetic records and indices. Checks the radix aggregation
// and priority sort against std::sort oracles, the child index against
// scans under adversarial insert and erase orders, and FlatChildTopology
// against a swap-with-last vector, and reports the pass counts and height
// bounds that the complexity tests check.
MreComplexityProbeResult mre_complexity_probe_for_tests(int item_count) {
  if (item_count < 2) {
    throw std::invalid_argument(
        "MRE complexity probe requires at least two items");
  }

  std::vector<SplitOccurrence256> occurrences;
  occurrences.reserve(static_cast<std::size_t>(item_count) * 2u);
  for (int i = 0; i < item_count; ++i) {
    const uint64_t index = static_cast<uint64_t>(i);
    SplitOccurrence256 occurrence;
    // Keys have the exact shape {size, {id, 0, 0, 0}}: few sizes, so ids
    // order most runs, and ids spread over both 16-bit digits. Multiplying
    // by an odd constant modulo 2^32 is a bijection, so every item has its
    // own key.
    occurrence.key.size = static_cast<uint32_t>(
        2 + stable_hash_u64(index % 23, index / 23, 0x082efa98ec4e6c89ULL) %
                61);
    occurrence.key.h[0] =
        static_cast<uint32_t>(index * uint64_t{0x9e3779b1});
    occurrence.gene = static_cast<uint32_t>((37LL * i + 11) % item_count);
    occurrence.edge = static_cast<uint32_t>(
        stable_hash_u64(index % 97, index / 97, 0xa4093822299f31d0ULL));
    occurrence.node = static_cast<uint32_t>(i);
    occurrence.one_taxon = static_cast<uint32_t>(item_count - i);
    occurrences.push_back(occurrence);
    if (i % 3 == 0) {
      SplitOccurrence256 duplicate = occurrence;
      duplicate.edge += 1;
      duplicate.node ^= 0x5a5a5a5aU;
      occurrences.push_back(duplicate);
    }
    if (i % 5 == 0) {
      SplitOccurrence256 second_gene = occurrence;
      second_gene.gene = static_cast<uint32_t>(
          (static_cast<uint64_t>(occurrence.gene) + 1u) %
          static_cast<uint64_t>(item_count));
      second_gene.edge ^= 0x80000000U;
      occurrences.push_back(second_gene);
    }
  }
  std::reverse(occurrences.begin(), occurrences.end());
  std::rotate(
      occurrences.begin(),
      occurrences.begin() +
          static_cast<std::ptrdiff_t>(occurrences.size() / 3u),
      occurrences.end());

  const auto occurrence_less_oracle = [](const SplitOccurrence256& a,
                                         const SplitOccurrence256& b) {
    if (a.key.size != b.key.size) return a.key.size < b.key.size;
    if (a.key.h != b.key.h) return a.key.h < b.key.h;
    if (a.gene != b.gene) return a.gene < b.gene;
    return a.edge < b.edge;
  };
  std::vector<SplitOccurrence256> ordered_occurrence_oracle = occurrences;
  std::sort(ordered_occurrence_oracle.begin(), ordered_occurrence_oracle.end(),
            occurrence_less_oracle);

  std::vector<MreEntry> aggregate_oracle;
  aggregate_oracle.reserve(static_cast<std::size_t>(item_count));
  std::size_t oracle_begin = 0;
  while (oracle_begin < ordered_occurrence_oracle.size()) {
    const SplitOccurrence256& first = ordered_occurrence_oracle[oracle_begin];
    MreEntry aggregate;
    aggregate.key = first.key;
    aggregate.split_id = static_cast<uint32_t>(aggregate_oracle.size());
    aggregate.first_tree = first.gene;
    aggregate.first_edge = first.edge;
    aggregate.rep_tree = first.gene;
    aggregate.rep_node = first.node;
    aggregate.one_taxon = first.one_taxon;
    uint32_t previous_gene = UINT32_MAX;
    std::size_t oracle_end = oracle_begin;
    while (oracle_end < ordered_occurrence_oracle.size() &&
           ordered_occurrence_oracle[oracle_end].key == first.key) {
      if (ordered_occurrence_oracle[oracle_end].gene != previous_gene) {
        previous_gene = ordered_occurrence_oracle[oracle_end].gene;
        ++aggregate.count;
      }
      ++oracle_end;
    }
    aggregate_oracle.push_back(aggregate);
    oracle_begin = oracle_end;
  }

  SplitRadixStats radix_stats;
  std::vector<MreEntry> entries =
      aggregate_split_occurrences256(occurrences, radix_stats);
  const auto occurrence_equal = [](const SplitOccurrence256& a,
                                   const SplitOccurrence256& b) {
    return a.key == b.key && a.gene == b.gene && a.edge == b.edge &&
           a.node == b.node && a.one_taxon == b.one_taxon;
  };
  const auto entry_equal = [](const MreEntry& a, const MreEntry& b) {
    return a.key == b.key && a.split_id == b.split_id &&
           a.count == b.count && a.first_tree == b.first_tree &&
           a.first_edge == b.first_edge && a.rep_tree == b.rep_tree &&
           a.rep_node == b.rep_node && a.one_taxon == b.one_taxon;
  };
  const bool occurrence_order_matches = std::equal(
      occurrences.begin(), occurrences.end(),
      ordered_occurrence_oracle.begin(), occurrence_equal);
  const bool aggregation_matches =
      entries.size() == aggregate_oracle.size() &&
      std::equal(entries.begin(), entries.end(), aggregate_oracle.begin(),
                 entry_equal);
  bool dense_ids_match = true;
  for (std::size_t i = 0; i < entries.size(); ++i) {
    dense_ids_match = dense_ids_match && entries[i].split_id == i;
  }

  std::vector<MreEntry> priority_oracle = aggregate_oracle;
  std::sort(priority_oracle.begin(), priority_oracle.end(), entry_less);
  radix_sort_split_priority256(entries, radix_stats);
  const bool priority_matches =
      entries.size() == priority_oracle.size() &&
      std::equal(entries.begin(), entries.end(), priority_oracle.begin(),
                 entry_equal);
  const bool radix_matches = occurrence_order_matches && aggregation_matches &&
                             dense_ids_match && priority_matches;

  StaticLcaIndex lca;
  lca.tin.resize(static_cast<std::size_t>(item_count + 1));
  for (int child = 0; child <= item_count; ++child) {
    lca.tin[static_cast<std::size_t>(child)] = child;
  }
  std::vector<StaticLcaIndex> lca_by_gene;
  lca_by_gene.push_back(std::move(lca));
  ChildRangeForest forest(item_count + 1, 1, lca_by_gene);
  std::vector<uint8_t> live(static_cast<std::size_t>(item_count + 1), 0);
  bool avl_valid = true;
  bool aggregates_valid = true;
  long long empty_complement_visits = 0;

  const auto leaf_mass = [](int child) {
    return static_cast<uint32_t>(1 + (child % 11));
  };
  const auto validate = [&]() {
    avl_valid = avl_valid && forest.validate_root_for_tests(0, 0);
  };
  const auto check_range = [&](int lo, int hi) {
    ChildAgg expected;
    std::vector<int> expected_inside;
    std::vector<int> expected_outside;
    for (int child = 1; child <= item_count; ++child) {
      if (live[static_cast<std::size_t>(child)] == 0) continue;
      if (child >= lo && child < hi) {
        ++expected.count;
        expected.leaf_sum += leaf_mass(child);
        if (expected.min_node < 0) expected.min_node = child;
        expected.max_node = child;
        expected_inside.push_back(child);
      } else {
        expected_outside.push_back(child);
      }
    }
    const ChildAgg observed = forest.range_query(0, 0, lo, hi);
    const std::vector<int> inside = forest.enumerate_range(0, 0, lo, hi);
    const long long before_complement = forest.enumeration_visit_count();
    const std::vector<int> outside =
        forest.enumerate_complement(0, 0, lo, hi);
    const long long complement_visits =
        forest.enumeration_visit_count() - before_complement;
    if (expected.count > 0 && expected_outside.empty()) {
      empty_complement_visits =
          std::max(empty_complement_visits, complement_visits);
    }
    aggregates_valid =
        aggregates_valid && observed.count == expected.count &&
        observed.leaf_sum == expected.leaf_sum &&
        observed.min_node == expected.min_node &&
        observed.max_node == expected.max_node && inside == expected_inside &&
        outside == expected_outside;
  };

  // Monotone insertion is adversarial for an unbalanced BST and for a treap
  // whose fixed priorities happen to be monotone in key order.
  for (int child = 1; child <= item_count; ++child) {
    forest.insert_child(0, 0, child, child, leaf_mass(child), child);
    live[static_cast<std::size_t>(child)] = 1;
  }
  validate();
  check_range(1, item_count + 1);
  check_range(item_count / 4, 3 * item_count / 4);
  check_range(item_count + 1, item_count + 2);

  for (int child = 1; child <= item_count; child += 2) {
    forest.erase_child(0, 0, child, child);
    live[static_cast<std::size_t>(child)] = 0;
  }
  validate();
  check_range(item_count / 5, 4 * item_count / 5);
  for (int child = item_count - (item_count & 1); child >= 2; child -= 2) {
    forest.erase_child(0, 0, child, child);
    live[static_cast<std::size_t>(child)] = 0;
  }
  validate();
  check_range(1, item_count + 1);

  for (int child = item_count; child >= 1; --child) {
    forest.insert_child(0, 0, child, child, leaf_mass(child), child);
    live[static_cast<std::size_t>(child)] = 1;
  }
  validate();
  for (int child = 3; child <= item_count; child += 3) {
    forest.erase_child(0, 0, child, child);
    live[static_cast<std::size_t>(child)] = 0;
  }
  validate();
  check_range(item_count / 7, 6 * item_count / 7);
  for (int child = 3; child <= item_count; child += 3) {
    forest.insert_child(0, 0, child, child, leaf_mass(child), child);
    live[static_cast<std::size_t>(child)] = 1;
  }
  validate();
  check_range(1, item_count + 1);

  FlatChildTopology topology(item_count + 2);
  std::vector<int> root_children;
  std::vector<int> nested_children;
  root_children.reserve(static_cast<std::size_t>(item_count + 1));
  nested_children.reserve(static_cast<std::size_t>(item_count));
  for (int child = 1; child <= item_count + 1; ++child) {
    topology.attach(0, child);
    root_children.push_back(child);
  }
  const auto oracle_detach = [](std::vector<int>& children, int child) {
    const auto found = std::find(children.begin(), children.end(), child);
    if (found == children.end()) return false;
    *found = children.back();
    children.pop_back();
    return true;
  };
  for (int child = 2; child <= item_count; child += 3) {
    topology.detach(0, child);
    if (!oracle_detach(root_children, child)) {
      throw std::logic_error("MRE topology probe oracle lost a child");
    }
    topology.attach(1, child);
    nested_children.push_back(child);
  }
  const std::vector<int> nested_snapshot = nested_children;
  for (std::size_t i = 1; i < nested_snapshot.size(); i += 5) {
    const int child = nested_snapshot[i];
    topology.detach(1, child);
    if (!oracle_detach(nested_children, child)) {
      throw std::logic_error("MRE topology probe oracle lost a move");
    }
    topology.attach(0, child);
    root_children.push_back(child);
  }
  const bool topology_matches = topology.validate_for_tests() &&
                                topology.children(0) == root_children &&
                                topology.children(1) == nested_children;

  MreComplexityProbeResult result;
  result.radix_matches_oracle = radix_matches;
  result.avl_invariants_hold =
      avl_valid && forest.height_violation_count() == 0;
  result.range_aggregates_match = aggregates_valid;
  result.radix_passes = static_cast<int>(radix_stats.total_passes());
  result.radix_items = radix_stats.total_items();
  result.occurrence_radix_passes =
      static_cast<int>(radix_stats.occurrence_passes);
  result.priority_radix_passes =
      static_cast<int>(radix_stats.priority_passes);
  result.occurrence_radix_items = radix_stats.occurrence_items;
  result.priority_radix_items = radix_stats.priority_items;
  result.emitted_occurrences = radix_stats.emitted_occurrences;
  result.deduplicated_gene_splits = radix_stats.deduplicated_gene_splits;
  result.distinct_splits = radix_stats.distinct_splits;
  result.dense_split_ids = dense_ids_match;
  result.max_avl_height = forest.maximum_height();
  result.avl_height_limit =
      ChildRangeForest::height_limit_for_count(item_count);
  result.empty_complement_visits = empty_complement_visits;
  result.empty_complement_visit_limit =
      2LL * result.avl_height_limit + 2;
  result.flat_topology_matches = topology_matches;
  result.topology_peak_slots = topology.peak_live_slots();
  result.topology_slot_limit = topology.slot_limit();
  return result;
}

// Builds the MRE starting tree; see the overview at the top of this file.
StartingTree build_mre_tree(const std::vector<GeneIndex>& genes,
                            int taxon_count,
                            const Options& options,
                            const GuideFilter* guide) {
  validate_mre_dimensions(genes, taxon_count);
  progress_mre_start(options, static_cast<int>(genes.size()),
                           taxon_count);
  const auto count_started = std::chrono::steady_clock::now();
  const std::vector<MreEntry> entries =
      count_mre_entries(genes, taxon_count);
  RunStats& stats = run_stats();
  stats_add_time(stats.mre_count_ns,
                      std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - count_started));
  stats.mre_distinct_splits.store(static_cast<long long>(entries.size()),
                                        std::memory_order_relaxed);
  // The count phase has just freed its occurrence records and radix scratch.
  release_retained_allocator_pages();

  const auto assembly_started = std::chrono::steady_clock::now();
  MreState st = initialize_state(genes, taxon_count);
  int accepted = 0;
  long long candidates_checked = 0;
  // A tree on n taxa has at most n - 3 nontrivial splits.
  const int max_accepted_splits = std::max(0, taxon_count - 3);
  progress_mre_assembly(0, static_cast<long long>(entries.size()), 0,
                              max_accepted_splits);
  MreGeneUpdateStream updates(st, max_accepted_splits);
  for (const MreEntry& entry : entries) {
    if (accepted >= max_accepted_splits) break;
    ++candidates_checked;
    const auto note_assembly_progress = [&]() {
      if ((candidates_checked & 1023ll) == 0 ||
          accepted >= max_accepted_splits) {
        progress_mre_assembly(
            candidates_checked, static_cast<long long>(entries.size()),
            accepted, max_accepted_splits);
      }
    };
    const int gene = static_cast<int>(entry.rep_tree);
    const int rep_node = static_cast<int>(entry.rep_node);
    const int taxon = static_cast<int>(entry.one_taxon);
    if (taxon < 0 || taxon >= taxon_count) {
      note_assembly_progress();
      continue;
    }
    const int leaf = st.taxon_to_leaf[static_cast<std::size_t>(taxon)];
    if (leaf < 0) {
      note_assembly_progress();
      continue;
    }
    // Only the representative gene's index is read, and the accept itself
    // changes only the shared topology; the genes' updates are streamed.
    std::unique_lock<std::mutex> gene_lock = updates.hold(gene);
    const int u0 = find_u0(st, leaf, gene, rep_node);
    const StaticLcaIndex& lca = st.static_lca[static_cast<std::size_t>(gene)];
    const int lo = lca.tin[static_cast<std::size_t>(rep_node)];
    const int hi = lca.tout[static_cast<std::size_t>(rep_node)];
    // The split's taxa are those of the gene below rep_node. It is compatible
    // with the tree exactly when the children of u0 whose gene LCA lies in
    // rep_node's subtree hold all of them (key.size taxa present in the gene)
    // and not all of u0's. Taxa the gene lacks are not counted, so they cannot
    // stand in for a split taxon outside u0. With a single such child the
    // split is already in the tree, restricted to the gene's taxa, and
    // accepting it would add a node with one child.
    const ChildAgg contained = st.child_index->range_query(u0, gene, lo, hi);
    if (contained.count < 2 ||
        contained.leaf_sum != static_cast<uint64_t>(entry.key.size) ||
        contained.leaf_sum == present_in_gene(st, u0, gene)) {
      note_assembly_progress();
      continue;
    }
    // With a guide tree: the cluster this split would add, the taxa of the
    // contained children, must not conflict with the guide.
    if (guide != nullptr) {
      std::vector<int> taxa;
      std::vector<int> stack = st.child_index->enumerate_range(u0, gene, lo, hi);
      while (!stack.empty()) {
        const int v = stack.back();
        stack.pop_back();
        if (st.taxon[static_cast<std::size_t>(v)] >= 0) {
          taxa.push_back(st.taxon[static_cast<std::size_t>(v)]);
        } else {
          for (int c : st.topology.children(v)) stack.push_back(c);
        }
      }
      if (guide->conflicts(taxa)) {
        stats_add(stats.guide_seed_rejected);
        note_assembly_progress();
        continue;
      }
    }
    // Move whichever child group is smaller.
    MreGeneUpdate update =
        static_cast<uint64_t>(contained.count) * 2u <=
                st.degree[static_cast<std::size_t>(u0)]
            ? accept_move_contained(
                  st, u0, st.child_index->enumerate_range(u0, gene, lo, hi))
            : accept_move_complement(
                  st, u0,
                  st.child_index->enumerate_complement(u0, gene, lo, hi));
    gene_lock.unlock();
    updates.issue(std::move(update));
    ++accepted;
    note_assembly_progress();
  }
  updates.finish();
  progress_mre_finish(accepted, max_accepted_splits);
  stats_add_time(stats.mre_assembly_ns,
                      std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - assembly_started));
  stats.mre_accepted_splits.store(accepted, std::memory_order_relaxed);
  const long long height_violations =
      st.child_index->height_violation_count();
  if (st.topology.live_slots() != st.node_count - 1 ||
      st.topology.allocated_slots() > st.topology.slot_limit() ||
      st.topology.peak_live_slots() > st.topology.slot_limit()) {
    stats.complexity_contract_satisfied.store(0, std::memory_order_relaxed);
    throw std::logic_error("MRE flat topology check failed");
  }
  if (height_violations != 0) {
    stats.complexity_contract_satisfied.store(0, std::memory_order_relaxed);
  }
  return materialize_mre_tree(st);
}

std::vector<SplitAggregate256> mre_split_entries_for_tests(
    const std::vector<GeneIndex>& genes,
    int taxon_count) {
  validate_mre_dimensions(genes, taxon_count);
  return count_mre_entries(genes, taxon_count);
}

}  // namespace qraft::internal
