#include "mre/exact_split_ids.hpp"

#include "internal.hpp"

#include <algorithm>
#include <bit>
#include <chrono>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

// Exact split identifiers for the MRE.
//
// A gene node is closed when it holds all or none of the copies of every
// taxon. Closed nodes with at least two distinct taxa, together with one leaf
// per taxon, form a single-copy forest whose clusters are the nodes' distinct
// taxon sets. Its clusters get exact identifiers by halving the taxon range:
// the identifier of a cluster on [lo, hi) is the rank of the pair of the
// identifiers of its parts in the two halves, which come from the forest
// restricted to each half. This is the halving labelling of Jansson, Sung,
// Tabatabaee and Yang, "A Faster Algorithm for Constructing the Frequency
// Difference Consensus Tree", STACS 2024 (Fast_Label_Trees, Algorithm 3),
// applied to the forest of closed gene clusters; it takes O(N + kn log n)
// time.
//
// Every forest here is a list of segments, one per gene, and no parent link
// crosses a segment. A call on many nodes runs its passes over groups of whole
// segments in parallel and sorts its pairs with a stable parallel radix sort.
// The recursion itself stays serial, so the live memory is that of one chain
// of calls, and every pass produces the arrays of the serial order, whatever
// the number of workers.

namespace qraft::internal {
namespace {

// Calls on fewer nodes run their passes serially.
constexpr std::size_t kParallelNodes = std::size_t{1} << 15;
// Every parallel part covers at least this many nodes.
constexpr std::size_t kNodesPerPart = std::size_t{1} << 13;
// The widest digit of the parallel radix sorts.
constexpr int kRadixBits = 11;

// An allocator whose resize leaves new trivial elements unset, so that large
// arrays are first written, and so first touched, inside the parallel passes.
// Every element is written before it is read.
template <class T>
struct DefaultInitAllocator : std::allocator<T> {
  template <class U>
  struct rebind {
    using other = DefaultInitAllocator<U>;
  };
  DefaultInitAllocator() noexcept = default;
  template <class U>
  DefaultInitAllocator(const DefaultInitAllocator<U>&) noexcept {}
  template <class U>
  void construct(U* p) noexcept(std::is_nothrow_default_constructible_v<U>) {
    ::new (static_cast<void*>(p)) U;
  }
  template <class U, class... Args>
  void construct(U* p, Args&&... args) {
    ::new (static_cast<void*>(p)) U(std::forward<Args>(args)...);
  }
};
template <class T>
using RawVec = std::vector<T, DefaultInitAllocator<T>>;

// The forest of one call, with unset growth.
struct Forest {
  RawVec<int32_t> parent;
  RawVec<int32_t> taxon;
};

// Stable counting sort of `in` by key(i) in [0, buckets) into `out`.
template <class Key>
void counting_sort(const RawVec<uint32_t>& in,
                   RawVec<uint32_t>& out,
                   std::size_t buckets,
                   std::vector<uint32_t>& count,
                   Key&& key) {
  count.assign(buckets + 1, 0);
  for (uint32_t i : in) ++count[static_cast<std::size_t>(key(i)) + 1];
  for (std::size_t b = 1; b <= buckets; ++b) count[b] += count[b - 1];
  out.resize(in.size());
  for (uint32_t i : in) out[count[static_cast<std::size_t>(key(i))]++] = i;
}

// How one call splits its work: `parts` equal position ranges, and chunks of
// whole segments, chunk c covering segments [chunk[c], chunk[c + 1]). A small
// call, or one with a single worker, has one part and one chunk.
struct CallPlan {
  int parts = 1;
  std::vector<uint32_t> chunk;
};

CallPlan plan_call(const std::vector<uint32_t>& seg, std::size_t m) {
  CallPlan plan;
  const std::size_t segments = seg.size() - 1;
  const int workers = current_inner_parallel_worker_limit();
  plan.chunk.push_back(0);
  if (workers > 1 && m >= kParallelNodes) {
    plan.parts = static_cast<int>(std::min<std::size_t>(
        static_cast<std::size_t>(workers), m / kNodesPerPart));
    // One scan: chunk c ends at the first segment end at or past c/parts of
    // the vertices.
    std::size_t s = 0;
    for (int c = 1; c < plan.parts; ++c) {
      const std::size_t target = m * static_cast<std::size_t>(c) /
                                 static_cast<std::size_t>(plan.parts);
      while (s < segments && seg[s] < target) ++s;
      if (s > plan.chunk.back() && s < segments) {
        plan.chunk.push_back(static_cast<uint32_t>(s));
      }
    }
  }
  plan.chunk.push_back(static_cast<uint32_t>(segments));
  return plan;
}

// Runs body(part, begin, end) over `parts` equal ranges of [0, n).
template <class Body>
void for_ranges(std::size_t n, int parts, Body&& body) {
  parallel_for_inner(0, parts, [&](int p) {
    const std::size_t pp = static_cast<std::size_t>(p);
    const std::size_t all = static_cast<std::size_t>(parts);
    body(p, n * pp / all, n * (pp + 1) / all);
  }, 1);
}

// Runs body(first segment, end segment) over the plan's chunks.
template <class Body>
void for_chunks(const CallPlan& plan, Body&& body) {
  const int chunks = static_cast<int>(plan.chunk.size()) - 1;
  parallel_for_inner(0, chunks, [&](int c) {
    body(plan.chunk[static_cast<std::size_t>(c)],
         plan.chunk[static_cast<std::size_t>(c) + 1]);
  }, 1);
}

// Sorts the items 0, 1, ..., m - 1 stably by key[i] < 2^bits, least
// significant digit first, over `parts` equal ranges. Each pass reads the
// (key, item) records in order and scatters them, so no pass gathers. Each
// range counts its digits, and its records follow those of the earlier ranges
// within each digit, so every pass equals the serial stable pass. On return
// `key` holds the sorted keys and the result the items in the same order.
RawVec<uint32_t> sort_records(RawVec<uint64_t>& key, int bits, int parts) {
  const std::size_t m = key.size();
  RawVec<uint32_t> item(m);
  if (bits <= 0) {
    for_ranges(m, parts, [&](int, std::size_t b, std::size_t e) {
      for (std::size_t i = b; i < e; ++i) item[i] = static_cast<uint32_t>(i);
    });
    return item;
  }
  const int passes = (bits + kRadixBits - 1) / kRadixBits;
  const int digit_bits = (bits + passes - 1) / passes;
  const std::size_t buckets = std::size_t{1} << digit_bits;
  const uint64_t digit_mask = buckets - 1;
  RawVec<uint64_t> key_out(m);
  RawVec<uint32_t> item_out(m);
  std::vector<uint32_t> counts(static_cast<std::size_t>(parts) * buckets);
  for (int pass = 0; pass < passes; ++pass) {
    const int shift = pass * digit_bits;
    std::fill(counts.begin(), counts.end(), 0u);
    for_ranges(m, parts, [&](int p, std::size_t b, std::size_t e) {
      uint32_t* const c = counts.data() + static_cast<std::size_t>(p) * buckets;
      const uint64_t* const k = key.data();
      for (std::size_t i = b; i < e; ++i) ++c[(k[i] >> shift) & digit_mask];
    });
    uint32_t running = 0;
    for (std::size_t d = 0; d < buckets; ++d) {
      for (int p = 0; p < parts; ++p) {
        uint32_t& c = counts[static_cast<std::size_t>(p) * buckets + d];
        const uint32_t here = c;
        c = running;
        running += here;
      }
    }
    for_ranges(m, parts, [&](int p, std::size_t b, std::size_t e) {
      uint32_t* const c = counts.data() + static_cast<std::size_t>(p) * buckets;
      const uint64_t* const k = key.data();
      const uint32_t* const it = item.data();
      uint64_t* const ko = key_out.data();
      uint32_t* const io = item_out.data();
      if (pass == 0) {
        // The items start in position order.
        for (std::size_t i = b; i < e; ++i) {
          const uint32_t at = c[(k[i] >> shift) & digit_mask]++;
          ko[at] = k[i];
          io[at] = static_cast<uint32_t>(i);
        }
      } else {
        for (std::size_t i = b; i < e; ++i) {
          const uint32_t at = c[(k[i] >> shift) & digit_mask]++;
          ko[at] = k[i];
          io[at] = it[i];
        }
      }
    });
    key.swap(key_out);
    item.swap(item_out);
  }
  return item;
}

// Numbers 1, 2, ... along sorted keys over `parts` ranges: equal keys share a
// number, and item[j] receives the number of key[j].
RawVec<uint32_t> rank_sorted(const RawVec<uint64_t>& key,
                             const RawVec<uint32_t>& item,
                             int parts) {
  const std::size_t m = key.size();
  RawVec<uint32_t> ids(m);
  std::vector<uint32_t> starts(static_cast<std::size_t>(parts) + 1, 0);
  for_ranges(m, parts, [&](int p, std::size_t b, std::size_t e) {
    uint32_t count = 0;
    for (std::size_t j = b; j < e; ++j) {
      count += j == 0 || key[j] != key[j - 1] ? 1u : 0u;
    }
    starts[static_cast<std::size_t>(p) + 1] = count;
  });
  for (int p = 0; p < parts; ++p) {
    starts[static_cast<std::size_t>(p) + 1] +=
        starts[static_cast<std::size_t>(p)];
  }
  for_ranges(m, parts, [&](int p, std::size_t b, std::size_t e) {
    uint32_t rank = starts[static_cast<std::size_t>(p)];
    const uint64_t* const k = key.data();
    const uint32_t* const it = item.data();
    uint32_t* const out = ids.data();
    for (std::size_t j = b; j < e; ++j) {
      if (j == 0 || k[j] != k[j - 1]) ++rank;
      out[it[j]] = rank;
    }
  });
  return ids;
}

// Serial ranking: the identity order sorted stably by the counting sorts of
// `sort(order, scratch)`, then numbered along it by key.
template <class Sort>
RawVec<uint32_t> rank_serial(const RawVec<uint64_t>& key, Sort&& sort) {
  const std::size_t m = key.size();
  RawVec<uint32_t> order(m);
  for (std::size_t i = 0; i < m; ++i) order[i] = static_cast<uint32_t>(i);
  RawVec<uint32_t> scratch;
  sort(order, scratch);
  RawVec<uint32_t> ids(m);
  uint32_t rank = 0;
  for (std::size_t j = 0; j < m; ++j) {
    const uint32_t i = order[j];
    if (j == 0 || key[i] != key[order[j - 1]]) ++rank;
    ids[i] = rank;
  }
  return ids;
}

// Identifiers of the clusters of `forest`, whose leaves carry taxa in
// [lo, hi), listed in segments `seg` that no parent link crosses. Each call is
// linear in its forest; the forests at one depth of the recursion have
// disjoint taxon ranges.
RawVec<uint32_t> halve(int lo,
                       int hi,
                       Forest&& forest,
                       std::vector<uint32_t>&& seg) {
  const std::size_t m = forest.parent.size();
  if (m == 0) return {};
  if (hi - lo == 1) return RawVec<uint32_t>(m, 1);
  const CallPlan plan = plan_call(seg, m);
  const bool parallel = plan.parts > 1;
  if (hi - lo <= 64) {
    // A mask is the exact set on this range. Lower taxa occupy higher bits,
    // so its numeric order is the same as the recursive pair order. At most
    // eight radix passes take O(m) time and words; replacing these bottom
    // levels preserves O(N + kn log n) time and O(N + kn) peak words.
    RawVec<uint64_t> masks(m);
    for_chunks(plan, [&](uint32_t s0, uint32_t s1) {
      const int top = hi - 1;
      const int32_t* const parent = forest.parent.data();
      const int32_t* const taxa = forest.taxon.data();
      uint64_t* const mask = masks.data();
      const std::size_t first = seg[s0];
      const std::size_t last = seg[s1];
      std::fill(mask + first, mask + last, uint64_t{0});
      for (std::size_t i = first; i < last; ++i) {
        if (taxa[i] >= 0) mask[i] = uint64_t{1} << (top - taxa[i]);
        if (parent[i] >= 0) {
          mask[static_cast<std::size_t>(parent[i])] |= mask[i];
        }
      }
    });
    RawVec<int32_t>().swap(forest.parent);
    RawVec<int32_t>().swap(forest.taxon);
    if (parallel) {
      const RawVec<uint32_t> item = sort_records(masks, hi - lo, plan.parts);
      return rank_sorted(masks, item, plan.parts);
    }
    return rank_serial(masks, [&](RawVec<uint32_t>& order,
                                  RawVec<uint32_t>& scratch) {
      std::vector<uint32_t> count;
      for (int shift = 0; shift < hi - lo; shift += 8) {
        counting_sort(order, scratch, 256, count, [&](uint32_t i) {
          return (masks[i] >> shift) & 255;
        });
        order.swap(scratch);
      }
    });
  }
  const int mid = lo + (hi - lo) / 2;
  const std::size_t segments = seg.size() - 1;

  // below_*[i]: the vertex below i that is the LCA of i's leaves in that half,
  // or -1 when there are none. Children come first, so a vertex's children
  // have reported to it before it is visited. A vertex is kept in a half when
  // it is its own LCA there; the same pass counts each segment's kept
  // vertices into start_*[s + 1].
  RawVec<int32_t> below_a(m);
  RawVec<int32_t> below_b(m);
  std::vector<uint32_t> start_a(segments + 1, 0);
  std::vector<uint32_t> start_b(segments + 1, 0);
  // The passes read bounds and arrays through locals: their int32 stores
  // could otherwise alias captured values and force reloads.
  for_chunks(plan, [&](uint32_t s0, uint32_t s1) {
    const int32_t split = mid;
    const int32_t* const parent = forest.parent.data();
    const int32_t* const taxa = forest.taxon.data();
    int32_t* const ba = below_a.data();
    int32_t* const bb = below_b.data();
    std::fill(ba + seg[s0], ba + seg[s1], -1);
    std::fill(bb + seg[s0], bb + seg[s1], -1);
    for (uint32_t s = s0; s < s1; ++s) {
      const std::size_t first = seg[s];
      const std::size_t last = seg[s + 1];
      uint32_t ca = 0;
      uint32_t cb = 0;
      for (std::size_t i = first; i < last; ++i) {
        const int32_t taxon = taxa[i];
        const int32_t self = static_cast<int32_t>(i);
        if (taxon >= 0) (taxon < split ? ba[i] : bb[i]) = self;
        // Every child has reported, so i's maps are final here.
        const int32_t a = ba[i];
        const int32_t b = bb[i];
        ca += a == self ? 1u : 0u;
        cb += b == self ? 1u : 0u;
        const int32_t p = parent[i];
        if (p < 0) continue;
        const std::size_t up = static_cast<std::size_t>(p);
        // The first contributing child supplies its LCA; two or more
        // children make the parent the LCA. The maps also serve as the
        // accumulators.
        if (a >= 0) ba[up] = ba[up] < 0 ? a : p;
        if (b >= 0) bb[up] = bb[up] < 0 ? b : p;
      }
      start_a[s + 1] = ca;
      start_b[s + 1] = cb;
    }
  });
  // start_*[s]: the index in that half of the first kept vertex of segment s.
  for (std::size_t s = 0; s < segments; ++s) {
    start_a[s + 1] += start_a[s];
    start_b[s + 1] += start_b[s];
  }

  // Restrict both halves in the same pass. A kept vertex's index stays fixed;
  // a suppressed vertex takes its nearest kept ancestor's index in the
  // parents-first sweep. Every nonempty below map names a kept vertex of its
  // own segment, so these same arrays then translate the maps. All scratch is
  // O(m) and released before descent: O(N + kn log n) time and O(N + kn) peak
  // words.
  Forest half_a;
  Forest half_b;
  half_a.parent.resize(start_a[segments]);
  half_a.taxon.resize(start_a[segments]);
  half_b.parent.resize(start_b[segments]);
  half_b.taxon.resize(start_b[segments]);
  {
    RawVec<int32_t> index_a(m);
    RawVec<int32_t> index_b(m);
    for_chunks(plan, [&](uint32_t s0, uint32_t s1) {
      const int32_t* const parent = forest.parent.data();
      const int32_t* const taxa = forest.taxon.data();
      int32_t* const ba = below_a.data();
      int32_t* const bb = below_b.data();
      int32_t* const ia = index_a.data();
      int32_t* const ib = index_b.data();
      int32_t* const parent_a = half_a.parent.data();
      int32_t* const taxon_a = half_a.taxon.data();
      int32_t* const parent_b = half_b.parent.data();
      int32_t* const taxon_b = half_b.taxon.data();
      const std::size_t chunk_first = seg[s0];
      const std::size_t chunk_last = seg[s1];
      for (uint32_t s = s0; s < s1; ++s) {
        const std::size_t first = seg[s];
        const std::size_t last = seg[s + 1];
        int32_t kept_a = static_cast<int32_t>(start_a[s]);
        int32_t kept_b = static_cast<int32_t>(start_b[s]);
        for (std::size_t i = first; i < last; ++i) {
          const int32_t self = static_cast<int32_t>(i);
          ia[i] = ba[i] == self ? kept_a++ : -1;
          ib[i] = bb[i] == self ? kept_b++ : -1;
        }
      }
      for (std::size_t i = chunk_last; i-- > chunk_first;) {
        const int32_t p = parent[i];
        const int32_t up_a = p < 0 ? -1 : ia[static_cast<std::size_t>(p)];
        const int32_t up_b = p < 0 ? -1 : ib[static_cast<std::size_t>(p)];
        const int32_t own_a = ia[i];
        const int32_t own_b = ib[i];
        if (own_a >= 0) {
          parent_a[static_cast<std::size_t>(own_a)] = up_a;
          taxon_a[static_cast<std::size_t>(own_a)] = taxa[i];
        } else {
          ia[i] = up_a;
        }
        if (own_b >= 0) {
          parent_b[static_cast<std::size_t>(own_b)] = up_b;
          taxon_b[static_cast<std::size_t>(own_b)] = taxa[i];
        } else {
          ib[i] = up_b;
        }
      }
      for (std::size_t i = chunk_first; i < chunk_last; ++i) {
        if (ba[i] >= 0) ba[i] = ia[static_cast<std::size_t>(ba[i])];
        if (bb[i] >= 0) bb[i] = ib[static_cast<std::size_t>(bb[i])];
      }
    });
  }
  // The halves keep the segments that kept a vertex.
  std::vector<uint32_t> seg_a{0};
  std::vector<uint32_t> seg_b{0};
  for (std::size_t s = 0; s < segments; ++s) {
    if (start_a[s + 1] > start_a[s]) seg_a.push_back(start_a[s + 1]);
    if (start_b[s + 1] > start_b[s]) seg_b.push_back(start_b[s + 1]);
  }
  std::vector<uint32_t>().swap(start_a);
  std::vector<uint32_t>().swap(start_b);
  std::vector<uint32_t>().swap(seg);
  RawVec<int32_t>().swap(forest.parent);
  RawVec<int32_t>().swap(forest.taxon);
  const std::size_t size_a = half_a.parent.size();
  const std::size_t size_b = half_b.parent.size();

  // The identifier of each vertex's part in a half, 0 for an empty part,
  // packed as the key part_a * 2^bits_b + part_b, whose order is the pair
  // order. Only the maps of this call stay live while a half is solved.
  const int bits_b = std::bit_width(size_b);
  RawVec<uint32_t> part_a;
  {
    const RawVec<uint32_t> ids_a =
        halve(lo, mid, std::move(half_a), std::move(seg_a));
    part_a.resize(m);
    for_ranges(m, plan.parts, [&](int, std::size_t b, std::size_t e) {
      const int32_t* const below = below_a.data();
      for (std::size_t i = b; i < e; ++i) {
        part_a[i] = below[i] >= 0 ? ids_a[static_cast<std::size_t>(below[i])]
                                  : 0;
      }
    });
    RawVec<int32_t>().swap(below_a);
  }
  RawVec<uint64_t> key;
  {
    const RawVec<uint32_t> ids_b =
        halve(mid, hi, std::move(half_b), std::move(seg_b));
    key.resize(m);
    for_ranges(m, plan.parts, [&](int, std::size_t b, std::size_t e) {
      const int32_t* const below = below_b.data();
      for (std::size_t i = b; i < e; ++i) {
        const uint32_t part_b =
            below[i] >= 0 ? ids_b[static_cast<std::size_t>(below[i])] : 0;
        key[i] = (static_cast<uint64_t>(part_a[i]) << bits_b) | part_b;
      }
    });
    RawVec<int32_t>().swap(below_b);
    RawVec<uint32_t>().swap(part_a);
  }

  // Rank the keys.
  if (parallel) {
    const RawVec<uint32_t> item =
        sort_records(key, std::bit_width(size_a) + bits_b, plan.parts);
    return rank_sorted(key, item, plan.parts);
  }
  // Counting sorts by part_b, then by part_a.
  return rank_serial(key, [&](RawVec<uint32_t>& order,
                              RawVec<uint32_t>& scratch) {
    const uint64_t low = (uint64_t{1} << bits_b) - 1;
    std::vector<uint32_t> count;
    counting_sort(order, scratch, size_b + 1, count,
                  [&](uint32_t i) { return key[i] & low; });
    counting_sort(scratch, order, size_a + 1, count,
                  [&](uint32_t i) { return key[i] >> bits_b; });
  });
}

// Checks the forest, then numbers its clusters; `seg` lists segments that no
// parent link crosses.
RawVec<uint32_t> cluster_ids(Forest&& forest,
                             std::vector<uint32_t>&& seg,
                             int taxon_count) {
  const std::size_t m = forest.parent.size();
  {
    // Children first, leaves childless, taxa in range, a leaf below every
    // vertex.
    std::vector<uint8_t> has_leaf(m, 0);
    for (std::size_t i = 0; i < m; ++i) {
      const int32_t taxon = forest.taxon[i];
      if (taxon >= taxon_count) {
        throw std::invalid_argument("cluster forest taxon out of range");
      }
      if (taxon >= 0) {
        if (has_leaf[i] != 0) {
          throw std::invalid_argument("cluster forest leaf has children");
        }
        has_leaf[i] = 1;
      } else if (has_leaf[i] == 0) {
        throw std::invalid_argument("cluster forest vertex has no leaf");
      }
      const int32_t p = forest.parent[i];
      if (p < -1 || (p >= 0 && (static_cast<std::size_t>(p) <= i ||
                                static_cast<std::size_t>(p) >= m))) {
        throw std::invalid_argument(
            "cluster forest is not listed children first");
      }
      if (p >= 0) has_leaf[static_cast<std::size_t>(p)] = 1;
    }
  }
  if (m == 0) return {};
  return halve(0, taxon_count, std::move(forest), std::move(seg));
}

// Scratch for building one gene's part of the forest.
struct GeneForestScratch {
  explicit GeneForestScratch(int taxon_count) : copies(taxon_count) {}
  GeneCopyScan copies;
  std::vector<int32_t> closed_above;
  std::vector<int32_t> vertex;
  std::vector<int32_t> parent_node;
  GeneTraversal traversal;
};

// The forest vertices of one gene, rooted at its reference root: their number
// when `forest` is null; otherwise they are also written at positions first,
// first + 1, ... and ids[v] is set to node v's vertex plus one.
std::size_t gene_forest(const GeneIndex& gene,
                        GeneForestScratch& s,
                        Forest* forest,
                        std::size_t first,
                        uint32_t* ids) {
  if (gene.nodes.empty()) return 0;
  const GeneTraversal& tr = gene_reference_traversal(gene, s.traversal);
  const std::size_t node_count = gene.nodes.size();
  const int root = gene.reference_root;
  // The MRE ignores genes with fewer than four distinct taxa.
  if (s.copies.scan_taxa(gene, tr) < 4) return 0;
  s.copies.scan_nodes(gene, tr);
  // A node is a forest vertex when it is closed and holds at least two
  // distinct taxa.
  const auto in_forest = [&](std::size_t v) {
    return s.copies.closed(static_cast<int>(v)) &&
           s.copies.distinct(static_cast<int>(v)) >= 2;
  };
  // closed_above[v]: v itself when it is a forest vertex, else its nearest
  // such ancestor, or -1.
  s.closed_above.assign(node_count, -1);
  for (int node : tr.preorder) {
    const std::size_t v = static_cast<std::size_t>(node);
    if (in_forest(v)) {
      s.closed_above[v] = node;
    } else if (node != root) {
      s.closed_above[v] =
          s.closed_above[static_cast<std::size_t>(tr.parent[v])];
    }
  }
  // Each taxon's leaf at its first copy, when a forest vertex holds it.
  const auto taxon_leaf = [&](std::size_t v, int taxon) {
    return taxon >= 0 &&
           s.copies.first(taxon) == s.copies.tin(static_cast<int>(v)) &&
           s.closed_above[v] >= 0;
  };
  if (forest == nullptr) {
    std::size_t count = 0;
    for (int node : tr.preorder) {
      const std::size_t v = static_cast<std::size_t>(node);
      if (taxon_leaf(v, gene.nodes[v].taxon())) ++count;
      if (in_forest(v)) ++count;
    }
    return count;
  }
  // Vertices in reverse preorder, so children come first: each taxon's leaf,
  // then the node's own vertex. A vertex's parent is the lowest forest vertex
  // strictly above it (for a leaf: holding its taxon), which comes later in
  // this order, so the gene node is recorded first and resolved after the
  // gene.
  s.vertex.assign(node_count, -1);
  s.parent_node.clear();
  std::size_t next = first;
  for (auto it = tr.preorder.rbegin(); it != tr.preorder.rend(); ++it) {
    const std::size_t v = static_cast<std::size_t>(*it);
    const int taxon = gene.nodes[v].taxon();
    if (taxon_leaf(v, taxon)) {
      s.parent_node.push_back(s.closed_above[v]);
      forest->taxon[next] = taxon;
      forest->parent[next] = -1;
      ++next;
    }
    if (in_forest(v)) {
      s.vertex[v] = static_cast<int32_t>(next);
      s.parent_node.push_back(
          s.closed_above[static_cast<std::size_t>(tr.parent[v])]);
      forest->taxon[next] = -1;
      forest->parent[next] = -1;
      ids[v] = static_cast<uint32_t>(next) + 1;
      ++next;
    }
  }
  for (std::size_t j = 0; j < s.parent_node.size(); ++j) {
    const int32_t node = s.parent_node[j];
    if (node >= 0) {
      forest->parent[first + j] = s.vertex[static_cast<std::size_t>(node)];
    }
  }
  return next - first;
}

}  // namespace

std::vector<uint32_t> exact_cluster_ids(ClusterForest&& forest,
                                        int taxon_count) {
  const std::size_t m = forest.parent.size();
  if (forest.taxon.size() != m) {
    throw std::invalid_argument("cluster forest arrays differ in size");
  }
  if (m >= static_cast<std::size_t>(std::numeric_limits<int32_t>::max())) {
    throw std::length_error("cluster forest exceeds its index domain");
  }
  // Segments end where every vertex so far has its parent among them. A
  // parent outside the forest leaves one segment, which the checks reject.
  std::vector<uint32_t> seg{0};
  int64_t reach = -1;
  for (std::size_t i = 0; i < m; ++i) {
    reach = std::max<int64_t>(reach, forest.parent[i]);
    if (reach <= static_cast<int64_t>(i)) {
      seg.push_back(static_cast<uint32_t>(i + 1));
    }
  }
  if (seg.back() != m) seg.push_back(static_cast<uint32_t>(m));
  Forest own;
  own.parent.assign(forest.parent.begin(), forest.parent.end());
  own.taxon.assign(forest.taxon.begin(), forest.taxon.end());
  std::vector<int32_t>().swap(forest.parent);
  std::vector<int32_t>().swap(forest.taxon);
  const RawVec<uint32_t> ids =
      cluster_ids(std::move(own), std::move(seg), taxon_count);
  return std::vector<uint32_t>(ids.begin(), ids.end());
}

GeneCopyScan::GeneCopyScan(int taxon_count)
    : taxon_count_(taxon_count),
      first_(static_cast<std::size_t>(std::max(0, taxon_count)), -1),
      last_(static_cast<std::size_t>(std::max(0, taxon_count)), -1) {
  if (taxon_count < 0) {
    throw std::invalid_argument("negative taxon count");
  }
}

int GeneCopyScan::scan_taxa(const GeneIndex& gene, const GeneTraversal& tr) {
  for (int32_t taxon : touched_) {
    first_[static_cast<std::size_t>(taxon)] = -1;
    last_[static_cast<std::size_t>(taxon)] = -1;
  }
  touched_.clear();
  root_ = gene.reference_root;
  tin_.assign(gene.nodes.size(), -1);
  for (std::size_t pos = 0; pos < tr.preorder.size(); ++pos) {
    const int node = tr.preorder[pos];
    tin_[static_cast<std::size_t>(node)] = static_cast<int32_t>(pos);
    const int taxon = gene.nodes[static_cast<std::size_t>(node)].taxon();
    if (taxon < 0) continue;
    if (taxon >= taxon_count_) {
      throw std::invalid_argument("gene taxon lies outside the split universe");
    }
    const std::size_t t = static_cast<std::size_t>(taxon);
    if (first_[t] < 0) {
      first_[t] = static_cast<int32_t>(pos);
      touched_.push_back(taxon);
    }
    last_[t] = static_cast<int32_t>(pos);
  }
  return static_cast<int>(touched_.size());
}

void GeneCopyScan::scan_nodes(const GeneIndex& gene, const GeneTraversal& tr) {
  // Subtree size, the smallest first and largest last copy position of the
  // taxa below, and the number of distinct taxa below (counted at each
  // taxon's first copy).
  const std::size_t node_count = gene.nodes.size();
  size_.assign(node_count, 1);
  min_first_.assign(node_count, std::numeric_limits<int32_t>::max());
  max_last_.assign(node_count, -1);
  distinct_.assign(node_count, 0);
  for (int node : tr.preorder) {
    const std::size_t v = static_cast<std::size_t>(node);
    const int taxon = gene.nodes[v].taxon();
    if (taxon < 0) continue;
    const std::size_t t = static_cast<std::size_t>(taxon);
    min_first_[v] = first_[t];
    max_last_[v] = last_[t];
    distinct_[v] = first_[t] == tin_[v] ? 1 : 0;
  }
  for (auto it = tr.preorder.rbegin(); it != tr.preorder.rend(); ++it) {
    const int node = *it;
    if (node == root_) continue;
    const std::size_t v = static_cast<std::size_t>(node);
    const std::size_t p = static_cast<std::size_t>(tr.parent[v]);
    size_[p] += size_[v];
    min_first_[p] = std::min(min_first_[p], min_first_[v]);
    max_last_[p] = std::max(max_last_[p], max_last_[v]);
    distinct_[p] += distinct_[v];
  }
}

bool GeneCopyScan::closed(int node) const {
  const std::size_t v = static_cast<std::size_t>(node);
  return node != root_ && min_first_[v] >= tin_[v] &&
         max_last_[v] < tin_[v] + size_[v];
}

GeneClusterIds exact_gene_cluster_ids(const std::vector<GeneIndex>& genes,
                                      int taxon_count) {
  if (taxon_count < 0) {
    throw std::invalid_argument("negative taxon count");
  }
  const auto started = std::chrono::steady_clock::now();
  GeneClusterIds out;
  out.offset.assign(genes.size() + 1, 0);
  for (std::size_t g = 0; g < genes.size(); ++g) {
    const std::size_t nodes = genes[g].nodes.size();
    if (nodes > std::numeric_limits<std::size_t>::max() - out.offset[g]) {
      throw std::length_error("MRE cluster ids exceed the addressable size");
    }
    out.offset[g + 1] = out.offset[g] + nodes;
  }
  out.id.assign(out.offset.back(), 0);

  // Each gene's vertices form one segment of the forest, in gene order;
  // first[g] is the gene's first vertex. With several workers each gene is
  // counted, then written at its offset; with one, written in one pass.
  const int gene_count = static_cast<int>(genes.size());
  const int workers =
      std::min(current_inner_parallel_worker_limit(), std::max(gene_count, 1));
  std::vector<std::size_t> first(genes.size() + 1, 0);
  Forest forest;
  // A gene adds at most two vertices per node: its own and a taxon leaf.
  const auto check_fits = [](std::size_t vertices) {
    if (vertices >= static_cast<std::size_t>(
                        std::numeric_limits<int32_t>::max())) {
      throw std::length_error("MRE cluster forest exceeds its index domain");
    }
  };
  if (workers <= 1) {
    // The arrays hold only written vertices between genes, so growth copies
    // no unset element.
    GeneForestScratch scratch(taxon_count);
    for (std::size_t g = 0; g < genes.size(); ++g) {
      const std::size_t at = first[g];
      const std::size_t bound = at + 2 * genes[g].nodes.size();
      check_fits(bound);
      if (forest.parent.capacity() < bound) {
        const std::size_t grown =
            std::max(bound, 2 * forest.parent.capacity());
        forest.parent.reserve(grown);
        forest.taxon.reserve(grown);
      }
      forest.parent.resize(bound);
      forest.taxon.resize(bound);
      first[g + 1] = at + gene_forest(genes[g], scratch, &forest, at,
                                      out.id.data() + out.offset[g]);
      forest.parent.resize(first[g + 1]);
      forest.taxon.resize(first[g + 1]);
    }
  } else {
    const auto gene_range = [&](int w) {
      return std::pair<int, int>(
          static_cast<int>(static_cast<int64_t>(gene_count) * w / workers),
          static_cast<int>(static_cast<int64_t>(gene_count) * (w + 1) /
                           workers));
    };
    parallel_for_inner(0, workers, [&](int w) {
      GeneForestScratch scratch(taxon_count);
      const auto [g0, g1] = gene_range(w);
      for (int g = g0; g < g1; ++g) {
        first[static_cast<std::size_t>(g) + 1] = gene_forest(
            genes[static_cast<std::size_t>(g)], scratch, nullptr, 0, nullptr);
      }
    }, 1);
    // The same bound as the serial build, so both fail on the same input.
    for (std::size_t g = 0; g < genes.size(); ++g) {
      check_fits(first[g] + 2 * genes[g].nodes.size());
      first[g + 1] += first[g];
    }
    forest.parent.resize(first.back());
    forest.taxon.resize(first.back());
    parallel_for_inner(0, workers, [&](int w) {
      GeneForestScratch scratch(taxon_count);
      const auto [g0, g1] = gene_range(w);
      for (int g = g0; g < g1; ++g) {
        const std::size_t gg = static_cast<std::size_t>(g);
        gene_forest(genes[gg], scratch, &forest, first[gg],
                    out.id.data() + out.offset[gg]);
      }
    }, 1);
  }
  std::vector<uint32_t> seg{0};
  for (std::size_t g = 0; g < genes.size(); ++g) {
    if (first[g + 1] > first[g]) {
      seg.push_back(static_cast<uint32_t>(first[g + 1]));
    }
  }
  stats_add_time(run_stats().mre_ids_forest_ns,
                 std::chrono::steady_clock::now() - started);

  // Replace each vertex index (stored plus one) by the vertex's identifier.
  const RawVec<uint32_t> ids =
      cluster_ids(std::move(forest), std::move(seg), taxon_count);
  const std::size_t entries = out.id.size();
  const int parts = static_cast<int>(std::clamp<std::size_t>(
      entries / kNodesPerPart, 1,
      static_cast<std::size_t>(current_inner_parallel_worker_limit())));
  for_ranges(entries, parts, [&](int, std::size_t b, std::size_t e) {
    for (std::size_t j = b; j < e; ++j) {
      uint32_t& id = out.id[j];
      if (id != 0) id = ids[static_cast<std::size_t>(id) - 1];
    }
  });
  stats_add_time(run_stats().mre_ids_ns,
                 std::chrono::steady_clock::now() - started);
  return out;
}

}  // namespace qraft::internal
