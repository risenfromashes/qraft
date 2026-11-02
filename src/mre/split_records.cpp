#include "mre/split_records.hpp"

#include "internal.hpp"

#include <algorithm>
#include <limits>
#include <span>
#include <stdexcept>
#include <tuple>

namespace qraft::internal {
namespace {

// Every pass is a stable counting sort on one 16-bit digit, O(items + 2^16)
// per part. A pass splits its input into up to one part per worker, each of at
// least kRadixItemsPerPart items, so the parts' bucket counts stay small next
// to the items they order.
constexpr int kRadixDigitBits = 16;
constexpr std::size_t kRadixBuckets = std::size_t{1} << kRadixDigitBits;
constexpr std::size_t kRadixItemsPerPart = std::size_t{1} << 18;
// Buckets per task when the offsets are computed in parallel.
constexpr std::size_t kRadixBucketBlock = 1024;

int radix_part_count(std::size_t items) {
  const std::size_t by_size =
      std::max<std::size_t>(1, items / kRadixItemsPerPart);
  return static_cast<int>(std::min<std::size_t>(
      by_size,
      static_cast<std::size_t>(current_inner_parallel_worker_limit())));
}

// Turns per-part bucket counts (part-major) into scatter offsets: within a
// bucket the parts' items are placed in part order, after every smaller
// bucket's items.
void radix_scatter_offsets(std::vector<std::size_t>& counts,
                           int parts,
                           std::vector<std::size_t>& totals) {
  if (parts == 1) {
    std::size_t prefix = 0;
    for (std::size_t& count : counts) {
      const std::size_t bucket_size = count;
      count = prefix;
      prefix += bucket_size;
    }
    return;
  }
  totals.resize(kRadixBuckets);
  const int blocks = static_cast<int>(kRadixBuckets / kRadixBucketBlock);
  parallel_for_inner(0, blocks, [&](int block) {
    const std::size_t first = static_cast<std::size_t>(block) * kRadixBucketBlock;
    for (std::size_t b = first; b < first + kRadixBucketBlock; ++b) {
      std::size_t running = 0;
      for (int part = 0; part < parts; ++part) {
        std::size_t& count =
            counts[static_cast<std::size_t>(part) * kRadixBuckets + b];
        const std::size_t part_size = count;
        count = running;
        running += part_size;
      }
      totals[b] = running;
    }
  }, 1);
  std::size_t prefix = 0;
  for (std::size_t& total : totals) {
    const std::size_t bucket_size = total;
    total = prefix;
    prefix += bucket_size;
  }
  parallel_for_inner(0, blocks, [&](int block) {
    const std::size_t first = static_cast<std::size_t>(block) * kRadixBucketBlock;
    for (int part = 0; part < parts; ++part) {
      std::size_t* offsets =
          counts.data() + static_cast<std::size_t>(part) * kRadixBuckets;
      for (std::size_t b = first; b < first + kRadixBucketBlock; ++b) {
        offsets[b] += totals[b];
      }
    }
  }, 1);
}

// The stable passes of an LSD radix sort. The first pass reads the input, part
// by part in order, and later passes alternate between the output vector and
// a scratch vector; finish() leaves the sorted items in the output. Every pass
// equals a serial stable pass over the concatenated input, whatever the part
// count, so the order does not depend on the workers.
template <class Value>
class RadixPasses {
 public:
  // Sorts `values` in place.
  RadixPasses(std::vector<Value>& values,
              long long& pass_count,
              long long& item_count)
      : output_(values),
        current_(&values),
        pass_count_(pass_count),
        item_count_(item_count) {}

  // Sorts the items of `chunks`, taken in order, into `output`. The chunks
  // are freed after the first pass.
  RadixPasses(std::vector<std::vector<Value>>& chunks,
              std::vector<Value>& output,
              long long& pass_count,
              long long& item_count)
      : output_(output),
        chunks_(&chunks),
        pass_count_(pass_count),
        item_count_(item_count) {}

  // Stable passes on a 32-bit field, low digit first: two passes.
  template <class Getter>
  void by_u32(Getter&& getter) {
    for (int shift = 0; shift < 32; shift += kRadixDigitBits) {
      pass([&](const Value& value) {
        return static_cast<std::uint16_t>(getter(value) >> shift);
      });
    }
  }

  void finish() {
    if (current_ != nullptr && current_ != &output_) output_.swap(*current_);
    std::vector<Value>().swap(scratch_);
  }

 private:
  template <class Digit>
  void pass(Digit&& digit) {
    // The input as spans in order, and each part's first span; a part of the
    // chunked input spans consecutive chunks.
    std::vector<std::span<const Value>> spans;
    std::vector<std::size_t> part_first;
    std::vector<Value>* destination = nullptr;
    if (chunks_ != nullptr) {
      std::size_t items = 0;
      for (const std::vector<Value>& chunk : *chunks_) items += chunk.size();
      const std::size_t part_count =
          static_cast<std::size_t>(radix_part_count(items));
      std::size_t seen = 0;
      for (const std::vector<Value>& chunk : *chunks_) {
        if (chunk.empty()) continue;
        // A new part starts at the first chunk past each part's share.
        if (part_first.empty() ||
            seen * part_count >= items * part_first.size()) {
          part_first.push_back(spans.size());
        }
        spans.emplace_back(chunk);
        seen += chunk.size();
      }
      // The allocator keeps pages freed before the sort resident until they
      // decay; return them before the output is allocated next to the chunks.
      release_retained_allocator_pages();
      output_.resize(items);
      destination = &output_;
    } else {
      const std::size_t items = current_->size();
      const int part_count = radix_part_count(items);
      for (int part = 0; part < part_count; ++part) {
        const std::size_t first = items * static_cast<std::size_t>(part) /
                                  static_cast<std::size_t>(part_count);
        const std::size_t last = items * static_cast<std::size_t>(part + 1) /
                                 static_cast<std::size_t>(part_count);
        part_first.push_back(spans.size());
        spans.emplace_back(current_->data() + first, last - first);
      }
      destination = current_ == &output_ ? &scratch_ : &output_;
      destination->resize(items);
    }
    const int part_count = static_cast<int>(part_first.size());
    part_first.push_back(spans.size());
    const auto for_each_value = [&](int part, auto&& visit) {
      for (std::size_t span = part_first[static_cast<std::size_t>(part)];
           span < part_first[static_cast<std::size_t>(part) + 1]; ++span) {
        for (const Value& value : spans[span]) visit(value);
      }
    };
    if (part_count > 0) {
      counts_.assign(static_cast<std::size_t>(part_count) * kRadixBuckets, 0);
      parallel_for_inner(0, part_count, [&](int part) {
        std::size_t* local =
            counts_.data() + static_cast<std::size_t>(part) * kRadixBuckets;
        for_each_value(part, [&](const Value& value) {
          ++local[static_cast<std::size_t>(digit(value))];
        });
      }, 1);
      radix_scatter_offsets(counts_, part_count, totals_);
      Value* out = destination->data();
      parallel_for_inner(0, part_count, [&](int part) {
        std::size_t* local =
            counts_.data() + static_cast<std::size_t>(part) * kRadixBuckets;
        for_each_value(part, [&](const Value& value) {
          out[local[static_cast<std::size_t>(digit(value))]++] = value;
        });
      }, 1);
    }
    if (chunks_ != nullptr) {
      for (std::vector<Value>& chunk : *chunks_) std::vector<Value>().swap(chunk);
      chunks_ = nullptr;
      // The allocator keeps freed chunks resident until they decay; return
      // them before the scratch of the same size is allocated.
      release_retained_allocator_pages();
    }
    current_ = destination;
    ++pass_count_;
    item_count_ += static_cast<long long>(destination->size());
  }

  std::vector<Value>& output_;
  std::vector<std::vector<Value>>* chunks_ = nullptr;
  std::vector<Value>* current_ = nullptr;
  std::vector<Value> scratch_;
  std::vector<std::size_t> counts_;
  std::vector<std::size_t> totals_;
  long long& pass_count_;
  long long& item_count_;
};

// Sorts the records of `chunks`, in order, by (key, gene, edge) into `sorted`.
// Keys hold a size and an id below 2^32 in h[0] (SplitKey256), so the stable
// LSD schedule (edge, gene, id, size) of 8 passes orders the full key; a key
// with another shape throws std::logic_error before any pass.
void radix_sort_occurrences(
    std::vector<std::vector<SplitOccurrence256>>& chunks,
    std::vector<SplitOccurrence256>& sorted,
    SplitRadixStats& stats) {
  sorted.clear();
  bool empty = true;
  for (const std::vector<SplitOccurrence256>& chunk : chunks) {
    empty = empty && chunk.empty();
    for (const SplitOccurrence256& value : chunk) {
      if (value.key.h[0] > std::numeric_limits<uint32_t>::max() ||
          value.key.h[1] != 0 || value.key.h[2] != 0 || value.key.h[3] != 0) {
        throw std::logic_error("split key is not a size and a 32-bit id");
      }
    }
  }
  if (empty) return;
  RadixPasses<SplitOccurrence256> passes(chunks, sorted,
                                         stats.occurrence_passes,
                                         stats.occurrence_items);
  passes.by_u32([](const SplitOccurrence256& value) { return value.edge; });
  passes.by_u32([](const SplitOccurrence256& value) { return value.gene; });
  passes.by_u32([](const SplitOccurrence256& value) {
    return static_cast<uint32_t>(value.key.h[0]);
  });
  passes.by_u32(
      [](const SplitOccurrence256& value) { return value.key.size; });
  passes.finish();

  if (stats.occurrence_passes != 8) {
    throw std::logic_error("split-occurrence radix schedule mismatch");
  }
}

// One aggregate per key run of the sorted records [first, last), which start
// and end at run boundaries, with split ids counted from zero.
struct AggregatedPart {
  std::vector<SplitAggregate256> aggregates;
  long long gene_splits = 0;
};

AggregatedPart aggregate_sorted_range(
    const std::vector<SplitOccurrence256>& occurrences,
    std::size_t first,
    std::size_t last) {
  AggregatedPart part;
  std::size_t begin = first;
  while (begin < last) {
    const SplitOccurrence256& head = occurrences[begin];
    SplitAggregate256 aggregate;
    aggregate.key = head.key;
    aggregate.split_id = static_cast<std::uint32_t>(part.aggregates.size());
    aggregate.first_tree = head.gene;
    aggregate.first_edge = head.edge;
    aggregate.rep_tree = head.gene;
    aggregate.rep_node = head.node;
    aggregate.one_taxon = head.one_taxon;

    // A key's run is sorted by gene, so each change of gene is one more
    // supporting gene.
    std::uint32_t previous_gene = std::numeric_limits<std::uint32_t>::max();
    std::size_t end = begin;
    while (end < last && occurrences[end].key == head.key) {
      const std::uint32_t gene = occurrences[end].gene;
      if (gene != previous_gene) {
        if (aggregate.count == std::numeric_limits<std::uint32_t>::max()) {
          throw std::length_error("split support exceeds its 32-bit domain");
        }
        ++aggregate.count;
        ++part.gene_splits;
        previous_gene = gene;
      }
      ++end;
    }
    part.aggregates.push_back(aggregate);
    begin = end;
  }
  return part;
}

// Aggregates the sorted records in parallel parts that begin at key runs, and
// numbers the splits by their global key rank.
std::vector<SplitAggregate256> aggregate_sorted_occurrences(
    const std::vector<SplitOccurrence256>& occurrences,
    SplitRadixStats& stats) {
  const std::size_t items = occurrences.size();
  const int part_count = radix_part_count(items);
  std::vector<std::size_t> bounds(static_cast<std::size_t>(part_count) + 1,
                                  items);
  bounds[0] = 0;
  for (int part = 1; part < part_count; ++part) {
    std::size_t bound = std::max(
        bounds[static_cast<std::size_t>(part) - 1],
        items * static_cast<std::size_t>(part) /
            static_cast<std::size_t>(part_count));
    while (bound > 0 && bound < items &&
           occurrences[bound].key == occurrences[bound - 1].key) {
      ++bound;
    }
    bounds[static_cast<std::size_t>(part)] = bound;
  }
  std::vector<AggregatedPart> parts(static_cast<std::size_t>(part_count));
  parallel_for_inner(0, part_count, [&](int part) {
    parts[static_cast<std::size_t>(part)] = aggregate_sorted_range(
        occurrences, bounds[static_cast<std::size_t>(part)],
        bounds[static_cast<std::size_t>(part) + 1]);
  }, 1);

  std::vector<std::size_t> offsets(static_cast<std::size_t>(part_count) + 1, 0);
  long long gene_splits = 0;
  for (int part = 0; part < part_count; ++part) {
    offsets[static_cast<std::size_t>(part) + 1] =
        offsets[static_cast<std::size_t>(part)] +
        parts[static_cast<std::size_t>(part)].aggregates.size();
    gene_splits += parts[static_cast<std::size_t>(part)].gene_splits;
  }
  const std::size_t distinct = offsets.back();
  if (distinct > static_cast<std::size_t>(
                     std::numeric_limits<std::uint32_t>::max()) + 1u) {
    throw std::length_error("dense split id exceeds its 32-bit domain");
  }
  std::vector<SplitAggregate256> aggregates(distinct);
  parallel_for_inner(0, part_count, [&](int part) {
    const std::size_t offset = offsets[static_cast<std::size_t>(part)];
    std::vector<SplitAggregate256>& local =
        parts[static_cast<std::size_t>(part)].aggregates;
    for (std::size_t i = 0; i < local.size(); ++i) {
      SplitAggregate256 aggregate = local[i];
      aggregate.split_id = static_cast<std::uint32_t>(offset + i);
      aggregates[offset + i] = aggregate;
    }
    std::vector<SplitAggregate256>().swap(local);
  }, 1);
  stats.deduplicated_gene_splits = gene_splits;
  stats.distinct_splits = static_cast<long long>(aggregates.size());
  return aggregates;
}

}  // namespace

bool split_key_less(const SplitKey256& left, const SplitKey256& right) {
  if (left.size != right.size) return left.size < right.size;
  return left.h < right.h;
}

bool split_priority_less_by_key(const SplitAggregate256& left,
                                const SplitAggregate256& right) {
  if (left.count != right.count) return left.count > right.count;
  if (left.first_tree != right.first_tree) {
    return left.first_tree < right.first_tree;
  }
  if (left.first_edge != right.first_edge) {
    return left.first_edge < right.first_edge;
  }
  return split_key_less(left.key, right.key);
}

std::vector<SplitAggregate256> aggregate_split_occurrences256(
    std::vector<std::vector<SplitOccurrence256>>& chunks,
    SplitRadixStats& stats) {
  std::size_t items = 0;
  for (const std::vector<SplitOccurrence256>& chunk : chunks) {
    if (chunk.size() > std::numeric_limits<std::size_t>::max() - items) {
      throw std::length_error(
          "split records exceed the addressable size domain");
    }
    items += chunk.size();
  }
  if (items > static_cast<std::size_t>(std::numeric_limits<long long>::max())) {
    throw std::length_error("split occurrence count exceeds telemetry domain");
  }
  stats.emitted_occurrences = static_cast<long long>(items);
  std::vector<SplitOccurrence256> sorted;
  radix_sort_occurrences(chunks, sorted, stats);
  std::vector<SplitAggregate256> aggregates =
      aggregate_sorted_occurrences(sorted, stats);
  return aggregates;
}

std::vector<SplitAggregate256> aggregate_split_occurrences256(
    std::vector<SplitOccurrence256>& occurrences,
    SplitRadixStats& stats) {
  std::vector<std::vector<SplitOccurrence256>> chunks(1);
  chunks.front().swap(occurrences);
  stats.emitted_occurrences = static_cast<long long>(chunks.front().size());
  radix_sort_occurrences(chunks, occurrences, stats);
  return aggregate_sorted_occurrences(occurrences, stats);
}

void radix_sort_split_priority256(
    std::vector<SplitAggregate256>& aggregates,
    SplitRadixStats& stats) {
  if (aggregates.empty()) return;
  RadixPasses<SplitAggregate256> passes(aggregates, stats.priority_passes,
                                        stats.priority_items);
  // Stable LSD schedule for (~count, first_tree, first_edge, split_id).
  // Because split_id is the key rank, this is the order of
  // split_priority_less_by_key.
  passes.by_u32(
      [](const SplitAggregate256& value) { return value.split_id; });
  passes.by_u32(
      [](const SplitAggregate256& value) { return value.first_edge; });
  passes.by_u32(
      [](const SplitAggregate256& value) { return value.first_tree; });
  passes.by_u32([](const SplitAggregate256& value) {
    return static_cast<std::uint32_t>(~value.count);
  });
  passes.finish();

  if (stats.priority_passes != kSplitPriorityRadixPasses) {
    throw std::logic_error("split-priority radix schedule mismatch");
  }
}

}  // namespace qraft::internal
