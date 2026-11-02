#include "solver/solver.hpp"
#include "core/selection.hpp"
#include "divide/bipartition_test.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <tuple>

namespace qraft::internal {
namespace {

constexpr int kMinSide = 1;
constexpr double kScoreTieEps = 1e-9;

// A nontrivial edge of the query tree (the starting tree induced on the
// active taxa), proposed as a bipartition. `split` is the side without active
// taxon 0, stored as at most two ranges of query leaf positions.
struct StartingTreeEdge {
  struct LeafSet {
    int first_begin = 0;
    int first_end = 0;
    int second_begin = 0;
    int second_end = 0;
    int count = 0;

    bool operator==(const LeafSet& other) const {
      return first_begin == other.first_begin && first_end == other.first_end &&
             second_begin == other.second_begin &&
             second_end == other.second_end && count == other.count;
    }
  };

  LeafSet split;
  // Set only for the trivial bipartition, which has no edge.
  Candidate materialized_candidate;
  int edge_child_node = -1;
  // Preorder index of edge_child_node.
  int representative_order = std::numeric_limits<int>::max();
  bool materialized = false;
  // Failed the bipartition test in this subproblem.
  bool rejected = false;
};

int effective_min_side(const SolverContext& ctx, int active_count) {
  (void)ctx;
  return std::min(kMinSide, std::max(1, active_count / 2));
}

int split_balance_from_count(int left, int active_count) {
  return std::min(left, active_count - left);
}

// Active taxon 0 versus the rest (1|m-1).
Candidate trivial_bipartition(int m) {
  Candidate candidate;
  candidate.split = Bitset(m);
  if (m > 0) candidate.split.set(0);
  candidate.representative_order = std::numeric_limits<int>::max();
  return candidate;
}

void add_leaf_interval(StartingTreeEdge::LeafSet& set,
                       int begin,
                       int end) {
  if (begin >= end) return;
  if (set.first_begin == set.first_end) {
    set.first_begin = begin;
    set.first_end = end;
  } else {
    set.second_begin = begin;
    set.second_end = end;
  }
  set.count += end - begin;
}

bool leaf_set_contains(const StartingTreeEdge::LeafSet& set,
                       int leaf_pos) {
  return (set.first_begin <= leaf_pos && leaf_pos < set.first_end) ||
         (set.second_begin <= leaf_pos && leaf_pos < set.second_end);
}

bool leaf_set_lex_less(const StartingTreeEdge::LeafSet& a,
                       const StartingTreeEdge::LeafSet& b,
                       const std::vector<int>& active_to_leaf_pos) {
  for (int active = 0; active < static_cast<int>(active_to_leaf_pos.size());
       ++active) {
    const int leaf_pos = active_to_leaf_pos[static_cast<std::size_t>(active)];
    if (leaf_pos < 0) continue;
    const bool av = leaf_set_contains(a, leaf_pos);
    const bool bv = leaf_set_contains(b, leaf_pos);
    if (av != bv) return !av && bv;
  }
  return false;
}

bool leaf_set_storage_less(const StartingTreeEdge::LeafSet& a,
                           const StartingTreeEdge::LeafSet& b) {
  return std::tie(a.first_begin, a.first_end, a.second_begin, a.second_end,
                  a.count) <
         std::tie(b.first_begin, b.first_end, b.second_begin, b.second_end,
                  b.count);
}

Bitset materialize_leaf_set(const QueryTree& query,
                            const StartingTreeEdge::LeafSet& set,
                            int m) {
  Bitset split(m);
  auto add_interval = [&](int begin, int end) {
    begin = std::clamp(begin, 0, static_cast<int>(query.leaf_active_order.size()));
    end = std::clamp(end, begin, static_cast<int>(query.leaf_active_order.size()));
    for (int pos = begin; pos < end; ++pos) {
      const int active = query.leaf_active_order[static_cast<std::size_t>(pos)];
      if (active >= 0 && active < m) split.set(active);
    }
  };
  add_interval(set.first_begin, set.first_end);
  add_interval(set.second_begin, set.second_end);
  return split;
}

// Every edge of `query` with at least two active taxa on each side, one per
// distinct bipartition (the two edges below a degree-two root induce the
// same one), ordered by preorder of the edge's child node.
std::vector<StartingTreeEdge> starting_tree_edges(
    const QueryTree& query) {
  const int m = query.leaf_count;
  if (m < 2 || query.root < 0) return {};

  std::vector<int> active_to_leaf_pos(static_cast<std::size_t>(m), -1);
  for (int pos = 0; pos < static_cast<int>(query.leaf_active_order.size());
       ++pos) {
    const int active = query.leaf_active_order[static_cast<std::size_t>(pos)];
    if (active >= 0 && active < m) {
      active_to_leaf_pos[static_cast<std::size_t>(active)] = pos;
    }
  }
  const int zero_leaf_pos = !active_to_leaf_pos.empty()
                                ? active_to_leaf_pos.front()
                                : -1;

  std::vector<StartingTreeEdge> candidates;
  for (int node : query.preorder) {
    if (node == query.root) continue;
    if (node < 0 || node >= static_cast<int>(query.subtree_leaf_first.size())) {
      continue;
    }
    const int begin = query.subtree_leaf_first[static_cast<std::size_t>(node)];
    const int count = query.subtree_leaf_count[static_cast<std::size_t>(node)];
    const int end = begin + count;
    if (count < 2 || count > m - 2) continue;
    StartingTreeEdge::LeafSet split;
    if (zero_leaf_pos >= begin && zero_leaf_pos < end) {
      add_leaf_interval(split, 0, begin);
      add_leaf_interval(split, end, m);
    } else {
      add_leaf_interval(split, begin, end);
    }

    const int order = query.tin[static_cast<std::size_t>(node)];
    StartingTreeEdge candidate;
    candidate.split = split;
    candidate.representative_order = order;
    candidate.edge_child_node = node;
    candidates.push_back(candidate);
  }

  // Merge equal bipartitions, keeping the edge earliest in preorder.
  std::sort(candidates.begin(), candidates.end(),
            [](const StartingTreeEdge& a,
               const StartingTreeEdge& b) {
              if (leaf_set_storage_less(a.split, b.split)) return true;
              if (leaf_set_storage_less(b.split, a.split)) return false;
              return a.representative_order < b.representative_order;
            });
  candidates.erase(
      std::unique(candidates.begin(), candidates.end(),
                  [](const StartingTreeEdge& a,
                     const StartingTreeEdge& b) {
                    return a.split == b.split;
                  }),
      candidates.end());
  std::sort(candidates.begin(), candidates.end(),
            [&](const StartingTreeEdge& a,
                const StartingTreeEdge& b) {
              if (a.representative_order != b.representative_order) {
                return a.representative_order < b.representative_order;
              }
              return leaf_set_lex_less(a.split, b.split, active_to_leaf_pos);
            });
  return candidates;
}

Candidate materialize_candidate(const QueryTree& query,
                                const StartingTreeEdge& item,
                                int m) {
  if (item.materialized) return item.materialized_candidate;
  Candidate candidate;
  candidate.split = materialize_leaf_set(query, item.split, m);
  candidate.representative_order = item.representative_order;
  return candidate;
}

// The balanced ablation: the unrejected edge with the largest smaller side,
// earliest in preorder among ties, or the trivial bipartition when none is
// eligible.
BipartitionChoice choose_balanced_bipartition(
    const SolverContext& ctx,
    const QueryTree& query,
    const std::vector<StartingTreeEdge>& edge_candidates,
    int active_count) {
  BipartitionChoice best;
  best.score = 0.0;
  best.representative_order = std::numeric_limits<int>::max();
  best.candidate_count = static_cast<int>(edge_candidates.size());
  int best_balance = -1;
  int best_index = -1;
  const int min_side = effective_min_side(ctx, active_count);
  for (int i = 0; i < static_cast<int>(edge_candidates.size()); ++i) {
    const StartingTreeEdge& candidate =
        edge_candidates[static_cast<std::size_t>(i)];
    const int left_count = candidate.split.count;
    const int right_count = active_count - left_count;
    if (candidate.rejected || left_count < min_side || right_count < min_side) {
      continue;
    }
    const int balance = split_balance_from_count(left_count, active_count);
    const int rep = candidate.representative_order;
    if (balance > best_balance ||
        (balance == best_balance && rep < best.representative_order)) {
      best_balance = balance;
      best.representative_order = rep;
      best_index = i;
    }
  }

  if (best_index >= 0) {
    best.candidate_index = best_index;
    best.candidate = materialize_candidate(
        query, edge_candidates[static_cast<std::size_t>(best_index)],
        active_count);
  } else {
    best.candidate = trivial_bipartition(active_count);
    best.representative_order = best.candidate.representative_order;
  }
  const std::vector<Side> assignment =
      assignment_from_split(best.candidate.split, active_count);
  best.left_count =
      static_cast<int>(std::count(assignment.begin(), assignment.end(),
                                  Side::Left));
  best.right_count = static_cast<int>(assignment.size()) - best.left_count;
  return best;
}

// The index of the unrejected edge with the highest score, or -1. Among the
// edges within kScoreTieEps of the exact maximum, the most balanced wins,
// then the earliest in preorder.
int select_scored_bipartition(
    const std::vector<StartingTreeEdge>& edge_candidates,
    const std::vector<double>& scores,
    int active_count,
    int min_side) {
  if (scores.size() != edge_candidates.size()) {
    throw std::logic_error("bipartition candidate/score size mismatch");
  }

  std::vector<int> eligible;
  eligible.reserve(edge_candidates.size());
  for (int i = 0; i < static_cast<int>(edge_candidates.size()); ++i) {
    const StartingTreeEdge& candidate =
        edge_candidates[static_cast<std::size_t>(i)];
    const int left_count = candidate.split.count;
    const int right_count = active_count - left_count;
    if (candidate.rejected || left_count < min_side || right_count < min_side) {
      continue;
    }
    eligible.push_back(i);
  }

  const auto selection = select_max_with_tie_band(
      eligible.begin(), eligible.end(), kScoreTieEps,
      [&](int index) { return scores[static_cast<std::size_t>(index)]; },
      [&](int index) {
        const StartingTreeEdge& candidate =
            edge_candidates[static_cast<std::size_t>(index)];
        const int balance = split_balance_from_count(candidate.split.count,
                                                     active_count);
        return std::tuple{-balance,
                          candidate.representative_order,
                          candidate.edge_child_node,
                          candidate.split.first_begin,
                          candidate.split.first_end,
                          candidate.split.second_begin,
                          candidate.split.second_end,
                          index};
      });
  return selection.selected == eligible.end() ? -1 : *selection.selected;
}

}  // namespace

int select_bipartition_for_tests(
    const std::vector<BipartitionProbe>& probes,
    int active_count,
    int min_side) {
  if (active_count < 0 || min_side < 0) {
    throw std::invalid_argument(
        "bipartition probe sizes must be non-negative");
  }

  std::vector<StartingTreeEdge> candidates;
  std::vector<double> scores;
  candidates.reserve(probes.size());
  scores.reserve(probes.size());
  for (const BipartitionProbe& probe : probes) {
    StartingTreeEdge candidate;
    candidate.split.first_begin = probe.leaf_intervals[0];
    candidate.split.first_end = probe.leaf_intervals[1];
    candidate.split.second_begin = probe.leaf_intervals[2];
    candidate.split.second_end = probe.leaf_intervals[3];
    candidate.split.count = probe.left_count;
    candidate.representative_order = probe.representative_order;
    candidate.edge_child_node = probe.edge_child_node;
    candidates.push_back(candidate);
    scores.push_back(probe.score);
  }
  return select_scored_bipartition(candidates, scores, active_count,
                                              min_side);
}

// The bipartition that divides `taxa`. Starting-tree edges are ranked by K=2
// quartet score (or by balance in the ablation), and the best unrejected edge
// is proposed to the bipartition test, up to `proposals` times. When no edge
// exists or none passes, the result is the trivial bipartition.
BipartitionChoice choose_bipartition(const SolverContext& ctx,
                                     const std::vector<int>& taxa,
                                     uint64_t path,
                                     int depth) {
  (void)depth;
  progress_note_phase(path, "choose bipartition");

  const ActiveTaxa active = make_active_taxa(taxa, context_taxa(ctx).size());
  const int active_count = static_cast<int>(taxa.size());
  const QueryTree query = induce_query_tree(
      context_starting_tree(ctx), active);
  std::vector<StartingTreeEdge> edge_candidates =
      starting_tree_edges(query);
  const bool trivial_completion = edge_candidates.empty();
  if (edge_candidates.empty()) {
    StartingTreeEdge completion_item;
    completion_item.materialized = true;
    completion_item.materialized_candidate =
        trivial_bipartition(active_count);
    completion_item.split.count =
        completion_item.materialized_candidate.split.count();
    completion_item.representative_order =
        completion_item.materialized_candidate.representative_order;
    edge_candidates.push_back(std::move(completion_item));
  }
  stats_add(run_stats().candidate_sets);
  stats_add(run_stats().candidates,
                        static_cast<long long>(edge_candidates.size()));

  if (trivial_completion) {
    BipartitionChoice completion;
    completion.candidate = materialize_candidate(
        query, edge_candidates.front(), active_count);
    completion.score = 0.0;
    completion.representative_order = completion.candidate.representative_order;
    completion.candidate_count = 1;
    completion.left_count = 1;
    completion.right_count = active_count - 1;
    return completion;
  }

  const bool balanced = ctx.options.divide_score != DivideScore::Quartet;
  std::vector<double> scores;
  if (!balanced) {
    stats_add(run_stats().bipartitions_scored,
                          static_cast<long long>(edge_candidates.size()));
    std::vector<int> edge_child_nodes;
    edge_child_nodes.reserve(edge_candidates.size());
    for (const StartingTreeEdge& item : edge_candidates) {
      edge_child_nodes.push_back(item.edge_child_node);
    }
    const auto score_started = std::chrono::steady_clock::now();
    scores = score_bipartitions(ctx, active, query, edge_child_nodes);
    stats_add_time(run_stats().bipartition_score_ns,
                        std::chrono::steady_clock::now() - score_started);
  }

  BipartitionChoice best;
  bool accepted = false;
  bool rejected_any = false;
  // Scores are computed once above; at most 16 proposals are tested,
  // independent of n. Each attempt rescans for the exact maximum, which keeps
  // the tie-band rule without an epsilon sort comparator. An unavailable test
  // accepts the first proposal.
  for (int attempt = 0; attempt < ctx.options.bipartition_test.proposals;
       ++attempt) {
    if (balanced) {
      best = choose_balanced_bipartition(
          ctx, query, edge_candidates, active_count);
    } else {
      best.candidate_index = select_scored_bipartition(
          edge_candidates, scores, active_count,
          effective_min_side(ctx, active_count));
      if (best.candidate_index >= 0) {
        const int i = best.candidate_index;
        best.score = scores[i];
        best.representative_order = edge_candidates[i].representative_order;
        best.candidate = materialize_candidate(
            query, edge_candidates[i], active_count);
      }
    }
    if (best.candidate_index < 0) break;
    const BipartitionTestResult test =
        test_bipartition(ctx, active.taxa, best.candidate.split);
    if (!test.available || test.accepted) {
      accepted = true;
      break;
    }
    edge_candidates[best.candidate_index].rejected = true;
    rejected_any = true;
  }
  best.candidate_count = static_cast<int>(edge_candidates.size());
  if (!accepted) {
    best.candidate = trivial_bipartition(active_count);
    best.candidate_index = -1;
    best.score = 0.0;
    best.representative_order = best.candidate.representative_order;
    if (rejected_any) {
      stats_add(run_stats().bipartition_test_singletons);
    }
  }
  const std::vector<Side> assignment =
      assignment_from_split(best.candidate.split, active_count);
  best.left_count =
      static_cast<int>(std::count(assignment.begin(), assignment.end(),
                                  Side::Left));
  best.right_count = static_cast<int>(assignment.size()) - best.left_count;

  return best;
}

}  // namespace qraft::internal
