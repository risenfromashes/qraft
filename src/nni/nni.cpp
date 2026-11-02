#include "solver/solver.hpp"
#include "core/selection.hpp"
#include "counting/count.hpp"
#include "genes/gene_lca.hpp"
#include "mre/exact_split_ids.hpp"
#include "counting/query_taxa.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>

namespace qraft::internal {
namespace {

// Gains within this of each other are tied, and a move must gain more.
constexpr double kNniGainEps = 1e-9;
// The gates: a move must gain at least kNniMinGainPerGene quartets per gene;
// the support gate (off under --no-nni-support-gate) stops NNI at a selected
// move that lowers the number of genes containing the replaced bipartition by
// more than kNniMaxSupportDropPerGene times the gene count; and at most
// kNniMaxMoves moves are applied (a constant, as the complexity bound
// requires; --no-nni-cap lifts it). Later moves on large trees are mostly
// neutral or wrong, so the cap is small; smaller trees rarely reach it.
constexpr double kNniMinGainPerGene = 5.0;
constexpr double kNniMaxSupportDropPerGene = 0.01;
constexpr int kNniMaxMoves = 20;

// Overflow-checked size arithmetic; each throws std::length_error(what).
std::size_t checked_size_product(std::size_t first,
                                         std::size_t second,
                                         const char* what) {
  if (first != 0 && second > std::numeric_limits<std::size_t>::max() / first) {
    throw std::length_error(what);
  }
  return first * second;
}

int checked_int_size(std::size_t size, const char* what) {
  if (size > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error(what);
  }
  return static_cast<int>(size);
}

std::size_t checked_doubled_capacity(std::size_t size,
                                              bool add_one,
                                              const char* what) {
  const std::size_t extra = add_one ? 1u : 0u;
  if (size > (std::numeric_limits<std::size_t>::max() - extra) / 2u) {
    throw std::length_error(what);
  }
  return size * 2u + extra;
}

long long checked_doubled_counter(std::size_t size,
                                           const char* what) {
  if (size > static_cast<std::size_t>(
                 std::numeric_limits<long long>::max() / 2ll)) {
    throw std::length_error(what);
  }
  return static_cast<long long>(size) * 2ll;
}

// An NNI around the query edge from vertex p to its child q. p's other two
// arms are p0 and p1; q's children are q0, the larger, and q1, the smaller
// (q_split_query). When p is counted, p0, p1 and q0 have the colours in
// p_colors and q_colors[0], and q1 is recoloured to 3. The alternatives swap
// p0 with q0 or with q1.
struct NniEdgeProposal {
  int p_query = -1;
  int q_query = -1;
  int q_split_query = -1;
  std::array<std::uint8_t, 2> p_colors{0, 1};
  std::array<std::uint8_t, 2> q_colors{2, 3};
  NniMove swap_p0_q0;
  NniMove swap_p0_q1;
};

// The key of an undirected tree edge.
uint64_t nni_edge_key(int a, int b) {
  const uint32_t lo = static_cast<uint32_t>(std::min(a, b));
  const uint32_t hi = static_cast<uint32_t>(std::max(a, b));
  return (static_cast<uint64_t>(lo) << 32) | static_cast<uint64_t>(hi);
}

bool nni_gain_meets_threshold(double gain, double min_gain) {
  return gain > kNniGainEps && gain + kNniGainEps >= min_gain;
}

double nni_gene_count(const SolverContext& ctx) {
  const int gene_count = checked_int_size(
      scoring_genes(ctx).size(), "NNI gene count exceeds its integer domain");
  return static_cast<double>(std::max(1, gene_count));
}

double nni_min_gain(const SolverContext& ctx) {
  return std::max(kNniGainEps, kNniMinGainPerGene * nni_gene_count(ctx));
}

std::vector<int> query_children(const QueryTree& query, int node) {
  std::vector<int> children;
  if (node < 0 || node >= static_cast<int>(query.nodes.size())) {
    return children;
  }
  const int first = query.nodes[static_cast<std::size_t>(node)].first_adj;
  const int degree = query.nodes[static_cast<std::size_t>(node)].degree;
  children.reserve(static_cast<std::size_t>(degree));
  for (int i = 0; i < degree; ++i) {
    const int nb = query.adjacency[static_cast<std::size_t>(first + i)];
    if (query.parent[static_cast<std::size_t>(nb)] == node) {
      children.push_back(nb);
    }
  }
  return children;
}

bool query_internal_binary(const QueryTree& query, int node) {
  return node >= 0 && node < static_cast<int>(query.nodes.size()) &&
         query.nodes[static_cast<std::size_t>(node)].taxon_active < 0 &&
         query.nodes[static_cast<std::size_t>(node)].degree == 3;
}

// The channel of the quartet topology ab|cd, for distinct colours a, b, c, d.
std::array<int, 3> sorted_pair_pair(int a, int b, int c, int d) {
  std::array<int, 2> left{std::min(a, b), std::max(a, b)};
  std::array<int, 2> right{std::min(c, d), std::max(c, d)};
  if (right < left) std::swap(left, right);
  return std::array<int, 3>{left[0], left[1], right[1]};
}

double pair_pair_score(const ChannelCounts& counts,
                       int a,
                       int b,
                       int c,
                       int d) {
  const std::array<int, 3> key = sorted_pair_pair(a, b, c, d);
  if (key == std::array<int, 3>{0, 1, 3}) return counts.q01_23;
  if (key == std::array<int, 3>{0, 2, 3}) return counts.q02_13;
  if (key == std::array<int, 3>{0, 3, 2}) return counts.q03_12;
  return 0.0;
}

int query_source_node(const QueryTree& query, int node) {
  if (node < 0 || node >= static_cast<int>(query.nodes.size())) return -1;
  return query.nodes[static_cast<std::size_t>(node)].source_node;
}

// One proposal per query edge whose ends are both binary internal vertices
// adjacent in `tree`, restricted to the sorted central edge keys
// `allowed_central_edges` when given.
std::vector<NniEdgeProposal> collect_nni_proposals(
    const SolverTree& tree,
    const QueryTree& query,
    const std::vector<uint64_t>* allowed_central_edges = nullptr) {
  std::vector<NniEdgeProposal> proposals;
  if (query.root < 0 || query.nodes.empty()) return proposals;
  for (int q : query.preorder) {
    if (q == query.root || !query_internal_binary(query, q)) continue;
    const int p = query.parent[static_cast<std::size_t>(q)];
    if (!query_internal_binary(query, p)) continue;

    const int edge_a = query_source_node(query, p);
    const int edge_b = query_source_node(query, q);
    if (edge_a < 0 || edge_b < 0 || !has_solver_edge(tree, edge_a, edge_b)) {
      continue;
    }
    if (allowed_central_edges != nullptr &&
        !std::binary_search(allowed_central_edges->begin(),
                            allowed_central_edges->end(),
                            nni_edge_key(edge_a, edge_b))) {
      continue;
    }

    const VertexArmColoring p_coloring =
        vertex_arm_coloring(query, p);
    const std::array<int, 3> p_branches = p_coloring.branch_nodes;
    const std::array<std::uint8_t, 3> p_branch_colors =
        p_coloring.branch_colors;
    int q_branch_index = -1;
    for (int i = 0; i < 3; ++i) {
      if (p_branches[static_cast<std::size_t>(i)] == q) {
        q_branch_index = i;
        break;
      }
    }
    if (q_branch_index < 0) continue;

    std::array<int, 2> p_side_query{-1, -1};
    std::array<std::uint8_t, 2> p_colors{0, 1};
    int p_side_count = 0;
    for (int i = 0; i < 3; ++i) {
      if (i == q_branch_index) continue;
      p_side_query[static_cast<std::size_t>(p_side_count)] =
          p_branches[static_cast<std::size_t>(i)];
      p_colors[static_cast<std::size_t>(p_side_count)] =
          p_branch_colors[static_cast<std::size_t>(i)];
      ++p_side_count;
    }
    if (p_side_count != 2) continue;

    const std::vector<int> q_children = query_children(query, q);
    if (q_children.size() != 2) continue;
    const int q0 = q_children[0];
    const int q1 = q_children[1];
    const int q0_size = query_leaf_count(query, q0);
    const int q1_size = query_leaf_count(query, q1);
    const int q_split = q0_size <= q1_size ? q0 : q1;
    const int q_base = q_split == q0 ? q1 : q0;

    const int p0_source = query_source_node(query, p_side_query[0]);
    const int p1_source = query_source_node(query, p_side_query[1]);
    const int q0_source = query_source_node(query, q_base);
    const int q1_source = query_source_node(query, q_split);
    if (p0_source < 0 || p1_source < 0 || q0_source < 0 || q1_source < 0 ||
        !has_solver_edge(tree, edge_a, p0_source) ||
        !has_solver_edge(tree, edge_a, p1_source) ||
        !has_solver_edge(tree, edge_b, q0_source) ||
        !has_solver_edge(tree, edge_b, q1_source)) {
      continue;
    }
    NniEdgeProposal proposal;
    proposal.p_query = p;
    proposal.q_query = q;
    proposal.q_split_query = q_split;
    proposal.p_colors = p_colors;
    proposal.q_colors = std::array<std::uint8_t, 2>{
        p_branch_colors[static_cast<std::size_t>(q_branch_index)], 3};

    proposal.swap_p0_q0.edge_a = edge_a;
    proposal.swap_p0_q0.edge_b = edge_b;
    proposal.swap_p0_q0.swap_a_neighbor = p0_source;
    proposal.swap_p0_q0.swap_b_neighbor = q0_source;
    proposal.swap_p0_q0.alternative = 1;

    proposal.swap_p0_q1.edge_a = edge_a;
    proposal.swap_p0_q1.edge_b = edge_b;
    proposal.swap_p0_q1.swap_a_neighbor = p0_source;
    proposal.swap_p0_q1.swap_b_neighbor = q1_source;
    proposal.swap_p0_q1.alternative = 2;

    proposals.push_back(proposal);
  }
  return proposals;
}

// Quartet support of the current topology around a central edge and of the
// two alternatives, for one gene or summed over genes.
struct NniProposalScores {
  double current = 0.0;
  double swap_p0_q0 = 0.0;
  double swap_p0_q1 = 0.0;
};

struct ScoredNniData {
  std::vector<NniEdgeProposal> proposals;
  std::vector<NniProposalScores> scores;
};

// The positions of q1's leaves in the query's leaf owner order.
struct NniOwnerRange {
  int begin = 0;
  int end = 0;
};

// Scores every proposal with K=4 counting, one score row per gene summed in
// gene order. With `per_gene_normalized`, each gene's three scores of a
// proposal are divided by their sum (when positive) before the sum over
// genes, as ASTER counts them for branch support (ASTRAL-III too on binary
// gene trees).
std::vector<NniProposalScores> score_nni_proposals_hdt(
    const SolverContext& ctx,
    const QueryTree& query,
    const std::vector<NniEdgeProposal>& proposals,
    const std::vector<int>& active_taxa,
    int taxon_count,
    bool per_gene_normalized = false) {
  std::vector<NniProposalScores> totals(proposals.size());
  if (proposals.empty()) return totals;
  const int proposal_count = checked_int_size(
      proposals.size(), "NNI proposal count exceeds its integer domain");
  (void)checked_int_size(
      query.nodes.size(), "NNI query exceeds its integer node domain");
  (void)checked_int_size(
      active_taxa.size(), "NNI active-taxon set exceeds its integer domain");

  std::vector<int> proposal_by_q(query.nodes.size(), -1);
  for (int proposal_index = 0; proposal_index < proposal_count;
       ++proposal_index) {
    const int q = proposals[static_cast<std::size_t>(proposal_index)].q_query;
    if (q >= 0 && q < static_cast<int>(proposal_by_q.size())) {
      proposal_by_q[static_cast<std::size_t>(q)] = proposal_index;
    }
  }

  std::vector<int> query_active_to_owner(query.leaf_count, -1);
  for (int active = 0; active < query.leaf_count; ++active) {
    if (active < static_cast<int>(query.active_to_node.size()) &&
        query.active_to_node[static_cast<std::size_t>(active)] >= 0 &&
        active < static_cast<int>(active_taxa.size())) {
      query_active_to_owner[static_cast<std::size_t>(active)] =
          active_taxa[static_cast<std::size_t>(active)];
    }
  }
  const std::vector<int> leaf_owner_order =
      query_leaf_owner_order(query, query_active_to_owner);
  std::vector<NniOwnerRange> q_split_ranges(proposals.size());
  for (int proposal_index = 0; proposal_index < proposal_count;
       ++proposal_index) {
    const int q_split =
        proposals[static_cast<std::size_t>(proposal_index)].q_split_query;
    if (q_split < 0 ||
        q_split >= static_cast<int>(query.subtree_leaf_first.size()) ||
        query.leaf_active_order.empty()) {
      continue;
    }
    const int raw_begin =
        query.subtree_leaf_first[static_cast<std::size_t>(q_split)];
    const int raw_end =
        raw_begin + query.subtree_leaf_count[static_cast<std::size_t>(q_split)];
    const int begin =
        std::clamp(raw_begin, 0,
                   static_cast<int>(leaf_owner_order.size()));
    const int end =
        std::clamp(raw_end, begin,
                   static_cast<int>(leaf_owner_order.size()));
    q_split_ranges[static_cast<std::size_t>(proposal_index)] =
        NniOwnerRange{begin, end};
  }

  // NNI scores the input genes, so every taxon is its own owner.
  OwnerProjection owners;
  owners.taxon_to_owner.resize(static_cast<std::size_t>(taxon_count));
  std::iota(owners.taxon_to_owner.begin(), owners.taxon_to_owner.end(), 0);
  owners.owner_count = taxon_count;
  owners.live_owner.assign(static_cast<std::size_t>(taxon_count), 0);
  owners.live_taxa.reserve(active_taxa.size());
  for (int taxon : active_taxa) {
    if (taxon >= 0 && taxon < taxon_count) {
      owners.live_owner[static_cast<std::size_t>(taxon)] = 1;
      owners.live_taxa.push_back(taxon);
    }
  }

  const std::vector<GeneIndex>& genes = scoring_genes(ctx);
  const int gene_count = checked_int_size(
      genes.size(), "NNI gene count exceeds its integer domain");
  const int chunks = gene_partition_count(gene_count);
  const std::size_t partial_count = checked_size_product(
      static_cast<std::size_t>(chunks), proposals.size(),
      "NNI partial score table exceeds the addressable size domain");
  RowTable<NniProposalScores> partial(static_cast<std::size_t>(chunks),
                                      proposals.size());
  const std::size_t partial_bytes = checked_size_product(
      partial_count, sizeof(NniProposalScores),
      "NNI partial score bytes exceed the addressable size domain");
  if (partial_bytes >
      static_cast<std::size_t>(std::numeric_limits<long long>::max())) {
    throw std::length_error(
        "NNI partial score bytes exceed the telemetry domain");
  }
  std::vector<GeneCountTelemetry> telemetry =
      gene_count_telemetry_slots(chunks);

  // Only vertices p of a proposal are read, and the traversal skips query
  // subtrees without one.
  std::vector<std::uint8_t> requested_vertices(query.nodes.size(), 0);
  for (int q = 0; q < static_cast<int>(proposal_by_q.size()); ++q) {
    if (proposal_by_q[static_cast<std::size_t>(q)] < 0) continue;
    const int vertex = query.parent[static_cast<std::size_t>(q)];
    if (vertex >= 0 &&
        static_cast<std::size_t>(vertex) < requested_vertices.size()) {
      requested_vertices[static_cast<std::size_t>(vertex)] = 1;
    }
  }

  // Genes differ widely in size, so hand them out one at a time; rows stay
  // one per gene and are reduced in gene order below.
  parallel_for_dynamic(0, chunks, [&](int chunk) {
    const int begin = chunk;
    const int end = begin + 1;
    NniProposalScores* row =
        partial.init_row(static_cast<std::size_t>(chunk));
    long long queries = 0;
    HdtTraversalStats hdt_stats;
    HdtExtractWorkspace hdt_workspace;
    // The workspace is local to the chunk, so worker threads do not keep
    // root-sized buffers between passes.
    for (int g = begin; g < end; ++g) {
      const GeneIndex& gene = genes[static_cast<std::size_t>(g)];
      if (!gene_has_owner_quartet(gene, owners)) continue;
      HdtColorPolicy colors;
      colors.default_current = 0;
      colors.default_background = 0;
      colors.current_by_owner.assign(
          static_cast<std::size_t>(owners.owner_count), 0);
      colors.background_by_owner.assign(
          static_cast<std::size_t>(owners.owner_count), 0);
      for (int owner : leaf_owner_order) {
        if (owner < 0 || owner >= owners.owner_count) continue;
        colors.current_by_owner[static_cast<std::size_t>(owner)] = 1;
      }
      HdtInitialStateResult initial = hdt_build_initial_state(
          build_base_component_graph(gene, owners, 4, colors));
      hdt_account_initial_state(initial, hdt_stats);
      count_vertices(
          std::move(initial.state), query, leaf_owner_order,
          hdt_workspace, hdt_stats,
          [&](int vertex, const std::array<int, 3>& branch_nodes,
              const std::array<std::uint8_t, 3>&,
              HdtState& current) {
            for (int branch_node : branch_nodes) {
              if (branch_node < 0 ||
                  branch_node >= static_cast<int>(proposal_by_q.size()) ||
                  query.parent[static_cast<std::size_t>(branch_node)] !=
                      vertex) {
                continue;
              }
              // Recolour q1 to 3, read the three topologies of p0, p1, q0
              // and q1, then restore q1 to q's arm colour.
              const int proposal_index =
                  proposal_by_q[static_cast<std::size_t>(branch_node)];
              if (proposal_index < 0) continue;
              const NniEdgeProposal& proposal =
                  proposals[static_cast<std::size_t>(proposal_index)];
              const NniOwnerRange range =
                  q_split_ranges[static_cast<std::size_t>(proposal_index)];
              hdt_recolor_owner_range(
                  current, leaf_owner_order, range.begin, range.end, 3);
              const std::array<double, 3> values =
                  hdt_root_counts(current);
              const ChannelCounts counts{values[0], values[1], values[2]};
              NniProposalScores& out =
                  row[static_cast<std::size_t>(proposal_index)];
              const int p0 = proposal.p_colors[0];
              const int p1 = proposal.p_colors[1];
              const int q0 = proposal.q_colors[0];
              const int q1 = proposal.q_colors[1];
              out.current += pair_pair_score(counts, p0, p1, q0, q1);
              out.swap_p0_q0 += pair_pair_score(counts, p1, q0, p0, q1);
              out.swap_p0_q1 += pair_pair_score(counts, p1, q1, p0, q0);
              hdt_recolor_owner_range(
                  current, leaf_owner_order, range.begin, range.end,
                  static_cast<std::uint8_t>(proposal.q_colors[0]));
              ++queries;
            }
          },
          4, &requested_vertices);
    }
    if (per_gene_normalized) {
      // One chunk holds one gene, so the row is that gene's scores.
      for (std::size_t j = 0; j < proposals.size(); ++j) {
        NniProposalScores& out = row[j];
        const double sum = out.current + out.swap_p0_q0 + out.swap_p0_q1;
        if (sum > 0.0) {
          out.current /= sum;
          out.swap_p0_q0 /= sum;
          out.swap_p0_q1 /= sum;
        }
      }
    }
    if (!telemetry.empty()) {
      telemetry[static_cast<std::size_t>(chunk)] = GeneCountTelemetry{
          queries, hdt_stats.recolor_changed_leaves,
          hdt_stats.recolor_marked_nodes, hdt_stats};
    }
  }, 1);

  add_rows_in_order(partial.data(), chunks, proposals.size(), totals.data(),
                    [](NniProposalScores& out, const NniProposalScores& in) {
                      out.current += in.current;
                      out.swap_p0_q0 += in.swap_p0_q0;
                      out.swap_p0_q1 += in.swap_p0_q1;
                    });
  merge_gene_count_telemetry(telemetry);
  stats_add(run_stats().hdt_builds,
                        static_cast<long long>(gene_count));
  return totals;
}

// The owners of a local pass: the query minus the proposals' central
// vertices falls into components, and every taxon of one component has the
// same colour in every colouring the pass reads, since each proposal's four
// arms are unions of components. owner_of_taxon maps a taxon to its
// component's owner, numbered by first leaf in query preorder, or -1.
struct NniLocalOwners {
  std::vector<int> owner_of_taxon;
  // Query leaf position of each owner's first leaf.
  std::vector<int> leaf_position;
  int owner_count = 0;
};

NniLocalOwners nni_local_owners(const QueryTree& query,
                                const std::vector<NniEdgeProposal>& proposals,
                                const std::vector<int>& active_taxa,
                                int taxon_count) {
  const std::size_t node_count = query.nodes.size();
  std::vector<std::uint8_t> central(node_count, 0);
  for (const NniEdgeProposal& proposal : proposals) {
    for (int node : {proposal.p_query, proposal.q_query}) {
      if (node >= 0 && static_cast<std::size_t>(node) < node_count) {
        central[static_cast<std::size_t>(node)] = 1;
      }
    }
  }
  // Preorder visits a parent first, so a non-central node joins its
  // parent's component unless the parent is central. The root is its own
  // parent in a query.
  std::vector<int> component(node_count, -1);
  std::vector<int> owner_of_component;
  NniLocalOwners owners;
  owners.owner_of_taxon.assign(static_cast<std::size_t>(taxon_count), -1);
  for (int node : query.preorder) {
    if (central[static_cast<std::size_t>(node)] != 0) continue;
    const int parent =
        node == query.root ? -1 : query.parent[static_cast<std::size_t>(node)];
    if (parent >= 0 && central[static_cast<std::size_t>(parent)] == 0) {
      component[static_cast<std::size_t>(node)] =
          component[static_cast<std::size_t>(parent)];
    } else {
      component[static_cast<std::size_t>(node)] =
          static_cast<int>(owner_of_component.size());
      owner_of_component.push_back(-1);
    }
    const int active = query.nodes[static_cast<std::size_t>(node)].taxon_active;
    if (active < 0 || active >= static_cast<int>(active_taxa.size())) continue;
    const int taxon = active_taxa[static_cast<std::size_t>(active)];
    if (taxon < 0 || taxon >= taxon_count) continue;
    int& owner = owner_of_component[static_cast<std::size_t>(
        component[static_cast<std::size_t>(node)])];
    if (owner < 0) {
      owner = owners.owner_count++;
      owners.leaf_position.push_back(
          query.subtree_leaf_first[static_cast<std::size_t>(node)]);
    }
    owners.owner_of_taxon[static_cast<std::size_t>(taxon)] = owner;
  }
  return owners;
}

// The colour of every owner at the read of `proposal`: its arm's colour
// from vertex_arm_coloring at p, or 3 inside q1, as count_vertices colours
// them.
std::vector<std::uint8_t> nni_local_owner_colors(
    const QueryTree& query,
    const NniEdgeProposal& proposal,
    const NniLocalOwners& owners) {
  const auto inside = [&](int node, int position) {
    const int first = query.subtree_leaf_first[static_cast<std::size_t>(node)];
    return position >= first &&
           position <
               first + query.subtree_leaf_count[static_cast<std::size_t>(node)];
  };
  const int p = proposal.p_query;
  const VertexArmColoring arms = vertex_arm_coloring(query, p);
  const int parent = query.parent[static_cast<std::size_t>(p)];
  std::vector<std::uint8_t> colors(static_cast<std::size_t>(owners.owner_count));
  for (int owner = 0; owner < owners.owner_count; ++owner) {
    const int position = owners.leaf_position[static_cast<std::size_t>(owner)];
    if (inside(proposal.q_split_query, position)) {
      colors[static_cast<std::size_t>(owner)] = 3;
      continue;
    }
    int arm = -1;
    for (int i = 0; i < 3; ++i) {
      const int branch = arms.branch_nodes[static_cast<std::size_t>(i)];
      if (branch < 0) continue;
      const bool member = p != query.root && branch == parent
                              ? !inside(p, position)
                              : inside(branch, position);
      if (member) {
        arm = i;
        break;
      }
    }
    if (arm < 0) {
      throw std::logic_error("NNI local owner lies in no arm of its proposal");
    }
    colors[static_cast<std::size_t>(owner)] =
        arms.branch_colors[static_cast<std::size_t>(arm)];
  }
  return colors;
}

// `gene` with every maximal subtree whose leaves all have one owner replaced
// by one leaf of that owner, weighted by the subtree's leaf mass. Counting
// reads such a subtree only through its mass: a quartet takes at most one
// leaf of a colour, and every leaf of the subtree meets the rest of the gene
// at the subtree's root. Retained vertices keep their preorder, which is the
// new numbering, and their incidence order. Returns false when every leaf
// has the same owner, so that no quartet has four colours.
bool collapse_gene_to_owners(const GeneIndex& gene,
                             const std::vector<int>& owner_of_taxon,
                             GeneIndex& out) {
  constexpr int kEmpty = -3;
  constexpr int kMixed = -2;
  const std::vector<int>& parent = gene.root_traversal.parent;
  const std::vector<int>& preorder = gene.root_traversal.preorder;
  const std::size_t node_count = gene.nodes.size();
  if (node_count == 0 || preorder.size() != node_count) return false;
  std::vector<int> subtree_owner(node_count, kEmpty);
  std::vector<double> mass(node_count, 0.0);
  std::vector<int> first_leaf(node_count, -1);
  for (auto it = preorder.rbegin(); it != preorder.rend(); ++it) {
    const int v = *it;
    const GeneIndex::Node& node = gene.nodes[static_cast<std::size_t>(v)];
    int& owner = subtree_owner[static_cast<std::size_t>(v)];
    if (node.key != -1) {
      const int taxon = node.taxon();
      const int leaf_owner =
          taxon >= 0 && taxon < static_cast<int>(owner_of_taxon.size())
              ? owner_of_taxon[static_cast<std::size_t>(taxon)]
              : -1;
      // A leaf without an owner is never merged.
      owner = leaf_owner >= 0 ? leaf_owner : kMixed;
      mass[static_cast<std::size_t>(v)] = gene_leaf_mass(gene, v);
      first_leaf[static_cast<std::size_t>(v)] = v;
    }
    for (int child : gene_neighbors(gene, v)) {
      if (parent[static_cast<std::size_t>(child)] != v) continue;
      const int child_owner = subtree_owner[static_cast<std::size_t>(child)];
      owner = owner == kEmpty ? child_owner
                              : (owner == child_owner ? owner : kMixed);
      mass[static_cast<std::size_t>(v)] += mass[static_cast<std::size_t>(child)];
      if (first_leaf[static_cast<std::size_t>(v)] < 0) {
        first_leaf[static_cast<std::size_t>(v)] =
            first_leaf[static_cast<std::size_t>(child)];
      }
    }
  }
  if (subtree_owner[static_cast<std::size_t>(gene.root)] >= 0) return false;

  // A vertex is kept when its parent is kept and mixed; a kept vertex with a
  // single owner becomes a leaf.
  std::vector<int> new_id(node_count, -1);
  std::vector<int> kept;
  kept.reserve(node_count);
  for (int v : preorder) {
    const int up = parent[static_cast<std::size_t>(v)];
    const bool keep =
        v == gene.root ||
        (up >= 0 && new_id[static_cast<std::size_t>(up)] >= 0 &&
         subtree_owner[static_cast<std::size_t>(up)] == kMixed);
    if (!keep) continue;
    new_id[static_cast<std::size_t>(v)] = static_cast<int>(kept.size());
    kept.push_back(v);
  }
  const int kept_count = static_cast<int>(kept.size());
  out = GeneIndex{};
  out.nodes.resize(static_cast<std::size_t>(kept_count));
  out.root = 0;
  out.reference_root = 0;
  out.root_traversal.parent.assign(static_cast<std::size_t>(kept_count), -1);
  out.root_traversal.preorder.resize(static_cast<std::size_t>(kept_count));
  out.adjacency_offsets.assign(static_cast<std::size_t>(kept_count) + 1u, 0);
  out.adjacency.reserve(static_cast<std::size_t>(kept_count) * 2u);
  std::vector<double> kept_mass(static_cast<std::size_t>(kept_count), 1.0);
  bool unit_masses = true;
  for (int id = 0; id < kept_count; ++id) {
    const int v = kept[static_cast<std::size_t>(id)];
    const bool collapsed =
        v != gene.root && subtree_owner[static_cast<std::size_t>(v)] >= 0;
    out.root_traversal.preorder[static_cast<std::size_t>(id)] = id;
    const int up = parent[static_cast<std::size_t>(v)];
    if (up >= 0) {
      out.root_traversal.parent[static_cast<std::size_t>(id)] =
          new_id[static_cast<std::size_t>(up)];
    }
    if (collapsed) {
      out.nodes[static_cast<std::size_t>(id)] =
          gene.nodes[static_cast<std::size_t>(
              first_leaf[static_cast<std::size_t>(v)])];
      kept_mass[static_cast<std::size_t>(id)] =
          mass[static_cast<std::size_t>(v)];
      out.adjacency.push_back(new_id[static_cast<std::size_t>(up)]);
    } else {
      out.nodes[static_cast<std::size_t>(id)] =
          gene.nodes[static_cast<std::size_t>(v)];
      if (gene.nodes[static_cast<std::size_t>(v)].key != -1) {
        kept_mass[static_cast<std::size_t>(id)] = gene_leaf_mass(gene, v);
      }
      for (int neighbor : gene_neighbors(gene, v)) {
        out.adjacency.push_back(new_id[static_cast<std::size_t>(neighbor)]);
      }
    }
    unit_masses = unit_masses && kept_mass[static_cast<std::size_t>(id)] == 1.0;
    out.adjacency_offsets[static_cast<std::size_t>(id) + 1u] =
        static_cast<int>(out.adjacency.size());
  }
  if (!unit_masses) out.leaf_mass = std::move(kept_mass);
  return true;
}

// Scores the proposals of a local pass, which share few central vertices,
// on genes collapsed to their owners (nni_local_owners). Each gene's state
// is built once and recoloured per proposal, so the pass costs a scan of
// each gene plus counting on its collapsed size. The scores equal
// score_nni_proposals_hdt's.
std::vector<NniProposalScores> score_nni_proposals_local(
    const SolverContext& ctx,
    const QueryTree& query,
    const std::vector<NniEdgeProposal>& proposals,
    const std::vector<int>& active_taxa,
    int taxon_count) {
  std::vector<NniProposalScores> totals(proposals.size());
  if (proposals.empty()) return totals;
  const NniLocalOwners local =
      nni_local_owners(query, proposals, active_taxa, taxon_count);
  std::vector<std::vector<std::uint8_t>> colors;
  colors.reserve(proposals.size());
  for (const NniEdgeProposal& proposal : proposals) {
    colors.push_back(nni_local_owner_colors(query, proposal, local));
  }
  OwnerProjection owners;
  owners.taxon_to_owner = local.owner_of_taxon;
  owners.owner_count = local.owner_count;

  const std::vector<GeneIndex>& genes = scoring_genes(ctx);
  const int gene_count = checked_int_size(
      genes.size(), "NNI gene count exceeds its integer domain");
  const int chunks = gene_partition_count(gene_count);
  RowTable<NniProposalScores> partial(static_cast<std::size_t>(chunks),
                                      proposals.size());
  std::vector<GeneCountTelemetry> telemetry =
      gene_count_telemetry_slots(chunks);
  parallel_for_dynamic(0, chunks, [&](int g) {
    NniProposalScores* row = partial.init_row(static_cast<std::size_t>(g));
    long long queries = 0;
    HdtTraversalStats hdt_stats;
    GeneIndex collapsed;
    if (collapse_gene_to_owners(genes[static_cast<std::size_t>(g)],
                                local.owner_of_taxon, collapsed) &&
        gene_has_owner_quartet(collapsed, owners)) {
      HdtColorPolicy policy;
      policy.default_current = 0;
      policy.default_background = 0;
      policy.current_by_owner = colors.front();
      policy.background_by_owner.assign(
          static_cast<std::size_t>(owners.owner_count), 0);
      HdtInitialStateResult initial = hdt_build_initial_state(
          build_base_component_graph(collapsed, owners, 4, policy));
      hdt_account_initial_state(initial, hdt_stats);
      HdtState& state = initial.state;
      std::array<std::vector<int>, 4> recolor;
      for (std::size_t j = 0; j < proposals.size(); ++j) {
        if (j > 0) {
          for (std::vector<int>& list : recolor) list.clear();
          for (int owner = 0; owner < owners.owner_count; ++owner) {
            const std::uint8_t color = colors[j][static_cast<std::size_t>(owner)];
            if (color != colors[j - 1][static_cast<std::size_t>(owner)]) {
              recolor[color].push_back(owner);
            }
          }
          for (std::uint8_t color = 0; color < 4; ++color) {
            if (!recolor[color].empty()) {
              hdt_recolor_owners(state, recolor[color], color);
            }
          }
        }
        const std::array<double, 3> values = hdt_root_counts(state);
        const ChannelCounts counts{values[0], values[1], values[2]};
        const NniEdgeProposal& proposal = proposals[j];
        const int p0 = proposal.p_colors[0];
        const int p1 = proposal.p_colors[1];
        const int q0 = proposal.q_colors[0];
        const int q1 = proposal.q_colors[1];
        row[j].current += pair_pair_score(counts, p0, p1, q0, q1);
        row[j].swap_p0_q0 += pair_pair_score(counts, p1, q0, p0, q1);
        row[j].swap_p0_q1 += pair_pair_score(counts, p1, q1, p0, q0);
        ++queries;
      }
    }
    if (!telemetry.empty()) {
      telemetry[static_cast<std::size_t>(g)] = GeneCountTelemetry{
          queries, hdt_stats.recolor_changed_leaves,
          hdt_stats.recolor_marked_nodes, hdt_stats};
    }
  }, 1);

  add_rows_in_order(partial.data(), chunks, proposals.size(), totals.data(),
                    [](NniProposalScores& out, const NniProposalScores& in) {
                      out.current += in.current;
                      out.swap_p0_q0 += in.swap_p0_q0;
                      out.swap_p0_q1 += in.swap_p0_q1;
                    });
  merge_gene_count_telemetry(telemetry);
  stats_add(run_stats().hdt_builds, static_cast<long long>(gene_count));
  return totals;
}

// Swaps swap_a_neighbor of edge_a with swap_b_neighbor of edge_b; false when
// the move does not apply to `tree`.
bool apply_nni_move(SolverTree& tree, const NniMove& move) {
  if (!has_solver_edge(tree, move.edge_a, move.edge_b) ||
      !has_solver_edge(tree, move.edge_a, move.swap_a_neighbor) ||
      !has_solver_edge(tree, move.edge_b, move.swap_b_neighbor)) {
    return false;
  }
  remove_solver_edge(tree, move.edge_a, move.swap_a_neighbor);
  remove_solver_edge(tree, move.edge_b, move.swap_b_neighbor);
  add_solver_edge(tree, move.edge_a, move.swap_b_neighbor);
  add_solver_edge(tree, move.edge_b, move.swap_a_neighbor);
  return true;
}

// Scores the proposals of `tree`, or only those whose central edge is in
// `allowed_central_edges`.
ScoredNniData score_nni_edges(
    const SolverContext& ctx,
    const SolverTree& tree,
    const std::vector<int>& active_taxa,
    const std::vector<uint64_t>* allowed_central_edges = nullptr) {
  const int taxon_count = static_cast<int>(context_taxa(ctx).size());
  const QueryTree query =
      query_tree_from_solver_tree(tree, active_taxa, taxon_count);
  ScoredNniData data;
  data.proposals = collect_nni_proposals(tree, query, allowed_central_edges);
  const long long scored_moves = checked_doubled_counter(
      data.proposals.size(), "NNI proposal counter exceeds its integer domain");
  stats_add(run_stats().nni_edges_scored,
                        scored_moves);
  if (allowed_central_edges == nullptr) {
    stats_add(run_stats().nni_full_score_passes);
  } else {
    stats_add(run_stats().nni_local_score_passes);
    stats_add(run_stats().nni_local_edges_scored,
                          scored_moves);
  }
  // A local pass scores a few edges around one move, so it counts on genes
  // collapsed to the owners of their neighbourhood.
  data.scores =
      allowed_central_edges == nullptr
          ? score_nni_proposals_hdt(ctx, query, data.proposals, active_taxa,
                                    taxon_count)
          : score_nni_proposals_local(ctx, query, data.proposals, active_taxa,
                                      taxon_count);
  return data;
}

struct CachedNniMove {
  NniMove move;
  uint64_t central_edge_key = 0;
};

// Caches every alternative whose gain passes the gain gate.
void add_cached_positive_moves(const ScoredNniData& scored,
                               std::vector<CachedNniMove>& cache,
                               double min_gain) {
  for (int i = 0; i < static_cast<int>(scored.proposals.size()); ++i) {
    const NniEdgeProposal& proposal =
        scored.proposals[static_cast<std::size_t>(i)];
    const NniProposalScores& scores =
        scored.scores[static_cast<std::size_t>(i)];
    const uint64_t central_key =
        nni_edge_key(proposal.swap_p0_q0.edge_a,
                              proposal.swap_p0_q0.edge_b);
    const double gain0 = scores.swap_p0_q0 - scores.current;
    if (nni_gain_meets_threshold(gain0, min_gain)) {
      NniMove move = proposal.swap_p0_q0;
      move.gain = gain0;
      cache.push_back(CachedNniMove{move, central_key});
    }
    const double gain1 = scores.swap_p0_q1 - scores.current;
    if (nni_gain_meets_threshold(gain1, min_gain)) {
      NniMove move = proposal.swap_p0_q1;
      move.gain = gain1;
      cache.push_back(CachedNniMove{move, central_key});
    }
  }
}

bool nni_move_applicable(const SolverTree& tree, const NniMove& move) {
  return has_solver_edge(tree, move.edge_a, move.edge_b) &&
         has_solver_edge(tree, move.edge_a, move.swap_a_neighbor) &&
         has_solver_edge(tree, move.edge_b, move.swap_b_neighbor);
}

// The tie-break key among moves of equal gain.
std::array<int, 4> nni_move_stable_key(const NniMove& move) {
  return std::array<int, 4>{move.edge_a, move.edge_b, move.swap_a_neighbor,
                            move.swap_b_neighbor};
}

// The applicable cached move of largest gain; the stable key picks among
// the moves within kNniGainEps of it.
NniMove best_cached_nni_move(const SolverTree& tree,
                                const std::vector<CachedNniMove>& cache) {
  const auto selection = select_max_with_tie_band(
      cache.begin(), cache.end(), kNniGainEps,
      [](const CachedNniMove& cached) { return cached.move.gain; },
      [](const CachedNniMove& cached) {
        return nni_move_stable_key(cached.move);
      },
      [&](const CachedNniMove& cached) {
        return cached.move.gain > kNniGainEps &&
               nni_move_applicable(tree, cached.move);
      });
  return selection.selected == cache.end() ? NniMove{}
                                            : selection.selected->move;
}

// The distinct nodes a move reconnects.
std::vector<int> nni_touched_nodes(const NniMove& move) {
  std::vector<int> seeds;
  seeds.reserve(4);
  const std::array<int, 4> nodes{move.edge_a, move.edge_b,
                                 move.swap_a_neighbor,
                                 move.swap_b_neighbor};
  for (int node : nodes) {
    if (node < 0) continue;
    if (std::find(seeds.begin(), seeds.end(), node) == seeds.end()) {
      seeds.push_back(node);
    }
  }
  return seeds;
}

bool unordered_edge_equal(int a, int b, int c, int d) {
  return (a == c && b == d) || (a == d && b == c);
}

// The taxa reachable from `start` without entering `blocked`, following the
// adjacency that `neighbors` reports, which may describe the tree after a
// move.
template <typename NeighborFn>
std::vector<int> component_taxa_with_neighbors(const SolverTree& tree,
                                               int start,
                                               int blocked,
                                               int taxon_count,
                                               NeighborFn&& neighbors) {
  std::vector<int> taxa;
  if (start < 0 || start >= static_cast<int>(tree.nodes.size())) return taxa;
  std::vector<uint8_t> seen(tree.nodes.size(), 0);
  std::vector<int> stack;
  std::vector<int> nbs;
  stack.push_back(start);
  if (blocked >= 0 && blocked < static_cast<int>(seen.size())) {
    seen[static_cast<std::size_t>(blocked)] = 1;
  }
  while (!stack.empty()) {
    const int node = stack.back();
    stack.pop_back();
    if (node < 0 || node >= static_cast<int>(tree.nodes.size()) ||
        seen[static_cast<std::size_t>(node)] ||
        !tree.nodes[static_cast<std::size_t>(node)].active) {
      continue;
    }
    seen[static_cast<std::size_t>(node)] = 1;
    const int taxon = tree.nodes[static_cast<std::size_t>(node)].taxon;
    if (taxon >= 0 && taxon < taxon_count) taxa.push_back(taxon);
    nbs.clear();
    nbs.reserve(tree.nodes[static_cast<std::size_t>(node)].neighbors.size() + 2);
    neighbors(node, nbs);
    for (int nb : nbs) {
      if (nb >= 0 && nb < static_cast<int>(seen.size()) &&
          !seen[static_cast<std::size_t>(nb)]) {
        stack.push_back(nb);
      }
    }
  }
  return taxa;
}

// The taxa of the b-side of the central edge (a, b) of `tree`.
std::vector<int> old_central_side_taxa(const SolverTree& tree,
                                       int a,
                                       int b,
                                       int taxon_count) {
  const auto neighbors = [&](int node, std::vector<int>& out) {
    for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
      out.push_back(nb);
    }
  };
  return component_taxa_with_neighbors(tree, b, a, taxon_count, neighbors);
}

// The taxa of the edge_b side of the central edge once `move` is applied,
// read through the moved adjacency without editing `tree`.
std::vector<int> new_central_side_taxa(const SolverTree& tree,
                                       const NniMove& move,
                                       int taxon_count) {
  const auto neighbors = [&](int node, std::vector<int>& out) {
    for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
      if (unordered_edge_equal(node, nb, move.edge_a, move.swap_a_neighbor) ||
          unordered_edge_equal(node, nb, move.edge_b, move.swap_b_neighbor)) {
        continue;
      }
      out.push_back(nb);
    }
    if (node == move.edge_a) {
      out.push_back(move.swap_b_neighbor);
    } else if (node == move.swap_b_neighbor) {
      out.push_back(move.edge_a);
    }
    if (node == move.edge_b) {
      out.push_back(move.swap_a_neighbor);
    } else if (node == move.swap_a_neighbor) {
      out.push_back(move.edge_b);
    }
  };
  return component_taxa_with_neighbors(tree, move.edge_b, move.edge_a,
                                       taxon_count, neighbors);
}

// The exact support gate. A gene supports the bipartition
// whose side without taxon 0 is W when one of its edges has W as the side
// away from the gene's smallest taxon, keeps all copies of every taxon on
// one side, and leaves at least two of the gene's taxa on each side, in a
// gene of at least four taxa. That holds exactly when, in the gene rooted at
// its reference root, the LCA v of all copies of W's taxa is not that root,
// the gene holds every taxon of W, v holds all or none of every taxon's
// copies and exactly |W| distinct taxa, and the gene has at least four taxa
// and at least |W| + 2. The LCA and the number of present taxa of a disjoint
// union come from its parts, so with the species tree rooted at taxon 0's
// leaf every node keeps both, per gene, for the taxa below it. A query is
// then one LCA per gene and arm, and an applied move changes only the taxa
// below its two central nodes.
class NniExactSupport {
 public:
  explicit NniExactSupport(const SolverContext& ctx) {
    const auto started = std::chrono::steady_clock::now();
    const std::vector<GeneIndex>& genes = scoring_genes(ctx);
    n_ = checked_int_size(context_taxa(ctx).size(),
                          "NNI taxon count exceeds its integer domain");
    k_ = checked_int_size(genes.size(),
                          "NNI gene count exceeds its integer domain");
    genes_.resize(static_cast<std::size_t>(k_));
    taxon_lca_.assign(checked_size_product(
                          static_cast<std::size_t>(n_),
                          static_cast<std::size_t>(k_),
                          "NNI support table exceeds the addressable size"),
                      -1);
    const int chunks =
        std::max(1, std::min(k_, hdt_workspace_worker_limit(ctx.options)));
    // Each gene writes only its own entries.
    parallel_for_inner(0, chunks, [&](int chunk) {
      const int begin = static_cast<int>(
          static_cast<int64_t>(k_) * chunk / chunks);
      const int end = static_cast<int>(
          static_cast<int64_t>(k_) * (chunk + 1) / chunks);
      GeneCopyScan copies(n_);
      GeneTraversal reference_scratch;
      for (int g = begin; g < end; ++g) {
        const GeneIndex& gene = genes[static_cast<std::size_t>(g)];
        Gene& out = genes_[static_cast<std::size_t>(g)];
        if (gene.nodes.empty()) continue;
        if (gene.nodes.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max())) {
          throw std::length_error(
              "NNI support gene exceeds its integer node domain");
        }
        out.lca = build_static_lca_index(gene, n_);
        // Only the climb of lca() is needed.
        std::vector<int>().swap(out.lca.tin);
        std::vector<int>().swap(out.lca.tout);
        out.root = gene.reference_root;
        const GeneTraversal& tr =
            gene_reference_traversal(gene, reference_scratch);
        out.taxa = copies.scan_taxa(gene, tr);
        copies.scan_nodes(gene, tr);
        out.closed_taxa.assign(gene.nodes.size(), 0);
        for (int node = 0; node < static_cast<int>(gene.nodes.size());
             ++node) {
          if (copies.closed(node)) {
            out.closed_taxa[static_cast<std::size_t>(node)] =
                copies.distinct(node);
          }
        }
        for (int32_t taxon : copies.taxa()) {
          taxon_lca_[slot(taxon, g)] = out.lca.lca(
              tr.preorder[static_cast<std::size_t>(copies.first(taxon))],
              tr.preorder[static_cast<std::size_t>(copies.last(taxon))]);
        }
      }
    }, 1);
    stats_add_time(run_stats().nni_split_support_ns,
                   std::chrono::duration_cast<std::chrono::nanoseconds>(
                       std::chrono::steady_clock::now() - started));
  }

  // The support of the bipartition `side_taxa | rest`, read on its side
  // without taxon 0; -1 when it is trivial. side_taxa holds distinct taxa.
  // O(k |side|) LCA steps.
  int support_of_side(const std::vector<int>& side_taxa) const {
    std::vector<int> side;
    const int size = normalized_side(side_taxa, side);
    return size < 0 ? -1 : support_of_taxa(side, size);
  }

  // The side without taxon 0 of `side_taxa | rest`, sorted, and its size; -1
  // when the bipartition is trivial.
  int normalized_side(const std::vector<int>& side_taxa,
                      std::vector<int>& side) const {
    side.clear();
    if (n_ < 4) return -1;
    const int side_size = static_cast<int>(side_taxa.size());
    const int other_size = n_ - side_size;
    if (side_size <= 1 || other_size <= 1) return -1;
    bool contains_ref = false;
    for (int taxon : side_taxa) {
      if (taxon < 0 || taxon >= n_) continue;
      contains_ref = contains_ref || taxon == 0;
      side.push_back(taxon);
    }
    std::sort(side.begin(), side.end());
    // A side that lists a taxon twice, or one outside [0, n), equals no gene
    // side (those are sets of the gene's taxa): it keeps its listed size for
    // the trivial-side test and supports nothing.
    const bool listed_once =
        static_cast<int>(side.size()) == side_size &&
        std::adjacent_find(side.begin(), side.end()) == side.end();
    if (contains_ref) {
      std::vector<uint8_t> in_side(static_cast<std::size_t>(n_), 0);
      for (int taxon : side) in_side[static_cast<std::size_t>(taxon)] = 1;
      side.clear();
      for (int taxon = 0; taxon < n_; ++taxon) {
        if (!in_side[static_cast<std::size_t>(taxon)]) side.push_back(taxon);
      }
    }
    if (!listed_once) side.clear();
    const int size = contains_ref ? other_size : side_size;
    return size < 2 || size > n_ - 2 ? -1 : size;
  }

  int support_of_taxa(const std::vector<int>& side, int size) const {
    return count(size, [&](int g) {
      int lca = -1;
      int present = 0;
      const StaticLcaIndex& index = genes_[static_cast<std::size_t>(g)].lca;
      for (int taxon : side) {
        const int node = taxon_lca_[slot(taxon, g)];
        if (node < 0) continue;
        lca = index.lca(lca, node);
        ++present;
      }
      return std::pair<int, int>{lca, present};
    });
  }

  // Roots `tree` at taxon 0's leaf and computes every node's aggregates. A
  // tree that does not hold each of the n taxa exactly once is answered
  // through support_of_side instead.
  void attach(const SolverTree& tree) {
    tree_complete_ = false;
    const int node_count = checked_int_size(
        tree.nodes.size(), "NNI support tree exceeds its integer domain");
    int root = -1;
    std::vector<uint8_t> seen_taxon(static_cast<std::size_t>(n_), 0);
    int taxa = 0;
    for (int v = 0; v < node_count; ++v) {
      const SolverTree::Node& node = tree.nodes[static_cast<std::size_t>(v)];
      if (!node.active || node.taxon < 0) continue;
      if (node.taxon >= n_ ||
          seen_taxon[static_cast<std::size_t>(node.taxon)] != 0) {
        return;
      }
      seen_taxon[static_cast<std::size_t>(node.taxon)] = 1;
      ++taxa;
      if (node.taxon == 0) root = v;
    }
    if (taxa != n_ || root < 0 || n_ < 4) return;

    parent_.assign(static_cast<std::size_t>(node_count), -1);
    below_.assign(static_cast<std::size_t>(node_count), 0);
    const std::size_t entries = checked_size_product(
        static_cast<std::size_t>(node_count), static_cast<std::size_t>(k_),
        "NNI support table exceeds the addressable size");
    node_lca_.assign(entries, -1);
    node_present_.assign(entries, 0);
    std::vector<int> order;
    order.reserve(static_cast<std::size_t>(node_count));
    std::vector<int> stack{root};
    std::vector<uint8_t> visited(static_cast<std::size_t>(node_count), 0);
    visited[static_cast<std::size_t>(root)] = 1;
    while (!stack.empty()) {
      const int v = stack.back();
      stack.pop_back();
      order.push_back(v);
      for (int nb : active_neighbors(tree, v)) {
        if (visited[static_cast<std::size_t>(nb)]) continue;
        visited[static_cast<std::size_t>(nb)] = 1;
        parent_[static_cast<std::size_t>(nb)] = v;
        stack.push_back(nb);
      }
    }
    // Children before parents; genes in parallel, each writing its own
    // entries.
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
      if (*it != root) set_below(tree, *it);
    }
    const int chunks = std::max(1, std::min(k_, current_inner_parallel_worker_limit()));
    parallel_for_inner(0, chunks, [&](int chunk) {
      const int begin = static_cast<int>(
          static_cast<int64_t>(k_) * chunk / chunks);
      const int end = static_cast<int>(
          static_cast<int64_t>(k_) * (chunk + 1) / chunks);
      for (auto it = order.rbegin(); it != order.rend(); ++it) {
        if (*it != root) recompute(tree, *it, begin, end);
      }
    }, 1);
    tree_complete_ = true;
  }

  // The support of the new central bipartition of the applicable `move` minus
  // that of the old one; 0 when either is trivial.
  int delta_for_move(const SolverTree& tree, const NniMove& move) const {
    if (n_ < 4) return 0;
    if (!tree_complete_) {
      // Equal sides count for neither.
      std::vector<int> old_side;
      std::vector<int> new_side;
      const int old_size = normalized_side(
          old_central_side_taxa(tree, move.edge_a, move.edge_b, n_), old_side);
      if (old_size < 0) return 0;
      const int new_size =
          normalized_side(new_central_side_taxa(tree, move, n_), new_side);
      if (new_size < 0 || new_side == old_side) return 0;
      return support_of_taxa(new_side, new_size) -
             support_of_taxa(old_side, old_size);
    }
    const int a = move.edge_a;
    const int b = move.edge_b;
    const int lower = parent_[static_cast<std::size_t>(b)] == a ? b : a;
    const int upper = lower == b ? a : b;
    if (parent_[static_cast<std::size_t>(lower)] != upper) {
      throw std::logic_error("NNI support tree lost its central edge");
    }
    const int root_arm = parent_[static_cast<std::size_t>(upper)];
    const int old_size = below_[static_cast<std::size_t>(lower)];
    if (std::min(old_size, n_ - old_size) < 2) return 0;
    // The arms of each central node after the move; the side without taxon
    // 0 is the one whose arms do not include the arm that holds it.
    std::vector<int> arms_a;
    std::vector<int> arms_b;
    for (int nb : active_neighbors(tree, a)) {
      if (nb != b && nb != move.swap_a_neighbor) arms_a.push_back(nb);
    }
    arms_a.push_back(move.swap_b_neighbor);
    for (int nb : active_neighbors(tree, b)) {
      if (nb != a && nb != move.swap_b_neighbor) arms_b.push_back(nb);
    }
    arms_b.push_back(move.swap_a_neighbor);
    const std::vector<int>& arms =
        std::find(arms_a.begin(), arms_a.end(), root_arm) != arms_a.end()
            ? arms_b
            : arms_a;
    int new_size = 0;
    for (int arm : arms) new_size += below_[static_cast<std::size_t>(arm)];
    if (std::min(new_size, n_ - new_size) < 2) return 0;
    return support_of_arms(arms, new_size) -
           support_of_arms({lower}, old_size);
  }

  // Brings the aggregates up to date after `move` was applied to `tree`.
  void applied(const SolverTree& tree, const NniMove& move) {
    if (!tree_complete_) return;
    const int a = move.edge_a;
    const int b = move.edge_b;
    const int upper = parent_[static_cast<std::size_t>(b)] == a ? a : b;
    const int root_arm = parent_[static_cast<std::size_t>(upper)];
    // The arm that holds taxon 0 now hangs from this central node.
    const int new_upper = root_arm == move.swap_a_neighbor   ? b
                          : root_arm == move.swap_b_neighbor ? a
                                                             : upper;
    const int new_lower = new_upper == a ? b : a;
    parent_[static_cast<std::size_t>(new_upper)] = root_arm;
    parent_[static_cast<std::size_t>(new_lower)] = new_upper;
    for (int center : {a, b}) {
      for (int nb : active_neighbors(tree, center)) {
        if (nb != a && nb != b && nb != root_arm) {
          parent_[static_cast<std::size_t>(nb)] = center;
        }
      }
    }
    for (int node : {new_lower, new_upper}) {
      set_below(tree, node);
      recompute(tree, node, 0, k_);
    }
  }

 private:
  struct Gene {
    StaticLcaIndex lca;
    // The number of distinct taxa below a node when it holds all or none of
    // every taxon's copies, else 0.
    std::vector<int32_t> closed_taxa;
    int root = -1;
    int taxa = 0;
  };

  std::size_t slot(int item, int gene) const {
    return static_cast<std::size_t>(item) * static_cast<std::size_t>(k_) +
           static_cast<std::size_t>(gene);
  }

  // The number of genes with an eligible side of `size` taxa equal to the
  // side whose (LCA, present taxa) in gene g is per_gene(g).
  template <class PerGene>
  int count(int size, PerGene&& per_gene) const {
    int supported = 0;
    for (int g = 0; g < k_; ++g) {
      const Gene& gene = genes_[static_cast<std::size_t>(g)];
      if (gene.taxa < 4 || size > gene.taxa - 2) continue;
      const auto [lca, present] = per_gene(g);
      if (present != size || lca < 0 || lca == gene.root) continue;
      if (gene.closed_taxa[static_cast<std::size_t>(lca)] == size) ++supported;
    }
    return supported;
  }

  int support_of_arms(const std::vector<int>& arms, int size) const {
    return count(size, [&](int g) {
      int lca = -1;
      int present = 0;
      const StaticLcaIndex& index = genes_[static_cast<std::size_t>(g)].lca;
      for (int arm : arms) {
        lca = index.lca(lca, node_lca_[slot(arm, g)]);
        present += node_present_[slot(arm, g)];
      }
      return std::pair<int, int>{lca, present};
    });
  }

  void set_below(const SolverTree& tree, int node) {
    const int taxon = tree.nodes[static_cast<std::size_t>(node)].taxon;
    if (taxon >= 0) {
      below_[static_cast<std::size_t>(node)] = 1;
      return;
    }
    int below = 0;
    for (int nb : active_neighbors(tree, node)) {
      if (nb != parent_[static_cast<std::size_t>(node)]) {
        below += below_[static_cast<std::size_t>(nb)];
      }
    }
    below_[static_cast<std::size_t>(node)] = below;
  }

  // Recomputes genes [begin, end) of `node` from its children.
  void recompute(const SolverTree& tree, int node, int begin, int end) {
    const int taxon = tree.nodes[static_cast<std::size_t>(node)].taxon;
    if (taxon >= 0) {
      for (int g = begin; g < end; ++g) {
        const int lca = taxon_lca_[slot(taxon, g)];
        node_lca_[slot(node, g)] = lca;
        node_present_[slot(node, g)] = lca >= 0 ? 1 : 0;
      }
      return;
    }
    const std::vector<int> neighbors = active_neighbors(tree, node);
    for (int g = begin; g < end; ++g) {
      const StaticLcaIndex& index = genes_[static_cast<std::size_t>(g)].lca;
      int lca = -1;
      int present = 0;
      for (int nb : neighbors) {
        if (nb == parent_[static_cast<std::size_t>(node)]) continue;
        lca = index.lca(lca, node_lca_[slot(nb, g)]);
        present += node_present_[slot(nb, g)];
      }
      node_lca_[slot(node, g)] = lca;
      node_present_[slot(node, g)] = present;
    }
  }

  int n_ = 0;
  int k_ = 0;
  std::vector<Gene> genes_;
  // [taxon * k + gene]: the LCA of the taxon's copies, -1 when absent.
  std::vector<int32_t> taxon_lca_;
  // The species tree rooted at taxon 0's leaf, when it holds every taxon
  // once: parents, the number of taxa below, and per (node, gene) the LCA of
  // the copies of those taxa and how many of them the gene holds.
  bool tree_complete_ = false;
  std::vector<int32_t> parent_;
  std::vector<int32_t> below_;
  std::vector<int32_t> node_lca_;
  std::vector<int32_t> node_present_;
};

// The sorted keys of every edge at a node within `radius` of a seed.
std::vector<uint64_t> collect_affected_nni_edge_keys(
    const SolverTree& tree,
    const std::vector<int>& seeds,
    int radius) {
  std::vector<uint64_t> affected;
  if (tree.nodes.empty() || seeds.empty()) return affected;

  std::vector<int> distance(tree.nodes.size(), -1);
  std::vector<int> queue;
  queue.reserve(tree.nodes.size());
  for (int seed : seeds) {
    if (seed < 0 || seed >= static_cast<int>(tree.nodes.size()) ||
        !tree.nodes[static_cast<std::size_t>(seed)].active ||
        distance[static_cast<std::size_t>(seed)] >= 0) {
      continue;
    }
    distance[static_cast<std::size_t>(seed)] = 0;
    queue.push_back(seed);
  }

  for (std::size_t head = 0; head < queue.size(); ++head) {
    const int node = queue[head];
    const int node_distance = distance[static_cast<std::size_t>(node)];
    for (int nb : active_neighbors(tree, node)) {
      affected.push_back(nni_edge_key(node, nb));
      if (node_distance < radius &&
          distance[static_cast<std::size_t>(nb)] < 0) {
        distance[static_cast<std::size_t>(nb)] = node_distance + 1;
        queue.push_back(nb);
      }
    }
  }
  std::sort(affected.begin(), affected.end());
  affected.erase(std::unique(affected.begin(), affected.end()), affected.end());
  return affected;
}

// Drops cached moves whose central edge is affected or that do not apply.
void erase_cached_nni_edges(std::vector<CachedNniMove>& cache,
                            const SolverTree& tree,
                            const std::vector<uint64_t>& affected) {
  cache.erase(
      std::remove_if(
          cache.begin(), cache.end(),
          [&](const CachedNniMove& cached) {
            return std::binary_search(affected.begin(), affected.end(),
                                      cached.central_edge_key) ||
                   !nni_move_applicable(tree, cached.move);
          }),
      cache.end());
}

// Sorted union.
void merge_edge_key_sets(std::vector<uint64_t>& into,
                         const std::vector<uint64_t>& from) {
  std::vector<uint64_t> merged;
  merged.reserve(into.size() + from.size());
  std::set_union(into.begin(), into.end(), from.begin(), from.end(),
                 std::back_inserter(merged));
  into = std::move(merged);
}

// One pass of --nni-batch: among the sets of cached moves whose gains add
// exactly, a set of largest total gain, sorted by stable key. The gains of a
// set add when no two of its central edges share a vertex and, of two central
// edges joined by one edge, both moves keep that edge. A move at central edge
// uw that swaps side subtree a at u with b at w gives the same tree as the
// swap of the other two side subtrees, so each move has two implementations,
// each keeping a different pair of side edges, and the set may use either.
// Moves that fail the support gate, when `split_support` is set, are left out.
//
// A dynamic program over the tree rooted at its first vertex of degree three
// finds the set in O(n log n). For a vertex v below parent p:
//   N:  pv is not central and no central edge ends at v;
//   S:  pv is not central, and a central edge at v, if any, goes to a child
//       and its move keeps pv;
//   A:  pv is not central, and a central edge at v, if any, goes to a child;
//   E0, E1: pv is central and its move keeps the edge to child 0 or 1 of v.
// A vertex without three edges is a barrier: no move has a central edge at
// it (proposals need binary ends), so it passes up the sum of its children's
// A values, and moves on either side of it never interact.
std::vector<NniMove> select_nni_batch(
    const SolverTree& tree,
    const std::vector<CachedNniMove>& cache,
    const NniExactSupport* split_support,
    double min_gain,
    double max_support_drop) {
  constexpr double kNegInf = -std::numeric_limits<double>::infinity();
  const auto slot = [](int node) { return static_cast<std::size_t>(node); };
  // The neighbour of `node` other than `skip_a` and `skip_b`.
  const auto other_neighbor = [&](int node, int skip_a, int skip_b) {
    for (int nb : active_neighbors(tree, node)) {
      if (nb != skip_a && nb != skip_b) return nb;
    }
    return -1;
  };

  // Both implementations of every qualifying move as (u, w, side cut at u,
  // side cut at w) with u < w; for each key, the largest gain comes first.
  struct Implementation {
    int u = -1;
    int w = -1;
    int cut_u = -1;
    int cut_w = -1;
    double gain = 0.0;
  };
  const auto key = [](const Implementation& x) {
    return std::array<int, 4>{x.u, x.w, x.cut_u, x.cut_w};
  };
  std::vector<Implementation> implementations;
  for (const CachedNniMove& cached : cache) {
    const NniMove& move = cached.move;
    if (!nni_gain_meets_threshold(move.gain, min_gain) ||
        !nni_move_applicable(tree, move)) {
      continue;
    }
    if (split_support != nullptr &&
        split_support->delta_for_move(tree, move) < -max_support_drop) {
      stats_add(run_stats().nni_support_stops);
      continue;
    }
    const int a = move.edge_a;
    const int b = move.edge_b;
    const std::array<std::array<int, 2>, 2> cuts{{
        {move.swap_a_neighbor, move.swap_b_neighbor},
        {other_neighbor(a, b, move.swap_a_neighbor),
         other_neighbor(b, a, move.swap_b_neighbor)},
    }};
    for (const std::array<int, 2>& cut : cuts) {
      if (cut[0] < 0 || cut[1] < 0) continue;
      implementations.push_back(
          a < b ? Implementation{a, b, cut[0], cut[1], move.gain}
                : Implementation{b, a, cut[1], cut[0], move.gain});
    }
  }
  std::sort(implementations.begin(), implementations.end(),
            [&](const Implementation& x, const Implementation& y) {
              if (key(x) != key(y)) return key(x) < key(y);
              return x.gain > y.gain;
            });
  std::vector<NniMove> batch;
  if (implementations.empty()) return batch;
  // The gain of the move at central edge vx that keeps side edge keep_v at v
  // and keep_x at x; -inf when no qualifying move does.
  const auto gain_of = [&](int v, int x, int keep_v, int keep_x) {
    const int cut_v = other_neighbor(v, x, keep_v);
    const int cut_x = other_neighbor(x, v, keep_x);
    const Implementation probe = v < x ? Implementation{v, x, cut_v, cut_x, 0.0}
                                       : Implementation{x, v, cut_x, cut_v, 0.0};
    const auto it = std::lower_bound(
        implementations.begin(), implementations.end(), probe,
        [&](const Implementation& entry, const Implementation& query) {
          return key(entry) < key(query);
        });
    return it != implementations.end() && key(*it) == key(probe) ? it->gain
                                                                 : kNegInf;
  };

  const int node_count = static_cast<int>(tree.nodes.size());
  int root = -1;
  for (int v = 0; v < node_count && root < 0; ++v) {
    if (tree.nodes[slot(v)].active && active_neighbors(tree, v).size() == 3) {
      root = v;
    }
  }
  for (int v = 0; v < node_count && root < 0; ++v) {
    if (tree.nodes[slot(v)].active && active_neighbors(tree, v).size() > 1) {
      root = v;
    }
  }
  if (root < 0) return batch;

  // Children in increasing node order: two for a binary vertex, three for a
  // binary root, any number below a barrier.
  std::vector<int> parent(slot(node_count), -1);
  std::vector<std::array<int, 2>> children(slot(node_count),
                                           std::array<int, 2>{-1, -1});
  std::vector<std::vector<int>> barrier_children(slot(node_count));
  std::array<int, 3> root_children{-1, -1, -1};
  std::vector<int> preorder;
  std::vector<char> seen(slot(node_count), 0);
  std::vector<int> stack{root};
  seen[slot(root)] = 1;
  while (!stack.empty()) {
    const int v = stack.back();
    stack.pop_back();
    preorder.push_back(v);
    std::vector<int> kids;
    for (int nb : active_neighbors(tree, v)) {
      if (seen[slot(nb)] == 0) kids.push_back(nb);
    }
    std::sort(kids.begin(), kids.end());
    if (v == root && kids.size() == 3) {
      root_children = {kids[0], kids[1], kids[2]};
    } else if (v != root && kids.size() == 2) {
      children[slot(v)] = {kids[0], kids[1]};
    } else if (!kids.empty()) {
      barrier_children[slot(v)] = kids;
    }
    for (int kid : kids) {
      seen[slot(kid)] = 1;
      parent[slot(kid)] = v;
      stack.push_back(kid);
    }
  }

  struct Values {
    double N = 0.0;
    double S = 0.0;
    double A = 0.0;
    double E0 = kNegInf;
    double E1 = kNegInf;
  };
  // How S and A were reached, with c0, c1 the children of v: 0, no central
  // edge at v; 1 + s, edge v-c0 keeping pv and child s of c0; 3 + s, edge
  // v-c1 keeping pv and child s of c1; 5 + s, edge v-c0 keeping c1 and child
  // s of c0; 7 + s, edge v-c1 keeping c0 and child s of c1.
  struct Choice {
    std::uint8_t S = 0;
    std::uint8_t A = 0;
    bool A_from_S = true;
  };
  std::vector<Values> values(slot(node_count));
  std::vector<Choice> choices(slot(node_count));
  // The best value with vx central, v keeping keep_v; *kept is the child of x
  // that the move keeps.
  const auto best_central = [&](int v, int x, int keep_v, int* kept) {
    const std::array<int, 2>& grandchildren = children[slot(x)];
    if (grandchildren[0] < 0) return kNegInf;
    const Values& below = values[slot(x)];
    const double keep0 = below.E0 + gain_of(v, x, keep_v, grandchildren[0]);
    const double keep1 = below.E1 + gain_of(v, x, keep_v, grandchildren[1]);
    *kept = keep1 > keep0 ? 1 : 0;
    return keep1 > keep0 ? keep1 : keep0;
  };
  for (auto it = preorder.rbegin(); it != preorder.rend(); ++it) {
    const int v = *it;
    if (v == root) continue;
    Values& out = values[slot(v)];
    if (!barrier_children[slot(v)].empty()) {
      double sum = 0.0;
      for (int kid : barrier_children[slot(v)]) sum += values[slot(kid)].A;
      out = Values{};
      out.N = out.S = out.A = sum;
      continue;
    }
    const std::array<int, 2>& c = children[slot(v)];
    if (c[0] < 0) {
      out = Values{};
      continue;
    }
    const Values& L = values[slot(c[0])];
    const Values& R = values[slot(c[1])];
    const int p = parent[slot(v)];
    int kept_0p = 0;
    int kept_1p = 0;
    int kept_01 = 0;
    int kept_10 = 0;
    const double c0_keep_parent = best_central(v, c[0], p, &kept_0p);
    const double c1_keep_parent = best_central(v, c[1], p, &kept_1p);
    const double c0_keep_sibling = best_central(v, c[0], c[1], &kept_01);
    const double c1_keep_sibling = best_central(v, c[1], c[0], &kept_10);
    Choice& choice = choices[slot(v)];
    out.N = L.A + R.A;
    out.S = out.N;
    choice.S = 0;
    if (c0_keep_parent + R.N > out.S) {
      out.S = c0_keep_parent + R.N;
      choice.S = static_cast<std::uint8_t>(1 + kept_0p);
    }
    if (c1_keep_parent + L.N > out.S) {
      out.S = c1_keep_parent + L.N;
      choice.S = static_cast<std::uint8_t>(3 + kept_1p);
    }
    out.A = out.S;
    choice.A_from_S = true;
    if (c0_keep_sibling + R.S > out.A) {
      out.A = c0_keep_sibling + R.S;
      choice.A = static_cast<std::uint8_t>(5 + kept_01);
      choice.A_from_S = false;
    }
    if (c1_keep_sibling + L.S > out.A) {
      out.A = c1_keep_sibling + L.S;
      choice.A = static_cast<std::uint8_t>(7 + kept_10);
      choice.A_from_S = false;
    }
    out.E0 = L.S + R.N;
    out.E1 = L.N + R.S;
  }

  // At the root: no central edge, or root edge i central, keeping root edge j
  // and cutting the third.
  const bool root_barrier = root_children[0] < 0;
  double best = 0.0;
  int best_i = -1;
  int best_j = -1;
  int best_kept = 0;
  if (root_barrier) {
    for (int kid : barrier_children[slot(root)]) best += values[slot(kid)].A;
  } else {
    for (int kid : root_children) best += values[slot(kid)].A;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        if (j == i) continue;
        const int ci = root_children[slot(i)];
        const int cj = root_children[slot(j)];
        const int ck = root_children[slot(3 - i - j)];
        int kept = 0;
        const double value = best_central(root, ci, cj, &kept) +
                             values[slot(cj)].S + values[slot(ck)].N;
        if (value > best) {
          best = value;
          best_i = i;
          best_j = j;
          best_kept = kept;
        }
      }
    }
  }

  // Rebuild the set top-down; states 0 N, 1 S, 2 A, 3 E0, 4 E1.
  const auto emit = [&](int v, int x, int keep_v, int keep_x) {
    NniMove move;
    move.edge_a = v;
    move.edge_b = x;
    move.swap_a_neighbor = other_neighbor(v, x, keep_v);
    move.swap_b_neighbor = other_neighbor(x, v, keep_x);
    move.gain = gain_of(v, x, keep_v, keep_x);
    batch.push_back(move);
  };
  std::vector<std::pair<int, int>> pending;
  if (root_barrier) {
    for (int kid : barrier_children[slot(root)]) pending.push_back({kid, 2});
  } else if (best_i < 0) {
    for (int kid : root_children) pending.push_back({kid, 2});
  } else {
    const int ci = root_children[slot(best_i)];
    const int cj = root_children[slot(best_j)];
    const int ck = root_children[slot(3 - best_i - best_j)];
    emit(root, ci, cj, children[slot(ci)][slot(best_kept)]);
    pending.push_back({ci, 3 + best_kept});
    pending.push_back({cj, 1});
    pending.push_back({ck, 0});
  }
  while (!pending.empty()) {
    const auto [v, state] = pending.back();
    pending.pop_back();
    if (!barrier_children[slot(v)].empty()) {
      for (int kid : barrier_children[slot(v)]) pending.push_back({kid, 2});
      continue;
    }
    const std::array<int, 2>& c = children[slot(v)];
    if (c[0] < 0) continue;
    if (state == 3) {
      pending.push_back({c[0], 1});
      pending.push_back({c[1], 0});
      continue;
    }
    if (state == 4) {
      pending.push_back({c[0], 0});
      pending.push_back({c[1], 1});
      continue;
    }
    const Choice& choice = choices[slot(v)];
    const int p = parent[slot(v)];
    int code = 0;
    if (state == 1) code = choice.S;
    if (state == 2) code = choice.A_from_S ? choice.S : choice.A;
    switch (code) {
      case 1:
      case 2:
        emit(v, c[0], p, children[slot(c[0])][slot(code - 1)]);
        pending.push_back({c[0], 3 + code - 1});
        pending.push_back({c[1], 0});
        break;
      case 3:
      case 4:
        emit(v, c[1], p, children[slot(c[1])][slot(code - 3)]);
        pending.push_back({c[1], 3 + code - 3});
        pending.push_back({c[0], 0});
        break;
      case 5:
      case 6:
        emit(v, c[0], c[1], children[slot(c[0])][slot(code - 5)]);
        pending.push_back({c[0], 3 + code - 5});
        pending.push_back({c[1], 1});
        break;
      case 7:
      case 8:
        emit(v, c[1], c[0], children[slot(c[1])][slot(code - 7)]);
        pending.push_back({c[1], 3 + code - 7});
        pending.push_back({c[0], 1});
        break;
      default:
        pending.push_back({c[0], 2});
        pending.push_back({c[1], 2});
        break;
    }
  }
  std::sort(batch.begin(), batch.end(), [](const NniMove& x, const NniMove& y) {
    return nni_move_stable_key(x) < nni_move_stable_key(y);
  });
  return batch;
}

// Keeps the `count` moves of largest gain, ties by stable key, in stable key
// order. A subset of a batch is a batch: its gains still add.
void keep_largest_nni_gains(std::vector<NniMove>& batch, std::size_t count) {
  if (batch.size() <= count) return;
  std::sort(batch.begin(), batch.end(), [](const NniMove& x, const NniMove& y) {
    if (x.gain != y.gain) return x.gain > y.gain;
    return nni_move_stable_key(x) < nni_move_stable_key(y);
  });
  batch.resize(count);
  std::sort(batch.begin(), batch.end(), [](const NniMove& x, const NniMove& y) {
    return nni_move_stable_key(x) < nni_move_stable_key(y);
  });
}

}  // namespace

int nni_split_support_for_tests(const SolverContext& ctx,
                                    const std::vector<int>& side_taxa) {
  return NniExactSupport(ctx).support_of_side(side_taxa);
}

std::vector<int> nni_support_gate_trace_for_tests(
    const SolverContext& ctx,
    SolverTree& tree,
    const std::vector<NniMove>& moves) {
  NniExactSupport gate(ctx);
  gate.attach(tree);
  std::vector<int> deltas;
  for (const NniMove& move : moves) {
    if (!nni_move_applicable(tree, move)) {
      deltas.push_back(std::numeric_limits<int>::min());
      continue;
    }
    deltas.push_back(gate.delta_for_move(tree, move));
    apply_nni_move(tree, move);
    gate.applied(tree, move);
  }
  return deltas;
}

namespace {

std::vector<NniMove> scored_moves(const ScoredNniData& scored) {
  std::vector<NniMove> moves;
  for (std::size_t i = 0; i < scored.proposals.size(); ++i) {
    const NniEdgeProposal& proposal = scored.proposals[i];
    const NniProposalScores& scores = scored.scores[i];
    moves.push_back(proposal.swap_p0_q0);
    moves.back().gain = scores.swap_p0_q0 - scores.current;
    moves.push_back(proposal.swap_p0_q1);
    moves.back().gain = scores.swap_p0_q1 - scores.current;
  }
  return moves;
}

}  // namespace

std::vector<NniMove> nni_scored_moves_for_tests(
    const SolverContext& ctx,
    const SolverTree& tree,
    const std::vector<int>& active_taxa) {
  return scored_moves(score_nni_edges(ctx, tree, active_taxa));
}

std::vector<NniMove> nni_local_scored_moves_for_tests(
    const SolverContext& ctx,
    const SolverTree& tree,
    const std::vector<int>& active_taxa,
    const std::vector<std::pair<int, int>>& central_edges) {
  std::vector<uint64_t> allowed;
  allowed.reserve(central_edges.size());
  for (const auto& [a, b] : central_edges) allowed.push_back(nni_edge_key(a, b));
  std::sort(allowed.begin(), allowed.end());
  allowed.erase(std::unique(allowed.begin(), allowed.end()), allowed.end());
  return scored_moves(score_nni_edges(ctx, tree, active_taxa, &allowed));
}

std::vector<NniMove> nni_batch_for_tests(const SolverContext& ctx,
                                         const SolverTree& tree,
                                         const std::vector<int>& active_taxa) {
  std::vector<CachedNniMove> cache;
  add_cached_positive_moves(score_nni_edges(ctx, tree, active_taxa), cache,
                            kNniGainEps);
  return select_nni_batch(tree, cache, nullptr, kNniGainEps, 0.0);
}

NniMove select_nni_move_for_tests(
    const SolverTree& tree,
    const std::vector<NniMove>& candidates,
    double min_gain) {
  std::vector<CachedNniMove> cache;
  cache.reserve(candidates.size());
  for (const NniMove& move : candidates) {
    if (!nni_gain_meets_threshold(move.gain, min_gain)) continue;
    cache.push_back(CachedNniMove{
        move, nni_edge_key(move.edge_a, move.edge_b)});
  }
  return best_cached_nni_move(tree, cache);
}

// Gated NNI on the final tree. Every move is scored once and cached; after
// moves are applied, only edges within distance two of them are rescored.
// One move at a time, the cached move of largest gain is applied until no
// move passes the gain gate, the support gate stops, or the cap is reached.
// With --nni-batch, each pass applies the select_nni_batch set, which leaves
// out moves that fail the support gate, until the set is empty or the cap is
// reached.
void refine_nni(SolverContext& ctx,
                        SolverTree& tree,
                        const std::vector<int>& active_taxa) {
  if (!ctx.options.nni) return;
  const double min_gain = nni_min_gain(ctx);
  const double max_support_drop =
      kNniMaxSupportDropPerGene * nni_gene_count(ctx);
  const int max_moves = ctx.options.nni_move_cap
                            ? kNniMaxMoves
                            : std::numeric_limits<int>::max();
  const int progress_cap = ctx.options.nni_move_cap ? kNniMaxMoves : 0;
  const auto started = std::chrono::steady_clock::now();
  progress_note_stage("nni");
  std::unique_ptr<NniExactSupport> support_gate;
  if (ctx.options.nni_support_gate) {
    support_gate = std::make_unique<NniExactSupport>(ctx);
    support_gate->attach(tree);
  }
  const auto progress_edges = []() {
    return run_stats().nni_edges_scored.load(
        std::memory_order_relaxed);
  };
  progress_note_nni("init", 0, progress_cap, progress_edges(), 0.0,
                        min_gain);

  ScoredNniData scored = score_nni_edges(ctx, tree, active_taxa);
  std::vector<CachedNniMove> cache;
  cache.reserve(checked_doubled_capacity(
      scored.proposals.size(), false,
      "NNI move cache exceeds the addressable size domain"));
  add_cached_positive_moves(scored, cache, min_gain);

  int moves = 0;
  double accepted_gain = 0.0;
  // The moves' integer gains, summed exactly over the run's NNI passes (the
  // SPR refinement runs NNI again after each round that moves).
  WideCount exact_gain =
      run_stats_enabled() ? run_stats().nni_gain : WideCount{};
  if (ctx.options.nni_batch) {
    while (moves < max_moves) {
      std::vector<NniMove> batch = select_nni_batch(
          tree, cache, support_gate.get(), min_gain, max_support_drop);
      if (batch.empty()) break;
      keep_largest_nni_gains(batch, static_cast<std::size_t>(max_moves - moves));
      std::vector<uint64_t> affected;
      for (const NniMove& move : batch) {
        const std::vector<int> seeds = nni_touched_nodes(move);
        merge_edge_key_sets(affected,
                            collect_affected_nni_edge_keys(tree, seeds, 2));
        if (!apply_nni_move(tree, move)) continue;
        if (support_gate != nullptr) support_gate->applied(tree, move);
        ++moves;
        accepted_gain += move.gain;
        wide_add(exact_gain, move.gain);
        stats_add(run_stats().nni_moves_accepted);
        merge_edge_key_sets(affected,
                            collect_affected_nni_edge_keys(tree, seeds, 2));
      }
      stats_add(run_stats().nni_batch_passes);
      // A local pass collapses each gene around every affected edge, so above
      // a quarter of the edges a full pass is cheaper; the scores, and so the
      // moves, are the same either way.
      std::size_t edges = 0;
      for (const SolverTree::Node& node : tree.nodes) edges += node.active ? 1 : 0;
      edges = edges > 0 ? edges - 1 : 0;
      if (affected.size() * 4 > edges) {
        cache.clear();
        if (moves < max_moves) {
          add_cached_positive_moves(score_nni_edges(ctx, tree, active_taxa),
                                    cache, min_gain);
        }
      } else {
        erase_cached_nni_edges(cache, tree, affected);
        if (moves < max_moves) {
          add_cached_positive_moves(
              score_nni_edges(ctx, tree, active_taxa, &affected), cache,
              min_gain);
        }
      }
      progress_note_nni("move", moves, progress_cap, progress_edges(),
                            accepted_gain, min_gain);
    }
  } else {
    while (moves < max_moves) {
      NniMove move = best_cached_nni_move(tree, cache);
      if (move.edge_a < 0 ||
          !nni_gain_meets_threshold(move.gain, min_gain)) {
        break;
      }
      if (support_gate != nullptr &&
          support_gate->delta_for_move(tree, move) < -max_support_drop) {
        stats_add(run_stats().nni_support_stops);
        break;
      }
      const std::vector<int> seeds = nni_touched_nodes(move);
      std::vector<uint64_t> affected =
          collect_affected_nni_edge_keys(tree, seeds, 2);
      if (!apply_nni_move(tree, move)) {
        // The move is stale: drop it and the cached moves near it.
        const uint64_t failed_key =
            nni_edge_key(move.edge_a, move.edge_b);
        const auto position =
            std::lower_bound(affected.begin(), affected.end(), failed_key);
        if (position == affected.end() || *position != failed_key) {
          affected.insert(position, failed_key);
        }
        erase_cached_nni_edges(cache, tree, affected);
        continue;
      }
      if (support_gate != nullptr) support_gate->applied(tree, move);
      ++moves;
      accepted_gain += move.gain;
      wide_add(exact_gain, move.gain);
      stats_add(run_stats().nni_moves_accepted);
      merge_edge_key_sets(affected, collect_affected_nni_edge_keys(tree, seeds, 2));
      erase_cached_nni_edges(cache, tree, affected);
      if (moves < max_moves) {
        add_cached_positive_moves(
            score_nni_edges(ctx, tree, active_taxa, &affected), cache,
            min_gain);
      }
      progress_note_nni("move", moves, progress_cap, progress_edges(),
                            accepted_gain, min_gain);
    }
  }
  progress_note_nni("done", moves, progress_cap, progress_edges(),
                        accepted_gain, min_gain);
  stats_set_nni_gain(exact_gain);
  stats_add_time(run_stats().nni_ns,
                      std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - started));
}

std::vector<BranchQuartetFrequencies> branch_quartet_frequencies(
    const SolverContext& ctx,
    const SolverTree& tree,
    const std::vector<int>& active_taxa,
    bool per_gene_normalized) {
  const int taxon_count = static_cast<int>(context_taxa(ctx).size());
  const QueryTree query =
      query_tree_from_solver_tree(tree, active_taxa, taxon_count);
  const std::vector<NniEdgeProposal> proposals =
      collect_nni_proposals(tree, query);
  const std::vector<NniProposalScores> scores = score_nni_proposals_hdt(
      ctx, query, proposals, active_taxa, taxon_count, per_gene_normalized);
  // The taxa on the far side of p's branch b: below b for a child, the rest
  // of the tree for p's parent.
  const auto arm_taxa = [&](int p, int b) -> double {
    if (p != query.root && b == query.parent[static_cast<std::size_t>(p)]) {
      return static_cast<double>(query.leaf_count - query_leaf_count(query, p));
    }
    return static_cast<double>(query_leaf_count(query, b));
  };
  std::vector<BranchQuartetFrequencies> out;
  out.reserve(proposals.size());
  for (std::size_t i = 0; i < proposals.size(); ++i) {
    const NniEdgeProposal& proposal = proposals[i];
    BranchQuartetFrequencies f;
    f.edge_a = proposal.swap_p0_q0.edge_a;
    f.edge_b = proposal.swap_p0_q0.edge_b;
    f.main = scores[i].current;
    f.alternative1 = scores[i].swap_p0_q0;
    f.alternative2 = scores[i].swap_p0_q1;
    const int p = proposal.p_query;
    const int q = proposal.q_query;
    double product = 1.0;
    for (int b : vertex_arm_coloring(query, p).branch_nodes) {
      if (b != q) product *= arm_taxa(p, b);
    }
    for (int c : query_children(query, q)) product *= arm_taxa(q, c);
    f.quartets_around = product;
    out.push_back(f);
  }
  return out;
}

}  // namespace qraft::internal
