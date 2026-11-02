// K=2 quartet scores of the starting-tree edges proposed as bipartitions.
// One HDT traversal of the query tree per gene colours the two sides of each
// scored edge 0 and 1 and reads W(00|11) - W(01|01), satisfied minus
// violated quartets.

#include "solver/solver.hpp"
#include "core/checked.hpp"
#include "counting/count.hpp"
#include "counting/query_taxa.hpp"

#include <algorithm>
#include <chrono>

namespace qraft::internal {
namespace {

void score_bipartitions_for_genes(
    const std::vector<GeneIndex>& genes,
    const OwnerProjection& owners,
    const QueryTree& query,
    const std::vector<unsigned char>& score_node,
    const std::vector<int>& node_to_edge,
    int edge_count,
    int begin_gene,
    int end_gene,
    bool induce_source,
    std::vector<double>& scores);


void score_bipartitions_for_genes(
    const std::vector<GeneIndex>& genes,
    const OwnerProjection& owners,
    const QueryTree& query,
    const std::vector<unsigned char>& score_node,
    const std::vector<int>& node_to_edge,
    int edge_count,
    int begin_gene,
    int end_gene,
    bool induce_source,
    std::vector<double>& scores) {
  begin_gene = std::max(0, begin_gene);
  end_gene = std::max(begin_gene,
                      std::min(end_gene, static_cast<int>(genes.size())));
  const int cpu_gene_count = end_gene - begin_gene;
  if (cpu_gene_count <= 0 || edge_count <= 0) return;

  const int chunks = gene_partition_count(cpu_gene_count);
  (void)checked_hdt_size_mul(
      static_cast<std::size_t>(chunks), static_cast<std::size_t>(edge_count),
      "bipartition partial score table exceeds the addressable size domain");
  RowTable<double> partial(static_cast<std::size_t>(chunks),
                           static_cast<std::size_t>(edge_count));
  std::vector<GeneCountTelemetry> telemetry =
      gene_count_telemetry_slots(chunks);
  const auto build_started = std::chrono::steady_clock::now();
  // Genes differ widely in size, so hand them out one at a time; rows stay
  // one per gene and are reduced in gene order below.
  parallel_for_dynamic(0, chunks, [&](int chunk) {
    const int begin = begin_gene + chunk;
    const int end = begin + 1;
    double* row = partial.init_row(static_cast<std::size_t>(chunk));
    long long visits = 0;
    HdtTraversalStats hdt_stats;
    // Scoped to the chunk: a root-level query can grow it to gene size, and
    // worker threads must not hold that across later subproblems.
    HdtExtractWorkspace hdt_workspace;
    for (int g = begin; g < end; ++g) {
      const GeneIndex& gene = genes[static_cast<std::size_t>(g)];
      if (!gene_has_owner_quartet(gene, owners)) continue;
      HdtColorPolicy colors;
      colors.default_current = 0;
      colors.default_background = 0;
      colors.current_by_owner.assign(
          static_cast<std::size_t>(owners.owner_count), 1);
      colors.background_by_owner.assign(
          static_cast<std::size_t>(owners.owner_count), 0);
      HdtInitialStateResult initial = hdt_build_initial_state(
          induce_source
              ? build_induced_base_component_graph(gene, owners, 2, colors)
              : build_base_component_graph(gene, owners, 2, colors));
      hdt_account_initial_state(initial, hdt_stats);
      count_bipartitions(
          std::move(initial.state), query, query.leaf_active_order, score_node,
          hdt_workspace, hdt_stats,
          [&](int query_node, HdtState& current) {
            const int edge =
                node_to_edge[static_cast<std::size_t>(query_node)];
            if (edge < 0) return;
            const std::array<double, 3> counts =
                hdt_root_counts(current);
            row[static_cast<std::size_t>(edge)] += counts[0] - counts[1];
            ++visits;
          },
          induce_source ? 16 : 4);
    }
    if (!telemetry.empty()) {
      telemetry[static_cast<std::size_t>(chunk)] = GeneCountTelemetry{
          visits, hdt_stats.recolor_changed_leaves,
          hdt_stats.recolor_marked_nodes, hdt_stats};
    }
  }, 1);

  add_rows_in_order(partial.data(), chunks,
                    static_cast<std::size_t>(edge_count), scores.data(),
                    [](double& total, double part) { total += part; });
  merge_gene_count_telemetry(telemetry);
  stats_add(run_stats().hdt_builds,
                        static_cast<long long>(cpu_gene_count));
  stats_add_time(
      run_stats().hdt_build_ns,
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now() - build_started));
}

}  // namespace

// The score of each edge, given by its child node in `query`, summed over the
// scoring genes in gene order. All scores are zero in the balanced ablation
// and with fewer than four active taxa.
std::vector<double> score_bipartitions(
    const SolverContext& ctx,
    const ActiveTaxa& active,
    const QueryTree& query,
    const std::vector<int>& edge_child_nodes) {
  const int edge_count = static_cast<int>(edge_child_nodes.size());
  std::vector<double> scores(static_cast<std::size_t>(edge_count), 0.0);
  if (ctx.options.divide_score == DivideScore::Balanced) {
    return scores;
  }
  if (edge_count == 0 || query.root < 0 ||
      static_cast<int>(active.taxa.size()) < 4) {
    return scores;
  }

  const std::vector<GeneIndex>& genes = scoring_genes(ctx);
  const int gene_count = static_cast<int>(genes.size());
  if (gene_count == 0) return scores;
  note_gene_set_use(genes);
  std::vector<unsigned char> score_node(query.nodes.size(), 0);
  std::vector<int> node_to_edge(query.nodes.size(), -1);
  int valid_edge_count = 0;
  for (int edge = 0; edge < edge_count; ++edge) {
    const int node = edge_child_nodes[static_cast<std::size_t>(edge)];
    if (node < 0 || node >= static_cast<int>(query.nodes.size())) continue;
    score_node[static_cast<std::size_t>(node)] = 1;
    node_to_edge[static_cast<std::size_t>(node)] = edge;
    ++valid_edge_count;
  }
  if (valid_edge_count == 0) return scores;

  const OwnerProjection owners = make_owner_projection(
      ctx, active.active_by_taxon);

  // Under the bipartition test, trivial bipartitions make heavy reuse common,
  // so a gene can hold many leaves outside the active taxa. Those leaves are
  // then pruned before the HDT is built, and the state-reuse bound factor of
  // the traversal rises from 4 to 16.
  score_bipartitions_for_genes(
      genes, owners, query, score_node, node_to_edge, edge_count, 0,
      gene_count,
      ctx.bipartition_test_index != nullptr, scores);
  return scores;
}

}  // namespace qraft::internal
