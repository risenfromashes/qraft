#include "solver/solver.hpp"
#include "regraft/regraft.hpp"
#include "core/checked.hpp"
#include "counting/count.hpp"
#include "counting/query_taxa.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <mutex>
#include <string>

namespace qraft::internal {
namespace {

// Edges within this of the best score are tied; the stable edge key breaks
// the tie.
constexpr double kRegraftTieEps = 1e-9;

// A solved side rooted for regrafting. tin/tout bound each node's subtree as
// preorder positions [tin, tout); subtree_taxa_count counts active taxa below
// a node; edge_child holds the child end of every edge; internal_nodes holds
// the vertices of active degree three.
struct RegraftView {
  std::vector<int> parent;
  std::vector<int> preorder;
  std::vector<int> postorder;
  std::vector<int> tin;
  std::vector<int> tout;
  std::vector<int> subtree_taxa_count;
  std::vector<int> edge_child;
  std::vector<int> internal_nodes;
  int root = -1;
};

int active_degree(const SolverTree& tree, int node) {
  int degree = 0;
  for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
    if (nb >= 0 && nb < static_cast<int>(tree.nodes.size()) &&
        tree.nodes[static_cast<std::size_t>(nb)].active) {
      ++degree;
    }
  }
  return degree;
}

// Roots the active part of `tree` at its first internal node, or at a leaf
// when it has none.
RegraftView build_regraft_view(const SolverTree& tree,
                                    const std::vector<int>& active_taxa,
                                    int taxon_count) {
  RegraftView view;
  const int n = static_cast<int>(tree.nodes.size());
  view.parent.assign(static_cast<std::size_t>(n), -1);
  view.tin.assign(static_cast<std::size_t>(n), 0);
  view.tout.assign(static_cast<std::size_t>(n), 0);
  view.subtree_taxa_count.assign(static_cast<std::size_t>(n), 0);
  if (n == 0) return view;

  const TaxonIndex active_by_taxon =
      make_taxon_index(active_taxa, taxon_count);
  std::vector<int> active_leaf_node(active_by_taxon.taxa.size(), -1);

  for (int i = 0; i < n; ++i) {
    if (!tree.nodes[static_cast<std::size_t>(i)].active) continue;
    const int taxon = tree.nodes[static_cast<std::size_t>(i)].taxon;
    const int active = taxon >= 0 ? active_by_taxon.find(taxon) : -1;
    if (active >= 0) active_leaf_node[static_cast<std::size_t>(active)] = i;
    if (view.root < 0 && taxon < 0 && active_degree(tree, i) >= 2) {
      view.root = i;
    }
  }
  if (view.root < 0) {
    for (int node : active_leaf_node) {
      if (node >= 0) {
        view.root = node;
        break;
      }
    }
  }
  if (view.root < 0) return view;

  std::vector<std::pair<int, int>> stack{{view.root, 0}};
  view.parent[static_cast<std::size_t>(view.root)] = view.root;
  while (!stack.empty()) {
    const int node = stack.back().first;
    int& next = stack.back().second;
    if (next == 0) {
      view.tin[static_cast<std::size_t>(node)] =
          static_cast<int>(view.preorder.size());
      view.preorder.push_back(node);
    }
    if (next < static_cast<int>(
                   tree.nodes[static_cast<std::size_t>(node)].neighbors.size())) {
      const int nb =
          tree.nodes[static_cast<std::size_t>(node)].neighbors
              [static_cast<std::size_t>(next++)];
      if (!tree.nodes[static_cast<std::size_t>(nb)].active ||
          nb == view.parent[static_cast<std::size_t>(node)]) {
        continue;
      }
      view.parent[static_cast<std::size_t>(nb)] = node;
      view.edge_child.push_back(nb);
      stack.push_back({nb, 0});
      continue;
    }
    view.tout[static_cast<std::size_t>(node)] =
        static_cast<int>(view.preorder.size());
    view.postorder.push_back(node);
    stack.pop_back();
  }

  for (int node : view.postorder) {
    int count = 0;
    const int taxon = tree.nodes[static_cast<std::size_t>(node)].taxon;
    if (taxon >= 0 && active_by_taxon.find(taxon) >= 0) count = 1;
    for (int nb : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
      if (nb >= 0 && nb < n && tree.nodes[static_cast<std::size_t>(nb)].active &&
          view.parent[static_cast<std::size_t>(nb)] == node) {
        count += view.subtree_taxa_count[static_cast<std::size_t>(nb)];
      }
    }
    view.subtree_taxa_count[static_cast<std::size_t>(node)] = count;
  }

  for (int node : view.preorder) {
    if (tree.nodes[static_cast<std::size_t>(node)].taxon < 0 &&
        active_degree(tree, node) == 3) {
      view.internal_nodes.push_back(node);
    }
  }
  return view;
}

// One plan per internal vertex. The arm sizes follow the vertex's active
// neighbors in adjacency order; the arm toward the root holds every taxon
// outside the vertex's subtree.
std::vector<VertexColorPlan> build_vertex_color_plans(
    const SolverTree& tree,
    const RegraftView& view,
    const std::vector<int>& active_taxa,
    int taxon_count) {
  (void)active_taxa;
  (void)taxon_count;
  std::vector<VertexColorPlan> plans;
  plans.reserve(view.internal_nodes.size());
  const int total_active =
      view.root >= 0 &&
              view.root < static_cast<int>(view.subtree_taxa_count.size())
          ? view.subtree_taxa_count[static_cast<std::size_t>(view.root)]
          : 0;
  for (int vertex : view.internal_nodes) {
    VertexColorPlan plan;
    plan.vertex_node = vertex;
    const int vertex_subtree =
        vertex >= 0 && vertex < static_cast<int>(view.subtree_taxa_count.size())
            ? view.subtree_taxa_count[static_cast<std::size_t>(vertex)]
            : 0;
    int direction = 0;
    for (int nb : tree.nodes[static_cast<std::size_t>(vertex)].neighbors) {
      if (!tree.nodes[static_cast<std::size_t>(nb)].active) continue;
      if (direction >= 3) break;
      int branch_size = 0;
      if (view.parent[static_cast<std::size_t>(nb)] == vertex) {
        branch_size = view.subtree_taxa_count[static_cast<std::size_t>(nb)];
      } else if (view.parent[static_cast<std::size_t>(vertex)] == nb) {
        branch_size = total_active - vertex_subtree;
      }
      plan.branch_size[static_cast<std::size_t>(direction)] = branch_size;
      ++direction;
    }
    plans.push_back(std::move(plan));
  }
  return plans;
}

int direction_for_neighbor(const SolverTree& tree, int vertex, int neighbor);

// Builds the query tree of one side and maps each scored query vertex's three
// branches to arm numbers. The owner fields hold taxon ids until
// remap_regraft_jobs_to_owners.
RegraftJob prepare_regraft_job(
    const SolverTree& child_tree,
    const std::vector<int>& active_taxa,
    const std::vector<int>& context_taxa_ids,
    int taxon_count,
    std::vector<VertexColorPlan> plans) {
  RegraftJob work;
  work.child_tree = &child_tree;
  work.active_taxa_ids = &active_taxa;
  work.context_taxa_ids = &context_taxa_ids;
  work.query =
      query_tree_from_solver_tree(child_tree, active_taxa, taxon_count);
  work.plans = std::move(plans);
  work.totals.assign(work.plans.size(), VertexScores{});
  work.source_to_plan.assign(child_tree.nodes.size(), -1);
  for (int index = 0; index < static_cast<int>(work.plans.size()); ++index) {
    const int source = work.plans[static_cast<std::size_t>(index)].vertex_node;
    if (source >= 0 && source < static_cast<int>(work.source_to_plan.size())) {
      work.source_to_plan[static_cast<std::size_t>(source)] = index;
    }
  }

  work.direction_by_query_vertex.assign(
      work.query.nodes.size(), std::array<int, 3>{-1, -1, -1});
  std::vector<int> vertices_of_plan(work.plans.size(), 0);
  for (int query_vertex = 0;
       query_vertex < static_cast<int>(work.query.nodes.size());
       ++query_vertex) {
    const int source =
        work.query.nodes[static_cast<std::size_t>(query_vertex)].source_node;
    if (source < 0 || source >= static_cast<int>(work.source_to_plan.size()) ||
        work.source_to_plan[static_cast<std::size_t>(source)] < 0) {
      continue;
    }
    std::array<int, 3> branch_nodes{-1, -1, -1};
    const int child_begin =
        work.query.child_first[static_cast<std::size_t>(query_vertex)];
    const int child_count =
        work.query.child_count[static_cast<std::size_t>(query_vertex)];
    if (query_vertex == work.query.root) {
      if (child_count != 3) continue;
      branch_nodes = std::array<int, 3>{
          work.query.child_nodes[static_cast<std::size_t>(child_begin)],
          work.query.child_nodes[static_cast<std::size_t>(child_begin + 1)],
          work.query.child_nodes[static_cast<std::size_t>(child_begin + 2)]};
    } else {
      if (child_count != 2) continue;
      branch_nodes =
          vertex_arm_coloring(work.query, query_vertex)
              .branch_nodes;
    }
    for (int i = 0; i < 3; ++i) {
      const int branch_node = branch_nodes[static_cast<std::size_t>(i)];
      const int branch_source =
          branch_node >= 0 &&
                  branch_node < static_cast<int>(work.query.nodes.size())
              ? work.query.nodes[static_cast<std::size_t>(branch_node)]
                    .source_node
              : -1;
      work.direction_by_query_vertex[static_cast<std::size_t>(query_vertex)]
                                    [static_cast<std::size_t>(i)] =
          direction_for_neighbor(child_tree, source, branch_source);
    }
    ++vertices_of_plan[static_cast<std::size_t>(
        work.source_to_plan[static_cast<std::size_t>(source)])];
  }
  for (int count : vertices_of_plan) {
    work.max_vertices_per_plan = std::max(work.max_vertices_per_plan, count);
  }

  work.query_active_to_owner.assign(active_taxa.size(), -1);
  for (int active = 0; active < static_cast<int>(active_taxa.size());
       ++active) {
    work.query_active_to_owner[static_cast<std::size_t>(active)] =
        active_taxa[static_cast<std::size_t>(active)];
  }
  work.query_leaf_owner_order =
      query_leaf_owner_order(work.query, work.query_active_to_owner);
  work.context_owner_ids = context_taxa_ids;
  return work;
}

// Adds `value` to the edges in the subtree of `node`, including the edge
// above it, through the difference array `diff` over preorder positions.
template <class Count>
void range_add(std::vector<Count>& diff,
               const RegraftView& view,
               int node,
               Count value) {
  if (node < 0 || node >= static_cast<int>(view.tin.size())) return;
  diff[static_cast<std::size_t>(view.tin[static_cast<std::size_t>(node)])] +=
      value;
  diff[static_cast<std::size_t>(view.tout[static_cast<std::size_t>(node)])] -=
      value;
}

// The arm number of `neighbor` at `vertex`, or -1 when it is not an active
// neighbor.
int direction_for_neighbor(const SolverTree& tree, int vertex, int neighbor) {
  int direction = 0;
  for (int nb : tree.nodes[static_cast<std::size_t>(vertex)].neighbors) {
    if (!tree.nodes[static_cast<std::size_t>(nb)].active) continue;
    if (nb == neighbor) return direction;
    ++direction;
  }
  return -1;
}

// Replaces taxon ids by the owners of `local_by_taxon`, the live taxa of the
// scoring pass. Context taxa without an owner are dropped.
void remap_regraft_jobs_to_owners(std::vector<RegraftJob>& jobs,
                                const TaxonIndex& local_by_taxon) {
  for (RegraftJob& job : jobs) {
    if (job.active_taxa_ids != nullptr) {
      const std::vector<int>& active_taxa = *job.active_taxa_ids;
      job.query_active_to_owner.assign(active_taxa.size(), -1);
      for (int active = 0; active < static_cast<int>(active_taxa.size());
           ++active) {
        job.query_active_to_owner[static_cast<std::size_t>(active)] =
            local_by_taxon.find(
                active_taxa[static_cast<std::size_t>(active)]);
      }
      job.query_leaf_owner_order =
          query_leaf_owner_order(job.query,
                                            job.query_active_to_owner);
    }
    job.context_owner_ids.clear();
    if (job.context_taxa_ids != nullptr) {
      job.context_owner_ids.reserve(job.context_taxa_ids->size());
      for (int taxon : *job.context_taxa_ids) {
        const int owner = local_by_taxon.find(taxon);
        if (owner >= 0) job.context_owner_ids.push_back(owner);
      }
    }
  }
}

void add_vertex_scores(VertexScores& total,
                      const VertexScores& part) {
  for (int direction = 0; direction < 3; ++direction) {
    total.raw[static_cast<std::size_t>(direction)] +=
        part.raw[static_cast<std::size_t>(direction)];
  }
}

// Scores every job on genes [begin_gene, end_gene) of the scoring genes and
// adds the results to totals_by_job.
void score_regraft_jobs_for_genes(
    const SolverContext& ctx,
    std::vector<RegraftJob>& jobs,
    const OwnerProjection& owners,
    int begin_gene,
    int end_gene,
    std::vector<std::vector<VertexScores>>& totals_by_job) {
  if (jobs.empty()) return;
  const std::vector<GeneIndex>& genes = scoring_genes(ctx);
  begin_gene = std::max(0, begin_gene);
  end_gene = std::max(begin_gene,
                      std::min(end_gene, static_cast<int>(genes.size())));
  const int gene_count = end_gene - begin_gene;
  if (gene_count <= 0) return;

  // One score row per gene. Genes are scheduled dynamically and the rows are
  // summed in gene order, so the result does not depend on the worker count.
  const int chunks = gene_partition_count(gene_count);

  std::vector<RowTable<VertexScores>> partials;
  partials.reserve(jobs.size());
  std::size_t partial_count = 0;
  for (const RegraftJob& job : jobs) {
    const std::size_t job_partial_count = checked_hdt_size_mul(
        static_cast<std::size_t>(chunks), job.plans.size(),
        "regraft partial score table exceeds the addressable size domain");
    partial_count = checked_hdt_size_add(
        partial_count, job_partial_count,
        "regraft partial score tables exceed the addressable size domain");
    partials.emplace_back(static_cast<std::size_t>(chunks), job.plans.size());
  }
  // For jobs that request exact totals, each gene's largest channel sum over
  // the scored vertices.
  std::vector<std::vector<double>> magnitudes(jobs.size());
  for (std::size_t j = 0; j < jobs.size(); ++j) {
    if (jobs[j].exact_requested) {
      magnitudes[j].assign(static_cast<std::size_t>(chunks), 0.0);
    }
  }
  std::vector<GeneCountTelemetry> telemetry =
      gene_count_telemetry_slots(chunks);
  parallel_for_dynamic(0, chunks, [&](int chunk) {
    for (RowTable<VertexScores>& table : partials) {
      table.init_row(static_cast<std::size_t>(chunk));
    }
    const int begin = begin_gene + chunk;
    const int end = begin + 1;
    long long queries = 0;
    HdtTraversalStats hdt_stats;
    HdtExtractWorkspace hdt_workspace;
    // The workspace is local to the chunk, so worker threads do not keep
    // root-sized buffers.
    for (int g = begin; g < end; ++g) {
      const GeneIndex& gene = genes[static_cast<std::size_t>(g)];
      if (!gene_has_owner_quartet(gene, owners)) continue;
      GeneIndex induced_gene;
      const GeneIndex* graph_gene = &gene;
      // With the bipartition test, rejected proposals leave 1|m-1 subproblems
      // whose large side reuses its parent's genes, so many gene leaves lie
      // outside this pass's taxa. Inducing the gene on the live taxa first,
      // once for both jobs, avoids building weight-zero leaves that Contract
      // would remove. Every duplicate leaf is kept.
      if (ctx.bipartition_test_index != nullptr) {
        bool changed = false;
        induced_gene =
            induce_gene_for_owner_projection(gene, owners, changed);
        if (changed) graph_gene = &induced_gene;
      }
      for (int job_index = 0; job_index < static_cast<int>(jobs.size());
           ++job_index) {
        RegraftJob& job = jobs[static_cast<std::size_t>(job_index)];
        VertexScores* row = partials[static_cast<std::size_t>(job_index)].row(
            static_cast<std::size_t>(chunk));
        auto record_vertex = [&](int query_vertex,
                                 const std::array<std::uint8_t, 3>&
                                     branch_colors,
                                 const ChannelCounts& counts) {
              // The arms have colours 0, 1, 2 and the other side colour 3.
              // Arm c's satisfied channel pairs c with 3: W(03|12), W(13|02)
              // or W(23|01). It scores that channel minus twice the other two.
              const int source =
                  job.query.nodes[static_cast<std::size_t>(query_vertex)]
                      .source_node;
              if (source < 0 ||
                  source >= static_cast<int>(job.source_to_plan.size())) {
                return;
              }
              const int plan =
                  job.source_to_plan[static_cast<std::size_t>(source)];
              if (plan < 0) return;
              const double q12_03 = counts.q03_12;
              if (job.exact_requested) {
                double& magnitude =
                    magnitudes[static_cast<std::size_t>(job_index)]
                              [static_cast<std::size_t>(chunk)];
                magnitude = std::max(
                    magnitude, q12_03 + counts.q02_13 + counts.q01_23);
              }
              const std::array<double, 3> raw_by_color{
                  q12_03 - 2.0 * (counts.q02_13 + counts.q01_23),
                  counts.q02_13 - 2.0 * (q12_03 + counts.q01_23),
                  counts.q01_23 - 2.0 * (q12_03 + counts.q02_13)};
              for (int i = 0; i < 3; ++i) {
                const int direction =
                    job.direction_by_query_vertex
                        [static_cast<std::size_t>(query_vertex)]
                        [static_cast<std::size_t>(i)];
                const int color = branch_colors[static_cast<std::size_t>(i)];
                if (direction >= 0 && direction < 3 && color >= 0 &&
                    color < 3) {
                  VertexScores& score =
                      row[static_cast<std::size_t>(plan)];
                  score.raw[static_cast<std::size_t>(direction)] +=
                      raw_by_color[static_cast<std::size_t>(color)];
                }
              }
              ++queries;
            };
        // The side's taxa start at colour 1, which the traversal recolours
        // arm by arm; the other side keeps colour 3 in both colourings.
        HdtColorPolicy colors;
        colors.default_current = 0;
        colors.default_background = 0;
        colors.current_by_owner.assign(
            static_cast<std::size_t>(owners.owner_count), 0);
        colors.background_by_owner.assign(
            static_cast<std::size_t>(owners.owner_count), 0);
        for (int owner : job.query_leaf_owner_order) {
          if (owner < 0 || owner >= owners.owner_count) continue;
          colors.current_by_owner[static_cast<std::size_t>(owner)] = 1;
        }
        for (int owner : job.context_owner_ids) {
          if (owner < 0 || owner >= owners.owner_count) continue;
          colors.current_by_owner[static_cast<std::size_t>(owner)] = 3;
          colors.background_by_owner[static_cast<std::size_t>(owner)] = 3;
        }
        HdtComponentGraph base =
            build_base_component_graph(*graph_gene, owners, 4, colors);
        // Weight-zero absent leaves need the initial Contract before an HDT
        // is built.
        HdtState entry_state;
        if (hdt_component_graph_has_absent_leaf(base)) {
          HdtInitialStateResult initial =
              hdt_build_initial_state(std::move(base));
          hdt_account_initial_state(initial, hdt_stats);
          entry_state = std::move(initial.state);
        } else {
          // The traversal entry either contracts the graph onto the side's
          // taxa or builds its HDT, so an HDT built here would be discarded.
          entry_state.components = std::move(base);
          entry_state.hdt_deferred = true;
        }
        count_vertices(
            std::move(entry_state), job.query, job.query_leaf_owner_order,
            hdt_workspace, hdt_stats,
            [&](int query_vertex, const std::array<int, 3>&,
                const std::array<std::uint8_t, 3>& branch_colors,
                HdtState& current) {
              const std::array<double, 3> values =
                  hdt_root_counts(current);
              record_vertex(
                  query_vertex, branch_colors,
                  ChannelCounts{values[0], values[1], values[2]});
            },
            // A state is reused for a child while its size is within this
            // many Contract bounds. Any constant keeps the complexity bounds;
            // the larger one saves Contract calls on the induced genes.
            ctx.bipartition_test_index != nullptr ? 8 : 4);
      }
    }
    if (!telemetry.empty()) {
      telemetry[static_cast<std::size_t>(chunk)] = GeneCountTelemetry{
          queries, hdt_stats.recolor_changed_leaves,
          hdt_stats.recolor_marked_nodes, hdt_stats};
    }
  }, 1);

  for (int job_index = 0; job_index < static_cast<int>(jobs.size());
       ++job_index) {
    RegraftJob& job = jobs[static_cast<std::size_t>(job_index)];
    if (job.exact_requested) {
      // The counters sum products of nonnegative counts, so a channel count
      // whose true value is below 2^53 is computed exactly. A vertex's terms
      // are at most twice its channel sum s, and a row adds one term from
      // each of at most max_vertices_per_plan vertices. If that bound,
      // 2 s max_vertices_per_plan, stays below 2^52 (half the range, for
      // the rounding of s itself), the gene's counts, terms and partial sums
      // are exact integers and its rows convert without loss; otherwise the
      // job's totals are inexact.
      const auto* rows = partials[static_cast<std::size_t>(job_index)].data();
      for (int gene = 0; gene < chunks; ++gene) {
        const double magnitude =
            2.0 * job.max_vertices_per_plan *
            magnitudes[static_cast<std::size_t>(job_index)]
                      [static_cast<std::size_t>(gene)];
        for (std::size_t plan = 0; plan < job.plans.size(); ++plan) {
          for (std::size_t arm = 0; arm < 3; ++arm) {
            const double value = rows[static_cast<std::size_t>(gene) *
                                      job.plans.size() + plan].raw[arm];
            if (!(magnitude < 0x1p52) || std::trunc(value) != value) {
              job.exact_valid = false;
            } else {
              job.exact_totals[plan][arm] += static_cast<RegraftCount>(value);
            }
          }
        }
      }
    }
    add_rows_in_order(
        partials[static_cast<std::size_t>(job_index)].data(), chunks,
        jobs[static_cast<std::size_t>(job_index)].plans.size(),
        totals_by_job[static_cast<std::size_t>(job_index)].data(),
        [](VertexScores& total, const VertexScores& part) {
          add_vertex_scores(total, part);
        });
  }
  merge_gene_count_telemetry(telemetry);
  stats_add(run_stats().hdt_builds,
                        static_cast<long long>(gene_count));
}

// Scores every job over all scoring genes, restricted to `live_taxa_ids`
// (both sides), and adds the results to each job's totals.
void score_regraft_jobs(
    const SolverContext& ctx,
    std::vector<RegraftJob>& jobs,
    const std::vector<int>& live_taxa_ids,
    int taxon_count) {
  if (jobs.empty()) return;

  const std::vector<GeneIndex>& genes = scoring_genes(ctx);
  note_gene_set_use(genes);
  const TaxonIndex local_by_taxon =
      make_taxon_index(live_taxa_ids, taxon_count);
  const OwnerProjection owners =
      make_owner_projection(ctx, local_by_taxon);
  remap_regraft_jobs_to_owners(jobs, local_by_taxon);
  std::vector<std::vector<VertexScores>> totals_by_job;
  totals_by_job.reserve(jobs.size());
  for (const RegraftJob& job : jobs) {
    totals_by_job.emplace_back(job.plans.size());
  }
  score_regraft_jobs_for_genes(
      ctx, jobs, owners, 0, static_cast<int>(genes.size()), totals_by_job);

  for (int job_index = 0; job_index < static_cast<int>(jobs.size());
       ++job_index) {
    RegraftJob& job = jobs[static_cast<std::size_t>(job_index)];
    for (int index = 0; index < static_cast<int>(job.plans.size()); ++index) {
      add_vertex_scores(job.totals[static_cast<std::size_t>(index)],
                       totals_by_job[static_cast<std::size_t>(job_index)]
                                    [static_cast<std::size_t>(index)]);
    }
  }
}

// A single-taxon side attaches at its only node.
RegraftEdge singleton_regraft_edge(const RegraftView& view) {
  RegraftEdge edge;
  edge.a = view.root;
  edge.b = -1;
  edge.ties = 1;
  edge.second_rank_score = edge.rank_score;
  edge.rank_margin = std::numeric_limits<double>::infinity();
  return edge;
}

struct RankedRegraftEdge {
  RegraftEdge edge;
  int input_order = -1;
};

// Score classes in decreasing preference: +Inf, finite, -Inf, NaN.
int regraft_score_class(double score) {
  if (std::isnan(score)) return 0;
  if (std::isinf(score)) return std::signbit(score) ? 1 : 3;
  return 2;
}

// The tie-break key: the unordered endpoints, then orientation, then input
// order.
bool regraft_stable_edge_less(const RankedRegraftEdge& a,
                            const RankedRegraftEdge& b) {
  const int a_first = std::min(a.edge.a, a.edge.b);
  const int a_second = std::max(a.edge.a, a.edge.b);
  const int b_first = std::min(b.edge.a, b.edge.b);
  const int b_second = std::max(b.edge.a, b.edge.b);
  if (a_first != b_first) return a_first < b_first;
  if (a_second != b_second) return a_second < b_second;
  if (a.edge.a != b.edge.a) return a.edge.a < b.edge.a;
  if (a.edge.b != b.edge.b) return a.edge.b < b.edge.b;
  return a.input_order < b.input_order;
}

// Best first by exact score, then by the stable key. It uses no epsilon, so
// it is a strict weak order for std::sort.
bool regraft_exact_rank_less(const RankedRegraftEdge& a,
                           const RankedRegraftEdge& b) {
  const int a_class = regraft_score_class(a.edge.rank_score);
  const int b_class = regraft_score_class(b.edge.rank_score);
  if (a_class != b_class) return a_class > b_class;
  if (a_class == 2 && a.edge.rank_score != b.edge.rank_score) {
    return a.edge.rank_score > b.edge.rank_score;
  }
  return regraft_stable_edge_less(a, b);
}

// Whether `score` ties the exact best: the same class, and within
// kRegraftTieEps when finite.
bool regraft_score_in_best_band(double score, double exact_best) {
  const int score_class = regraft_score_class(score);
  const int best_class = regraft_score_class(exact_best);
  if (score_class != best_class) return false;
  if (best_class != 2) return true;
  return std::fabs(score - exact_best) <= kRegraftTieEps;
}

// How far the best score leads the next edge outside its tie band.
double regraft_score_margin(double best, double next) {
  const int best_class = regraft_score_class(best);
  const int next_class = regraft_score_class(next);
  if (best_class != next_class) return std::numeric_limits<double>::infinity();
  if (best_class != 2) return 0.0;
  const double margin = best - next;
  return std::isfinite(margin) ? margin
                               : std::numeric_limits<double>::infinity();
}

// Turns vertex arm scores into edge scores and selects the best edge. An
// edge's score sums, over the scored vertices, the score of the arm that
// contains the edge. With `allowed`, only edges with both ends marked compete;
// every edge is still reported in `edge_scores`. A child arm adds to its subtree; the arm toward the root
// adds to every edge and is subtracted from the vertex's child subtrees. One
// prefix sum over preorder then gives every edge in O(m). The exact maximum
// is found first, and the stable key picks among the edges tied with it.
RegraftEdge finish_regraft_edge(
    const SolverContext& ctx,
    const SolverTree& child_tree,
    const RegraftView& view,
    const std::vector<VertexColorPlan>& plans,
    const std::vector<VertexScores>& scores,
    std::size_t context_taxa_count,
    std::vector<RegraftEdge>* edge_scores = nullptr,
    bool force_raw = false,
    const std::vector<std::array<RegraftCount, 3>>* exact = nullptr,
    const std::vector<char>* allowed = nullptr) {
  std::vector<std::array<double, 3>> normalized_coefficients(plans.size());
  for (int index = 0; index < static_cast<int>(plans.size()); ++index) {
    std::array<double, 3> normalized =
        scores[static_cast<std::size_t>(index)].raw;
    const VertexColorPlan& plan = plans[static_cast<std::size_t>(index)];
    const double z =
        static_cast<double>(plan.branch_size[0]) *
        static_cast<double>(plan.branch_size[1]) *
        static_cast<double>(plan.branch_size[2]) *
        static_cast<double>(context_taxa_count);
    if (z > 0.0) {
      const double scale = 1.0 / z;
      for (double& value : normalized) value *= scale;
    }
    normalized_coefficients[static_cast<std::size_t>(index)] = normalized;
  }

  std::vector<double> raw_diff(view.preorder.size() + 1, 0.0);
  std::vector<double> normalized_diff(view.preorder.size() + 1, 0.0);
  std::vector<RegraftCount> exact_diff(exact ? view.preorder.size() + 1 : 0);
  RegraftCount exact_global = 0;
  double raw_global = 0.0;
  double normalized_global = 0.0;
  for (int index = 0; index < static_cast<int>(plans.size()); ++index) {
    const int vertex = plans[static_cast<std::size_t>(index)].vertex_node;
    int direction = 0;
    for (int nb : child_tree.nodes[static_cast<std::size_t>(vertex)].neighbors) {
      if (!child_tree.nodes[static_cast<std::size_t>(nb)].active) continue;
      if (direction >= 3) break;
      const double raw_value =
          scores[static_cast<std::size_t>(index)]
              .raw[static_cast<std::size_t>(direction)];
      const double normalized_value =
          normalized_coefficients[static_cast<std::size_t>(index)]
                                 [static_cast<std::size_t>(direction)];
      const RegraftCount exact_value = exact
          ? (*exact)[static_cast<std::size_t>(index)][static_cast<std::size_t>(direction)]
          : 0;
      ++direction;
      if (view.parent[static_cast<std::size_t>(nb)] == vertex) {
        if (exact) range_add(exact_diff, view, nb, exact_value);
        range_add(raw_diff, view, nb, raw_value);
        range_add(normalized_diff, view, nb, normalized_value);
      } else if (view.parent[static_cast<std::size_t>(vertex)] == nb) {
        if (exact) exact_global += exact_value;
        raw_global += raw_value;
        normalized_global += normalized_value;
        for (int child :
             child_tree.nodes[static_cast<std::size_t>(vertex)].neighbors) {
          if (child_tree.nodes[static_cast<std::size_t>(child)].active &&
              view.parent[static_cast<std::size_t>(child)] == vertex) {
            if (exact) range_add(exact_diff, view, child, -exact_value);
            range_add(raw_diff, view, child, -raw_value);
            range_add(normalized_diff, view, child, -normalized_value);
          }
        }
      }
    }
  }

  std::vector<double> raw_prefix(view.preorder.size(), 0.0);
  std::vector<double> normalized_prefix(view.preorder.size(), 0.0);
  RegraftCount exact_running = 0;
  double raw_running = 0.0;
  double normalized_running = 0.0;
  for (int i = 0; i < static_cast<int>(view.preorder.size()); ++i) {
    if (exact) {
      exact_running += exact_diff[static_cast<std::size_t>(i)];
      exact_diff[static_cast<std::size_t>(i)] = exact_running + exact_global;
    }
    raw_running += raw_diff[static_cast<std::size_t>(i)];
    normalized_running += normalized_diff[static_cast<std::size_t>(i)];
    raw_prefix[static_cast<std::size_t>(i)] = raw_running + raw_global;
    normalized_prefix[static_cast<std::size_t>(i)] =
        normalized_running + normalized_global;
  }

  std::vector<RankedRegraftEdge> ranked;
  ranked.reserve(view.edge_child.size());
  for (int child : view.edge_child) {
    RegraftEdge edge;
    edge.a = view.parent[static_cast<std::size_t>(child)];
    edge.b = child;
    edge.raw_score =
        raw_prefix[static_cast<std::size_t>(
            view.tin[static_cast<std::size_t>(child)])];
    if (exact) {
      edge.exact_raw_score = exact_diff[static_cast<std::size_t>(
          view.tin[static_cast<std::size_t>(child)])];
      edge.exact_raw_score_valid = true;
    }
    edge.normalized_score =
        normalized_prefix[static_cast<std::size_t>(
            view.tin[static_cast<std::size_t>(child)])];
    switch (force_raw ? RegraftScore::Raw : ctx.options.regraft_score) {
      case RegraftScore::Raw:
        edge.rank_score = edge.raw_score;
        break;
      case RegraftScore::Normalized:
        edge.rank_score = edge.normalized_score;
        break;
    }
    if (edge_scores) edge_scores->push_back(edge);
    if (allowed != nullptr && ((*allowed)[static_cast<std::size_t>(edge.a)] == 0 ||
                               (*allowed)[static_cast<std::size_t>(edge.b)] == 0)) {
      continue;
    }
    ranked.push_back(
        RankedRegraftEdge{edge, static_cast<int>(ranked.size())});
  }
  std::sort(ranked.begin(), ranked.end(), regraft_exact_rank_less);
  if (ranked.empty()) return RegraftEdge{};
  const double exact_best = ranked.front().edge.rank_score;
  int best_ties = 0;
  while (best_ties < static_cast<int>(ranked.size()) &&
         regraft_score_in_best_band(
             ranked[static_cast<std::size_t>(best_ties)].edge.rank_score,
             exact_best)) {
    ++best_ties;
  }
  if (best_ties > 1) {
    stats_add(run_stats().regraft_edge_ties,
                          best_ties - 1);
  }
  const auto selected = std::min_element(
      ranked.begin(), ranked.begin() + best_ties, regraft_stable_edge_less);
  RegraftEdge best = selected->edge;
  best.rank = 0;
  best.ties = best_ties;
  if (best_ties < static_cast<int>(ranked.size())) {
    best.second_rank_score =
        ranked[static_cast<std::size_t>(best_ties)].edge.rank_score;
    best.rank_margin =
        regraft_score_margin(exact_best, best.second_rank_score);
  } else {
    best.second_rank_score = best.rank_score;
    best.rank_margin = std::numeric_limits<double>::infinity();
  }
  return best;
}

}  // namespace

// The edge of `child_tree` where `context_taxa_ids` attaches. A side without
// an internal vertex attaches at its first edge, or at its only node.
RegraftEdge choose_regraft_edge(
    const SolverContext& ctx,
    const SolverTree& child_tree,
    const std::vector<int>& active_taxa,
    const std::vector<int>& context_taxa_ids) {
  const auto total_started = std::chrono::steady_clock::now();
  stats_add(run_stats().regraft_calls);

  const int taxon_count = static_cast<int>(context_taxa(ctx).size());
  RegraftView view =
      build_regraft_view(child_tree, active_taxa, taxon_count);
  if (view.root < 0) return RegraftEdge{};
  if (view.edge_child.empty()) {
    stats_add_time(run_stats().regraft_ns,
                        std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - total_started));
    return singleton_regraft_edge(view);
  }
  if (view.internal_nodes.empty()) {
    RegraftEdge edge;
    const int child = view.edge_child.front();
    edge.a = view.parent[static_cast<std::size_t>(child)];
    edge.b = child;
    edge.ties = 1;
    edge.second_rank_score = edge.rank_score;
    edge.rank_margin = std::numeric_limits<double>::infinity();
    stats_add_time(run_stats().regraft_ns,
                        std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::steady_clock::now() - total_started));
    return edge;
  }

  std::vector<VertexColorPlan> plans =
      build_vertex_color_plans(child_tree, view, active_taxa, taxon_count);
  stats_add(run_stats().regraft_vertices,
                        static_cast<long long>(plans.size()));

  std::vector<RegraftJob> jobs;
  jobs.push_back(prepare_regraft_job(child_tree, active_taxa,
                                          context_taxa_ids, taxon_count,
                                          std::move(plans)));
  std::vector<int> live_taxa;
  live_taxa.reserve(active_taxa.size() + context_taxa_ids.size());
  live_taxa.insert(live_taxa.end(), active_taxa.begin(), active_taxa.end());
  live_taxa.insert(live_taxa.end(), context_taxa_ids.begin(),
                   context_taxa_ids.end());
  score_regraft_jobs(ctx, jobs, live_taxa, taxon_count);

  const RegraftJob& job = jobs.front();
  RegraftEdge best = finish_regraft_edge(
      ctx, child_tree, view, job.plans, job.totals, context_taxa_ids.size());
  stats_add_time(run_stats().regraft_ns,
                      std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - total_started));
  return best;
}

// For the two solved sides of a bipartition, the edge of each where the
// other side attaches. Both sides are scored in one pass over the genes.
std::pair<RegraftEdge, RegraftEdge> choose_regraft_edges(
    const SolverContext& ctx,
    const SolverTree& left_tree,
    const SolverTree& right_tree,
    const std::vector<int>& left_taxa,
    const std::vector<int>& right_taxa,
    std::vector<RegraftEdge>* left_scores,
    std::vector<RegraftEdge>* right_scores,
    bool force_raw,
    bool exact_scores,
    const std::vector<char>* left_allowed,
    const std::vector<char>* right_allowed) {
  if (left_scores) left_scores->clear();
  if (right_scores) right_scores->clear();
  const auto total_started = std::chrono::steady_clock::now();
  stats_add(run_stats().regraft_calls, 2);

  const int taxon_count = static_cast<int>(context_taxa(ctx).size());
  RegraftView left_view =
      build_regraft_view(left_tree, left_taxa, taxon_count);
  RegraftView right_view =
      build_regraft_view(right_tree, right_taxa, taxon_count);

  RegraftEdge left_edge;
  RegraftEdge right_edge;
  bool left_needs_score = false;
  bool right_needs_score = false;

  auto initialize_edge = [](const SolverTree& tree,
                            const RegraftView& view,
                            RegraftEdge& edge) -> bool {
    (void)tree;
    if (view.root < 0) {
      edge = RegraftEdge{};
      return false;
    }
    if (view.edge_child.empty()) {
      edge = singleton_regraft_edge(view);
      return false;
    }
    if (view.internal_nodes.empty()) {
      const int child = view.edge_child.front();
      edge.a = view.parent[static_cast<std::size_t>(child)];
      edge.b = child;
      edge.ties = 1;
      edge.second_rank_score = edge.rank_score;
      edge.rank_margin = std::numeric_limits<double>::infinity();
      return false;
    }
    return true;
  };

  left_needs_score = initialize_edge(left_tree, left_view, left_edge);
  right_needs_score = initialize_edge(right_tree, right_view, right_edge);

  std::vector<RegraftJob> jobs;
  jobs.reserve(2);
  int left_job = -1;
  int right_job = -1;
  if (left_needs_score) {
    std::vector<VertexColorPlan> plans =
        build_vertex_color_plans(left_tree, left_view, left_taxa, taxon_count);
    stats_add(run_stats().regraft_vertices,
                          static_cast<long long>(plans.size()));
    left_job = static_cast<int>(jobs.size());
    jobs.push_back(prepare_regraft_job(left_tree, left_taxa, right_taxa,
                                            taxon_count, std::move(plans)));
  }
  if (right_needs_score) {
    std::vector<VertexColorPlan> plans =
        build_vertex_color_plans(right_tree, right_view, right_taxa,
                                 taxon_count);
    stats_add(run_stats().regraft_vertices,
                          static_cast<long long>(plans.size()));
    right_job = static_cast<int>(jobs.size());
    jobs.push_back(prepare_regraft_job(right_tree, right_taxa, left_taxa,
                                            taxon_count, std::move(plans)));
  }

  for (auto& job : jobs) {
    job.exact_requested = exact_scores;
    job.exact_valid = exact_scores;
    if (exact_scores) job.exact_totals.resize(job.plans.size());
  }
  if (!left_needs_score) left_edge.exact_raw_score_valid = exact_scores;
  if (!right_needs_score) right_edge.exact_raw_score_valid = exact_scores;

  if (!jobs.empty()) {
    long long vertex_count = 0;
    for (const RegraftJob& job : jobs) {
      vertex_count += static_cast<long long>(job.plans.size());
    }
    progress_note_stage("regraft L=" + std::to_string(left_taxa.size()) +
                        " R=" + std::to_string(right_taxa.size()) +
                        " vertices=" + std::to_string(vertex_count));
    std::vector<int> live_taxa;
    live_taxa.reserve(left_taxa.size() + right_taxa.size());
    live_taxa.insert(live_taxa.end(), left_taxa.begin(), left_taxa.end());
    live_taxa.insert(live_taxa.end(), right_taxa.begin(), right_taxa.end());
    score_regraft_jobs(ctx, jobs, live_taxa, taxon_count);
  }

  if (left_job >= 0) {
    const RegraftJob& job = jobs[static_cast<std::size_t>(left_job)];
    left_edge = finish_regraft_edge(
        ctx, left_tree, left_view, job.plans, job.totals, right_taxa.size(),
        left_scores, force_raw, job.exact_valid ? &job.exact_totals : nullptr,
        left_allowed);
  }
  if (right_job >= 0) {
    const RegraftJob& job = jobs[static_cast<std::size_t>(right_job)];
    right_edge = finish_regraft_edge(
        ctx, right_tree, right_view, job.plans, job.totals, left_taxa.size(),
        right_scores, force_raw, job.exact_valid ? &job.exact_totals : nullptr,
        right_allowed);
  }

  stats_add_time(run_stats().regraft_ns,
                      std::chrono::duration_cast<std::chrono::nanoseconds>(
                          std::chrono::steady_clock::now() - total_started));
  // A two-taxon side can contain degree-two rooting vertices. All its
  // edges have zero sigma, and callers inspecting an old join need them all.
  const auto trivial_scores = [](const RegraftView& view, const RegraftEdge& edge,
                                 std::vector<RegraftEdge>* output) {
    if (!output) return;
    if (view.edge_child.empty()) {
      output->push_back(edge);
    } else {
      for (int child : view.edge_child) {
        RegraftEdge value = edge;
        value.a = view.parent[static_cast<std::size_t>(child)];
        value.b = child;
        output->push_back(value);
      }
    }
  };
  if (!left_needs_score) trivial_scores(left_view, left_edge, left_scores);
  if (!right_needs_score) trivial_scores(right_view, right_edge, right_scores);
  return {left_edge, right_edge};
}

}  // namespace qraft::internal
