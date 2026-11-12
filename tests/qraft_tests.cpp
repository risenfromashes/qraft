#include "qraft/qraft.hpp"
#include "support/rf.hpp"
#include "qraft/tree.hpp"

#include "internal.hpp"
#include "solver/solver.hpp"
#include "core/selection.hpp"
#include "support/small_stack.hpp"
#include "divide/bipartition_test.hpp"
#include "counting/counters.hpp"
#include "counting/hdt.hpp"
#include "counting/extract_contract.hpp"
#include "counting/count.hpp"
#include "counting/units.hpp"
#include "counting/query_taxa.hpp"
#include "mre/exact_split_ids.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

[[noreturn]] void fail(const std::string& message) {
  throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
  if (!condition) fail(message);
}

// The bipartition-test index distance between the leaves of two taxa.
int test_index_distance(const qraft::internal::GeneDistanceIndex& index,
                        const qraft::internal::GeneIndex& gene, int a, int b) {
  return index.distance(index.leaf_rank(a, gene.taxon_to_node[a]),
                        index.leaf_rank(b, gene.taxon_to_node[b]));
}

bool finite_close(double actual,
                  double expected,
                  double absolute_tolerance = 1e-9,
                  double relative_tolerance = 0.0) {
  return std::isfinite(actual) && std::isfinite(expected) &&
         std::fabs(actual - expected) <=
             absolute_tolerance + relative_tolerance * std::fabs(expected);
}

void require_finite(double value, const std::string& message) {
  if (!std::isfinite(value)) {
    fail(message + ": non-finite value=" + std::to_string(value));
  }
}

void require_close(double actual,
                   double expected,
                   const std::string& message,
                   double absolute_tolerance = 1e-9,
                   double relative_tolerance = 0.0) {
  if (!finite_close(actual, expected, absolute_tolerance,
                    relative_tolerance)) {
    fail(message + ": actual=" + std::to_string(actual) +
         " expected=" + std::to_string(expected));
  }
}

std::string read_file(const std::filesystem::path& path) {
  std::ifstream in(path);
  require(static_cast<bool>(in), "failed to open " + path.string());
  return std::string((std::istreambuf_iterator<char>(in)),
                     std::istreambuf_iterator<char>());
}

std::vector<qraft::GeneTree> parse_gene_tree_file(
    const std::filesystem::path& path) {
  std::ifstream in(path);
  require(static_cast<bool>(in), "failed to open " + path.string());
  std::vector<qraft::GeneTree> genes;
  std::string line;
  while (std::getline(in, line)) {
    if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
    genes.push_back(qraft::parse_newick(line));
  }
  require(!genes.empty(), "gene file contained no trees: " + path.string());
  return genes;
}

std::vector<std::string> sorted_labels(const qraft::Tree& tree) {
  std::vector<std::string> out = tree.leaf_labels();
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<qraft::GeneTree> small_gene_set() {
  return {
      qraft::parse_newick("((A,B),(C,D));"),
      qraft::parse_newick("((A,B),(C,D));"),
      qraft::parse_newick("((A,C),(B,D));"),
      qraft::parse_newick("((A,B),(E,F));"),
      qraft::parse_newick("((C,D),(E,F));"),
      qraft::parse_newick("((A,E),(B,F));"),
  };
}

void add_adjacency_edge(std::vector<std::vector<int>>& adjacency, int a, int b) {
  const std::size_t need = static_cast<std::size_t>(std::max(a, b)) + 1u;
  if (adjacency.size() < need) adjacency.resize(need);
  adjacency[static_cast<std::size_t>(a)].push_back(b);
  adjacency[static_cast<std::size_t>(b)].push_back(a);
}

void finish_gene_adjacency(qraft::internal::GeneIndex& gene,
                           std::vector<std::vector<int>>& adjacency) {
  adjacency.resize(gene.nodes.size());
  qraft::internal::set_gene_adjacency(gene, adjacency);
}

qraft::internal::GeneIndex small_internal_gene_index() {
  qraft::internal::GeneIndex gene;
  gene.nodes.resize(7);
  std::vector<std::vector<int>> gene_adjacency;
  auto add_edge = [&](int a, int b) {
    add_adjacency_edge(gene_adjacency, a, b);
  };
  add_edge(0, 1);
  add_edge(0, 2);
  add_edge(1, 3);
  add_edge(1, 4);
  add_edge(2, 5);
  add_edge(2, 6);
  for (int taxon = 0; taxon < 4; ++taxon) {
    gene.nodes[static_cast<std::size_t>(taxon + 3)].set_taxon(taxon);
  }
  gene.root = 0;
  gene.reference_root = 3;
  gene.taxon_to_node = {3, 4, 5, 6};
  finish_gene_adjacency(gene, gene_adjacency);
  gene.root_traversal = qraft::internal::traversal(gene, gene.root);
  return gene;
}

qraft::internal::GeneIndex three_leaf_star_gene_index() {
  qraft::internal::GeneIndex gene;
  gene.nodes.resize(4);
  std::vector<std::vector<int>> gene_adjacency;
  auto add_edge = [&](int a, int b) {
    add_adjacency_edge(gene_adjacency, a, b);
  };
  add_edge(0, 1);
  add_edge(0, 2);
  add_edge(0, 3);
  for (int taxon = 0; taxon < 3; ++taxon) {
    gene.nodes[static_cast<std::size_t>(taxon + 1)].set_taxon(taxon);
  }
  gene.root = 0;
  gene.reference_root = 1;
  gene.taxon_to_node = {1, 2, 3};
  finish_gene_adjacency(gene, gene_adjacency);
  gene.root_traversal = qraft::internal::traversal(gene, gene.root);
  return gene;
}

qraft::internal::GeneIndex star_gene_index(int leaf_count) {
  qraft::internal::GeneIndex gene;
  gene.nodes.resize(static_cast<std::size_t>(leaf_count + 1));
  std::vector<std::vector<int>> gene_adjacency;
  auto add_edge = [&](int a, int b) {
    add_adjacency_edge(gene_adjacency, a, b);
  };
  for (int taxon = 0; taxon < leaf_count; ++taxon) {
    add_edge(0, taxon + 1);
    gene.nodes[static_cast<std::size_t>(taxon + 1)].set_taxon(taxon);
  }
  gene.root = 0;
  gene.reference_root = leaf_count > 0 ? 1 : 0;
  gene.taxon_to_node.resize(static_cast<std::size_t>(leaf_count));
  for (int taxon = 0; taxon < leaf_count; ++taxon) {
    gene.taxon_to_node[static_cast<std::size_t>(taxon)] = taxon + 1;
  }
  finish_gene_adjacency(gene, gene_adjacency);
  gene.root_traversal = qraft::internal::traversal(gene, gene.root);
  return gene;
}

qraft::internal::GeneIndex caterpillar_gene_index(int leaf_count) {
  qraft::internal::GeneIndex gene;
  if (leaf_count <= 2) return star_gene_index(std::max(leaf_count, 1));
  const int internal_count = leaf_count - 2;
  gene.nodes.resize(static_cast<std::size_t>(internal_count + leaf_count));
  auto leaf_node = [&](int taxon) { return internal_count + taxon; };
  std::vector<std::vector<int>> gene_adjacency;
  auto add_edge = [&](int a, int b) {
    add_adjacency_edge(gene_adjacency, a, b);
  };
  add_edge(0, leaf_node(0));
  add_edge(0, leaf_node(1));
  if (internal_count == 1) {
    add_edge(0, leaf_node(2));
  } else {
    add_edge(0, 1);
    for (int i = 1; i < internal_count - 1; ++i) {
      add_edge(i, leaf_node(i + 1));
      add_edge(i, i + 1);
    }
    add_edge(internal_count - 1, leaf_node(leaf_count - 2));
    add_edge(internal_count - 1, leaf_node(leaf_count - 1));
  }
  for (int taxon = 0; taxon < leaf_count; ++taxon) {
    gene.nodes[static_cast<std::size_t>(leaf_node(taxon))].set_taxon(taxon);
  }
  gene.root = 0;
  gene.reference_root = leaf_node(0);
  gene.taxon_to_node.resize(static_cast<std::size_t>(leaf_count));
  for (int taxon = 0; taxon < leaf_count; ++taxon) {
    gene.taxon_to_node[static_cast<std::size_t>(taxon)] = leaf_node(taxon);
  }
  finish_gene_adjacency(gene, gene_adjacency);
  gene.root_traversal = qraft::internal::traversal(gene, gene.root);
  return gene;
}

// The "broom": a root of degree `arms` whose arms are resolved cherries. The
// HDT must fold the high-degree root into a binary decomposition exactly.
qraft::internal::GeneIndex broom_gene_index(int arms) {
  qraft::internal::GeneIndex gene;
  const int leaf_count = 2 * arms;
  const int node_count = 1 + arms + leaf_count;
  gene.nodes.resize(static_cast<std::size_t>(node_count));
  std::vector<std::vector<int>> gene_adjacency;
  auto add_edge = [&](int a, int b) {
    add_adjacency_edge(gene_adjacency, a, b);
  };
  auto leaf_node = [&](int taxon) { return 1 + arms + taxon; };
  for (int i = 0; i < arms; ++i) {
    const int center = 1 + i;
    add_edge(0, center);
    add_edge(center, leaf_node(2 * i));
    add_edge(center, leaf_node(2 * i + 1));
  }
  for (int taxon = 0; taxon < leaf_count; ++taxon) {
    gene.nodes[static_cast<std::size_t>(leaf_node(taxon))].set_taxon(taxon);
  }
  gene.root = 0;
  gene.reference_root = leaf_node(0);
  gene.taxon_to_node.resize(static_cast<std::size_t>(leaf_count));
  for (int taxon = 0; taxon < leaf_count; ++taxon) {
    gene.taxon_to_node[static_cast<std::size_t>(taxon)] = leaf_node(taxon);
  }
  finish_gene_adjacency(gene, gene_adjacency);
  gene.root_traversal = qraft::internal::traversal(gene, gene.root);
  return gene;
}

// Random rooted gene tree with `leaf_count` leaves that mixes resolved vertices
// and polytomies: the root has three children and every other internal node
// 2..5. Deterministic in `seed`.
qraft::internal::GeneIndex random_arbitrary_gene_index(int leaf_count,
                                                      unsigned seed) {
  qraft::internal::GeneIndex gene;
  gene.nodes.clear();
  auto add_node = [&]() {
    gene.nodes.emplace_back();
    return static_cast<int>(gene.nodes.size()) - 1;
  };
  std::vector<std::vector<int>> gene_adjacency;
  auto add_edge = [&](int a, int b) {
    add_adjacency_edge(gene_adjacency, a, b);
  };
  gene.taxon_to_node.resize(static_cast<std::size_t>(leaf_count));
  std::vector<int> roots;
  for (int taxon = 0; taxon < leaf_count; ++taxon) {
    const int node = add_node();
    gene.nodes[static_cast<std::size_t>(node)].set_taxon(taxon);
    gene.taxon_to_node[static_cast<std::size_t>(taxon)] = node;
    roots.push_back(node);
  }
  std::mt19937 rng(seed);
  // Keep >=3 roots until the final grouping so the tree's root ends up degree
  // >=3 (real gene trees have no degree-2 nodes; index_gene suppresses them).
  while (roots.size() > 3) {
    std::shuffle(roots.begin(), roots.end(), rng);
    const int max_k =
        std::min<int>(5, static_cast<int>(roots.size()) - 2);
    std::uniform_int_distribution<int> pick(2, std::max(2, max_k));
    const int k = std::min<int>(pick(rng), static_cast<int>(roots.size()) - 2);
    const int parent = add_node();
    for (int i = 0; i < k; ++i) {
      add_edge(parent, roots.back());
      roots.pop_back();
    }
    roots.push_back(parent);
  }
  if (roots.size() > 1) {
    std::shuffle(roots.begin(), roots.end(), rng);
    const int parent = add_node();
    while (!roots.empty()) {
      add_edge(parent, roots.back());
      roots.pop_back();
    }
    roots.push_back(parent);
  }
  gene.root = roots.empty() ? 0 : roots.front();
  gene.reference_root = gene.taxon_to_node.empty()
                            ? gene.root
                            : gene.taxon_to_node[0];
  finish_gene_adjacency(gene, gene_adjacency);
  gene.root_traversal = qraft::internal::traversal(gene, gene.root);
  return gene;
}

qraft::internal::GeneIndex nested_polytomy_gene_index() {
  qraft::internal::GeneIndex gene;
  gene.nodes.resize(13);
  std::vector<std::vector<int>> gene_adjacency;
  auto add_edge = [&](int a, int b) {
    add_adjacency_edge(gene_adjacency, a, b);
  };
  // Adjacent degree-four centers 0/1, plus a nested degree-four center 3
  // reached through resolved degree-three vertex 2.
  add_edge(0, 1);
  add_edge(0, 2);
  add_edge(0, 4);
  add_edge(0, 5);
  add_edge(1, 6);
  add_edge(1, 7);
  add_edge(1, 8);
  add_edge(2, 3);
  add_edge(2, 9);
  add_edge(3, 10);
  add_edge(3, 11);
  add_edge(3, 12);
  gene.taxon_to_node.resize(9);
  for (int taxon = 0; taxon < 9; ++taxon) {
    const int node = taxon + 4;
    gene.nodes[static_cast<std::size_t>(node)].set_taxon(taxon);
    gene.taxon_to_node[static_cast<std::size_t>(taxon)] = node;
  }
  gene.root = 0;
  gene.reference_root = 4;
  finish_gene_adjacency(gene, gene_adjacency);
  gene.root_traversal = qraft::internal::traversal(gene, gene.root);
  return gene;
}

qraft::internal::GeneIndex alternating_corridor_gene_index(int segments) {
  using qraft::internal::GeneIndex;
  require(segments >= 2, "alternating corridor needs at least two segments");
  GeneIndex gene;
  auto add_node = [&]() {
    gene.nodes.emplace_back();
    return static_cast<int>(gene.nodes.size()) - 1;
  };
  std::vector<std::vector<int>> gene_adjacency;
  auto add_edge = [&](int a, int b) {
    add_adjacency_edge(gene_adjacency, a, b);
  };
  std::vector<int> spine(static_cast<std::size_t>(segments));
  for (int& node : spine) node = add_node();
  for (int i = 1; i < segments; ++i) {
    add_edge(spine[static_cast<std::size_t>(i - 1)],
             spine[static_cast<std::size_t>(i)]);
  }
  auto leaf = [&]() {
    const int node = add_node();
    const int taxon = static_cast<int>(gene.taxon_to_node.size());
    gene.nodes[static_cast<std::size_t>(node)].set_taxon(taxon);
    gene.taxon_to_node.push_back(node);
    return node;
  };
  add_edge(spine.front(), leaf());
  add_edge(spine.front(), leaf());
  add_edge(spine.back(), leaf());
  add_edge(spine.back(), leaf());
  for (int i = 1; i + 1 < segments; ++i) {
    if ((i & 1) == 0) {
      add_edge(spine[static_cast<std::size_t>(i)], leaf());
    } else {
      const int cherry = add_node();
      add_edge(spine[static_cast<std::size_t>(i)], cherry);
      add_edge(cherry, leaf());
      add_edge(cherry, leaf());
    }
  }
  gene.root = spine.front();
  gene.reference_root = gene.taxon_to_node.front();
  finish_gene_adjacency(gene, gene_adjacency);
  gene.root_traversal = qraft::internal::traversal(gene, gene.root);
  return gene;
}

qraft::internal::GeneIndex balanced_gene_index(int leaf_count) {
  using qraft::internal::GeneIndex;
  require(leaf_count >= 3, "balanced gene needs at least three leaves");
  GeneIndex gene;
  auto add_node = [&]() {
    gene.nodes.emplace_back();
    return static_cast<int>(gene.nodes.size()) - 1;
  };
  std::vector<std::vector<int>> gene_adjacency;
  auto add_edge = [&](int a, int b) {
    add_adjacency_edge(gene_adjacency, a, b);
  };
  auto build = [&](auto&& self, int count) -> int {
    if (count == 1) {
      const int leaf = add_node();
      const int taxon = static_cast<int>(gene.taxon_to_node.size());
      gene.nodes[static_cast<std::size_t>(leaf)].set_taxon(taxon);
      gene.taxon_to_node.push_back(leaf);
      return leaf;
    }
    const int parent = add_node();
    const int left = self(self, count / 2);
    const int right = self(self, count - count / 2);
    add_edge(parent, left);
    add_edge(parent, right);
    return parent;
  };
  const int root = add_node();
  const int first = leaf_count / 4;
  const int second = leaf_count / 4;
  const int third = leaf_count - first - second;
  for (int count : {first, second, third}) {
    add_edge(root, build(build, count));
  }
  gene.root = root;
  gene.reference_root = gene.taxon_to_node.front();
  finish_gene_adjacency(gene, gene_adjacency);
  gene.root_traversal = qraft::internal::traversal(gene, gene.root);
  return gene;
}

void bit_reverse_gene_taxa(qraft::internal::GeneIndex& gene) {
  const int leaf_count = static_cast<int>(gene.taxon_to_node.size());
  require(leaf_count > 0 && (leaf_count & (leaf_count - 1)) == 0,
          "bit-reversal fixture requires a power-of-two leaf count");
  int bits = 0;
  while ((1 << bits) < leaf_count) ++bits;
  std::vector<int> remapped(static_cast<std::size_t>(leaf_count), -1);
  for (int taxon = 0; taxon < leaf_count; ++taxon) {
    int value = taxon;
    int reversed = 0;
    for (int bit = 0; bit < bits; ++bit) {
      reversed = (reversed << 1) | (value & 1);
      value >>= 1;
    }
    const int node = gene.taxon_to_node[static_cast<std::size_t>(taxon)];
    gene.nodes[static_cast<std::size_t>(node)].set_taxon(reversed);
    remapped[static_cast<std::size_t>(reversed)] = node;
  }
  gene.taxon_to_node = std::move(remapped);
}

qraft::internal::GeneIndex random_linear_gene_index(int leaf_count,
                                                   unsigned seed) {
  using qraft::internal::GeneIndex;
  require(leaf_count >= 3, "random linear gene needs at least three leaves");
  GeneIndex gene;
  auto add_node = [&]() {
    gene.nodes.emplace_back();
    return static_cast<int>(gene.nodes.size()) - 1;
  };
  std::vector<std::vector<int>> gene_adjacency;
  auto add_edge = [&](int a, int b) {
    add_adjacency_edge(gene_adjacency, a, b);
  };
  std::vector<int> frontier;
  frontier.reserve(static_cast<std::size_t>(leaf_count));
  gene.taxon_to_node.resize(static_cast<std::size_t>(leaf_count), -1);
  for (int taxon = 0; taxon < leaf_count; ++taxon) {
    const int leaf = add_node();
    gene.nodes[static_cast<std::size_t>(leaf)].set_taxon(taxon);
    gene.taxon_to_node[static_cast<std::size_t>(taxon)] = leaf;
    frontier.push_back(leaf);
  }
  std::mt19937 rng(seed);
  std::shuffle(frontier.begin(), frontier.end(), rng);
  while (frontier.size() > 10u) {
    std::vector<int> next;
    next.reserve((frontier.size() + 1u) / 2u);
    std::size_t position = 0;
    while (position < frontier.size()) {
      const std::size_t remaining = frontier.size() - position;
      int degree = remaining <= 5u
                       ? static_cast<int>(remaining)
                       : 2 + static_cast<int>(rng() % 4u);
      if (remaining > static_cast<std::size_t>(degree) &&
          remaining - static_cast<std::size_t>(degree) == 1u) {
        --degree;
      }
      require(degree >= 2, "random linear grouping produced a unary node");
      const int parent = add_node();
      for (int i = 0; i < degree; ++i) {
        add_edge(parent, frontier[position++]);
      }
      next.push_back(parent);
    }
    std::shuffle(next.begin(), next.end(), rng);
    frontier = std::move(next);
  }
  const int root = add_node();
  for (int child : frontier) add_edge(root, child);
  gene.root = root;
  gene.reference_root = gene.taxon_to_node.front();
  finish_gene_adjacency(gene, gene_adjacency);
  gene.root_traversal = qraft::internal::traversal(gene, gene.root);
  return gene;
}

int independent_quartet_channel(int a, int b, int c, int d, int palette) {
  std::array<int, 2> left{std::min(a, b), std::max(a, b)};
  std::array<int, 2> right{std::min(c, d), std::max(c, d)};
  if (right < left) std::swap(left, right);
  if (palette == 2) {
    if (left == std::array<int, 2>{0, 0} &&
        right == std::array<int, 2>{1, 1}) {
      return 0;
    }
    if (left == std::array<int, 2>{0, 1} &&
        right == std::array<int, 2>{0, 1}) {
      return 1;
    }
    return -1;
  }
  if (left == std::array<int, 2>{0, 1} &&
      right == std::array<int, 2>{2, 3}) {
    return 0;
  }
  if (left == std::array<int, 2>{0, 2} &&
      right == std::array<int, 2>{1, 3}) {
    return 1;
  }
  if (left == std::array<int, 2>{0, 3} &&
      right == std::array<int, 2>{1, 2}) {
    return 2;
  }
  return -1;
}

std::array<double, 3> independent_gene_quartet_counts(
    const qraft::internal::GeneIndex& gene,
    const qraft::internal::OwnerProjection& owners,
    int palette,
    const std::vector<std::uint8_t>& owner_colors) {
  struct Leaf {
    int node = -1;
    int owner = -1;
    double weight = 0.0;
  };
  std::vector<Leaf> leaves;
  for (int node = 0; node < static_cast<int>(gene.nodes.size()); ++node) {
    const auto& item = gene.nodes[static_cast<std::size_t>(node)];
    const int owner = qraft::internal::owner_for_gene_leaf(item, owners);
    if (owner < 0 || owner >= owners.owner_count ||
        owner >= static_cast<int>(owner_colors.size()) ||
        (!owners.live_owner.empty() &&
         (owner >= static_cast<int>(owners.live_owner.size()) ||
          owners.live_owner[static_cast<std::size_t>(owner)] == 0))) {
      continue;
    }
    require_finite(qraft::internal::gene_leaf_mass(gene, node), "independent quartet oracle leaf weight");
    leaves.push_back(Leaf{node, owner, qraft::internal::gene_leaf_mass(gene, node)});
  }

  const int node_count = static_cast<int>(gene.nodes.size());
  const int leaf_count = static_cast<int>(leaves.size());
  const std::size_t path_words =
      (static_cast<std::size_t>(node_count) + 63u) / 64u;
  const std::size_t path_stride = path_words;
  std::vector<std::uint64_t> pair_paths(
      static_cast<std::size_t>(leaf_count) *
          static_cast<std::size_t>(leaf_count) * path_stride,
      0);
  std::vector<int> parent(static_cast<std::size_t>(node_count), -2);
  std::vector<int> pending;
  pending.reserve(static_cast<std::size_t>(node_count));
  for (int first = 0; first < leaf_count; ++first) {
    for (int second = first + 1; second < leaf_count; ++second) {
      std::fill(parent.begin(), parent.end(), -2);
      pending.clear();
      const int source = leaves[static_cast<std::size_t>(first)].node;
      const int target = leaves[static_cast<std::size_t>(second)].node;
      parent[static_cast<std::size_t>(source)] = -1;
      pending.push_back(source);
      for (std::size_t head = 0;
           head < pending.size() &&
           parent[static_cast<std::size_t>(target)] == -2;
           ++head) {
        const int node = pending[head];
        for (int next : qraft::internal::gene_neighbors(gene, node)) {
          if (parent[static_cast<std::size_t>(next)] != -2) continue;
          parent[static_cast<std::size_t>(next)] = node;
          pending.push_back(next);
        }
      }
      require(parent[static_cast<std::size_t>(target)] != -2,
              "independent quartet oracle found a disconnected gene tree");
      std::uint64_t* bits =
          pair_paths.data() +
          (static_cast<std::size_t>(first) *
               static_cast<std::size_t>(leaf_count) +
           static_cast<std::size_t>(second)) *
              path_stride;
      for (int node = target; node >= 0;
           node = parent[static_cast<std::size_t>(node)]) {
        bits[static_cast<std::size_t>(node) / 64u] |=
            std::uint64_t{1} << (static_cast<unsigned>(node) & 63u);
      }
    }
  }

  const auto paths_are_disjoint = [&](int first_a,
                                      int first_b,
                                      int second_a,
                                      int second_b) {
    if (first_a > first_b) std::swap(first_a, first_b);
    if (second_a > second_b) std::swap(second_a, second_b);
    const std::uint64_t* first =
        pair_paths.data() +
        (static_cast<std::size_t>(first_a) *
             static_cast<std::size_t>(leaf_count) +
         static_cast<std::size_t>(first_b)) *
            path_stride;
    const std::uint64_t* second =
        pair_paths.data() +
        (static_cast<std::size_t>(second_a) *
             static_cast<std::size_t>(leaf_count) +
         static_cast<std::size_t>(second_b)) *
            path_stride;
    for (std::size_t word = 0; word < path_words; ++word) {
      if ((first[word] & second[word]) != 0) return false;
    }
    return true;
  };

  std::array<double, 3> result{0.0, 0.0, 0.0};
  for (int i = 0; i < static_cast<int>(leaves.size()); ++i) {
    for (int j = i + 1; j < static_cast<int>(leaves.size()); ++j) {
      for (int k = j + 1; k < static_cast<int>(leaves.size()); ++k) {
        for (int l = k + 1; l < static_cast<int>(leaves.size()); ++l) {
          const std::array<int, 4> indices{i, j, k, l};
          const std::array<std::array<int, 4>, 3> pairings{{
              {0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}}};
          double weight = 1.0;
          for (int index : indices) {
            weight *= leaves[static_cast<std::size_t>(index)].weight;
          }
          require_finite(weight, "independent quartet oracle quartet weight");
          for (const auto& pairing : pairings) {
            const Leaf& a = leaves[static_cast<std::size_t>(
                indices[static_cast<std::size_t>(pairing[0])])];
            const Leaf& b = leaves[static_cast<std::size_t>(
                indices[static_cast<std::size_t>(pairing[1])])];
            const Leaf& c = leaves[static_cast<std::size_t>(
                indices[static_cast<std::size_t>(pairing[2])])];
            const Leaf& d = leaves[static_cast<std::size_t>(
                indices[static_cast<std::size_t>(pairing[3])])];
            const int a_index = indices[static_cast<std::size_t>(pairing[0])];
            const int b_index = indices[static_cast<std::size_t>(pairing[1])];
            const int c_index = indices[static_cast<std::size_t>(pairing[2])];
            const int d_index = indices[static_cast<std::size_t>(pairing[3])];
            if (!paths_are_disjoint(a_index, b_index, c_index, d_index)) {
              continue;
            }
            const int channel = independent_quartet_channel(
                owner_colors[static_cast<std::size_t>(a.owner)],
                owner_colors[static_cast<std::size_t>(b.owner)],
                owner_colors[static_cast<std::size_t>(c.owner)],
                owner_colors[static_cast<std::size_t>(d.owner)], palette);
            if (channel >= 0) {
              result[static_cast<std::size_t>(channel)] += weight;
              require_finite(result[static_cast<std::size_t>(channel)],
                             "independent quartet oracle channel sum");
            }
            break;
          }
        }
      }
    }
  }
  return result;
}

void require_quartet_counts(
    const std::array<double, 3>& actual,
    const std::array<double, 3>& expected,
    int palette,
    const std::string& label) {
  const int channels = palette == 2 ? 2 : 3;
  for (int channel = 0; channel < channels; ++channel) {
    require_close(actual[static_cast<std::size_t>(channel)],
                  expected[static_cast<std::size_t>(channel)],
                  label + " channel=" + std::to_string(channel), 1e-10,
                  1e-10);
  }
}

qraft::internal::OwnerProjection identity_owner_map(int count) {
  qraft::internal::OwnerProjection owners;
  owners.taxon_to_owner.resize(static_cast<std::size_t>(count));
  std::iota(owners.taxon_to_owner.begin(), owners.taxon_to_owner.end(), 0);
  owners.owner_count = count;
  owners.live_owner.assign(static_cast<std::size_t>(count), 1);
  owners.live_taxa.resize(static_cast<std::size_t>(count));
  std::iota(owners.live_taxa.begin(), owners.live_taxa.end(), 0);
  return owners;
}

qraft::internal::ActiveTaxa active_set(std::vector<int> taxa,
                                        int taxon_count) {
  qraft::internal::ActiveTaxa active;
  active.active_by_taxon =
      qraft::internal::make_taxon_index(taxa, taxon_count);
  active.taxa = active.active_by_taxon.taxa;
  return active;
}

qraft::internal::StartingTree small_starting_tree() {
  qraft::internal::StartingTree starting_tree;
  starting_tree.nodes.resize(7);
  starting_tree.root = 0;
  starting_tree.nodes[0].children = {1, 2};
  starting_tree.nodes[1].parent = 0;
  starting_tree.nodes[2].parent = 0;
  starting_tree.nodes[1].children = {3, 4};
  starting_tree.nodes[3].parent = 1;
  starting_tree.nodes[4].parent = 1;
  starting_tree.nodes[2].children = {5, 6};
  starting_tree.nodes[5].parent = 2;
  starting_tree.nodes[6].parent = 2;
  starting_tree.nodes[3].taxon = 0;
  starting_tree.nodes[4].taxon = 1;
  starting_tree.nodes[5].taxon = 2;
  starting_tree.nodes[6].taxon = 3;
  starting_tree.taxon_to_leaf = {3, 4, 5, 6};
  return starting_tree;
}

qraft::internal::StartingTree balanced_starting_tree(int leaf_count) {
  using qraft::internal::StartingTree;
  StartingTree starting_tree;
  starting_tree.taxon_to_leaf.assign(static_cast<std::size_t>(leaf_count), -1);
  auto add_node = [&]() {
    starting_tree.nodes.emplace_back();
    return static_cast<int>(starting_tree.nodes.size()) - 1;
  };
  auto build = [&](auto&& self, int begin, int end) -> int {
    const int node = add_node();
    if (end - begin == 1) {
      starting_tree.nodes[static_cast<std::size_t>(node)].taxon = begin;
      starting_tree.taxon_to_leaf[static_cast<std::size_t>(begin)] = node;
      return node;
    }
    const int middle = begin + (end - begin) / 2;
    const int left = self(self, begin, middle);
    const int right = self(self, middle, end);
    starting_tree.nodes[static_cast<std::size_t>(node)].children = {left, right};
    starting_tree.nodes[static_cast<std::size_t>(left)].parent = node;
    starting_tree.nodes[static_cast<std::size_t>(right)].parent = node;
    return node;
  };
  starting_tree.root = leaf_count > 0 ? build(build, 0, leaf_count) : -1;
  return starting_tree;
}

qraft::internal::StartingTree caterpillar_starting_tree(int leaf_count) {
  using qraft::internal::StartingTree;
  StartingTree starting_tree;
  starting_tree.taxon_to_leaf.assign(static_cast<std::size_t>(leaf_count), -1);
  auto leaf = [&](int taxon) {
    const int node = static_cast<int>(starting_tree.nodes.size());
    starting_tree.nodes.emplace_back();
    starting_tree.nodes.back().taxon = taxon;
    starting_tree.taxon_to_leaf[static_cast<std::size_t>(taxon)] = node;
    return node;
  };
  if (leaf_count <= 0) return starting_tree;
  int root = leaf(leaf_count - 1);
  for (int taxon = leaf_count - 2; taxon >= 0; --taxon) {
    const int next_leaf = leaf(taxon);
    const int parent = static_cast<int>(starting_tree.nodes.size());
    starting_tree.nodes.emplace_back();
    starting_tree.nodes.back().children = {next_leaf, root};
    starting_tree.nodes[static_cast<std::size_t>(next_leaf)].parent = parent;
    starting_tree.nodes[static_cast<std::size_t>(root)].parent = parent;
    root = parent;
  }
  starting_tree.root = root;
  return starting_tree;
}

qraft::Options test_options(int threads = 2) {
  qraft::Options options;
  options.threads = threads;
  options.progress = qraft::ProgressMode::Off;
  options.divide_score = qraft::DivideScore::Quartet;
  options.regraft_score = qraft::RegraftScore::Raw;
  return options;
}

void test_tie_band_selection_permutations() {
  using namespace qraft::internal;
  constexpr double kEpsilon = 1e-9;
  constexpr double kBaseGain = 2.0;

  struct Candidate {
    int stable_key = -1;
    double score = 0.0;
    bool eligible = true;
  };

  const std::array<Candidate, 3> canonical{{
      {0, kBaseGain, true},
      {1, kBaseGain + 0.75 * kEpsilon, true},
      {2, kBaseGain + 1.5 * kEpsilon, true},
  }};
  for (const Candidate& candidate : canonical) {
    require_finite(candidate.score,
                   "shared score-selector adversarial input");
  }
  require_close(canonical[1].score - canonical[0].score,
                0.75 * kEpsilon, "shared selector A-to-B spacing", 1e-15,
                1e-6);
  require_close(canonical[2].score - canonical[0].score,
                1.5 * kEpsilon, "shared selector A-to-C spacing", 1e-15,
                1e-6);
  require(canonical[2].score - canonical[1].score <= kEpsilon &&
              canonical[2].score - canonical[0].score > kEpsilon,
          "shared selector adversarial scores do not define the intended "
          "maximum band");

  SolverTree applicability_tree;
  applicability_tree.nodes.resize(6);
  for (const std::array<int, 2>& edge :
       {std::array<int, 2>{0, 1}, std::array<int, 2>{0, 2},
        std::array<int, 2>{0, 3}, std::array<int, 2>{0, 4},
        std::array<int, 2>{1, 5}}) {
    add_solver_edge(applicability_tree, edge[0], edge[1]);
  }

  std::array<NniMove, 3> nni_moves{};
  for (int i = 0; i < static_cast<int>(nni_moves.size()); ++i) {
    NniMove& move = nni_moves[static_cast<std::size_t>(i)];
    move.edge_a = 0;
    move.edge_b = 1;
    move.swap_a_neighbor = i + 2;
    move.swap_b_neighbor = 5;
    move.alternative = i + 1;
    move.gain = canonical[static_cast<std::size_t>(i)].score;
  }
  const auto nni_key = [](const NniMove& move) {
    return std::array<int, 4>{move.edge_a, move.edge_b,
                              move.swap_a_neighbor, move.swap_b_neighbor};
  };
  require(nni_key(nni_moves[0]) < nni_key(nni_moves[1]) &&
              nni_key(nni_moves[1]) < nni_key(nni_moves[2]),
          "NNI adversarial keys are not ordered A < B < C");

  std::array<int, 3> permutation{0, 1, 2};
  do {
    std::vector<Candidate> candidates;
    std::vector<NniMove> permuted_moves;
    candidates.reserve(permutation.size());
    permuted_moves.reserve(permutation.size());
    for (int index : permutation) {
      candidates.push_back(canonical[static_cast<std::size_t>(index)]);
      permuted_moves.push_back(nni_moves[static_cast<std::size_t>(index)]);
    }

    const auto selection = select_max_with_tie_band(
        candidates.begin(), candidates.end(), kEpsilon,
        [](const Candidate& candidate) { return candidate.score; },
        [](const Candidate& candidate) { return candidate.stable_key; },
        [](const Candidate& candidate) { return candidate.eligible; });
    require(selection.exact_maximum != candidates.end() &&
                selection.exact_maximum->stable_key == 2,
            "shared selector did not find exact numerical maximum C");
    require(selection.selected != candidates.end() &&
                selection.selected->stable_key == 1,
            "shared selector did not select stable-key B from C's band");
    require(selection.band_size == 2,
            "shared selector maximum band should contain only B and C");
    require_finite(selection.selected->score,
                   "shared selector adversarial result");

    const NniMove cached = select_nni_move_for_tests(
        applicability_tree, permuted_moves, kEpsilon);
    const auto require_nni_result = [&](const char* path,
                                        const NniMove& result) {
      require(nni_key(result) == nni_key(nni_moves[1]),
              std::string("NNI ") + path +
                  " selector did not choose stable-key B");
      require_finite(result.gain,
                     std::string("NNI ") + path +
                         " adversarial gain");
      require_close(result.gain, canonical[1].score,
                    std::string("NNI ") + path +
                        " selected gain",
                    1e-15, 1e-12);
    };
    require_nni_result("cached", cached);
  } while (std::next_permutation(permutation.begin(), permutation.end()));

  const std::array<Candidate, 2> outside_band{{
      {0, kBaseGain, true},
      {1, kBaseGain + 1.5 * kEpsilon, true},
  }};
  const auto outside_selection = select_max_with_tie_band(
      outside_band.begin(), outside_band.end(), kEpsilon,
      [](const Candidate& candidate) { return candidate.score; },
      [](const Candidate& candidate) { return candidate.stable_key; });
  require(outside_selection.selected == outside_band.begin() + 1 &&
              outside_selection.band_size == 1,
          "stable key overrode an exact score outside the epsilon band");

  const std::array<Candidate, 2> exact_tie{{
      {2, kBaseGain, true},
      {1, kBaseGain, true},
  }};
  const auto tie_selection = select_max_with_tie_band(
      exact_tie.begin(), exact_tie.end(), kEpsilon,
      [](const Candidate& candidate) { return candidate.score; },
      [](const Candidate& candidate) { return candidate.stable_key; });
  require(tie_selection.selected == exact_tie.begin() + 1 &&
              tie_selection.band_size == 2,
          "exact score tie did not preserve stable-key semantics");

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double negative_infinity =
      -std::numeric_limits<double>::infinity();
  const double positive_infinity = std::numeric_limits<double>::infinity();
  require(score_class(nan) == ScoreClass::NotANumber &&
              score_class(negative_infinity) ==
                  ScoreClass::NegativeInfinity &&
              score_class(kBaseGain) == ScoreClass::Finite &&
              score_class(positive_infinity) ==
                  ScoreClass::PositiveInfinity,
          "shared selector score classes diverged from the declared order");
  require(exact_score_less(nan, negative_infinity) &&
              exact_score_less(negative_infinity, kBaseGain) &&
              exact_score_less(kBaseGain, positive_infinity),
          "shared selector non-finite class order is not total");

  const std::array<Candidate, 5> classified{{
      {0, nan, true},
      {1, negative_infinity, true},
      {2, kBaseGain, true},
      {4, positive_infinity, true},
      {3, positive_infinity, true},
  }};
  std::array<int, 5> classified_permutation{0, 1, 2, 3, 4};
  do {
    std::vector<Candidate> candidates;
    candidates.reserve(classified_permutation.size());
    for (int index : classified_permutation) {
      candidates.push_back(classified[static_cast<std::size_t>(index)]);
    }
    const auto selection = select_max_with_tie_band(
        candidates.begin(), candidates.end(), kEpsilon,
        [](const Candidate& candidate) { return candidate.score; },
        [](const Candidate& candidate) { return candidate.stable_key; });
    require(selection.exact_maximum != candidates.end() &&
                score_class(selection.exact_maximum->score) ==
                    ScoreClass::PositiveInfinity &&
                selection.selected != candidates.end() &&
                selection.selected->stable_key == 3 &&
                selection.band_size == 2,
            "non-finite score classes changed selection under permutation");
  } while (std::next_permutation(classified_permutation.begin(),
                                 classified_permutation.end()));
}

void test_bipartition_near_ties() {
  using namespace qraft::internal;
  constexpr double kEpsilon = 1e-9;
  constexpr double kBaseScore = 3.0;

  const std::array<BipartitionProbe, 3> canonical{{
      {kBaseScore, 4, 0, 10, {0, 2, 4, 6}},
      {kBaseScore + 0.75 * kEpsilon, 4, 1, 11, {0, 2, 5, 7}},
      {kBaseScore + 1.5 * kEpsilon, 4, 2, 12, {0, 3, 5, 6}},
  }};
  std::array<int, 3> permutation{0, 1, 2};
  for (int workers : {1, 2, 8, 16}) {
    qraft::Options options = test_options(workers);
    configure_threads(options);
    permutation = {0, 1, 2};
    do {
      std::vector<BipartitionProbe> probes;
      probes.reserve(permutation.size());
      for (int index : permutation) {
        probes.push_back(canonical[static_cast<std::size_t>(index)]);
      }
      const int selected =
          select_bipartition_for_tests(probes, 8, 1);
      require(selected >= 0 &&
                  probes[static_cast<std::size_t>(selected)]
                          .representative_order == 1,
              "bipartition near-tie selection changed under permutation or "
              "worker count");
    } while (std::next_permutation(permutation.begin(), permutation.end()));
  }

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const std::vector<BipartitionProbe> classified{
      {nan, 4, 0, 10, {0, 2, 4, 6}},
      {-inf, 4, 1, 11, {0, 2, 5, 7}},
      {kBaseScore, 4, 2, 12, {0, 3, 5, 6}},
      {inf, 4, 4, 14, {0, 1, 4, 7}},
      {inf, 4, 3, 13, {0, 1, 5, 8}},
  };
  const int selected =
      select_bipartition_for_tests(classified, 8, 1);
  require(selected >= 0 &&
              classified[static_cast<std::size_t>(selected)]
                      .representative_order == 3,
          "bipartition call site diverged from the declared non-finite "
          "score-class order");

  const std::vector<BipartitionProbe> ineligible{
      {inf, 0, 0, 10, {0, 0, 0, 0}},
      {kBaseScore, 4, 1, 11, {0, 2, 5, 7}},
  };
  require(select_bipartition_for_tests(ineligible, 8, 1) == 1,
          "bipartition eligibility filter admitted an empty side");
}

void test_finite_score_guards() {
  // Gene indexing yields finite leaf masses and the API has no contract for
  // non-finite ones, so the semantic tests assume finite inputs and
  // intermediates. The comparison helpers reject non-finite values so that a
  // NaN cannot turn a semantic test into a false pass.
  for (double value : {std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::infinity(),
                       -std::numeric_limits<double>::infinity()}) {
    bool rejected = false;
    try {
      require_close(value, 0.0, "injected non-finite comparison");
    } catch (const std::runtime_error&) {
      rejected = true;
    }
    require(rejected, "finite-aware oracle accepted a non-finite actual value");

    rejected = false;
    try {
      require_close(0.0, value, "injected non-finite reference");
    } catch (const std::runtime_error&) {
      rejected = true;
    }
    require(rejected,
            "finite-aware oracle accepted a non-finite reference value");
  }
  std::cout << "finite-intermediate theorem assumption: semantic tests "
               "inject finite masses only; non-finite values are not part "
               "of the API contract\n";
}

void test_hdt_weighted_fixtures() {
  using namespace qraft::internal;
  struct Fixture {
    const char* name;
    GeneIndex gene;
    std::vector<std::uint8_t> colors2;
    std::vector<std::uint8_t> colors4;
  };

  GeneIndex unresolved = star_gene_index(4);
  const std::array<double, 4> unresolved_weights{0.5, 2.0 / 3.0, 0.75,
                                                 1.25};
  for (int taxon = 0; taxon < 4; ++taxon) {
    qraft::internal::set_gene_leaf_mass(
        unresolved, unresolved.taxon_to_node[static_cast<std::size_t>(taxon)],
        unresolved_weights[static_cast<std::size_t>(taxon)]);
  }

  GeneIndex localized = broom_gene_index(2);
  const std::array<double, 4> localized_weights{0.5, 1.5, 2.5, 3.5};
  for (int taxon = 0; taxon < 4; ++taxon) {
    qraft::internal::set_gene_leaf_mass(
        localized, localized.taxon_to_node[static_cast<std::size_t>(taxon)],
        localized_weights[static_cast<std::size_t>(taxon)]);
  }

  GeneIndex nested = nested_polytomy_gene_index();
  for (int taxon = 0; taxon < 9; ++taxon) {
    qraft::internal::set_gene_leaf_mass(
        nested, nested.taxon_to_node[static_cast<std::size_t>(taxon)],
        static_cast<double>(taxon + 2) / static_cast<double>(taxon + 3));
  }

  std::vector<Fixture> fixtures;
  fixtures.push_back(
      {"degree-four unresolved star", std::move(unresolved), {0, 0, 1, 1},
       {0, 1, 2, 3}});
  fixtures.push_back(
      {"pair localized in resolved broom arms", std::move(localized),
       {0, 0, 1, 1}, {0, 1, 2, 3}});
  fixtures.push_back(
      {"nested and adjacent arbitrary polytomies", std::move(nested),
       {0, 1, 0, 1, 1, 0, 1, 0, 1},
       {0, 1, 2, 3, 0, 3, 2, 1, 0}});

  for (const Fixture& fixture : fixtures) {
    const int leaves = static_cast<int>(fixture.gene.taxon_to_node.size());
    const OwnerProjection owners = identity_owner_map(leaves);
    for (int palette : {2, 4}) {
      const std::vector<std::uint8_t>& colors =
          palette == 2 ? fixture.colors2 : fixture.colors4;
      HdtColorPolicy policy;
      policy.current_by_owner = colors;
      policy.background_by_owner.assign(static_cast<std::size_t>(leaves), 0);
      const HdtState state = build_hdt_state(
          build_base_component_graph(fixture.gene, owners, palette, policy));
      const std::array<double, 3> expected =
          independent_gene_quartet_counts(fixture.gene, owners, palette,
                                           colors);
      require_quartet_counts(hdt_root_counts(state), expected, palette,
                             fixture.name);
    }
  }

  const std::array<double, 3> unresolved_expected =
      independent_gene_quartet_counts(fixtures[0].gene, identity_owner_map(4),
                                       4, fixtures[0].colors4);
  require_close(unresolved_expected[0], 0.0,
                "four leaves in four star arms must be unresolved");
  require_close(unresolved_expected[1], 0.0,
                "four leaves in four star arms must be unresolved");
  require_close(unresolved_expected[2], 0.0,
                "four leaves in four star arms must be unresolved");
  require_close(independent_gene_quartet_counts(
                    fixtures[1].gene, identity_owner_map(4), 4,
                    fixtures[1].colors4)[0],
                0.5 * 1.5 * 2.5 * 3.5,
                "localized broom pairs should resolve q01_23", 1e-12,
                1e-12);
}

void test_hdt_duplicate_owner_contexts() {
  using namespace qraft::internal;
  GeneIndex gene = nested_polytomy_gene_index();
  for (int taxon = 0; taxon < 9; ++taxon) {
    set_gene_leaf_mass(
        gene, gene.taxon_to_node[static_cast<std::size_t>(taxon)],
        static_cast<double>(2 * taxon + 1) / static_cast<double>(taxon + 4));
  }

  OwnerProjection owners;
  owners.taxon_to_owner = {0, 0, 1, 2, 2, 3, 4, -1, 5};
  owners.owner_count = 6;
  owners.live_owner.assign(6, 1);
  owners.live_taxa = {0, 1, 2, 3, 4, 5, 6, 8};

  for (int palette : {2, 4}) {
    HdtColorPolicy policy;
    policy.current_by_owner = palette == 2
                                  ? std::vector<std::uint8_t>{1, 1, 0, 1, 0, 1}
                                  : std::vector<std::uint8_t>{1, 2, 0, 3, 2, 1};
    policy.background_by_owner =
        palette == 2 ? std::vector<std::uint8_t>(6, 0)
                     : std::vector<std::uint8_t>{0, 3, 0, 3, 0, 3};

    auto make_state = [&] {
      return build_hdt_state(
          build_base_component_graph(gene, owners, palette, policy));
    };
    const std::array<double, 3> base_expected =
        independent_gene_quartet_counts(gene, owners, palette,
                                         policy.current_by_owner);
    HdtState base = make_state();
    require_quartet_counts(hdt_root_counts(base), base_expected, palette,
                           "duplicate-owner weighted base");
    HdtInitialStateResult induced = hdt_build_initial_state(
        build_induced_base_component_graph(gene, owners, palette, policy));
    require(!induced.contracted,
            "duplicate-owner induction retained an absent source leaf");
    require_quartet_counts(hdt_root_counts(induced.state), base_expected,
                           palette, "duplicate-owner direct induction");
    require(hdt_recolor_owner(
                induced.state, 0,
                static_cast<std::uint8_t>(palette == 2 ? 0 : 2)) == 2,
            "direct induction lost a duplicate physical owner occurrence");

    const int changed = hdt_recolor_owner(
        base, 0, static_cast<std::uint8_t>(palette == 2 ? 0 : 2));
    require(changed == 2,
            "recoloring a duplicate owner must update both physical leaves");
    std::vector<std::uint8_t> duplicate_recolored = policy.current_by_owner;
    duplicate_recolored[0] = static_cast<std::uint8_t>(palette == 2 ? 0 : 2);
    require_quartet_counts(
        hdt_root_counts(base),
        independent_gene_quartet_counts(gene, owners, palette,
                                         duplicate_recolored),
        palette, "multiple leaves per owner recolor");

    for (const std::vector<int>& selected :
         {std::vector<int>{}, std::vector<int>{1}, std::vector<int>{0},
          std::vector<int>{0, 2, 4}}) {
      const HdtLiveSelector selector =
          hdt_owner_selector(owners.owner_count, selected);
      std::vector<std::uint8_t> expected_colors = policy.current_by_owner;
      for (int owner = 0; owner < owners.owner_count; ++owner) {
        if (!selector.contains_owner(owner)) {
          expected_colors[static_cast<std::size_t>(owner)] =
              policy.background_by_owner[static_cast<std::size_t>(owner)];
        }
      }
      const std::array<double, 3> expected =
          independent_gene_quartet_counts(gene, owners, palette,
                                           expected_colors);

      const HdtContractResult contracted = hdt_contract(
          make_state(), selector, HdtContractMode::LightCopy);
      int expected_selected_atoms = 0;
      for (int mapped_owner : owners.taxon_to_owner) {
        if (mapped_owner >= 0 && selector.contains_owner(mapped_owner)) {
          ++expected_selected_atoms;
        }
      }
      require(contracted.stats.selected_live_atoms == expected_selected_atoms,
              "Contract selected-atom count lost duplicate-owner leaves");
      require_quartet_counts(hdt_root_counts(contracted.state), expected,
                             palette,
                             "s=0/s=1 duplicate-owner Contract semantics");

      const HdtExtractContractResult extracted =
          hdt_extract_contract(make_state(), selector);
      require_quartet_counts(hdt_root_counts(extracted.state), expected,
                             palette,
                             "s=0/s=1 duplicate-owner Extract semantics");
    }

    const std::vector<int> selected{0, 2, 4};
    const HdtLiveSelector selector =
        hdt_owner_selector(owners.owner_count, selected);
    HdtContractResult contracted = hdt_contract(
        make_state(), selector, HdtContractMode::LightCopy);
    HdtExtractContractResult extracted =
        hdt_extract_contract(make_state(), selector);
    const std::uint8_t next_color =
        static_cast<std::uint8_t>(palette == 2 ? 0 : 3);
    hdt_recolor_owner(contracted.state, 0, next_color);
    hdt_recolor_owner(extracted.state, 0, next_color);
    std::vector<std::uint8_t> contextual = policy.background_by_owner;
    for (int owner : selected) {
      contextual[static_cast<std::size_t>(owner)] =
          policy.current_by_owner[static_cast<std::size_t>(owner)];
    }
    contextual[0] = next_color;
    const std::array<double, 3> contextual_expected =
        independent_gene_quartet_counts(gene, owners, palette, contextual);
    require_quartet_counts(hdt_root_counts(contracted.state),
                           contextual_expected, palette,
                           "Contract contextual recolor equivalence");
    require_quartet_counts(hdt_root_counts(extracted.state),
                           contextual_expected, palette,
                           "ExtractContract contextual recolor equivalence");

    const HdtLiveSelector nested_selector =
        hdt_owner_selector(owners.owner_count, {0});
    const HdtContractResult nested = hdt_contract(
        extracted.state, nested_selector, HdtContractMode::HeavyReplace);
    std::vector<std::uint8_t> nested_colors = policy.background_by_owner;
    nested_colors[0] = next_color;
    require_quartet_counts(
        hdt_root_counts(nested.state),
        independent_gene_quartet_counts(gene, owners, palette, nested_colors),
        palette, "ExtractContract contextual nested Contract equivalence");
  }
}

qraft::internal::GeneIndex regraft_support_gene(int supported_taxon,
                                             double quartet_weight) {
  using qraft::internal::GeneIndex;
  require(supported_taxon >= 0 && supported_taxon < 3,
          "regraft support taxon outside the attachment star");
  GeneIndex gene;
  gene.nodes.resize(6);
  std::vector<std::vector<int>> gene_adjacency;
  auto add_edge = [&](int a, int b) {
    add_adjacency_edge(gene_adjacency, a, b);
  };
  std::array<int, 2> other{};
  int position = 0;
  for (int taxon = 0; taxon < 3; ++taxon) {
    if (taxon != supported_taxon) other[static_cast<std::size_t>(position++)] = taxon;
  }
  // (supported, context=3) | (other[0], other[1]).
  add_edge(0, 1);
  add_edge(0, supported_taxon + 2);
  add_edge(0, 5);
  add_edge(1, other[0] + 2);
  add_edge(1, other[1] + 2);
  gene.taxon_to_node.resize(4);
  for (int taxon = 0; taxon < 4; ++taxon) {
    const int node = taxon + 2;
    gene.nodes[static_cast<std::size_t>(node)].set_taxon(taxon);
    qraft::internal::set_gene_leaf_mass(gene, node,
                                       taxon == 3 ? quartet_weight : 1.0);
    gene.taxon_to_node[static_cast<std::size_t>(taxon)] = node;
  }
  gene.root = 0;
  gene.reference_root = 2;
  finish_gene_adjacency(gene, gene_adjacency);
  gene.root_traversal = qraft::internal::traversal(gene, gene.root);
  return gene;
}

qraft::internal::SolverTree regraft_star(
    const std::array<int, 3>& neighbor_order) {
  using qraft::internal::SolverTree;
  SolverTree tree;
  tree.nodes.resize(4);
  tree.nodes[0].taxon = -1;
  for (int taxon = 0; taxon < 3; ++taxon) {
    tree.nodes[static_cast<std::size_t>(taxon + 1)].taxon = taxon;
  }
  for (int taxon : neighbor_order) {
    const int leaf = taxon + 1;
    tree.nodes[0].neighbors.push_back(leaf);
    tree.nodes[static_cast<std::size_t>(leaf)].neighbors.push_back(0);
  }
  return tree;
}

void test_regraft_near_ties() {
  using namespace qraft::internal;
  constexpr double kTieEpsilon = 1e-9;
  const std::array<double, 3> weights{
      1.0, 1.0 + 0.25 * kTieEpsilon, 1.0 + 0.5 * kTieEpsilon};
  const std::array<double, 3> expected_scores{
      weights[0] - 2.0 * (weights[1] + weights[2]),
      weights[1] - 2.0 * (weights[0] + weights[2]),
      weights[2] - 2.0 * (weights[0] + weights[1])};
  for (double score : expected_scores) {
    require_finite(score, "independent regraft near-tie score");
  }
  require(expected_scores[2] > expected_scores[1] &&
              expected_scores[1] > expected_scores[0],
          "independent near-tie oracle did not produce a strict triple");
  require_close(expected_scores[1] - expected_scores[0],
                0.75 * kTieEpsilon, "near-tie lower spacing", 1e-15,
                1e-6);
  require_close(expected_scores[2] - expected_scores[0],
                1.5 * kTieEpsilon, "near-tie outer spacing", 1e-15, 1e-6);

  SolverContext context;
  context.taxa = qraft::TaxonRegistry::from_ordered_unique(
      {"active-0", "active-1", "active-2", "context"});
  for (int supported = 0; supported < 3; ++supported) {
    context.genes.push_back(regraft_support_gene(
        supported, weights[static_cast<std::size_t>(supported)]));
  }
  context.options = test_options(1);
  context.options.regraft_score = qraft::RegraftScore::Raw;
  const std::vector<int> active_taxa{0, 1, 2};
  const std::vector<int> context_taxa{3};
  std::array<int, 3> permutation{0, 1, 2};

  for (int workers : {1, 2, 8, 16}) {
    context.options.threads = workers;
    context.workers = configure_threads(context.options);
    permutation = {0, 1, 2};
    do {
      const SolverTree tree = regraft_star(permutation);
      const RegraftEdge selected = choose_regraft_edge(
          context, tree, active_taxa, context_taxa);
      // Exact ranking identifies edge 2 as the numeric maximum. The separate
      // decision band contains edges 1 and 2, so the stable edge key selects
      // edge 1; edge 0 is outside the epsilon band.
      require(selected.a == 0 && selected.b == 2,
              "regraft near-tie ranking changed under permutation/workers");
      require_close(selected.rank_score, expected_scores[1],
                    "regraft near-tie selected score", 1e-12, 1e-12);
      require_close(selected.second_rank_score, expected_scores[0],
                    "regraft near-tie first score outside decision band", 1e-12,
                    1e-12);
      require(selected.ties == 2,
              "regraft decision epsilon should tie only the top two edges");
      require_close(selected.rank_margin,
                    expected_scores[2] - expected_scores[0],
                    "regraft near-tie rank margin", 1e-15, 1e-6);
    } while (std::next_permutation(permutation.begin(), permutation.end()));
  }
}

void test_newick_parser() {
  const qraft::Tree tree =
      qraft::parse_newick("((A:0.1,B:0.2),(C:0.3,D:0.4));");
  require(sorted_labels(tree) == std::vector<std::string>({"A", "B", "C", "D"}),
          "Newick parser should retain leaf labels");
  const std::string serialized = qraft::to_newick(tree);
  require(!serialized.empty() && serialized.back() == ';',
          "Newick serializer should emit a terminated tree");
}

// The parser's node and adjacency order, which index_gene and the MRE read,
// its sorted output, and its errors.
void test_newick_structure() {
  const qraft::Tree tree =
      qraft::parse_newick(" ((A:1,B)x:2.5,(C,'D d')E)[comment]; ");
  const std::vector<std::vector<int>> neighbors = {
      {1, 4}, {2, 3, 0}, {1}, {1}, {5, 6, 0}, {4}, {4}};
  const std::vector<std::string> labels = {"", "", "A", "B", "", "C", "D d"};
  require(tree.root == 0 && tree.nodes.size() == neighbors.size(),
          "Newick parser changed the node count or root");
  for (std::size_t i = 0; i < neighbors.size(); ++i) {
    require(tree.nodes[i].neighbors == neighbors[i] &&
                tree.nodes[i].label == labels[i],
            "Newick parser changed the order of node " + std::to_string(i));
  }
  require(qraft::to_newick(tree) == "(('D d',C),(A,B));",
          "Newick writer changed its sorted output");

  const auto require_error = [](const std::string& text,
                                const std::string& message) {
    std::string actual;
    try {
      (void)qraft::parse_newick(text);
    } catch (const std::invalid_argument& error) {
      actual = error.what();
    }
    require(actual == message,
            "Newick parser gave '" + actual + "' for '" + text + "'");
  };
  require_error("(A,B", "expected ',' or ')' in Newick subtree");
  require_error("(A B)", "expected ',' or ')' in Newick subtree");
  require_error("(A,)", "leaf without label in Newick tree");
  require_error("()", "leaf without label in Newick tree");
  require_error("(A,B)C);", "unexpected trailing Newick text");
  require_error("(A,[x", "unterminated Newick comment");
  require_error("('A", "unexpected end of Newick text");
}

// A caterpillar nests as deep as it has taxa, and so does a chain of
// one-child nodes. The parser and the writer keep their open nodes on the
// heap, so both run here on a thread with a 256 KiB stack.
void test_newick_deep_trees() {
  const auto caterpillar = [](int taxa) {
    std::string text(static_cast<std::size_t>(taxa - 1), '(');
    text += "t0";
    for (int i = 1; i < taxa; ++i) text += ",t" + std::to_string(i) + ")";
    return text + ";";
  };
  // Writing copies each subtree's text into its parent's, which is quadratic
  // on a caterpillar, so the written one is smaller.
  constexpr int kParsedTaxa = 200000;
  constexpr int kWrittenTaxa = 20000;
  constexpr int kChainDepth = 20000;
  const std::string parsed_text = caterpillar(kParsedTaxa);
  const std::string written_text = caterpillar(kWrittenTaxa);
  std::string chain(static_cast<std::size_t>(kChainDepth), '(');
  chain += "a";
  chain.append(static_cast<std::size_t>(kChainDepth), ')');
  chain += ";";

  qraft::test::run_with_stack(256 * 1024, [&] {
    const qraft::Tree parsed = qraft::parse_newick(parsed_text);
    require(parsed.root == 0 &&
                parsed.nodes.size() ==
                    static_cast<std::size_t>(2 * kParsedTaxa - 1) &&
                parsed.leaf_labels().size() ==
                    static_cast<std::size_t>(kParsedTaxa),
            "deep caterpillar parsed into the wrong tree");
    require(qraft::to_newick(qraft::parse_newick(written_text)) ==
                written_text,
            "deep caterpillar did not round-trip");
    const qraft::Tree unary = qraft::parse_newick(chain);
    require(unary.nodes.size() == static_cast<std::size_t>(kChainDepth + 1) &&
                unary.nodes[static_cast<std::size_t>(kChainDepth)].label ==
                    "a",
            "deep one-child chain parsed into the wrong tree");
    require(qraft::to_newick(unary) == chain,
            "deep one-child chain did not round-trip");
  });
}

void test_gene_structure_validation() {
  const qraft::TaxonRegistry taxa =
      qraft::TaxonRegistry::from_ordered_unique({"A", "B", "C"});
  auto rejected = [&](const qraft::GeneTree& gene, const std::string& context) {
    bool threw = false;
    try {
      (void)qraft::internal::index_gene(gene, taxa);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    require(threw, "public gene indexing accepted " + context);
  };

  qraft::GeneTree labeled_internal;
  labeled_internal.root = 0;
  labeled_internal.nodes.resize(3);
  labeled_internal.nodes[0].label = "A";
  labeled_internal.nodes[0].neighbors = {1, 2};
  labeled_internal.nodes[1].label = "B";
  labeled_internal.nodes[1].neighbors = {0};
  labeled_internal.nodes[2].label = "C";
  labeled_internal.nodes[2].neighbors = {0};
  rejected(labeled_internal, "a labeled internal node");

  qraft::GeneTree unlabeled_leaf;
  unlabeled_leaf.root = 0;
  unlabeled_leaf.nodes.resize(3);
  unlabeled_leaf.nodes[0].neighbors = {1, 2};
  unlabeled_leaf.nodes[1].label = "A";
  unlabeled_leaf.nodes[1].neighbors = {0};
  unlabeled_leaf.nodes[2].neighbors = {0};
  rejected(unlabeled_leaf, "an unlabeled non-root leaf");

  qraft::GeneTree invalid_neighbor = qraft::parse_newick("(A,B);");
  invalid_neighbor.nodes[static_cast<std::size_t>(invalid_neighbor.root)]
      .neighbors[0] = static_cast<int>(invalid_neighbor.nodes.size());
  rejected(invalid_neighbor, "an out-of-domain adjacency endpoint");

  qraft::GeneTree asymmetric = qraft::parse_newick("(A,B);");
  asymmetric.nodes[1].neighbors.clear();
  rejected(asymmetric, "an asymmetric adjacency");

  qraft::GeneTree cyclic;
  cyclic.root = 0;
  cyclic.nodes.resize(3);
  for (int node = 0; node < 3; ++node) {
    cyclic.nodes[static_cast<std::size_t>(node)].neighbors = {
        (node + 1) % 3, (node + 2) % 3};
  }
  rejected(cyclic, "a cyclic adjacency");

  qraft::GeneTree disconnected;
  disconnected.root = 0;
  disconnected.nodes.resize(4);
  disconnected.nodes[0].label = "A";
  disconnected.nodes[0].neighbors = {1};
  disconnected.nodes[1].label = "B";
  disconnected.nodes[1].neighbors = {0};
  disconnected.nodes[2].label = "B";
  disconnected.nodes[2].neighbors = {3};
  disconnected.nodes[3].label = "C";
  disconnected.nodes[3].neighbors = {2};
  rejected(disconnected, "a disconnected adjacency");

  for (const std::string& valid : {"A;", "(A);", "(A,B);"}) {
    const qraft::GeneTree parsed = qraft::parse_newick(valid);
    const qraft::internal::GeneIndex indexed =
        qraft::internal::index_gene(parsed, taxa);
    require(!indexed.nodes.empty(),
            "parser-compatible singleton/unary tree was rejected: " + valid);
  }
  // A unary parse root is a dangling vertex: indexing removes it, so the
  // same unrooted tree indexes the same way however it is written.
  for (const auto& [plain, unary] :
       {std::pair<std::string, std::string>{"A;", "(A);"},
        {"(A,B);", "((A,B));"},
        {"((A,B),C);", "(((A,B),C));"}}) {
    const auto a = qraft::internal::index_gene(qraft::parse_newick(plain), taxa);
    const auto b = qraft::internal::index_gene(qraft::parse_newick(unary), taxa);
    std::vector<int> degrees_a;
    std::vector<int> degrees_b;
    for (int v = 0; v < static_cast<int>(a.nodes.size()); ++v) degrees_a.push_back(qraft::internal::gene_degree(a, v));
    for (int v = 0; v < static_cast<int>(b.nodes.size()); ++v) degrees_b.push_back(qraft::internal::gene_degree(b, v));
    std::sort(degrees_a.begin(), degrees_a.end());
    std::sort(degrees_b.begin(), degrees_b.end());
    require(degrees_a == degrees_b,
            "a unary parse root left a dangling vertex: " + unary);
  }
  {
    // The test index needs four taxa.
    const qraft::TaxonRegistry four({"A", "B", "C", "D"});
    const auto first = qraft::internal::index_gene(qraft::parse_newick("(A,C);"), four);
    const auto second = qraft::internal::index_gene(qraft::parse_newick("((A,C));"), four);
    const auto index = qraft::internal::make_bipartition_test_index({first, second}, four.size());
    require(index->single_copy && index->genes.size() == 2,
            "two-taxon genes did not get a test index");
    const int a = four.id_for("A");
    const int c = four.id_for("C");
    require(test_index_distance(index->genes[0], first, a, c) == 1 &&
                test_index_distance(index->genes[1], second, a, c) == 1,
            "a unary parse root lengthened a gene-tree distance");
  }
}

void test_checked_arithmetic_guards() {
  using namespace qraft::internal;
  require(checked_hdt_size_add(2u, 3u, "unexpected") == 5u &&
              checked_hdt_size_mul(3u, 4u, "unexpected") == 12u &&
              checked_hdt_long_long_size(13u, "unexpected") == 13LL &&
              checked_hdt_int_add(4, 5, "unexpected") == 9 &&
              checked_hdt_int_mul(6, 7, "unexpected") == 42 &&
              checked_hdt_int_mul_add(5, 4, 3, "unexpected") == 23 &&
              hdt_contract_size_bound(5) == 23,
          "checked arithmetic changed an in-domain result");

  std::size_t byte_total = 5u;
  checked_hdt_accumulate_bytes(byte_total, 3u, 4u, "unexpected");
  require(byte_total == 17u,
          "checked byte accumulation changed an in-domain result");

  int rejected = 0;
  try {
    (void)checked_hdt_size_add(std::numeric_limits<std::size_t>::max(), 1u,
                               "expected size addition overflow");
  } catch (const std::length_error&) {
    ++rejected;
  }
  try {
    (void)checked_hdt_size_mul(std::numeric_limits<std::size_t>::max(), 2u,
                               "expected size multiplication overflow");
  } catch (const std::length_error&) {
    ++rejected;
  }
  try {
    std::size_t total = std::numeric_limits<std::size_t>::max() - 1u;
    checked_hdt_accumulate_bytes(total, 1u, 2u,
                                 "expected byte accumulation overflow");
  } catch (const std::length_error&) {
    ++rejected;
  }
  try {
    (void)checked_hdt_int_size(
        static_cast<std::size_t>(std::numeric_limits<int>::max()) + 1u,
        "expected integer narrowing overflow");
  } catch (const std::length_error&) {
    ++rejected;
  }
  try {
    (void)checked_hdt_int_add(std::numeric_limits<int>::max(), 1,
                              "expected integer addition overflow");
  } catch (const std::length_error&) {
    ++rejected;
  }
  try {
    (void)checked_hdt_int_mul(std::numeric_limits<int>::max(), 2,
                              "expected integer multiplication overflow");
  } catch (const std::length_error&) {
    ++rejected;
  }
  try {
    (void)hdt_contract_size_bound(std::numeric_limits<int>::max());
  } catch (const std::length_error&) {
    ++rejected;
  }
  int expected_rejections = 7;
  if (std::numeric_limits<std::size_t>::max() >
      static_cast<std::size_t>(std::numeric_limits<long long>::max())) {
    ++expected_rejections;
    try {
      (void)checked_hdt_long_long_size(
          std::numeric_limits<std::size_t>::max(),
          "expected telemetry counter narrowing overflow");
    } catch (const std::length_error&) {
      ++rejected;
    }
  }
  require(rejected == expected_rejections,
          "arithmetic guard missed an overflow boundary");
}

void test_hdt_fail_closed_guards() {
  using namespace qraft::internal;

  HdtColorPolicy colors;
  colors.current_by_owner = {0, 1, 0, 1};
  colors.background_by_owner = {0, 0, 0, 0};
  HdtInitialStateResult initial;
  initial.state = build_hdt_state(build_base_component_graph(
      star_gene_index(4), identity_owner_map(4), 2, colors));
  initial.contract.size_bound_ok = true;
  initial.contracted = true;
  initial.state.hdt.stats.node_count_ok = false;
  HdtTraversalStats traversal_stats;
  bool rejected_bad_hdt = false;
  try {
    hdt_account_initial_state(initial, traversal_stats);
  } catch (const std::logic_error&) {
    rejected_bad_hdt = true;
  }
  require(rejected_bad_hdt &&
              !traversal_stats.complexity_contract_satisfied,
          "traversal did not fail closed on invalid HDT statistics");

  auto singleton = [](int taxon) {
    SolverTree tree;
    tree.nodes.push_back(SolverTree::Node{{}, taxon, true});
    return tree;
  };
  auto edge_tree = [] {
    SolverTree tree;
    tree.nodes.resize(2);
    tree.nodes[0].neighbors = {1};
    tree.nodes[1].neighbors = {0};
    tree.nodes[0].taxon = 0;
    tree.nodes[1].taxon = 1;
    return tree;
  };
  auto rejects_connection = [](SolverTree left,
                         SolverTree right,
                         RegraftEdge left_edge,
                         RegraftEdge right_edge) {
    try {
      (void)connect_at_edges(std::move(left), std::move(right), left_edge, right_edge);
    } catch (const std::logic_error&) {
      return true;
    }
    return false;
  };

  require(rejects_connection(singleton(0), singleton(1), RegraftEdge{7, -1},
                       RegraftEdge{0, -1}),
          "connection accepted an invalid left singleton endpoint");
  require(rejects_connection(singleton(0), singleton(1), RegraftEdge{0, -1},
                       RegraftEdge{-1, -1}),
          "connection accepted a right endpoint whose sentinel is offset");
  require(rejects_connection(edge_tree(), singleton(2), RegraftEdge{0, 0},
                       RegraftEdge{0, -1}),
          "connection accepted a non-edge attachment pair");

  const SolverTree connected = connect_at_edges(singleton(0), singleton(1),
                                    RegraftEdge{0, -1},
                                    RegraftEdge{0, -1});
  require(connected.nodes.size() == 2u && has_solver_edge(connected, 0, 1) &&
              connected.nodes[0].taxon == 0 &&
              connected.nodes[1].taxon == 1,
          "valid singleton connection did not join the sole vertices directly");

  auto require_binary_tree = [](const SolverTree& tree,
                                const std::string& context) {
    require(!tree.nodes.empty(), context + ": empty tree");
    long long directed_edges = 0;
    std::vector<unsigned char> seen(tree.nodes.size(), 0);
    std::vector<int> stack{0};
    seen[0] = 1;
    while (!stack.empty()) {
      const int node = stack.back();
      stack.pop_back();
      const SolverTree::Node& item =
          tree.nodes[static_cast<std::size_t>(node)];
      require(item.active, context + ": inactive node");
      const std::size_t expected_degree =
          item.taxon >= 0 ? 1u : 3u;
      require(item.neighbors.size() == expected_degree,
              context + ": non-binary degree");
      directed_edges += static_cast<long long>(item.neighbors.size());
      for (int nb : item.neighbors) {
        require(nb >= 0 && nb < static_cast<int>(tree.nodes.size()) &&
                    nb != node && has_solver_edge(tree, nb, node),
                context + ": invalid or asymmetric adjacency");
        if (seen[static_cast<std::size_t>(nb)] == 0) {
          seen[static_cast<std::size_t>(nb)] = 1;
          stack.push_back(nb);
        }
      }
    }
    require(std::all_of(seen.begin(), seen.end(), [](unsigned char value) {
              return value != 0;
            }),
            context + ": disconnected tree");
    require(directed_edges ==
                2LL * (static_cast<long long>(tree.nodes.size()) - 1LL),
            context + ": edge count is not a tree");
  };

  SolverTree three_leaf;
  three_leaf.nodes.resize(4);
  three_leaf.nodes[0].taxon = -1;
  for (int leaf = 1; leaf <= 3; ++leaf) {
    three_leaf.nodes[0].neighbors.push_back(leaf);
    three_leaf.nodes[static_cast<std::size_t>(leaf)].neighbors = {0};
    three_leaf.nodes[static_cast<std::size_t>(leaf)].taxon = leaf;
  }
  const SolverTree singleton_left_join = connect_at_edges(
      singleton(0), three_leaf, RegraftEdge{0, -1},
      RegraftEdge{0, 1});
  require(singleton_left_join.nodes.size() == 6u &&
              has_solver_edge(singleton_left_join, 0, 5),
          "singleton-left connection retained an extra connector");
  require_binary_tree(singleton_left_join, "singleton-left connection");

  const SolverTree singleton_right_join = connect_at_edges(
      three_leaf, singleton(0), RegraftEdge{0, 1},
      RegraftEdge{0, -1});
  require(singleton_right_join.nodes.size() == 6u &&
              has_solver_edge(singleton_right_join, 4, 5),
          "singleton-right connection retained an extra connector");
  require_binary_tree(singleton_right_join, "singleton-right connection");
}

void test_taxon_registry_scaling() {
  constexpr int kGenes = 257;
  std::vector<qraft::GeneTree> genes;
  genes.reserve(kGenes);
  for (int gene = 0; gene < kGenes; ++gene) {
    genes.push_back(qraft::parse_newick("((D,C),(B,A));"));
  }

  qraft::Options options;
  options.divide = false;
  options.nni = false;
  options.threads = 1;
  options.progress = qraft::ProgressMode::Off;
  qraft::internal::reset_run_stats();
  const qraft::SpeciesTree inferred =
      qraft::infer_species_tree(std::move(genes), options);
  require(sorted_labels(inferred) ==
              std::vector<std::string>({"A", "B", "C", "D"}),
          "deduplicated registry changed sorted taxon semantics");

  const qraft::internal::RunStats& stats =
      qraft::internal::run_stats();
  require(stats.taxon_registry_physical_labels.load(
              std::memory_order_relaxed) == 4LL * kGenes,
          "registry did not count every physical label");
  require(stats.taxon_registry_unique_labels.load(std::memory_order_relaxed) ==
              4 &&
              stats.taxon_registry_sort_items.load(std::memory_order_relaxed) ==
                  4,
          "registry comparison-sorted physical duplicate labels");
}

void test_mre_duplicate_and_incomplete_genes() {
  auto infer_with_starting_tree = [](std::vector<qraft::GeneTree> genes, int workers) {
    qraft::Options options;
    options.divide = false;
    options.nni = false;
    options.threads = workers;
    options.progress = qraft::ProgressMode::Off;
    return qraft::infer_species_tree(std::move(genes), options);
  };

  const qraft::Tree duplicate_expected =
      qraft::parse_newick("((A,B),(C,D));");
  std::string duplicate_output;
  for (int workers : {1, 4}) {
    std::vector<qraft::GeneTree> genes{
        qraft::parse_newick("((A,B),((C,C),D));")};
    const qraft::SpeciesTree inferred =
        infer_with_starting_tree(std::move(genes), workers);
    const auto [rf, nrf] = qraft::rf_distance(inferred, duplicate_expected);
    (void)nrf;
    require(rf == 0,
            "MRE discarded a valid split containing duplicate physical leaves");
    const std::string serialized = qraft::to_newick(inferred);
    if (duplicate_output.empty()) {
      duplicate_output = serialized;
    } else {
      require(serialized == duplicate_output,
              "duplicate-label MRE output changed with worker count");
    }
  }

  const qraft::Tree incomplete_expected =
      qraft::parse_newick("((D,E),(A,B,C));");
  std::string canonical_output;
  for (int workers : {1, 4}) {
    for (const std::string& rooted_gene : {
             std::string("(B,C,(D,E));"),
             std::string("(D,E,(B,C));"),
         }) {
      std::vector<qraft::GeneTree> genes{
          qraft::parse_newick("(A,B,C,D,E);"),
          qraft::parse_newick(rooted_gene),
      };
      const qraft::SpeciesTree inferred =
          infer_with_starting_tree(std::move(genes), workers);
      const auto [rf, nrf] = qraft::rf_distance(inferred, incomplete_expected);
      (void)nrf;
      require(rf == 0,
              "MRE split depended on the parse root of an incomplete gene");
      const std::string serialized = qraft::to_newick(inferred);
      if (canonical_output.empty()) {
        canonical_output = serialized;
      } else {
        require(serialized == canonical_output,
                "equivalent incomplete genes produced different MRE trees");
      }
    }
  }
}

// Exact MRE split identifiers (src/mre/exact_split_ids.cpp) against explicit
// taxon sets.

// A random tree on `labels` (repeats allowed) with polytomies: groups of two
// or three items merge until at most three remain.
std::string random_multifurcating_newick(std::mt19937& rng,
                                         std::vector<std::string> items) {
  std::shuffle(items.begin(), items.end(), rng);
  while (items.size() > 3) {
    const std::size_t take =
        std::min<std::size_t>(items.size() - 1, 2 + rng() % 2);
    std::string group = "(";
    for (std::size_t i = 0; i < take; ++i) {
      const std::size_t pick = rng() % items.size();
      group += (i == 0 ? "" : ",") + items[pick];
      items.erase(items.begin() + static_cast<std::ptrdiff_t>(pick));
    }
    items.push_back(group + ")");
  }
  std::string out = "(";
  for (std::size_t i = 0; i < items.size(); ++i) {
    out += (i == 0 ? "" : ",") + items[i];
  }
  return out + ");";
}

std::string exact_test_label(int taxon) {
  std::string label = std::to_string(taxon);
  return "t" + std::string(3 - label.size(), '0') + label;
}

// Every non-root node of a gene in its reference traversal, with its edge
// ordinal and its explicit taxon data.
struct ExplicitGeneSide {
  int node = -1;
  uint32_t edge = 0;
  std::vector<int> taxa;  // distinct, sorted
  bool closed = false;    // all or none of every taxon's copies
  // The taxon the MRE reports for the side: the node's own taxon when the
  // node is its taxon's last copy, else the first child's, in incidence order.
  int one_taxon = -1;
};

std::vector<ExplicitGeneSide> explicit_gene_sides(
    const qraft::internal::GeneIndex& gene, int& distinct_taxa) {
  using namespace qraft::internal;
  GeneTraversal scratch;
  const GeneTraversal& tr = gene_reference_traversal(gene, scratch);
  const int node_count = static_cast<int>(gene.nodes.size());
  std::map<int, int> copies;
  for (int v = 0; v < node_count; ++v) {
    const int taxon = gene.nodes[static_cast<std::size_t>(v)].taxon();
    if (taxon >= 0) ++copies[taxon];
  }
  distinct_taxa = static_cast<int>(copies.size());
  std::vector<std::vector<int>> children(static_cast<std::size_t>(node_count));
  for (int v : tr.preorder) {
    if (v != gene.reference_root) {
      children[static_cast<std::size_t>(tr.parent[static_cast<std::size_t>(v)])]
          .push_back(v);
    }
  }
  std::vector<int> last_copy(static_cast<std::size_t>(
                                 copies.empty() ? 0 : copies.rbegin()->first + 1),
                             -1);
  for (int v = 0; v < node_count; ++v) {
    const int taxon = gene.nodes[static_cast<std::size_t>(v)].taxon();
    if (taxon >= 0) last_copy[static_cast<std::size_t>(taxon)] = v;
  }
  std::vector<int> first_taxon(static_cast<std::size_t>(node_count), -1);
  for (auto it = tr.preorder.rbegin(); it != tr.preorder.rend(); ++it) {
    const int v = *it;
    const int taxon = gene.nodes[static_cast<std::size_t>(v)].taxon();
    if (taxon >= 0 && last_copy[static_cast<std::size_t>(taxon)] == v) {
      first_taxon[static_cast<std::size_t>(v)] = taxon;
      continue;
    }
    for (int c : gene_neighbors(gene, v)) {
      if (c == gene.reference_root ||
          tr.parent[static_cast<std::size_t>(c)] != v) {
        continue;
      }
      if (first_taxon[static_cast<std::size_t>(c)] >= 0) {
        first_taxon[static_cast<std::size_t>(v)] =
            first_taxon[static_cast<std::size_t>(c)];
        break;
      }
    }
  }
  std::vector<ExplicitGeneSide> sides;
  uint32_t edge = 0;
  for (int v : tr.preorder) {
    if (v == gene.reference_root) continue;
    std::map<int, int> inside;
    std::vector<int> stack{v};
    while (!stack.empty()) {
      const int x = stack.back();
      stack.pop_back();
      const int taxon = gene.nodes[static_cast<std::size_t>(x)].taxon();
      if (taxon >= 0) ++inside[taxon];
      for (int c : children[static_cast<std::size_t>(x)]) stack.push_back(c);
    }
    ExplicitGeneSide side;
    side.node = v;
    side.edge = edge++;
    side.closed = true;
    side.one_taxon = first_taxon[static_cast<std::size_t>(v)];
    for (const auto& [taxon, count] : inside) {
      side.taxa.push_back(taxon);
      side.closed = side.closed && count == copies[taxon];
    }
    sides.push_back(std::move(side));
  }
  return sides;
}

// The MRE split entries from explicit sets: (count desc, first gene, first
// edge), with the representative node and side size.
struct ExplicitSplit {
  std::vector<int> taxa;
  uint32_t count = 0;
  uint32_t first_tree = 0;
  uint32_t first_edge = 0;
  uint32_t rep_node = 0;
  int one_taxon = -1;
};

std::vector<ExplicitSplit> explicit_mre_splits(
    const std::vector<qraft::internal::GeneIndex>& genes) {
  std::map<std::vector<int>, ExplicitSplit> by_set;
  for (std::size_t g = 0; g < genes.size(); ++g) {
    int distinct = 0;
    const std::vector<ExplicitGeneSide> sides =
        explicit_gene_sides(genes[g], distinct);
    if (distinct < 4) continue;
    std::set<std::vector<int>> seen;
    for (const ExplicitGeneSide& side : sides) {
      const int size = static_cast<int>(side.taxa.size());
      if (!side.closed || size < 2 || size > distinct - 2) continue;
      auto [it, inserted] = by_set.try_emplace(side.taxa);
      ExplicitSplit& split = it->second;
      if (inserted) {
        split.taxa = side.taxa;
        split.first_tree = static_cast<uint32_t>(g);
        split.first_edge = side.edge;
        split.rep_node = static_cast<uint32_t>(side.node);
        split.one_taxon = side.one_taxon;
      }
      if (seen.insert(side.taxa).second) ++split.count;
    }
  }
  std::vector<ExplicitSplit> out;
  for (auto& [taxa, split] : by_set) out.push_back(std::move(split));
  std::sort(out.begin(), out.end(),
            [](const ExplicitSplit& a, const ExplicitSplit& b) {
              if (a.count != b.count) return a.count > b.count;
              if (a.first_tree != b.first_tree) return a.first_tree < b.first_tree;
              return a.first_edge < b.first_edge;
            });
  return out;
}

// Checks the exact identifiers and the split entries of `genes` against
// explicit sets; the assembly reads nothing else. With `complete`, every gene
// holds every taxon once and the MRE tree is also checked against a naive
// greedy consensus.
void check_exact_mre(const std::vector<qraft::internal::GeneIndex>& genes,
                     const qraft::TaxonRegistry& taxa,
                     bool complete,
                     const std::string& context) {
  using namespace qraft::internal;
  const int n = taxa.size();

  // Identifiers: nonzero exactly on closed sides of at least two taxa in
  // genes of at least four, and equal exactly when the sets are equal.
  const GeneClusterIds ids = exact_gene_cluster_ids(genes, n);
  std::map<std::vector<int>, uint32_t> id_of_set;
  std::map<uint32_t, std::vector<int>> set_of_id;
  for (std::size_t g = 0; g < genes.size(); ++g) {
    int distinct = 0;
    const std::vector<ExplicitGeneSide> sides =
        explicit_gene_sides(genes[g], distinct);
    require(ids.at(static_cast<int>(g), genes[g].reference_root) == 0,
            "exact id on a reference root: " + context);
    for (const ExplicitGeneSide& side : sides) {
      const uint32_t id = ids.at(static_cast<int>(g), side.node);
      const bool expected = distinct >= 4 && side.closed && side.taxa.size() >= 2;
      require((id != 0) == expected,
              "exact id present on the wrong gene node: " + context);
      if (id == 0) continue;
      const auto [set_it, set_new] = id_of_set.try_emplace(side.taxa, id);
      require(set_it->second == id, "equal sets got different exact ids: " + context);
      const auto [id_it, id_new] = set_of_id.try_emplace(id, side.taxa);
      require(id_it->second == side.taxa,
              "different sets got one exact id: " + context);
      (void)set_new;
      (void)id_new;
    }
  }

  // Entries: every field and the priority order against enumeration.
  const std::vector<ExplicitSplit> expected = explicit_mre_splits(genes);
  const std::vector<SplitAggregate256> exact =
      mre_split_entries_for_tests(genes, n);
  require(exact.size() == expected.size(),
          "split entry counts differ from enumeration: " + context);
  for (std::size_t i = 0; i < expected.size(); ++i) {
    const ExplicitSplit& want = expected[i];
    const SplitAggregate256& got = exact[i];
    require(got.count == want.count && got.first_tree == want.first_tree &&
                got.first_edge == want.first_edge &&
                got.rep_tree == want.first_tree &&
                got.rep_node == want.rep_node &&
                got.key.size == want.taxa.size() &&
                static_cast<int>(got.one_taxon) == want.one_taxon,
            "split entry differs from enumeration at rank " +
                std::to_string(i) + ": " + context);
    require(std::binary_search(want.taxa.begin(), want.taxa.end(),
                               want.one_taxon),
            "split entry's taxon lies outside its side: " + context);
    const auto id = id_of_set.find(want.taxa);
    require(id != id_of_set.end() &&
                exact[i].key.h == std::array<uint64_t, 4>{id->second, 0, 0, 0},
            "exact entry key is not its side's identifier: " + context);
  }

  if (!complete) return;
  qraft::Options options;
  options.progress = qraft::ProgressMode::Off;
  const StartingTree exact_tree = build_mre_tree(genes, n, options);

  // Naive greedy consensus on explicit sets, up to n - 3 clusters.
  std::set<std::vector<int>> accepted;
  for (const ExplicitSplit& split : expected) {
    if (static_cast<int>(accepted.size()) >= std::max(0, n - 3)) break;
    bool compatible = true;
    for (const std::vector<int>& other : accepted) {
      std::vector<int> common;
      std::set_intersection(split.taxa.begin(), split.taxa.end(), other.begin(),
                            other.end(), std::back_inserter(common));
      compatible = compatible &&
                   (common.empty() || common.size() == split.taxa.size() ||
                    common.size() == other.size());
    }
    if (compatible) accepted.insert(split.taxa);
  }
  std::set<std::vector<int>> built;
  std::vector<std::vector<int>> below(exact_tree.nodes.size());
  for (auto it = exact_tree.preorder.rbegin(); it != exact_tree.preorder.rend();
       ++it) {
    const StartingTree::Node& node = exact_tree.nodes[static_cast<std::size_t>(*it)];
    std::vector<int>& set = below[static_cast<std::size_t>(*it)];
    if (node.taxon >= 0) set.push_back(node.taxon);
    for (int c : node.children) {
      const std::vector<int>& part = below[static_cast<std::size_t>(c)];
      set.insert(set.end(), part.begin(), part.end());
    }
    std::sort(set.begin(), set.end());
    if (*it != exact_tree.root && node.taxon < 0) built.insert(set);
  }
  require(built == accepted,
          "exact MRE tree differs from the naive greedy consensus: " + context);
}

std::vector<qraft::internal::GeneIndex> index_exact_test_genes(
    const std::vector<std::string>& newicks, qraft::TaxonRegistry& taxa) {
  std::vector<qraft::Tree> trees;
  std::vector<std::string> labels;
  for (const std::string& newick : newicks) {
    trees.push_back(qraft::parse_newick(newick));
    for (const std::string& label : trees.back().leaf_labels()) {
      labels.push_back(label);
    }
  }
  taxa = qraft::TaxonRegistry(labels);
  std::vector<qraft::internal::GeneIndex> genes;
  for (const qraft::Tree& tree : trees) {
    genes.push_back(qraft::internal::index_gene(tree, taxa));
  }
  return genes;
}

void test_mre_exact_cluster_ids() {
  using namespace qraft::internal;
  // Random forests of random trees, each on a random subset of the taxa, with
  // single-vertex trees and unary-free polytomies; ids against explicit sets.
  std::mt19937 rng(20261030u);
  for (int round = 0; round < 400; ++round) {
    const int boundaries[] = {63, 64, 65, 127, 128, 129, 257};
    const int n = round < 300
                      ? 1 + static_cast<int>(rng() % (round < 100 ? 6 : 40))
                      : boundaries[round % 7];
    const int trees = 1 + static_cast<int>(rng() % 6);
    ClusterForest forest;
    std::vector<std::vector<int>> sets;
    for (int t = 0; t < trees; ++t) {
      std::vector<int> taxa;
      for (int x = 0; x < n; ++x) {
        if (rng() % 4 != 0) taxa.push_back(x);
      }
      if (taxa.empty()) taxa.push_back(static_cast<int>(rng() % n));
      std::shuffle(taxa.begin(), taxa.end(), rng);
      // Items are (vertex, set); merge two or three until one remains. A
      // merged vertex comes after its children.
      std::vector<std::pair<int, std::vector<int>>> items;
      for (int x : taxa) {
        items.push_back({static_cast<int>(forest.parent.size()), {x}});
        forest.parent.push_back(-1);
        forest.taxon.push_back(x);
        sets.push_back({x});
      }
      while (items.size() > 1) {
        const std::size_t take =
            std::min<std::size_t>(items.size(), 2 + rng() % 2);
        const int vertex = static_cast<int>(forest.parent.size());
        forest.parent.push_back(-1);
        forest.taxon.push_back(-1);
        std::vector<int> set;
        for (std::size_t i = 0; i < take; ++i) {
          const std::size_t pick = rng() % items.size();
          forest.parent[static_cast<std::size_t>(items[pick].first)] = vertex;
          set.insert(set.end(), items[pick].second.begin(),
                     items[pick].second.end());
          items.erase(items.begin() + static_cast<std::ptrdiff_t>(pick));
        }
        std::sort(set.begin(), set.end());
        sets.push_back(set);
        items.push_back({vertex, set});
      }
    }
    const std::size_t m = forest.parent.size();
    const std::vector<uint32_t> ids = exact_cluster_ids(std::move(forest), n);
    require(ids.size() == m, "exact cluster ids have the wrong size");
    // The recursive pair order is lexicographic on taxon presence (absent
    // before present). Check dense ranks as well as equality across the
    // machine-word base case and uneven recursive boundaries.
    std::map<std::vector<int>, uint32_t> ranks;
    std::vector<std::vector<int>> presence(m, std::vector<int>(n, 0));
    for (std::size_t i = 0; i < m; ++i) {
      for (int taxon : sets[i]) presence[i][taxon] = 1;
      ranks.emplace(presence[i], 0);
    }
    uint32_t rank = 0;
    for (auto& entry : ranks) entry.second = ++rank;
    for (std::size_t i = 0; i < m; ++i) {
      require(ids[i] == ranks.at(presence[i]),
              "exact cluster rank differs from taxon presence order");
      require(ids[i] >= 1 && ids[i] <= m, "exact cluster id out of range");
      for (std::size_t j = 0; j < i; ++j) {
        require((ids[i] == ids[j]) == (sets[i] == sets[j]),
                "exact cluster ids disagree with explicit sets, n=" +
                    std::to_string(n));
      }
    }
  }
  // Malformed forests fail closed.
  bool threw = false;
  try {
    exact_cluster_ids(ClusterForest{{-1, -1}, {0, -1}}, 2);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  require(threw, "a vertex without a leaf was accepted");
}

// The numbering runs a call's passes in parallel when it has at least 2^15
// vertices and more than one worker. On forests and gene sets that large, the
// identifiers are the same on every worker count: the taxon-presence ranks for
// the forests, and equal exactly on equal taxon sets for the genes.
void test_mre_exact_cluster_ids_parallel() {
  using namespace qraft::internal;
  std::mt19937 rng(20261102u);
  const int worker_counts[] = {1, 2, 3, 4, 8};
  // Presence order: at the first taxon where two sorted sets differ, the set
  // holding it is the larger; a proper prefix is the smaller.
  struct PresenceLess {
    bool operator()(const std::vector<int>& a,
                    const std::vector<int>& b) const {
      const auto [ia, ib] = std::mismatch(a.begin(), a.end(), b.begin(), b.end());
      if (ib == b.end()) return false;
      if (ia == a.end()) return true;
      return *ia > *ib;
    }
  };
  // n = 60 and 64 reach the mask base case in the top call (64 fills the
  // whole word); n = 700 halves in parallel before it.
  for (const auto& [n, trees] :
       {std::pair<int, int>{60, 1200}, std::pair<int, int>{64, 1100},
        std::pair<int, int>{700, 90}}) {
    ClusterForest forest;
    std::vector<std::vector<int>> sets;
    for (int t = 0; t < trees; ++t) {
      std::vector<std::pair<int, std::vector<int>>> items;
      for (int x = 0; x < n; ++x) {
        if (rng() % 4 == 0) continue;
        items.push_back({static_cast<int>(forest.parent.size()), {x}});
        forest.parent.push_back(-1);
        forest.taxon.push_back(x);
        sets.push_back({x});
      }
      while (items.size() > 1) {
        const std::size_t take =
            std::min<std::size_t>(items.size(), 2 + rng() % 2);
        const int vertex = static_cast<int>(forest.parent.size());
        forest.parent.push_back(-1);
        forest.taxon.push_back(-1);
        std::vector<int> set;
        for (std::size_t i = 0; i < take; ++i) {
          const std::size_t pick = rng() % items.size();
          forest.parent[static_cast<std::size_t>(items[pick].first)] = vertex;
          set.insert(set.end(), items[pick].second.begin(),
                     items[pick].second.end());
          items[pick] = std::move(items.back());
          items.pop_back();
        }
        std::sort(set.begin(), set.end());
        sets.push_back(set);
        items.push_back({vertex, std::move(set)});
      }
    }
    const std::size_t m = forest.parent.size();
    require(m >= (std::size_t{1} << 15) * 2,
            "parallel numbering test forest is too small");
    std::map<std::vector<int>, uint32_t, PresenceLess> ranks;
    for (const std::vector<int>& set : sets) ranks.emplace(set, 0);
    uint32_t rank = 0;
    for (auto& entry : ranks) entry.second = ++rank;
    std::vector<uint32_t> first_ids;
    for (int workers : worker_counts) {
      require(configure_threads(test_options(workers)) == workers,
              "failed to configure the numbering workers");
      const std::vector<uint32_t> ids =
          exact_cluster_ids(ClusterForest(forest), n);
      if (first_ids.empty()) {
        require(ids.size() == m, "exact cluster ids have the wrong size");
        for (std::size_t i = 0; i < m; ++i) {
          require(ids[i] == ranks.at(sets[i]),
                  "exact cluster rank differs from taxon presence order, n=" +
                      std::to_string(n));
        }
        first_ids = ids;
      }
      require(ids == first_ids,
              "exact cluster ids depend on the worker count, n=" +
                  std::to_string(n) + " workers=" + std::to_string(workers));
    }
  }

  // Genes that miss taxa and repeat some, numbered in one forest.
  std::vector<std::string> newicks;
  for (int g = 0; g < 200; ++g) {
    std::vector<std::string> items;
    for (int x = 0; x < 300; ++x) {
      if (rng() % 8 == 0) continue;
      items.push_back(exact_test_label(x));
      if (rng() % 16 == 0) items.push_back(exact_test_label(x));
    }
    newicks.push_back(random_multifurcating_newick(rng, items));
  }
  qraft::TaxonRegistry taxa;
  const std::vector<GeneIndex> genes = index_exact_test_genes(newicks, taxa);
  GeneClusterIds first;
  for (int workers : worker_counts) {
    require(configure_threads(test_options(workers)) == workers,
            "failed to configure the numbering workers");
    GeneClusterIds ids = exact_gene_cluster_ids(genes, taxa.size());
    if (first.id.empty()) {
      std::map<std::vector<int>, uint32_t> id_of_set;
      std::map<uint32_t, std::vector<int>> set_of_id;
      // The forest holds every numbered node and one leaf per taxon below
      // one.
      std::size_t forest_size = 0;
      for (std::size_t g = 0; g < genes.size(); ++g) {
        int distinct = 0;
        std::set<int> leaves;
        for (const ExplicitGeneSide& side :
             explicit_gene_sides(genes[g], distinct)) {
          const uint32_t id = ids.at(static_cast<int>(g), side.node);
          require((id != 0) ==
                      (distinct >= 4 && side.closed && side.taxa.size() >= 2),
                  "parallel exact id present on the wrong gene node");
          if (id == 0) continue;
          ++forest_size;
          leaves.insert(side.taxa.begin(), side.taxa.end());
          require(id_of_set.try_emplace(side.taxa, id).first->second == id,
                  "equal sets got different parallel exact ids");
          require(set_of_id.try_emplace(id, side.taxa).first->second ==
                      side.taxa,
                  "different sets got one parallel exact id");
        }
        forest_size += leaves.size();
      }
      require(forest_size >= (std::size_t{1} << 15) * 2,
              "parallel numbering test genes are too small");
      first = std::move(ids);
      continue;
    }
    require(ids.offset == first.offset && ids.id == first.id,
            "exact gene cluster ids depend on the worker count, workers=" +
                std::to_string(workers));
  }
  require(configure_threads(test_options(1)) == 1,
          "failed to restore one worker");
}

void test_mre_exact_split_entries() {
  {
    using namespace qraft::internal;
    // A key that is not a size and a 32-bit id fails closed.
    for (const std::array<uint64_t, 4>& h :
         {std::array<uint64_t, 4>{uint64_t{1} << 32, 0, 0, 0},
          std::array<uint64_t, 4>{1, 0, 0, 1}}) {
      std::vector<std::vector<SplitOccurrence256>> chunks(1);
      SplitKey256 key;
      key.size = 2;
      key.h = h;
      chunks[0].push_back({key, 0, 0, 0, 0});
      SplitRadixStats stats;
      bool threw = false;
      try {
        aggregate_split_occurrences256(chunks, stats);
      } catch (const std::logic_error&) {
        threw = true;
      }
      require(threw, "a split key with high bits was accepted");
    }
  }
  // Fixed cases: a gene of one taxon, two, three; a caterpillar and a star;
  // repeats of the smallest taxon and of others; genes missing the smallest
  // taxon.
  const std::vector<std::vector<std::string>> fixed = {
      {"((A,B),(C,D),E);", "(A);", "(C);"},
      {"(A,B);"},
      {"(A,B,C);"},
      {"(A,B,C,D,E,F);"},
      {"(A,(B,(C,(D,(E,(F,G))))));", "(A,(B,(C,(D,(E,(F,G))))));"},
      {"((A,B),((C,C),D));"},
      {"((A,A),(B,C),(D,E));", "(A,(B,C),(D,(E,A)));", "((A,B),C,(D,E));"},
      {"((A,B),(A,C),(D,E),F);", "(((B,C),A),(D,E),F);"},
      {"(B,C,(D,E));", "(A,B,C,D,E);", "(D,E,(B,C));"},
      {"((B,C),(D,E),F);", "((A,B),(C,D),(E,F));", "((C,D),(B,E),G);"},
      {"(((A,B),(C,D)),((E,F),(G,H)));", "(((A,C),(B,D)),((E,G),(F,H)));",
       "(((A,B),(C,D)),((E,G),(F,H)));"},
  };
  for (std::size_t c = 0; c < fixed.size(); ++c) {
    qraft::TaxonRegistry taxa;
    const auto genes = index_exact_test_genes(fixed[c], taxa);
    // Complete: every gene holds every taxon exactly once.
    bool complete = true;
    for (const auto& gene : genes) {
      int leaves = 0;
      for (const auto& node : gene.nodes) leaves += node.taxon() >= 0 ? 1 : 0;
      int distinct = 0;
      explicit_gene_sides(gene, distinct);
      complete = complete && leaves == taxa.size() && distinct == taxa.size();
    }
    check_exact_mre(genes, taxa, complete, "fixed case " + std::to_string(c));
  }

  // Exhaustive: every ordered pair of the 26 trees on five taxa (star, ten
  // with one cherry, fifteen binary) and every ordered triple of the binary
  // ones.
  const std::vector<std::string> five = {"A", "B", "C", "D", "E"};
  std::vector<std::string> all_trees{"(A,B,C,D,E);"};
  std::vector<std::string> binary_trees;
  for (int x = 0; x < 5; ++x) {
    for (int y = x + 1; y < 5; ++y) {
      std::vector<std::string> rest;
      for (int z = 0; z < 5; ++z) {
        if (z != x && z != y) rest.push_back(five[z]);
      }
      all_trees.push_back("((" + five[x] + "," + five[y] + ")," + rest[0] +
                          "," + rest[1] + "," + rest[2] + ");");
      // Binary trees with cherries {x,y} and {u,v}, the fifth leaf joining
      // them; each counted once with {x,y} before {u,v}.
      for (int u = 0; u < 5; ++u) {
        for (int v = u + 1; v < 5; ++v) {
          if (u == x || u == y || v == x || v == y || u < x) continue;
          int w = 0;
          while (w == x || w == y || w == u || w == v) ++w;
          binary_trees.push_back("((" + five[x] + "," + five[y] + ")," +
                                 five[w] + ",(" + five[u] + "," + five[v] +
                                 "));");
        }
      }
    }
  }
  require(binary_trees.size() == 15, "five-taxon binary trees miscounted");
  all_trees.insert(all_trees.end(), binary_trees.begin(), binary_trees.end());
  require(all_trees.size() == 26, "five-taxon trees miscounted");
  {
  // One worker for the exhaustive cases; the fuzz below uses the helpers.
  const qraft::internal::ScopedInnerParallelBudgetSuspend one_worker;
  for (const std::string& a : all_trees) {
    for (const std::string& b : all_trees) {
      qraft::TaxonRegistry taxa;
      const auto genes = index_exact_test_genes({a, b}, taxa);
      check_exact_mre(genes, taxa, true, a + " " + b);
    }
  }
  for (const std::string& a : binary_trees) {
    for (const std::string& b : binary_trees) {
      for (const std::string& c : binary_trees) {
        qraft::TaxonRegistry taxa;
        const auto genes = index_exact_test_genes({a, b, c}, taxa);
        check_exact_mre(genes, taxa, true, a + " " + b + " " + c);
      }
    }
  }
  }

  // Seeded fuzz: complete single-copy collections with repeated topologies
  // (ties and supports above one), then collections with missing taxa
  // (including the smallest) and repeated taxa (including the smallest).
  std::mt19937 rng(20261031u);
  for (int round = 0; round < 300; ++round) {
    const int n = 4 + static_cast<int>(rng() % 21);
    const int k = 1 + static_cast<int>(rng() % 8);
    const int mode = round % 3;  // 0 complete, 1 missing, 2 missing + copies
    std::vector<std::string> labels;
    for (int x = 0; x < n; ++x) labels.push_back(exact_test_label(x));
    std::vector<std::string> newicks;
    for (int g = 0; g < k; ++g) {
      if (g > 0 && rng() % 3 == 0) {
        newicks.push_back(newicks[rng() % newicks.size()]);
        continue;
      }
      std::vector<std::string> items;
      for (int x = 0; x < n; ++x) {
        // Keep taxon 0 in the first gene so the registry holds every taxon.
        const bool keep = mode == 0 || (g == 0) || rng() % 5 != 0;
        if (!keep) continue;
        items.push_back(labels[static_cast<std::size_t>(x)]);
        if (mode == 2 && rng() % 6 == 0) {
          items.push_back(labels[static_cast<std::size_t>(x)]);
        }
      }
      while (items.size() < 2) items.push_back(labels[rng() % labels.size()]);
      newicks.push_back(random_multifurcating_newick(rng, items));
    }
    qraft::TaxonRegistry taxa;
    const auto genes = index_exact_test_genes(newicks, taxa);
    check_exact_mre(genes, taxa, mode == 0,
                    "fuzz round " + std::to_string(round));
  }
}

// The O(1) leaf-distance index of the bipartition test against breadth-first
// distances, on random trees of every shape the index handles: binary,
// multifurcating, caterpillar and star; rooted at an internal node or at a
// leaf; with ranks keyed by taxon or by node; and with leaf counts that cross
// the block (32), superblock (1,024) and sparse-table boundaries. The
// sources are random, and every leaf is a target, so the ranges end anywhere.
void test_bipartition_test_distance_index() {
  using namespace qraft::internal;
  std::mt19937 rng(20261102);
  enum class Shape { binary, multifurcating, caterpillar, star };
  int checked = 0;
  // 1024 and 1026 leaves give 1023 and 1025 LCA depths, around one
  // superblock; 9000 leaves give nine superblocks, so queries reach sparse
  // level 2.
  for (const int leaves : {1, 2, 3, 4, 5, 31, 32, 33, 34, 65, 100, 1024, 1025,
                           1026, 1057, 2081, 4200, 9000}) {
    for (const Shape shape : {Shape::binary, Shape::multifurcating,
                              Shape::caterpillar, Shape::star}) {
      for (const bool leaf_root : {false, true}) {
        for (const bool padded : {false, true}) {
          // An unrooted tree: leaves are added by subdividing an edge or, for
          // multifurcations and stars, by joining an internal node.
          qraft::Tree tree;
          const auto add_node = [&](const std::string& label) {
            tree.nodes.push_back({{}, label});
            return static_cast<int>(tree.nodes.size()) - 1;
          };
          const auto link = [&](int u, int v) {
            tree.nodes[u].neighbors.push_back(v);
            tree.nodes[v].neighbors.push_back(u);
          };
          const auto unlink = [&](int u, int v) {
            auto& nu = tree.nodes[u].neighbors;
            auto& nv = tree.nodes[v].neighbors;
            nu.erase(std::find(nu.begin(), nu.end(), v));
            nv.erase(std::find(nv.begin(), nv.end(), u));
          };
          std::vector<std::pair<int, int>> edges;
          std::vector<int> internal;
          int last_leaf = add_node("t0");
          if (leaves >= 2) {
            const int second = add_node("t1");
            link(last_leaf, second);
            edges.push_back({last_leaf, second});
            last_leaf = second;
          }
          for (int t = 2; t < leaves; ++t) {
            const std::string label = "t" + std::to_string(t);
            const bool join = !internal.empty() &&
                ((shape == Shape::star) ||
                 (shape == Shape::multifurcating && rng() % 3 == 0));
            if (join) {
              const int target = shape == Shape::star
                  ? internal.front()
                  : internal[rng() % internal.size()];
              const int leaf = add_node(label);
              link(target, leaf);
              edges.push_back({target, leaf});
              last_leaf = leaf;
              continue;
            }
            std::size_t e = rng() % edges.size();
            if (shape == Shape::caterpillar) {
              // The edge to the last leaf, so the tree grows along one path.
              e = edges.size() - 1;
            }
            const auto [u, v] = edges[e];
            const int middle = add_node("");
            internal.push_back(middle);
            unlink(u, v);
            link(u, middle);
            link(middle, v);
            edges[e] = {u, middle};
            edges.push_back({middle, v});
            const int leaf = add_node(label);
            link(middle, leaf);
            edges.push_back({middle, leaf});
            last_leaf = leaf;
          }
          std::vector<int> leaf_nodes;
          for (int u = 0; u < static_cast<int>(tree.nodes.size()); ++u) {
            if (!tree.nodes[u].label.empty()) leaf_nodes.push_back(u);
          }
          tree.root = leaf_root || internal.empty()
              ? leaf_nodes[rng() % leaf_nodes.size()]
              : internal[rng() % internal.size()];
          std::vector<std::string> labels = tree.leaf_labels();
          // A padded registry holds more taxa than the gene has nodes, so the
          // ranks are keyed by node.
          const int registry_size = padded
              ? static_cast<int>(tree.nodes.size()) + 4
              : std::max(leaves, 4);
          for (int t = leaves; t < registry_size; ++t) {
            labels.push_back("p" + std::to_string(t));
          }
          std::shuffle(labels.begin(), labels.end(), rng);
          const qraft::TaxonRegistry taxa(labels);
          const GeneIndex gene = index_gene(tree, taxa);
          const auto index = make_bipartition_test_index({gene}, taxa.size());
          require(index->single_copy && index->genes.size() == 1,
                  "random gene did not get a test index");
          const GeneDistanceIndex& distances = index->genes.front();
          require(distances.rank_by_taxon ==
                      (gene.taxon_to_node.size() <= gene.nodes.size()),
                  "distance index chose the larger rank key");
          std::vector<int> present;
          for (int t = 0; t < taxa.size(); ++t) {
            if (gene.taxon_to_node[t] >= 0) present.push_back(t);
          }
          require(static_cast<int>(present.size()) == leaves,
                  "random gene lost a leaf");
          // The range minimum against a running minimum over lca_depth, from
          // every start rank (or 64 random ones), to every later rank.
          const int ranks = leaves;
          require(static_cast<int>(distances.lca_depth.size()) == std::max(0, ranks - 1),
                  "distance index has the wrong number of LCA depths");
          std::vector<int> starts(ranks);
          std::iota(starts.begin(), starts.end(), 0);
          std::shuffle(starts.begin(), starts.end(), rng);
          if (ranks > 1100) starts.resize(64);
          for (const int r : starts) {
            int running = std::numeric_limits<int>::max();
            for (int s = r + 1; s < ranks; ++s) {
              running = std::min(running, distances.lca_depth[s - 1]);
              require(distances.common_depth(r, s) == running &&
                          distances.common_depth(s, r) == running,
                      "range minimum disagrees with a running minimum at " +
                          std::to_string(leaves) + " leaves");
              ++checked;
            }
          }
          std::vector<int> sources = present;
          std::shuffle(sources.begin(), sources.end(), rng);
          if (sources.size() > 24) sources.resize(24);
          for (const int a : sources) {
            std::vector<int> row(gene.nodes.size(), -1);
            std::queue<int> pending;
            row[gene.taxon_to_node[a]] = 0;
            pending.push(gene.taxon_to_node[a]);
            while (!pending.empty()) {
              const int u = pending.front();
              pending.pop();
              for (const int v : gene_neighbors(gene, u)) {
                if (row[v] >= 0) continue;
                row[v] = row[u] + 1;
                pending.push(v);
              }
            }
            for (const int b : present) {
              const int forward = test_index_distance(distances, gene, a, b);
              require(forward == row[gene.taxon_to_node[b]] &&
                          test_index_distance(distances, gene, b, a) == forward &&
                          test_index_distance(distances, gene, a, b) == forward,
                      "distance index disagrees with BFS at " +
                          std::to_string(leaves) + " leaves");
              ++checked;
            }
          }
        }
      }
    }
  }
  require(checked > 100000, "distance index test checked too few pairs");
}

void test_bipartition_test() {
  using namespace qraft::internal;
  // Independent breadth-first path distances, including a multifurcation and
  // a rooted Newick whose degree-two root must be suppressed before indexing.
  for (const std::string topology : {
           "((A,B),(C,D));", "((A,B),(C,D),(E,F));",
           "(A,B,C,D,E,F);", "(A,(B,(C,(D,(E,F)))));"}) {
    const auto tree = qraft::parse_newick(topology);
    SolverContext ctx;
    ctx.options.bipartition_test.anchor_pairs = 1;
    ctx.options.bipartition_test.crossover = 0.0;
    ctx.options.bipartition_test.coefficient = 0.01;
    ctx.taxa = qraft::TaxonRegistry(tree.leaf_labels());
    ctx.genes.push_back(index_gene(tree, ctx.taxa));
    const int n = ctx.taxa.size();
    ctx.bipartition_test_index = make_bipartition_test_index(ctx.genes, n);
    require(ctx.bipartition_test_index->complete_single_copy, "single-copy test index unavailable");
    const auto& gene = ctx.genes.front();
    std::vector<std::vector<int>> distances(n);
    std::set<unsigned> splits;
    for (int a = 0; a < n; ++a) {
      auto& row = distances[a];
      row.assign(gene.nodes.size(), -1);
      std::queue<int> pending;
      const int start = gene.taxon_to_node[a];
      row[start] = 0;
      pending.push(start);
      while (!pending.empty()) {
        const int u = pending.front();
        pending.pop();
        for (int v : qraft::internal::gene_neighbors(gene, u)) {
          if (row[v] >= 0) continue;
          row[v] = row[u] + 1;
          pending.push(v);
        }
      }
      for (int b = 0; b < n; ++b) {
        require(test_index_distance(ctx.bipartition_test_index->genes[0], gene, a, b) ==
                    row[gene.taxon_to_node[b]],
                "bipartition-test distance disagrees with BFS");
      }
    }
    if (n == 4) {
      require(gene.nodes.size() == 6,
              "indexing retained the rooted quartet's degree-two root");
      require(*std::max_element(distances[0].begin(), distances[0].end()) == 3,
              "bipartition test must use suppressed unrooted edge counts");
    }
    for (int u = 0; u < static_cast<int>(gene.nodes.size()); ++u) {
      for (int v : qraft::internal::gene_neighbors(gene, u)) {
        if (u > v) continue;
        unsigned side = 0;
        std::vector<std::pair<int, int>> pending{{u, v}};
        while (!pending.empty()) {
          const auto [node, parent] = pending.back();
          pending.pop_back();
          if (gene.nodes[node].taxon() >= 0) side |= 1u << gene.nodes[node].taxon();
          for (int next : qraft::internal::gene_neighbors(gene, node)) {
            if (next != parent) pending.push_back({next, node});
          }
        }
        splits.insert(side);
        splits.insert(((1u << n) - 1) ^ side);
      }
    }
    // Every active subset, candidate, and anchor pair. This checks restriction
    // of original distances, which generally differs from recomputing unit
    // distances after suppressing inactive corridors.
    for (unsigned subset = 1; subset < (1u << n); ++subset) {
      if (std::popcount(subset) < 4) continue;
      for (unsigned mask = (subset - 1) & subset; mask;
           mask = (mask - 1) & subset) {
        bool displayed = false;
        for (unsigned full_split : splits) {
          displayed = displayed || (full_split & subset) == mask;
        }
        for (int a = 0; a < n; ++a) {
          if (!(mask & (1u << a))) continue;
          for (int b = 0; b < n; ++b) {
            if (!(subset & (1u << b)) || (mask & (1u << b))) continue;
            std::vector<int> order{a, b};
            for (int x = 0; x < n; ++x) {
              if (x != a && x != b && (subset & (1u << x))) order.push_back(x);
            }
            Bitset candidate(static_cast<int>(order.size()));
            int max_a = std::numeric_limits<int>::min();
            int min_b = std::numeric_limits<int>::max();
            for (int i = 0; i < static_cast<int>(order.size()); ++i) {
              const int x = order[i];
              const int z = distances[a][gene.taxon_to_node[x]] -
                            distances[b][gene.taxon_to_node[x]];
              if (mask & (1u << x)) {
                candidate.set(i);
                max_a = std::max(max_a, z);
              } else {
                min_b = std::min(min_b, z);
              }
            }
            const auto result = test_bipartition(ctx, order, candidate);
            require_close(result.gap, (min_b - max_a) / 2.0,
                          "bipartition-test gap disagrees with BFS");
            require(result.accepted == displayed,
                    "projection gap misclassified additive split");
          }
        }
      }
    }
  }

  // Equal quartet frequencies give a false bipartition a zero gap, which the
  // positive threshold rejects; thirty copies of each topology keep it zero.
  SolverContext ctx;
  ctx.options.bipartition_test.crossover = 0.0;
  ctx.taxa = qraft::TaxonRegistry({"A", "B", "C", "D"});
  ctx.options.bipartition_test.coefficient = 0.01;
  const std::vector<int> active{0, 1, 2, 3};
  Bitset ab(4);
  ab.set(0);
  ab.set(1);
  for (int i = 0; i < 30; ++i) {
    for (const auto* tree : {"((A,B),(C,D));", "((A,C),(B,D));", "((A,D),(B,C));"}) {
      ctx.genes.push_back(index_gene(qraft::parse_newick(tree), ctx.taxa));
    }
  }
  ctx.bipartition_test_index = make_bipartition_test_index(ctx.genes, 4);
  const auto tied = test_bipartition(ctx, active, ab);
  require(tied.available && !tied.accepted && tied.gap == 0.0, "zero-gap split was accepted");
  require(tied.threshold > 0.0, "bipartition-test threshold must be positive");
  require(bipartition_test_threshold(10000, 4, 0.05, 0.0) <
              bipartition_test_threshold(100, 4, 0.05, 0.0),
          "positive threshold tail does not decrease");
  require(bipartition_test_threshold(999, 4, 0.05, 1000.0) < 0.0 &&
              bipartition_test_threshold(1000, 4, 0.05, 1000.0) == 0.0 &&
              bipartition_test_threshold(1001, 4, 0.05, 1000.0) > 0.0,
          "shifted threshold has the wrong crossover signs");
  const double tail_a = bipartition_test_threshold(1000000, 4, 0.05, 10000.0);
  const double tail_b = bipartition_test_threshold(1000000000, 4, 0.05, 10000.0);
  require(tail_b > 0.0 && tail_b < tail_a &&
              tail_b * std::sqrt(1000000000.0) > tail_a * 1000.0,
          "threshold must decay while growing in standard-error units");
  const double near_k0 = std::nextafter(1000.0, 0.0);
  require(bipartition_test_threshold(1000, 4, 0.05, near_k0) > 0.0,
          "near-crossover cancellation erased the threshold sign");

  // Ten more AB|CD genes give topology frequencies 40/30/30: AB|CD has gap
  // 0.1 and the two alternatives nonpositive gaps. Checks the deterministic
  // parallel integer reductions for 1, 2 and 4 anchor pairs.
  const auto uniform_genes = ctx.genes;
  const auto uniform_index = ctx.bipartition_test_index;
  for (int i = 0; i < 10; ++i) {
    ctx.genes.push_back(index_gene(qraft::parse_newick("((A,B),(C,D));"), ctx.taxa));
  }
  ctx.bipartition_test_index = make_bipartition_test_index(ctx.genes, 4);
  for (int anchors : {1, 2, 4}) {
    ctx.options.bipartition_test.anchor_pairs = anchors;
    const auto positive = test_bipartition(ctx, active, ab);
    require(positive.accepted, "positive mixed-gene split was rejected");
    require_close(positive.gap, 0.1, "mixed-gene distance mean is incorrect");
    const int a = 0, c = 2;
    Bitset wrong(4);
    wrong.set(a);
    wrong.set(c);
    require(!test_bipartition(ctx, active, wrong).accepted,
            "nonpositive mixed-gene split was accepted");
  }

  // A nonadditive empirical mean can separate a bipartition for only some
  // anchor pairs. The test accepts when any tested pair separates it; here the
  // first two pairs do not.
  ctx.genes.clear();
  for (const auto& [count, tree] :
       std::vector<std::pair<int, std::string>>{
           {35, "((A,B),(C,D));"}, {25, "((A,C),(B,D));"},
           {40, "((A,D),(B,C));"}}) {
    for (int i = 0; i < count; ++i) {
      ctx.genes.push_back(index_gene(qraft::parse_newick(tree), ctx.taxa));
    }
  }
  ctx.bipartition_test_index = make_bipartition_test_index(ctx.genes, 4);
  ctx.options.bipartition_test.anchor_pairs = 2;
  require(!test_bipartition(ctx, active, ab).accepted,
          "first two anchor pairs should reject this empirical split");
  ctx.options.bipartition_test.anchor_pairs = 4;
  const auto rescued = test_bipartition(ctx, active, ab);
  require(rescued.accepted, "later anchor pair was not considered");
  require_close(rescued.gap, 0.1, "later anchor pair has the wrong gap");
  ctx.options.bipartition_test.anchor_pairs = 1;
  ctx.options.bipartition_test.coefficient = 0.1;
  ctx.options.bipartition_test.crossover = 1000000.0;
  const auto permissive = test_bipartition(ctx, active, ab);
  require(permissive.accepted && permissive.gap < 0.0 &&
              permissive.gap > permissive.threshold,
          "finite-sample offset did not accept a slightly negative gap");
  require_close(permissive.gap, -0.05, "negative empirical gap is incorrect");
  ctx.options.bipartition_test.crossover = 100.0;
  require(!test_bipartition(ctx, active, ab).accepted,
          "negative empirical gap passed at the crossover");
  ctx.options.bipartition_test.coefficient = 0.01;
  ctx.options.bipartition_test.crossover = 0.0;
  ctx.genes = uniform_genes;
  ctx.bipartition_test_index = uniform_index;

  // The bipartition test must use the original genes' distances even when
  // recursion installs a different active gene set.
  const std::vector<GeneIndex> unrelated{
      index_gene(qraft::parse_newick("((A,B),(C,D));"), ctx.taxa)};
  ctx.active_genes = &unrelated;
  require(test_bipartition(ctx, active, ab).gap == tied.gap,
          "bipartition test read restricted gene-set distances");
  for (const auto* tree : {"((A,A),(B,(C,D)));", "(A,B,C);"}) {
    const auto index = make_bipartition_test_index(
        {index_gene(qraft::parse_newick(tree), ctx.taxa)}, 4);
    require(!index->complete_single_copy, "test index silently used duplicate or missing taxa");
  }

  // A wrong, fully resolved starting tree cannot force a wrong quartet.
  const auto start_path = std::filesystem::temp_directory_path() / "qraft-bipartition-start.tre";
  { std::ofstream out(start_path); out << "((A,C),(B,D));\n"; }
  for (const auto regraft : {qraft::RegraftScore::Raw, qraft::RegraftScore::Normalized}) {
    for (const auto mode : {qraft::DivideScore::Quartet, qraft::DivideScore::Balanced}) {
      std::string first;
      for (int workers : {1, 4}) {
        qraft::Options options;
        options.threads = workers;
        options.starting_tree_path = start_path.string();
        options.regraft_score = regraft;
        options.divide_score = mode;
        // With 50 genes above the crossover the threshold is positive; below
        // the crossover a zero-gap false bipartition is accepted by design.
        options.bipartition_test.crossover = 10.0;
        options.nni = false;
        std::vector<qraft::GeneTree> genes(50, qraft::parse_newick("((A,B),(C,D));"));
        const auto inferred = qraft::infer_species_tree(genes, options);
        require(qraft::rf_distance(inferred, genes[0]).first == 0,
                "bipartition test plus singleton attachment failed on a wrong "
                "starting tree");
        const auto text = qraft::to_newick(inferred);
        if (first.empty()) first = text;
        else require(text == first, "wrong-start output changed with workers");

        // A completely wrong six-taxon starting tree must also be overcome
        // after gene-set restriction and regrafting on both sides.
        { std::ofstream out(start_path); out << "((A,C),(B,E),(D,F));\n"; }
        const auto truth = qraft::parse_newick("((A,B),(C,D),(E,F));");
        genes.assign(50, truth);
        const auto six = qraft::infer_species_tree(genes, options);
        require(qraft::rf_distance(six, truth).first == 0,
                "recursive bipartition test and regrafting failed on six taxa");
        options.nni = true;
        require(qraft::rf_distance(qraft::infer_species_tree(genes, options), truth).first == 0,
                "NNI moved the correct tree");
        { std::ofstream out(start_path); out << "((A,C),(B,D));\n"; }
      }
    }
  }
  std::filesystem::remove(start_path);
  for (double invalid : {0.0, -1.0, std::numeric_limits<double>::denorm_min(),
                         std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
    bool threw = false;
    try { (void)bipartition_test_threshold(10000, 4, invalid, 0.0); }
    catch (const std::invalid_argument&) { threw = true; }
    require(threw, "invalid or underflowed threshold coefficient accepted");
  }
  for (double invalid : {-1.0, std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
    bool threw = false;
    try { (void)bipartition_test_threshold(1000, 4, 0.01, invalid); }
    catch (const std::invalid_argument&) { threw = true; }
    require(threw, "invalid crossover accepted");
  }
}

// Random binary unrooted Newick over labels t00..t(n-1).
std::string random_binary_newick(std::mt19937& rng, int leaf_count) {
  std::vector<std::string> parts;
  for (int i = 0; i < leaf_count; ++i) {
    parts.push_back((i < 10 ? "t0" : "t") + std::to_string(i));
  }
  while (parts.size() > 3) {
    std::uniform_int_distribution<std::size_t> pick(0, parts.size() - 1);
    const std::size_t a = pick(rng);
    std::size_t b = pick(rng);
    while (b == a) b = pick(rng);
    const std::string joined = "(" + parts[a] + "," + parts[b] + ")";
    parts.erase(parts.begin() + static_cast<std::ptrdiff_t>(std::max(a, b)));
    parts.erase(parts.begin() + static_cast<std::ptrdiff_t>(std::min(a, b)));
    parts.push_back(joined);
  }
  return "(" + parts[0] + "," + parts[1] + "," + parts[2] + ");";
}

// Every edge bipartition of `tree`, as the sorted side without taxon 0.
std::set<std::vector<int>> brute_force_gene_splits(
    const qraft::Tree& tree, const qraft::TaxonRegistry& taxa) {
  std::set<std::vector<int>> splits;
  const int n = static_cast<int>(tree.nodes.size());
  for (int u = 0; u < n; ++u) {
    for (int v : tree.nodes[static_cast<std::size_t>(u)].neighbors) {
      std::vector<uint8_t> seen(static_cast<std::size_t>(n), 0);
      seen[static_cast<std::size_t>(u)] = 1;
      std::vector<int> stack{v};
      std::vector<int> side;
      while (!stack.empty()) {
        const int node = stack.back();
        stack.pop_back();
        if (seen[static_cast<std::size_t>(node)]) continue;
        seen[static_cast<std::size_t>(node)] = 1;
        const std::string& label = tree.nodes[static_cast<std::size_t>(node)].label;
        if (!label.empty()) side.push_back(taxa.id_for(label));
        for (int next : tree.nodes[static_cast<std::size_t>(node)].neighbors) {
          stack.push_back(next);
        }
      }
      if (std::find(side.begin(), side.end(), 0) != side.end()) continue;
      std::sort(side.begin(), side.end());
      if (side.size() >= 2 &&
          side.size() + 2 <= static_cast<std::size_t>(taxa.size())) {
        splits.insert(side);
      }
    }
  }
  return splits;
}

void test_nni_split_support() {
  using namespace qraft::internal;
  // The NNI gate reads gene-tree support for either orientation of a
  // bipartition. Compare it with explicit enumeration of gene splits,
  // including sides without taxon 0 that hold more than half the taxa.
  std::mt19937 rng(20260925u);
  for (int leaf_count : {6, 7, 9, 12}) {
    SolverContext ctx;
    std::vector<qraft::Tree> trees;
    for (int g = 0; g < 24; ++g) {
      // Repeat topologies so that supports above one occur.
      if (g % 3 == 2) {
        trees.push_back(trees.back());
      } else {
        trees.push_back(qraft::parse_newick(random_binary_newick(rng, leaf_count)));
      }
    }
    ctx.taxa = qraft::TaxonRegistry(trees.front().leaf_labels());
    std::map<std::vector<int>, int> expected;
    for (const qraft::Tree& tree : trees) {
      ctx.genes.push_back(index_gene(tree, ctx.taxa));
      for (const std::vector<int>& split : brute_force_gene_splits(tree, ctx.taxa)) {
        ++expected[split];
      }
    }
    int larger_sides = 0;
    for (const auto& [side, support] : expected) {
      std::vector<int> other;
      for (int taxon = 0; taxon < ctx.taxa.size(); ++taxon) {
        if (!std::binary_search(side.begin(), side.end(), taxon)) {
          other.push_back(taxon);
        }
      }
      larger_sides += side.size() > other.size() ? 1 : 0;
      const std::string context = "n=" + std::to_string(leaf_count) +
                                  " side size " + std::to_string(side.size());
      require(nni_split_support_for_tests(ctx, side) == support,
              "NNI split support disagrees with enumeration: " + context);
      require(nni_split_support_for_tests(ctx, other) == support,
              "NNI split support depends on orientation: " + context);
    }
    require(larger_sides > 0, "fixture lacks a side larger than its complement");
  }
}

// ---------------------------------------------------------------------------
// Direct quartet enumeration over physical leaves, for solver-level checks.
// ---------------------------------------------------------------------------

// Random unrooted tree over `labels` (duplicates allowed); with `polytomies`,
// some joins take three parts.
std::string random_newick_over(std::mt19937& rng,
                               const std::vector<std::string>& labels,
                               bool polytomies) {
  std::vector<std::string> parts = labels;
  while (parts.size() > 3) {
    std::shuffle(parts.begin(), parts.end(), rng);
    const std::size_t take =
        polytomies && parts.size() > 4 && rng() % 4 == 0 ? 3 : 2;
    std::string joined = "(" + parts[0];
    for (std::size_t i = 1; i < take; ++i) joined += "," + parts[i];
    joined += ")";
    parts.erase(parts.begin(), parts.begin() + static_cast<std::ptrdiff_t>(take));
    parts.push_back(joined);
  }
  std::string out = "(" + parts[0];
  for (std::size_t i = 1; i < parts.size(); ++i) out += "," + parts[i];
  return out + ");";
}

std::vector<std::string> numbered_labels(int count) {
  std::vector<std::string> labels;
  for (int i = 0; i < count; ++i) {
    labels.push_back((i < 10 ? "t0" : "t") + std::to_string(i));
  }
  return labels;
}

// The nontrivial unrooted splits of a tree with distinct labels, each stored
// as the side without the smallest label.
std::set<std::set<std::string>> label_splits(const qraft::Tree& tree) {
  std::set<std::string> all;
  for (const qraft::Tree::Node& node : tree.nodes) {
    if (!node.label.empty()) all.insert(node.label);
  }
  std::set<std::set<std::string>> splits;
  if (tree.nodes.empty() || all.size() < 4) return splits;
  // Leaf sets below every node, rooted at node 0.
  std::vector<std::set<std::string>> below(tree.nodes.size());
  std::vector<int> parent(tree.nodes.size(), -1);
  std::vector<int> order{0};
  parent[0] = 0;
  for (std::size_t i = 0; i < order.size(); ++i) {
    for (int nb : tree.nodes[static_cast<std::size_t>(order[i])].neighbors) {
      if (parent[static_cast<std::size_t>(nb)] >= 0) continue;
      parent[static_cast<std::size_t>(nb)] = order[i];
      order.push_back(nb);
    }
  }
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    const int node = *it;
    const std::string& label = tree.nodes[static_cast<std::size_t>(node)].label;
    if (!label.empty()) below[static_cast<std::size_t>(node)].insert(label);
    if (node == 0) continue;
    const int up = parent[static_cast<std::size_t>(node)];
    below[static_cast<std::size_t>(up)].insert(
        below[static_cast<std::size_t>(node)].begin(),
        below[static_cast<std::size_t>(node)].end());
    std::set<std::string> side = below[static_cast<std::size_t>(node)];
    if (side.size() < 2 || all.size() - side.size() < 2) continue;
    if (side.count(*all.begin()) != 0) {
      std::set<std::string> other;
      std::set_difference(all.begin(), all.end(), side.begin(), side.end(),
                          std::inserter(other, other.begin()));
      side = std::move(other);
    }
    splits.insert(std::move(side));
  }
  return splits;
}

// Newick of `tree` induced on the leaves in `keep`, unary nodes suppressed.
std::string restricted_newick(const qraft::Tree& tree,
                              const std::set<std::string>& keep) {
  std::function<std::string(int, int)> emit = [&](int node, int from) {
    const qraft::Tree::Node& current = tree.nodes[static_cast<std::size_t>(node)];
    if (!current.label.empty()) {
      return keep.count(current.label) != 0 ? current.label : std::string();
    }
    std::vector<std::string> parts;
    for (int nb : current.neighbors) {
      if (nb == from) continue;
      std::string part = emit(nb, node);
      if (!part.empty()) parts.push_back(std::move(part));
    }
    if (parts.empty()) return std::string();
    if (parts.size() == 1) return parts.front();
    std::string out = "(" + parts[0];
    for (std::size_t i = 1; i < parts.size(); ++i) out += "," + parts[i];
    return out + ")";
  };
  int start = 0;
  while (!tree.nodes[static_cast<std::size_t>(start)].label.empty()) ++start;
  return emit(start, -1) + ";";
}

// Gene trees that miss taxa. A split proposed by such a gene is judged by the
// taxa the gene has: the starting tree may neither gain a node with one child
// nor a cluster that no gene shows. The whole pipeline runs on such genes,
// including genes that keep at most two of a subproblem's taxa, and its
// output does not depend on the worker count.
// Starts and stops the MRE's gene update stream many times on four workers.
// A helper that checked for the stop before reading the wake count could miss
// the stop and sleep for ever, so the run never returned.
void test_mre_update_stream_stop() {
  using namespace qraft::internal;
  qraft::Options options = test_options(4);
  require(configure_threads(options) == 4, "failed to configure four workers");
  const std::vector<std::string> newicks = {
      "((a,(c,d,e,f)),(((t,y),w,z),b));", "((a,(c,d,e,f)),(((t,y),w,z),b));",
      "((a,(c,d,e,f)),((y,w,b),z));", "((a,(c,d,e,f)),(((y,w),b),z));"};
  std::vector<qraft::Tree> trees;
  std::set<std::string> labels;
  for (const std::string& text : newicks) {
    trees.push_back(qraft::parse_newick(text));
    for (const std::string& label : trees.back().leaf_labels()) {
      labels.insert(label);
    }
  }
  const qraft::TaxonRegistry taxa(
      std::vector<std::string>(labels.begin(), labels.end()));
  std::vector<GeneIndex> genes;
  for (const qraft::Tree& tree : trees) genes.push_back(index_gene(tree, taxa));
  for (int round = 0; round < 5000; ++round) {
    const StartingTree tree = build_mre_tree(genes, taxa.size(), options);
    require(tree.root >= 0, "the MRE returned no starting tree");
  }
  options.threads = 1;
  require(configure_threads(options) == 1, "failed to restore one worker");
}

void test_incomplete_gene_trees() {
  using namespace qraft::internal;
  auto starting_tree_of = [](const std::vector<std::string>& newicks,
                             int workers) {
    std::vector<qraft::Tree> trees;
    std::set<std::string> labels;
    for (const std::string& text : newicks) {
      trees.push_back(qraft::parse_newick(text));
      for (const std::string& label : trees.back().leaf_labels()) {
        labels.insert(label);
      }
    }
    SolverContext ctx;
    ctx.options = test_options(workers);
    ctx.workers = configure_threads(ctx.options);
    ctx.taxa = qraft::TaxonRegistry(
        std::vector<std::string>(labels.begin(), labels.end()));
    for (const qraft::Tree& tree : trees) {
      ctx.genes.push_back(index_gene(tree, ctx.taxa));
    }
    return std::pair{build_mre_tree(ctx.genes, ctx.taxa.size(), ctx.options),
                     ctx.taxa};
  };
  // Checks the tree's shape and that every cluster, restricted to the taxa of
  // some gene, is a nontrivial split of that gene.
  auto require_supported = [](const StartingTree& tree,
                              const qraft::TaxonRegistry& taxa,
                              const std::vector<std::string>& newicks,
                              const std::string& label) {
    std::vector<std::set<std::string>> gene_taxa;
    std::vector<std::set<std::set<std::string>>> gene_splits;
    for (const std::string& text : newicks) {
      const qraft::Tree gene = qraft::parse_newick(text);
      const std::vector<std::string> leaves = gene.leaf_labels();
      gene_taxa.emplace_back(leaves.begin(), leaves.end());
      gene_splits.push_back(label_splits(gene));
    }
    std::vector<std::set<std::string>> cluster(tree.nodes.size());
    for (auto it = tree.preorder.rbegin(); it != tree.preorder.rend(); ++it) {
      const StartingTree::Node& node = tree.nodes[static_cast<std::size_t>(*it)];
      if (node.taxon >= 0) {
        cluster[static_cast<std::size_t>(*it)].insert(taxa.label_for(node.taxon));
        continue;
      }
      require(node.children.size() >= 2,
              label + ": starting tree has a node with one child");
      for (int child : node.children) {
        cluster[static_cast<std::size_t>(*it)].insert(
            cluster[static_cast<std::size_t>(child)].begin(),
            cluster[static_cast<std::size_t>(child)].end());
      }
      if (*it == tree.root) continue;
      bool shown = false;
      for (std::size_t g = 0; g < newicks.size() && !shown; ++g) {
        std::set<std::string> side;
        std::set_intersection(
            cluster[static_cast<std::size_t>(*it)].begin(),
            cluster[static_cast<std::size_t>(*it)].end(), gene_taxa[g].begin(),
            gene_taxa[g].end(), std::inserter(side, side.begin()));
        if (side.size() < 2 || gene_taxa[g].size() - side.size() < 2) continue;
        if (side.count(*gene_taxa[g].begin()) != 0) {
          std::set<std::string> other;
          std::set_difference(gene_taxa[g].begin(), gene_taxa[g].end(),
                              side.begin(), side.end(),
                              std::inserter(other, other.begin()));
          side = std::move(other);
        }
        shown = gene_splits[g].count(side) != 0;
      }
      require(shown, label + ": starting-tree cluster is shown by no gene");
    }
  };

  const std::string complete = "((a,(c,d,e,f)),(((t,y),w,z),b));";
  const std::vector<std::string> base(3, complete);
  // Split {w,b,y} of a gene without t conflicts with {t,w,y,z}: {t,y} must not
  // stand in for the missing b.
  {
    std::vector<std::string> genes = base;
    genes.push_back("((a,(c,d,e,f)),((y,w,b),z));");
    const auto [tree, taxa] = starting_tree_of(genes, 1);
    require_supported(tree, taxa, genes, "conflicting incomplete split");
    std::vector<qraft::GeneTree> input;
    for (const std::string& text : genes) input.push_back(qraft::parse_newick(text));
    qraft::Options options = test_options(1);
    options.divide = false;
    options.nni = false;
    const auto [rf, nrf] = qraft::rf_distance(
        qraft::infer_species_tree(std::move(input), options),
        qraft::parse_newick(complete));
    (void)nrf;
    require(rf == 0, "a conflicting incomplete split entered the starting tree");
  }
  // Split {w,y} of a gene without t is compatible and adds {t,w,y}.
  {
    std::vector<std::string> genes = base;
    genes.push_back("((a,(c,d,e,f)),(((y,w),b),z));");
    const auto [tree, taxa] = starting_tree_of(genes, 1);
    require_supported(tree, taxa, genes, "compatible incomplete split");
    std::vector<qraft::GeneTree> input;
    for (const std::string& text : genes) input.push_back(qraft::parse_newick(text));
    qraft::Options options = test_options(1);
    options.divide = false;
    options.nni = false;
    const auto [rf, nrf] = qraft::rf_distance(
        qraft::infer_species_tree(std::move(input), options),
        qraft::parse_newick("((a,(c,d,e,f)),((((t,y),w),z),b));"));
    (void)nrf;
    require(rf == 0, "a compatible incomplete split was not accepted");
  }
  // Split {b,y} of a gene without t equals {t,y} there; it must not be added
  // again above {t,y}.
  {
    std::vector<std::string> genes(3, "((a,(c,d,e,f)),(((t,y),z),b));");
    genes.push_back("((a,(c,d,e,f)),((y,b),z));");
    const auto [tree, taxa] = starting_tree_of(genes, 1);
    require_supported(tree, taxa, genes, "redundant incomplete split");
  }

  // Random genes, some complete and some missing taxa, drawn around a common
  // tree so that splits recur.
  std::mt19937 rng(20260925u);
  for (int trial = 0; trial < 300; ++trial) {
    const int n = 5 + static_cast<int>(rng() % 6);
    std::vector<std::string> labels;
    for (int i = 0; i < n; ++i) labels.push_back("x" + std::to_string(i));
    const qraft::Tree species =
        qraft::parse_newick(random_newick_over(rng, labels, true));
    std::vector<std::string> genes;
    const int k = 2 + static_cast<int>(rng() % 7);
    for (int g = 0; g < k; ++g) {
      std::set<std::string> keep;
      for (const std::string& label : labels) {
        if (rng() % 4 != 0) keep.insert(label);
      }
      if (keep.size() < 4) keep.insert(labels.begin(), labels.end());
      if (rng() % 3 == 0) {
        genes.push_back(random_newick_over(
            rng, std::vector<std::string>(keep.begin(), keep.end()), true));
      } else {
        genes.push_back(restricted_newick(species, keep));
      }
    }
    std::string first_output;
    for (int workers : {1, 4}) {
      const auto [tree, taxa] = starting_tree_of(genes, workers);
      require_supported(tree, taxa, genes,
                        "random incomplete genes, trial " +
                            std::to_string(trial));
      std::vector<qraft::GeneTree> input;
      for (const std::string& text : genes) {
        input.push_back(qraft::parse_newick(text));
      }
      qraft::Options options = test_options(workers);
      options.divide = false;
      options.nni = false;
      const std::string output =
          qraft::to_newick(qraft::infer_species_tree(std::move(input), options));
      if (first_output.empty()) {
        first_output = output;
      } else {
        require(output == first_output,
                "MRE on incomplete genes changed with the worker count");
      }
    }
  }

  // The default pipeline and its ablations on 40 taxa, with genes that keep
  // anything from all taxa down to two.
  for (int trial = 0; trial < 6; ++trial) {
    std::vector<std::string> labels;
    for (int i = 0; i < 40; ++i) labels.push_back("y" + std::to_string(i));
    const qraft::Tree species =
        qraft::parse_newick(random_newick_over(rng, labels, false));
    std::vector<std::string> genes;
    for (int g = 0; g < 30; ++g) {
      std::vector<std::string> kept = labels;
      std::shuffle(kept.begin(), kept.end(), rng);
      const std::size_t sizes[] = {40, 30, 12, 6, 3, 2};
      kept.resize(sizes[g % 6]);
      genes.push_back(restricted_newick(
          species, std::set<std::string>(kept.begin(), kept.end())));
    }
    for (const auto& [divide, nni, raw] :
         {std::tuple{true, true, false}, std::tuple{true, false, true},
          std::tuple{false, true, false}}) {
      std::string first_output;
      for (int workers : {1, 4}) {
        std::vector<qraft::GeneTree> input;
        for (const std::string& text : genes) {
          input.push_back(qraft::parse_newick(text));
        }
        qraft::Options options = test_options(workers);
        options.divide = divide;
        options.nni = nni;
        options.regraft_score =
            raw ? qraft::RegraftScore::Raw : qraft::RegraftScore::Normalized;
        const qraft::SpeciesTree inferred =
            qraft::infer_species_tree(std::move(input), options);
        require(inferred.leaf_labels().size() == labels.size(),
                "pipeline on incomplete genes lost a taxon");
        const std::string output = qraft::to_newick(inferred);
        if (first_output.empty()) {
          first_output = output;
        } else {
          require(output == first_output,
                  "pipeline on incomplete genes changed with the worker count");
        }
      }
    }
  }
}

// Physical leaves of a tree with their taxa and unit-edge path lengths.
struct PhysicalLeaves {
  std::vector<int> taxon;
  std::vector<std::vector<int>> distance;
};

template <typename Neighbors, typename TaxonOf>
PhysicalLeaves physical_leaves_of(int node_count,
                                  Neighbors&& neighbors,
                                  TaxonOf&& taxon_of) {
  PhysicalLeaves out;
  std::vector<int> leaf_nodes;
  for (int node = 0; node < node_count; ++node) {
    const int taxon = taxon_of(node);
    if (taxon >= 0) {
      leaf_nodes.push_back(node);
      out.taxon.push_back(taxon);
    }
  }
  for (int source : leaf_nodes) {
    std::vector<int> depth(static_cast<std::size_t>(node_count), -1);
    std::queue<int> pending;
    depth[static_cast<std::size_t>(source)] = 0;
    pending.push(source);
    while (!pending.empty()) {
      const int node = pending.front();
      pending.pop();
      for (int next : neighbors(node)) {
        if (depth[static_cast<std::size_t>(next)] >= 0) continue;
        depth[static_cast<std::size_t>(next)] =
            depth[static_cast<std::size_t>(node)] + 1;
        pending.push(next);
      }
    }
    std::vector<int> row;
    for (int target : leaf_nodes) row.push_back(depth[static_cast<std::size_t>(target)]);
    out.distance.push_back(std::move(row));
  }
  return out;
}

PhysicalLeaves physical_leaves(const qraft::Tree& tree,
                               const qraft::TaxonRegistry& taxa) {
  return physical_leaves_of(
      static_cast<int>(tree.nodes.size()),
      [&](int node) { return tree.nodes[static_cast<std::size_t>(node)].neighbors; },
      [&](int node) {
        const std::string& label = tree.nodes[static_cast<std::size_t>(node)].label;
        return label.empty() ? -1 : taxa.id_for(label);
      });
}

PhysicalLeaves physical_leaves(const qraft::internal::SolverTree& tree) {
  return physical_leaves_of(
      static_cast<int>(tree.nodes.size()),
      [&](int node) { return tree.nodes[static_cast<std::size_t>(node)].neighbors; },
      [&](int node) {
        return static_cast<int>(tree.nodes[static_cast<std::size_t>(node)].taxon);
      });
}

// Resolution of leaves (a,b,c,d): 0 for ab|cd, 1 for ac|bd, 2 for ad|bc, and
// -1 for a star, by the four-point condition on path lengths.
// The resolution of leaves a, b, c, d: 0 for ab|cd, 1 for ac|bd, 2 for
// ad|bc, and -1 for a star or when a leaf is missing (-1).
int quartet_resolution(const PhysicalLeaves& leaves, int a, int b, int c, int d) {
  if (a < 0 || b < 0 || c < 0 || d < 0) return -1;
  const auto& D = leaves.distance;
  const int s0 = D[a][b] + D[c][d];
  const int s1 = D[a][c] + D[b][d];
  const int s2 = D[a][d] + D[b][c];
  if (s0 < s1 && s0 < s2) return 0;
  if (s1 < s0 && s1 < s2) return 1;
  if (s2 < s0 && s2 < s1) return 2;
  return -1;
}

std::vector<int> iota_taxa(int taxon_count) {
  std::vector<int> taxa(static_cast<std::size_t>(taxon_count));
  std::iota(taxa.begin(), taxa.end(), 0);
  return taxa;
}

// Physical leaf of each taxon in a single-copy tree, or -1 when it is missing.
std::vector<int> leaf_of_taxon(const PhysicalLeaves& leaves, int taxon_count) {
  std::vector<int> out(static_cast<std::size_t>(taxon_count), -1);
  for (int leaf = 0; leaf < static_cast<int>(leaves.taxon.size()); ++leaf) {
    out[static_cast<std::size_t>(leaves.taxon[static_cast<std::size_t>(leaf)])] = leaf;
  }
  return out;
}

qraft::internal::SolverTree solver_tree_from(const qraft::Tree& tree,
                                            const qraft::TaxonRegistry& taxa) {
  qraft::internal::SolverTree out;
  out.nodes.resize(tree.nodes.size());
  for (std::size_t node = 0; node < tree.nodes.size(); ++node) {
    out.nodes[node].neighbors = tree.nodes[node].neighbors;
    out.nodes[node].taxon =
        tree.nodes[node].label.empty() ? -1 : taxa.id_for(tree.nodes[node].label);
  }
  return out;
}

// A random gene over all of `labels` (one gene in three), over about three
// quarters of them, or over only two to four, so that genes can miss most
// of a subproblem's taxa.
std::string random_incomplete_gene(std::mt19937& rng,
                                   const std::vector<std::string>& labels,
                                   int g) {
  std::vector<std::string> kept = labels;
  std::shuffle(kept.begin(), kept.end(), rng);
  if (g % 3 == 1) {
    kept.resize(std::max<std::size_t>(4, labels.size() * 3 / 4));
  } else if (g % 3 == 2) {
    kept.resize(2 + rng() % 3);
  }
  return random_newick_over(rng, kept, g % 5 == 0);
}

// The NNI support gate (NniExactSupport) against explicit enumeration, for
// arbitrary sides and for sequences of applied moves.

// The b-side taxa of the edge (a, b) of `tree`, sorted.
std::vector<int> solver_side_taxa(const qraft::internal::SolverTree& tree,
                                  int a,
                                  int b) {
  std::vector<int> taxa;
  std::vector<uint8_t> seen(tree.nodes.size(), 0);
  seen[static_cast<std::size_t>(a)] = 1;
  std::vector<int> stack{b};
  while (!stack.empty()) {
    const int v = stack.back();
    stack.pop_back();
    if (seen[static_cast<std::size_t>(v)]) continue;
    seen[static_cast<std::size_t>(v)] = 1;
    if (tree.nodes[static_cast<std::size_t>(v)].taxon >= 0) {
      taxa.push_back(tree.nodes[static_cast<std::size_t>(v)].taxon);
    }
    for (int nb : tree.nodes[static_cast<std::size_t>(v)].neighbors) {
      stack.push_back(nb);
    }
  }
  std::sort(taxa.begin(), taxa.end());
  return taxa;
}

// The gate's support of `side | rest` by enumeration: -1 when trivial, else
// the number of genes with an eligible side equal to the side without taxon 0.
int enumerated_nni_support(const std::vector<std::set<std::vector<int>>>& sides,
                           int n,
                           std::vector<int> side) {
  if (n < 4) return -1;
  const int side_size = static_cast<int>(side.size());
  if (side_size <= 1 || n - side_size <= 1) return -1;
  // A side listing a taxon twice, or one outside [0, n), equals no gene side.
  std::vector<int> sorted = side;
  std::sort(sorted.begin(), sorted.end());
  const bool listed_once =
      std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end() &&
      (sorted.empty() || (sorted.front() >= 0 && sorted.back() < n));
  const bool contains_ref = std::find(side.begin(), side.end(), 0) != side.end();
  // The trivial-side test counts the listed taxa.
  const int size = contains_ref ? n - side_size : side_size;
  if (size < 2 || size > n - 2) return -1;
  if (!listed_once) return 0;
  if (contains_ref) {
    std::vector<int> rest;
    for (int taxon = 0; taxon < n; ++taxon) {
      if (std::find(side.begin(), side.end(), taxon) == side.end()) {
        rest.push_back(taxon);
      }
    }
    side = rest;
  }
  std::sort(side.begin(), side.end());
  int support = 0;
  for (const std::set<std::vector<int>>& gene : sides) {
    support += gene.count(side) != 0 ? 1 : 0;
  }
  return support;
}

void test_nni_exact_split_support() {
  using namespace qraft::internal;
  std::mt19937 rng(20261031u);
  int traced_moves = 0;
  int root_arm_moves = 0;
  int unlabeled_leaf_trees = 0;
  int fallback_trees = 0;
  for (int round = 0; round < 120; ++round) {
    const int n = 4 + static_cast<int>(rng() % 9);
    std::vector<std::string> labels;
    for (int x = 0; x < n; ++x) labels.push_back(exact_test_label(x));
    // Genes: complete, incomplete (the smallest taxon too), with repeated
    // taxa (the smallest too), with polytomies; repeated topologies.
    std::vector<std::string> newicks;
    newicks.push_back(random_newick_over(rng, labels, false));
    for (int g = 1; g < 10; ++g) {
      if (rng() % 4 == 0) {
        newicks.push_back(newicks[rng() % newicks.size()]);
      } else if (g % 3 == 0) {
        std::vector<std::string> items;
        for (const std::string& label : labels) {
          if (rng() % 6 == 0) continue;
          items.push_back(label);
          if (rng() % 5 == 0) items.push_back(label);
        }
        while (items.size() < 2) items.push_back(labels[rng() % labels.size()]);
        newicks.push_back(random_multifurcating_newick(rng, items));
      } else {
        newicks.push_back(random_incomplete_gene(rng, labels, g));
      }
    }
    SolverContext ctx;
    ctx.taxa = qraft::TaxonRegistry(labels);
    std::vector<std::set<std::vector<int>>> eligible;
    for (const std::string& newick : newicks) {
      ctx.genes.push_back(index_gene(qraft::parse_newick(newick), ctx.taxa));
      int distinct = 0;
      std::set<std::vector<int>> gene_sides;
      for (const ExplicitGeneSide& side :
           explicit_gene_sides(ctx.genes.back(), distinct)) {
        const int size = static_cast<int>(side.taxa.size());
        if (distinct >= 4 && side.closed && size >= 2 && size <= distinct - 2) {
          gene_sides.insert(side.taxa);
        }
      }
      eligible.push_back(gene_sides);
    }
    const std::string context = "round " + std::to_string(round);

    // Arbitrary sides: the gene sides, their complements and random subsets.
    std::vector<std::vector<int>> queries;
    for (const auto& gene : eligible) {
      for (const std::vector<int>& side : gene) {
        queries.push_back(side);
        std::vector<int> rest;
        for (int taxon = 0; taxon < n; ++taxon) {
          if (!std::binary_search(side.begin(), side.end(), taxon)) {
            rest.push_back(taxon);
          }
        }
        queries.push_back(rest);
      }
    }
    for (int q = 0; q < 20; ++q) {
      std::vector<int> side;
      for (int taxon = 0; taxon < n; ++taxon) {
        if (rng() % 2 == 0) side.push_back(taxon);
      }
      queries.push_back(side);
    }
    // Sides that list a taxon outside the range or twice support nothing.
    queries.push_back({1, 2, n});
    queries.push_back({0, 1, n});
    queries.push_back({1, 1, 2});
    for (const std::vector<int>& side : queries) {
      const int want = enumerated_nni_support(eligible, n, side);
      require(nni_split_support_for_tests(ctx, side) == want,
              "NNI support differs from enumeration: " + context);
    }

    // Move traces: moves picked one at a time from the scored alternatives of
    // the current tree, so that every move applies after the earlier ones.
    std::string species = random_newick_over(rng, labels, false);
    if (round % 4 == 3) {
      // An outer unary root leaves an active unlabeled leaf next to a
      // central node.
      species = "(" + species.substr(0, species.size() - 1) + ");";
      ++unlabeled_leaf_trees;
    }
    const SolverTree initial =
        solver_tree_from(qraft::parse_newick(species), ctx.taxa);
    std::vector<int> taxa(static_cast<std::size_t>(n));
    std::iota(taxa.begin(), taxa.end(), 0);
    SolverTree current = initial;
    std::vector<NniMove> sequence;
    std::vector<int> want;
    for (int step = 0; step < 12; ++step) {
      const std::vector<NniMove> moves =
          nni_scored_moves_for_tests(ctx, current, taxa);
      if (moves.empty()) break;
      const NniMove move = moves[rng() % moves.size()];
      sequence.push_back(move);
      const std::vector<int> old_side =
          solver_side_taxa(current, move.edge_a, move.edge_b);
      // The arm holding taxon 0 moves when it is a swapped neighbour's side.
      for (int arm : {move.swap_a_neighbor, move.swap_b_neighbor}) {
        const int center = arm == move.swap_a_neighbor ? move.edge_a : move.edge_b;
        const std::vector<int> arm_side = solver_side_taxa(current, center, arm);
        if (std::binary_search(arm_side.begin(), arm_side.end(), 0)) {
          ++root_arm_moves;
        }
      }
      remove_solver_edge(current, move.edge_a, move.swap_a_neighbor);
      remove_solver_edge(current, move.edge_b, move.swap_b_neighbor);
      add_solver_edge(current, move.edge_a, move.swap_b_neighbor);
      add_solver_edge(current, move.edge_b, move.swap_a_neighbor);
      const std::vector<int> new_side =
          solver_side_taxa(current, move.edge_a, move.edge_b);
      const int old_support = enumerated_nni_support(eligible, n, old_side);
      const int new_support = enumerated_nni_support(eligible, n, new_side);
      want.push_back(old_support < 0 || new_support < 0
                         ? 0
                         : new_support - old_support);
    }
    traced_moves += static_cast<int>(sequence.size());
    SolverTree traced = initial;
    require(nni_support_gate_trace_for_tests(ctx, traced, sequence) == want,
            "NNI gate deltas differ from enumeration: " + context);

    // The same moves on a tree that repeats a taxon or misses one (the gate
    // then reads the taxa of each side directly): relabel one leaf.
    SolverTree variant = initial;
    std::vector<int> leaves;
    for (int v = 0; v < static_cast<int>(variant.nodes.size()); ++v) {
      if (variant.nodes[static_cast<std::size_t>(v)].taxon >= 0) leaves.push_back(v);
    }
    const int relabeled = leaves[rng() % leaves.size()];
    variant.nodes[static_cast<std::size_t>(relabeled)].taxon =
        round % 2 == 0 ? static_cast<int>(rng() % n) : -1;
    ++fallback_trees;
    SolverTree variant_current = variant;
    std::vector<int> variant_want;
    for (const NniMove& move : sequence) {
      const std::vector<int> old_side =
          solver_side_taxa(variant_current, move.edge_a, move.edge_b);
      remove_solver_edge(variant_current, move.edge_a, move.swap_a_neighbor);
      remove_solver_edge(variant_current, move.edge_b, move.swap_b_neighbor);
      add_solver_edge(variant_current, move.edge_a, move.swap_b_neighbor);
      add_solver_edge(variant_current, move.edge_b, move.swap_a_neighbor);
      const std::vector<int> new_side =
          solver_side_taxa(variant_current, move.edge_a, move.edge_b);
      const int old_support = enumerated_nni_support(eligible, n, old_side);
      const int new_support = enumerated_nni_support(eligible, n, new_side);
      variant_want.push_back(old_support < 0 || new_support < 0 ||
                                     old_side == new_side
                                 ? 0
                                 : new_support - old_support);
    }
    require(nni_support_gate_trace_for_tests(ctx, variant, sequence) ==
                variant_want,
            "NNI gate deltas differ from enumeration on a tree with a "
            "repeated or missing taxon: " + context);
  }
  require(traced_moves > 500 && root_arm_moves > 20 && unlabeled_leaf_trees > 10 &&
              fallback_trees > 50,
          "NNI gate fixture lacks moves of the arm holding taxon 0");
}

void test_bipartition_scores_match_enumeration() {
  using namespace qraft::internal;
  // Every candidate bipartition's score is the satisfied-minus-violated count
  // over physical 2+2 leaf sets, including duplicated and missing taxa, and
  // with taxa outside the subproblem.
  std::mt19937 rng(20260926u);
  for (int taxon_count : {6, 8, 10, 12}) {
    const std::vector<std::string> labels = numbered_labels(taxon_count);
    const std::string base = random_newick_over(rng, labels, false);
    std::vector<qraft::Tree> trees;
    for (int g = 0; g < 18; ++g) {
      std::vector<std::string> gene_labels = labels;
      if (g % 3 == 1) gene_labels.push_back(labels[rng() % labels.size()]);
      if (g % 4 == 3) gene_labels.erase(gene_labels.begin() + (rng() % labels.size()));
      const bool copy_base = g % 2 == 0 && gene_labels.size() == labels.size();
      trees.push_back(qraft::parse_newick(
          copy_base ? base : random_newick_over(rng, gene_labels, g % 5 == 0)));
    }
    for (int g = 0; g < 9; ++g) {
      trees.push_back(qraft::parse_newick(random_incomplete_gene(rng, labels, g)));
    }
    SolverContext ctx;
    ctx.options = test_options(2);
    ctx.workers = configure_threads(ctx.options);
    ctx.taxa = qraft::TaxonRegistry(labels);
    for (const qraft::Tree& tree : trees) ctx.genes.push_back(index_gene(tree, ctx.taxa));
    ctx.starting_tree = build_mre_tree(ctx.genes, taxon_count, ctx.options);
    // From 10 taxa on, the subproblem leaves out two taxa.
    std::vector<int> active_taxa = iota_taxa(taxon_count);
    if (taxon_count >= 10) active_taxa.resize(active_taxa.size() - 2);
    std::vector<std::uint8_t> is_active(static_cast<std::size_t>(taxon_count), 0);
    for (int taxon : active_taxa) is_active[static_cast<std::size_t>(taxon)] = 1;
    const int active_count = static_cast<int>(active_taxa.size());
    const ActiveTaxa active = make_active_taxa(active_taxa, taxon_count);
    const QueryTree query =
        induce_query_tree(ctx.starting_tree, active);
    std::vector<int> edges;
    for (int node : query.preorder) {
      const int size = query.subtree_leaf_count[static_cast<std::size_t>(node)];
      if (node != query.root && size >= 2 && size <= active_count - 2) {
        edges.push_back(node);
      }
    }
    require(!edges.empty(), "divide fixture has no candidate bipartition");
    const std::vector<double> scores =
        score_bipartitions(ctx, active, query, edges);
    std::vector<PhysicalLeaves> genes;
    for (const qraft::Tree& tree : trees) genes.push_back(physical_leaves(tree, ctx.taxa));
    for (std::size_t e = 0; e < edges.size(); ++e) {
      std::vector<std::uint8_t> left(static_cast<std::size_t>(taxon_count), 0);
      const int first = query.subtree_leaf_first[static_cast<std::size_t>(edges[e])];
      const int count = query.subtree_leaf_count[static_cast<std::size_t>(edges[e])];
      for (int i = first; i < first + count; ++i) {
        left[static_cast<std::size_t>(active.taxa[static_cast<std::size_t>(
            query.leaf_active_order[static_cast<std::size_t>(i)])])] = 1;
      }
      double expected = 0.0;
      for (const PhysicalLeaves& gene : genes) {
        const int n = static_cast<int>(gene.taxon.size());
        for (int a = 0; a < n; ++a)
          for (int b = a + 1; b < n; ++b)
            for (int c = b + 1; c < n; ++c)
              for (int d = c + 1; d < n; ++d) {
                std::array<int, 4> leaves{a, b, c, d};
                if (std::any_of(leaves.begin(), leaves.end(), [&](int leaf) {
                      return !is_active[static_cast<std::size_t>(
                          gene.taxon[static_cast<std::size_t>(leaf)])];
                    })) {
                  continue;
                }
                std::array<int, 4> in{}, out{};
                int ins = 0, outs = 0;
                for (int leaf : leaves) {
                  if (left[static_cast<std::size_t>(gene.taxon[static_cast<std::size_t>(leaf)])]) {
                    in[static_cast<std::size_t>(ins++)] = leaf;
                  } else {
                    out[static_cast<std::size_t>(outs++)] = leaf;
                  }
                }
                if (ins != 2) continue;
                const int resolution =
                    quartet_resolution(gene, in[0], in[1], out[0], out[1]);
                expected += resolution == 0 ? 1.0 : resolution > 0 ? -1.0 : 0.0;
              }
      }
      require_close(scores[e], expected,
                    "divide score disagrees with enumeration (n=" +
                        std::to_string(taxon_count) + ")");
    }
  }
}

// The contextual-regrafting objective of attaching `context` at every edge of
// `side`, by direct enumeration: at each vertex u with arms R0, R1, R2, q_r(u)
// counts gene quartets x|y|z|b (one taxon per arm, b in the context) that put
// b with the arm-r taxon; the edge score is sum_u c_u (3 q_r(u) - 2 sum q(u))
// with r the arm holding the edge, and c_u = 1 or 1/(|R0||R1||R2||context|).
std::map<std::pair<int, int>, double> enumerated_regraft_scores(
    const qraft::internal::SolverTree& side,
    const std::vector<int>& context,
    const std::vector<PhysicalLeaves>& genes,
    int taxon_count,
    bool normalized) {
  const int nodes = static_cast<int>(side.nodes.size());
  const auto arm_taxa = [&](int from, int into) {
    std::vector<int> taxa;
    std::vector<int> stack{into};
    std::vector<std::uint8_t> seen(static_cast<std::size_t>(nodes), 0);
    seen[static_cast<std::size_t>(from)] = 1;
    while (!stack.empty()) {
      const int node = stack.back();
      stack.pop_back();
      if (seen[static_cast<std::size_t>(node)]) continue;
      seen[static_cast<std::size_t>(node)] = 1;
      const auto id = side.nodes[static_cast<std::size_t>(node)].taxon;
      if (id >= 0) taxa.push_back(static_cast<int>(id));
      for (int next : side.nodes[static_cast<std::size_t>(node)].neighbors) {
        stack.push_back(next);
      }
    }
    return taxa;
  };
  std::vector<std::vector<int>> gene_leaf;
  for (const PhysicalLeaves& gene : genes) {
    gene_leaf.push_back(leaf_of_taxon(gene, taxon_count));
  }
  struct ScoredVertex {
    int node;
    std::array<int, 3> arm_root;
    std::array<double, 3> q;
    double weight;
  };
  std::vector<ScoredVertex> vertices;
  for (int node = 0; node < nodes; ++node) {
    const auto& item = side.nodes[static_cast<std::size_t>(node)];
    if (item.taxon >= 0 || item.neighbors.size() != 3) continue;
    ScoredVertex vertex{node, {item.neighbors[0], item.neighbors[1], item.neighbors[2]},
                  {0.0, 0.0, 0.0}, 1.0};
    std::array<std::vector<int>, 3> arms;
    for (int r = 0; r < 3; ++r) {
      arms[static_cast<std::size_t>(r)] =
          arm_taxa(node, vertex.arm_root[static_cast<std::size_t>(r)]);
    }
    for (std::size_t g = 0; g < genes.size(); ++g) {
      const auto& leaf = gene_leaf[g];
      for (int x : arms[0])
        for (int y : arms[1])
          for (int z : arms[2])
            for (int b : context) {
              const int resolution = quartet_resolution(
                  genes[g], leaf[static_cast<std::size_t>(b)],
                  leaf[static_cast<std::size_t>(x)], leaf[static_cast<std::size_t>(y)],
                  leaf[static_cast<std::size_t>(z)]);
              if (resolution >= 0) vertex.q[static_cast<std::size_t>(resolution)] += 1.0;
            }
    }
    if (normalized) {
      vertex.weight = 1.0 / (static_cast<double>(context.size()) *
                             static_cast<double>(arms[0].size()) *
                             static_cast<double>(arms[1].size()) *
                             static_cast<double>(arms[2].size()));
    }
    vertices.push_back(vertex);
  }
  std::map<std::pair<int, int>, double> scores;
  for (int u = 0; u < nodes; ++u) {
    for (int v : side.nodes[static_cast<std::size_t>(u)].neighbors) {
      if (u > v) continue;
      double score = 0.0;
      for (const ScoredVertex& vertex : vertices) {
        // The arm of the vertex that holds edge (u, v).
        int arm = -1;
        for (int r = 0; r < 3 && arm < 0; ++r) {
          const int root = vertex.arm_root[static_cast<std::size_t>(r)];
          if (u == vertex.node) {
            if (v == root) arm = r;
            continue;
          }
          if (v == vertex.node) {
            if (u == root) arm = r;
            continue;
          }
          // Both endpoints lie in the same arm; test u.
          std::vector<int> stack{root};
          std::vector<std::uint8_t> seen(static_cast<std::size_t>(nodes), 0);
          seen[static_cast<std::size_t>(vertex.node)] = 1;
          while (!stack.empty() && arm < 0) {
            const int node = stack.back();
            stack.pop_back();
            if (seen[static_cast<std::size_t>(node)]) continue;
            seen[static_cast<std::size_t>(node)] = 1;
            if (node == u) arm = r;
            for (int next : side.nodes[static_cast<std::size_t>(node)].neighbors) {
              stack.push_back(next);
            }
          }
        }
        const double total = vertex.q[0] + vertex.q[1] + vertex.q[2];
        score += vertex.weight *
                 (3.0 * vertex.q[static_cast<std::size_t>(arm)] - 2.0 * total);
      }
      scores[{u, v}] = score;
    }
  }
  return scores;
}

void test_regraft_scores_match_enumeration() {
  using namespace qraft::internal;
  // Each side attaches at an edge maximizing the enumerated objective, and the
  // reported score equals the enumerated score of that edge. Genes may miss
  // taxa, and from the fourth trial on two taxa lie outside both sides.
  std::mt19937 rng(20260927u);
  for (qraft::RegraftScore mode :
       {qraft::RegraftScore::Raw, qraft::RegraftScore::Normalized}) {
    for (int trial = 0; trial < 6; ++trial) {
      const int taxon_count = 8 + trial % 3 + (trial >= 3 ? 2 : 0);
      const std::vector<std::string> labels = numbered_labels(taxon_count);
      SolverContext ctx;
      ctx.options = test_options(2);
      ctx.options.regraft_score = mode;
      ctx.workers = configure_threads(ctx.options);
      ctx.taxa = qraft::TaxonRegistry(labels);
      std::vector<PhysicalLeaves> genes;
      for (int g = 0; g < 16; ++g) {
        const qraft::Tree tree = qraft::parse_newick(
            g < 8 ? random_newick_over(rng, labels, false)
                  : random_incomplete_gene(rng, labels, g));
        ctx.genes.push_back(index_gene(tree, ctx.taxa));
        genes.push_back(physical_leaves(tree, ctx.taxa));
      }
      std::vector<int> order(static_cast<std::size_t>(taxon_count));
      std::iota(order.begin(), order.end(), 0);
      std::shuffle(order.begin(), order.end(), rng);
      const int left_size = 3 + trial % 3;
      const int live = taxon_count - (trial >= 3 ? 2 : 0);
      std::vector<int> left_taxa(order.begin(), order.begin() + left_size);
      std::vector<int> right_taxa(order.begin() + left_size,
                                  order.begin() + live);
      std::sort(left_taxa.begin(), left_taxa.end());
      std::sort(right_taxa.begin(), right_taxa.end());
      const auto side_tree = [&](const std::vector<int>& taxa) {
        std::vector<std::string> side_labels;
        for (int taxon : taxa) side_labels.push_back(ctx.taxa.label_for(taxon));
        return solver_tree_from(
            qraft::parse_newick(random_newick_over(rng, side_labels, false)), ctx.taxa);
      };
      const SolverTree left = side_tree(left_taxa);
      const SolverTree right = side_tree(right_taxa);
      const auto [left_edge, right_edge] = choose_regraft_edges(
          ctx, left, right, left_taxa, right_taxa);
      const bool normalized = mode == qraft::RegraftScore::Normalized;
      const auto check = [&](const SolverTree& side, const std::vector<int>& context,
                             const RegraftEdge& chosen) {
        const auto scores = enumerated_regraft_scores(side, context, genes,
                                                      taxon_count, normalized);
        const auto found = scores.find({std::min(chosen.a, chosen.b),
                                        std::max(chosen.a, chosen.b)});
        require(found != scores.end(), "regraft edge is not an edge of its side");
        double best = -std::numeric_limits<double>::infinity();
        for (const auto& [edge, score] : scores) best = std::max(best, score);
        const double tolerance = normalized ? 1e-12 : 1e-9;
        require(found->second >= best - 1e-9,
                "regraft edge does not maximize the enumerated objective");
        require_close(chosen.rank_score, found->second,
                      "regraft score disagrees with enumeration", tolerance,
                      1e-12);
      };
      check(left, right_taxa, left_edge);
      check(right, left_taxa, right_edge);
    }
  }
}

double enumerated_quartet_score(const qraft::internal::SolverTree& tree,
                                const std::vector<PhysicalLeaves>& genes,
                                int taxon_count) {
  const PhysicalLeaves species = physical_leaves(tree);
  const std::vector<int> species_leaf = leaf_of_taxon(species, taxon_count);
  double score = 0.0;
  for (const PhysicalLeaves& gene : genes) {
    const std::vector<int> leaf = leaf_of_taxon(gene, taxon_count);
    for (int a = 0; a < taxon_count; ++a)
      for (int b = a + 1; b < taxon_count; ++b)
        for (int c = b + 1; c < taxon_count; ++c)
          for (int d = c + 1; d < taxon_count; ++d) {
            const auto at = [](const std::vector<int>& map, int taxon) {
              return map[static_cast<std::size_t>(taxon)];
            };
            const int in_species =
                quartet_resolution(species, at(species_leaf, a), at(species_leaf, b),
                                   at(species_leaf, c), at(species_leaf, d));
            const int in_gene = quartet_resolution(gene, at(leaf, a), at(leaf, b),
                                                   at(leaf, c), at(leaf, d));
            score += in_species >= 0 && in_species == in_gene ? 1.0 : 0.0;
          }
  }
  return score;
}

namespace skeleton_test {

// Distances between all active nodes of a small tree.
std::vector<std::vector<int>> node_distances(const qraft::internal::SolverTree& tree) {
  const int n = static_cast<int>(tree.nodes.size());
  std::vector<std::vector<int>> dist(static_cast<std::size_t>(n),
                                     std::vector<int>(static_cast<std::size_t>(n), -1));
  for (int s = 0; s < n; ++s) {
    if (!tree.nodes[static_cast<std::size_t>(s)].active) continue;
    std::queue<int> pending;
    dist[static_cast<std::size_t>(s)][static_cast<std::size_t>(s)] = 0;
    pending.push(s);
    while (!pending.empty()) {
      const int u = pending.front();
      pending.pop();
      for (int v : tree.nodes[static_cast<std::size_t>(u)].neighbors) {
        if (!tree.nodes[static_cast<std::size_t>(v)].active) continue;
        if (dist[static_cast<std::size_t>(s)][static_cast<std::size_t>(v)] >= 0) continue;
        dist[static_cast<std::size_t>(s)][static_cast<std::size_t>(v)] =
            dist[static_cast<std::size_t>(s)][static_cast<std::size_t>(u)] + 1;
        pending.push(v);
      }
    }
  }
  return dist;
}

// The side of edge (a, b) containing a, as taxon flags.
std::vector<char> side_taxa(const qraft::internal::SolverTree& tree, int a, int b,
                            int taxon_count) {
  std::vector<char> in(static_cast<std::size_t>(taxon_count), 0);
  std::vector<int> stack{a};
  std::vector<char> seen(tree.nodes.size(), 0);
  seen[static_cast<std::size_t>(a)] = seen[static_cast<std::size_t>(b)] = 1;
  while (!stack.empty()) {
    const int u = stack.back();
    stack.pop_back();
    const int taxon = tree.nodes[static_cast<std::size_t>(u)].taxon;
    if (taxon >= 0) in[static_cast<std::size_t>(taxon)] = 1;
    for (int v : tree.nodes[static_cast<std::size_t>(u)].neighbors) {
      if (!tree.nodes[static_cast<std::size_t>(v)].active || seen[static_cast<std::size_t>(v)]) {
        continue;
      }
      seen[static_cast<std::size_t>(v)] = 1;
      stack.push_back(v);
    }
  }
  return in;
}

// The tree with the side of (a, b) containing a cut off and re-attached on
// the edge (u, v) of the other side; b is suppressed when it becomes a
// degree-two vertex (then (u, v) may be the merged edge).
qraft::internal::SolverTree reattach(const qraft::internal::SolverTree& tree, int a, int b,
                                     int u, int v) {
  using namespace qraft::internal;
  SolverTree t = tree;
  remove_solver_edge(t, a, b);
  int du = u, dv = v;
  if (t.nodes[static_cast<std::size_t>(b)].taxon < 0) {
    const auto rest = active_neighbors(t, b);
    if (rest.size() == 2) {
      remove_solver_edge(t, b, rest[0]);
      remove_solver_edge(t, b, rest[1]);
      add_solver_edge(t, rest[0], rest[1]);
      t.nodes[static_cast<std::size_t>(b)].active = false;
      if (du == b || dv == b) {
        du = rest[0];
        dv = rest[1];
      }
    }
  }
  const int m = static_cast<int>(t.nodes.size());
  t.nodes.push_back(SolverTree::Node{});
  remove_solver_edge(t, du, dv);
  add_solver_edge(t, du, m);
  add_solver_edge(t, m, dv);
  add_solver_edge(t, m, a);
  return t;
}

}  // namespace skeleton_test

// The skeleton pass against enumeration on random small inputs, genes that
// miss taxa included: the core is the finest that fits the read budget;
// every channel equals the enumerated count of its quartets; every
// destination's gain equals the enumerated change of the move over the
// quartets whose outside triple is centred at a core vertex, and the whole
// change when the moved side's edge is heavy or every vertex is read; each
// side's best destination is the largest; and the channels do not depend on
// the worker count or on when Extract/Contract replaces in-place reuse.
void test_skeleton_gains_exactness() {
  using namespace qraft::internal;
  using namespace skeleton_test;
  std::mt19937 rng(20261020u);
  int channels_checked = 0, gains_checked = 0, restricted_checked = 0, heavy_gains = 0,
      tight_budgets = 0;
  for (int trial = 0; trial < 60; ++trial) {
    const int n = 8 + trial % 6;
    const auto labels = numbered_labels(n);
    SolverContext ctx;
    ctx.options = test_options(1 + trial % 2 * 2);
    ctx.workers = configure_threads(ctx.options);
    ctx.taxa = qraft::TaxonRegistry(labels);
    std::vector<PhysicalLeaves> genes;
    for (int g = 0; g < 4 + trial % 7; ++g) {
      const auto gene = qraft::parse_newick(g % 2
          ? random_incomplete_gene(rng, labels, g)
          : random_newick_over(rng, labels, false));
      ctx.genes.push_back(index_gene(gene, ctx.taxa));
      genes.push_back(physical_leaves(gene, ctx.taxa));
    }
    const SolverTree tree = solver_tree_from(
        qraft::parse_newick(random_newick_over(rng, labels, false)), ctx.taxa);
    const std::vector<int> taxa = iota_taxa(n);
    // Every fourth trial reads every vertex; the others a budget of 1-6.
    const int budget = trial % 4 == 0 ? 64 : 1 + trial % 6;
    const int min_clade = 1 + trial % 2;
    const SkeletonPass pass = skeleton_pass(ctx, tree, taxa, budget, min_clade);
    const SolverTree& t = pass.tree;
    const int resolution = pass.resolution;
    const int core_count = static_cast<int>(pass.core.size());
    const int edge_count = static_cast<int>(pass.edge_nodes.size());
    require(core_count > 0 && edge_count > 0, "skeleton pass found no core or no edge");
    // The core: the vertices whose largest smaller side reaches the
    // resolution, which is the finest fitting the budget.
    int finer = 0, max_h = 0;
    for (int u = 0; u < static_cast<int>(t.nodes.size()); ++u) {
      const auto& node = t.nodes[static_cast<std::size_t>(u)];
      if (!node.active || node.taxon >= 0) continue;
      int h = 0;
      for (int v : active_neighbors(t, u)) {
        const std::vector<char> side = side_taxa(t, v, u, n);
        int count = 0;
        for (char c : side) count += c;
        h = std::max(h, std::min(count, n - count));
      }
      require((pass.core_index_of[static_cast<std::size_t>(u)] >= 0) == (h >= resolution),
              "core vertices are not the heavy-edge endpoints");
      finer += h >= resolution - 1 ? 1 : 0;
      max_h = std::max(max_h, h);
    }
    require(resolution >= 2, "core resolution below 2");
    require(core_count <= budget || resolution == std::max(2, max_h),
            "core exceeds its read budget without a tie at the top");
    require(resolution == 2 || finer > budget, "core is not the finest fitting its budget");
    tight_budgets += resolution > 2 ? 1 : 0;
    // 1. Channels against enumeration.
    for (int i = 0; i < core_count; ++i) {
      const int u = pass.core[static_cast<std::size_t>(i)];
      std::array<std::vector<char>, 3> arm;
      for (int c = 0; c < 3; ++c) {
        arm[static_cast<std::size_t>(c)] = side_taxa(
            t, pass.arm_neighbor[static_cast<std::size_t>(i)][static_cast<std::size_t>(c)], u, n);
      }
      for (int e = 0; e < edge_count; ++e) {
        const int a = pass.edge_child_source[static_cast<std::size_t>(e)];
        const int b = pass.edge_parent_source[static_cast<std::size_t>(e)];
        const std::size_t cell = static_cast<std::size_t>(i) * edge_count + e;
        const std::vector<char> away =
            pass.away_is_child[cell] != 0 ? side_taxa(t, a, b, n) : side_taxa(t, b, a, n);
        int arm_count[3] = {0, 0, 0};
        for (int x = 0; x < n; ++x) {
          if (!away[static_cast<std::size_t>(x)]) continue;
          for (int c = 0; c < 3; ++c) arm_count[c] += arm[static_cast<std::size_t>(c)][static_cast<std::size_t>(x)];
        }
        for (int c = 0; c < 3; ++c) {
          require((arm_count[c] > 0) == (c == pass.arm_of[cell]),
                  "skeleton pass: the away side is not in the recorded arm");
        }
        long long expect[3] = {0, 0, 0};
        for (const PhysicalLeaves& gene : genes) {
          const std::vector<int> leaf = leaf_of_taxon(gene, n);
          const auto in_arm = [&](int c, int y) {
            return !away[static_cast<std::size_t>(y)] &&
                   arm[static_cast<std::size_t>(c)][static_cast<std::size_t>(y)] &&
                   leaf[static_cast<std::size_t>(y)] >= 0;
          };
          for (int x = 0; x < n; ++x) {
            if (!away[static_cast<std::size_t>(x)] || leaf[static_cast<std::size_t>(x)] < 0) continue;
            for (int y0 = 0; y0 < n; ++y0) {
              if (!in_arm(0, y0)) continue;
              for (int y1 = 0; y1 < n; ++y1) {
                if (!in_arm(1, y1)) continue;
                for (int y2 = 0; y2 < n; ++y2) {
                  if (!in_arm(2, y2)) continue;
                  // 0: x y0 | y1 y2, 1: x y1 | y0 y2, 2: x y2 | y0 y1.
                  const int r = quartet_resolution(
                      gene, leaf[static_cast<std::size_t>(x)], leaf[static_cast<std::size_t>(y0)],
                      leaf[static_cast<std::size_t>(y1)], leaf[static_cast<std::size_t>(y2)]);
                  if (r >= 0) ++expect[r];
                }
              }
            }
          }
        }
        for (int c = 0; c < 3; ++c) {
          require(pass.channels[cell * 3 + static_cast<std::size_t>(c)] ==
                      static_cast<SkeletonCount>(expect[c]),
                  "skeleton channel differs from enumeration");
        }
        ++channels_checked;
      }
    }
    // 2. Gains against enumeration and the restricted identity.
    const PhysicalLeaves species = physical_leaves(t);
    const std::vector<int> species_leaf = leaf_of_taxon(species, n);
    const auto dist = node_distances(t);
    std::vector<int> leaf_node(static_cast<std::size_t>(n), -1);
    for (int u = 0; u < static_cast<int>(t.nodes.size()); ++u) {
      if (t.nodes[static_cast<std::size_t>(u)].active && t.nodes[static_cast<std::size_t>(u)].taxon >= 0) {
        leaf_node[static_cast<std::size_t>(t.nodes[static_cast<std::size_t>(u)].taxon)] = u;
      }
    }
    const double before = enumerated_quartet_score(t, genes, n);
    for (int e = 0; e < edge_count; ++e) {
      std::vector<SkeletonDestination> child_profile, parent_profile;
      const SkeletonEdgeGains g = skeleton_edge_gains(pass, e, &child_profile, &parent_profile);
      const int a = pass.edge_child_source[static_cast<std::size_t>(e)];
      const int b = pass.edge_parent_source[static_cast<std::size_t>(e)];
      for (int side = 0; side < 2; ++side) {
        const bool has = side == 0 ? g.child_has_data : g.parent_has_data;
        const auto& profile = side == 0 ? child_profile : parent_profile;
        const auto& best = side == 0 ? g.child_best : g.parent_best;
        if (!has) {
          require(profile.empty() && best.gain == 0, "a side without data has a profile");
          continue;
        }
        const int moved = side == 0 ? a : b;
        const int other = side == 0 ? b : a;
        const std::vector<char> moved_taxa = side_taxa(t, moved, other, n);
        SkeletonCount best_seen = 0;
        for (const SkeletonDestination& d : profile) {
          const SolverTree after_tree = reattach(t, moved, other, d.vertex, d.neighbor);
          const double after = enumerated_quartet_score(after_tree, genes, n);
          const PhysicalLeaves species_after = physical_leaves(after_tree);
          const std::vector<int> species_after_leaf = leaf_of_taxon(species_after, n);
          // The change over the quartets with one leaf in the moved side
          // whose outside triple is centred at a core vertex.
          long long restricted = 0;
          bool all_core_centred = true;
          for (const PhysicalLeaves& gene : genes) {
            const std::vector<int> leaf = leaf_of_taxon(gene, n);
            const auto at = [](const std::vector<int>& map, int taxon) {
              return map[static_cast<std::size_t>(taxon)];
            };
            for (int x = 0; x < n; ++x) {
              if (!moved_taxa[static_cast<std::size_t>(x)] || at(leaf, x) < 0) continue;
              for (int y0 = 0; y0 < n; ++y0)
                for (int y1 = y0 + 1; y1 < n; ++y1)
                  for (int y2 = y1 + 1; y2 < n; ++y2) {
                    if (moved_taxa[static_cast<std::size_t>(y0)] || moved_taxa[static_cast<std::size_t>(y1)] ||
                        moved_taxa[static_cast<std::size_t>(y2)]) {
                      continue;
                    }
                    if (at(leaf, y0) < 0 || at(leaf, y1) < 0 || at(leaf, y2) < 0) continue;
                    const int in_gene = quartet_resolution(gene, at(leaf, x), at(leaf, y0),
                                                           at(leaf, y1), at(leaf, y2));
                    const int was = quartet_resolution(
                        species, at(species_leaf, x), at(species_leaf, y0),
                        at(species_leaf, y1), at(species_leaf, y2));
                    const int now = quartet_resolution(
                        species_after, at(species_after_leaf, x), at(species_after_leaf, y0),
                        at(species_after_leaf, y1), at(species_after_leaf, y2));
                    const int change = (now >= 0 && now == in_gene ? 1 : 0) -
                                       (was >= 0 && was == in_gene ? 1 : 0);
                    if (change == 0) continue;
                    // The median of the outside triple in the original tree.
                    const int l0 = leaf_node[static_cast<std::size_t>(y0)];
                    const int l1 = leaf_node[static_cast<std::size_t>(y1)];
                    const int l2 = leaf_node[static_cast<std::size_t>(y2)];
                    const auto d2 = [&](int p, int q) {
                      return dist[static_cast<std::size_t>(p)][static_cast<std::size_t>(q)];
                    };
                    const int target = (d2(l0, l1) + d2(l0, l2) + d2(l1, l2)) / 2;
                    int centre = -1;
                    for (int m = 0; m < static_cast<int>(t.nodes.size()) && centre < 0; ++m) {
                      if (!t.nodes[static_cast<std::size_t>(m)].active) continue;
                      if (d2(l0, m) + d2(l1, m) + d2(l2, m) == target) centre = m;
                    }
                    require(centre >= 0, "no centre for an outside triple");
                    if (pass.core_index_of[static_cast<std::size_t>(centre)] >= 0) {
                      restricted += change;
                    } else {
                      all_core_centred = false;
                    }
                  }
            }
          }
          require(static_cast<SkeletonCount>(restricted) == d.gain,
                  "skeleton gain differs from the enumerated change over core-centred quartets");
          ++restricted_checked;
          if (pass.edge_heavy[static_cast<std::size_t>(e)] != 0 || resolution == 2) {
            require_close(after - before, static_cast<double>(d.gain),
                          "skeleton gain differs from the enumerated change of the move",
                          1e-9, 1e-6);
            require(all_core_centred,
                    "a heavy edge's move changed a quartet the core does not count");
            ++gains_checked;
            heavy_gains += pass.edge_heavy[static_cast<std::size_t>(e)] != 0 ? 1 : 0;
          }
          best_seen = std::max<SkeletonCount>(best_seen, d.gain);
        }
        require(best.gain == best_seen, "best skeleton destination is not the largest");
      }
    }
    // 3. The channels at another worker count and with Extract/Contract
    // forced on small inputs (reuse bound factor 2).
    if (trial % 3 == 0) {
      ctx.options.threads = trial % 2 == 0 ? 3 : 1;
      ctx.workers = configure_threads(ctx.options);
      const SkeletonPass again = skeleton_pass(ctx, tree, taxa, budget, min_clade, 2);
      require(again.core == pass.core && again.channels == pass.channels &&
                  again.arm_of == pass.arm_of && again.away_is_child == pass.away_is_child,
              "skeleton pass depends on the worker count or on state reuse");
    }
  }
  require(channels_checked > 500 && restricted_checked > 1000 && heavy_gains > 100 &&
              tight_budgets > 10,
          "skeleton campaign exercised too little");
  std::cout << "skeleton gains: channels " << channels_checked << ", restricted identities "
            << restricted_checked << ", exact gains " << gains_checked << " (heavy "
            << heavy_gains << "), trials with a binding budget " << tight_budgets << "\n";
}

// The SPR refinement on random small inputs, NNI off: every applied move
// raises the enumerated quartet score by exactly the oracle's gain, within
// the candidate budget, and the tree and gains are the same at one and four
// workers. A planted misplaced clade is moved home with NNI on; a non-binary
// tree is refused; the stage is off by default.
void test_spr_gain_exactness() {
  using namespace qraft::internal;
  SolverContext ctx;
  ctx.options = test_options(1);
  ctx.options.nni = false;
  ctx.options.spr_refine = true;
  ctx.options.spr_min_clade = 2;
  ctx.options.spr_candidates = 1;
  ctx.options.spr_rounds = 1;
  ctx.workers = configure_threads(ctx.options);
  // The planted clade move from test_spr_exactness, with three taxa per tip.
  const std::string ab = "((A1,A2),A3),(B1,(B2,B3))";
  const std::string cd = "((C1,C2),C3),(D1,(D2,D3))";
  const std::string ef = "((E1,E2),E3),(F1,(F2,F3))";
  const std::string gh = "((G1,G2),G3),(H1,(H2,H3))";
  std::vector<std::string> labels;
  for (char c = 'A'; c <= 'H'; ++c) {
    for (int i = 1; i <= 3; ++i) labels.push_back(std::string(1, c) + std::to_string(i));
  }
  ctx.taxa = qraft::TaxonRegistry(labels);
  const auto good = qraft::parse_newick("(((" + ab + "),(" + cd + ")),((" +
                                        ef + "),(" + gh + ")));" );
  SolverTree tree = solver_tree_from(qraft::parse_newick(
      "((" + cd + "),(((" + ab + "),(" + ef + ")),(" + gh + ")));"), ctx.taxa);
  std::vector<PhysicalLeaves> genes;
  for (int g = 0; g < 3; ++g) {
    ctx.genes.push_back(index_gene(good, ctx.taxa));
    genes.push_back(physical_leaves(good, ctx.taxa));
  }
  reset_run_stats();
  set_run_stats_enabled(true);
  RegraftCount expected = 0;
  for (int call = 0; call < 5; ++call) {
    const auto before = static_cast<RegraftCount>(enumerated_quartet_score(tree, genes, 24));
    const bool moved = refine_spr(ctx, tree, iota_taxa(24));
    const auto after = static_cast<RegraftCount>(enumerated_quartet_score(tree, genes, 24));
    expected += after - before;
    require(run_stats().spr_gain_exact, "SPR exact gain unavailable");
    require(run_stats().spr_gain == expected, "SPR integer gain differs from enumeration");
    if (!moved) break;
  }
  require(run_stats().spr_moves.load() > 0, "SPR exact gain fixture applied no moves");
  set_run_stats_enabled(false);
  reset_run_stats();
  require(run_stats().spr_gain == 0 && run_stats().spr_gain_exact,
          "SPR exact gain statistics did not reset");

  // One more gene, with 8000 copies each of A1, C1, E1 and G1: a vertex
  // whose three arms and other side hold one of them each has a channel sum
  // near 8000^4 > 2^51, which exceeds the exactness bound, so the planted
  // move's gain is reported inexact.
  std::string heavy_text = "(((" + ab + "),(" + cd + ")),((" + ef + "),(" + gh + ")));";
  for (const std::string taxon : {"A1", "C1", "E1", "G1"}) {
    std::string copies = "(";
    for (int i = 0; i < 8000; ++i) copies += (i == 0 ? "" : ",") + taxon;
    heavy_text.replace(heavy_text.find(taxon), taxon.size(), copies + ")");
  }
  ctx.genes.push_back(index_gene(qraft::parse_newick(heavy_text), ctx.taxa));
  SolverTree heavy_tree = solver_tree_from(qraft::parse_newick(
      "((" + cd + "),(((" + ab + "),(" + ef + ")),(" + gh + ")));"), ctx.taxa);
  set_run_stats_enabled(true);
  (void)refine_spr(ctx, heavy_tree, iota_taxa(24));
  require(run_stats().spr_moves.load() > 0 && !run_stats().spr_gain_exact,
          "SPR gain past the exact range was not reported inexact");
  set_run_stats_enabled(false);
  reset_run_stats();
}

void test_spr_exactness() {
  using namespace qraft::internal;
  std::mt19937 rng(20261019u);
  int accepted = 0;
  int single_moves = 0;
  for (int trial = 0; trial < 210; ++trial) {
    const int n = 8 + trial % 7;
    const auto labels = numbered_labels(n);
    SolverContext ctx;
    ctx.options = test_options(1);
    ctx.options.nni = false;  // isolate the re-joins
    ctx.options.spr_refine = true;
    ctx.options.spr_rounds = 1;  // one round per call: inspect each separately
    ctx.options.spr_min_clade = 2;
    ctx.options.spr_reads = trial % 3 == 0 ? 2 + trial % 4 : 128;
    ctx.options.spr_candidates = trial % 5 == 0 ? 1 : 150;
    ctx.options.regraft_score = trial % 2 ? qraft::RegraftScore::Normalized
                                          : qraft::RegraftScore::Raw;
    ctx.taxa = qraft::TaxonRegistry(labels);
    std::vector<PhysicalLeaves> genes;
    for (int g = 0; g < 5 + trial % 16; ++g) {
      const auto gene = qraft::parse_newick(g % 2
          ? random_incomplete_gene(rng, labels, g)
          : random_newick_over(rng, labels, false));
      ctx.genes.push_back(index_gene(gene, ctx.taxa));
      genes.push_back(physical_leaves(gene, ctx.taxa));
    }
    SolverTree initial = solver_tree_from(
        qraft::parse_newick(random_newick_over(rng, labels, false)), ctx.taxa);
    // Degree-two Newick roots, including roots on pendant edges.
    if (trial % 3 == 0) {
      const int a = static_cast<int>(rng() % initial.nodes.size());
      const int b = initial.nodes[static_cast<std::size_t>(a)].neighbors.front();
      const int root = static_cast<int>(initial.nodes.size());
      initial.nodes.push_back(SolverTree::Node{});
      remove_solver_edge(initial, a, b);
      add_solver_edge(initial, a, root);
      add_solver_edge(initial, root, b);
    }
    const std::vector<int> taxa = iota_taxa(n);
    SolverTree reference;
    RegraftCount reference_gain = 0;
    for (int workers : {1, 4}) {
      ctx.options.threads = workers;
      ctx.workers = configure_threads(ctx.options);
      reset_run_stats();
      set_run_stats_enabled(true);
      SolverTree tree = initial;
      double before = enumerated_quartet_score(tree, genes, n);
      // With one candidate a call applies at most one move, so each applied
      // move is checked against enumeration on its own.
      const int calls = ctx.options.spr_candidates == 1 ? 12 : 3;
      for (int call = 0; call < calls; ++call) {
        const RegraftCount previous_gain = run_stats().spr_gain;
        const long long previous_scored = run_stats().spr_candidates_scored.load();
        const bool moved = refine_spr(ctx, tree, taxa);
        const double after = enumerated_quartet_score(tree, genes, n);
        require(after >= before, "SPR decreased the exact quartet score");
        require(run_stats().spr_gain_exact &&
                    static_cast<RegraftCount>(after - before) ==
                        run_stats().spr_gain - previous_gain,
                "SPR gain disagrees with enumeration");
        require(run_stats().spr_candidates_scored.load() - previous_scored <=
                    ctx.options.spr_candidates,
                "SPR exceeded its candidate budget");
        int active_nodes = 0;
        for (const auto& node : tree.nodes) active_nodes += node.active ? 1 : 0;
        require(active_nodes == 2 * n - 2,
                "SPR re-join left a degree-two vertex or lost one");
        if (!moved) break;
        require(after > before, "SPR applied a move that does not improve");
        if (workers == 1) {
          ++accepted;
          if (ctx.options.spr_candidates == 1) ++single_moves;
        }
        before = after;
      }
      if (workers == 1) {
        reference = tree;
        reference_gain = run_stats().spr_gain;
      } else {
        require(reference.nodes.size() == tree.nodes.size(), "SPR worker size mismatch");
        for (std::size_t i = 0; i < tree.nodes.size(); ++i) {
          require(reference.nodes[i].neighbors == tree.nodes[i].neighbors &&
                      reference.nodes[i].taxon == tree.nodes[i].taxon,
                  "SPR differs between one and four workers");
        }
        require(reference_gain == run_stats().spr_gain, "SPR worker gain mismatch");
      }
      set_run_stats_enabled(false);
    }
  }
  require(accepted > 100, "SPR random campaign applied too few moves");
  require(single_moves > 30, "SPR random campaign checked too few single moves");

  // Every round ranks clades at the same floor: three rounds in one call
  // equal three one-round calls. The trials where a later round at floor 10
  // would move differently tell this apart from a floor lowered after the
  // first round.
  const auto splits = [](const SolverTree& tree, int n) {
    std::set<std::vector<char>> out;
    for (int a = 0; a < static_cast<int>(tree.nodes.size()); ++a) {
      if (!tree.nodes[static_cast<std::size_t>(a)].active) continue;
      for (int b : tree.nodes[static_cast<std::size_t>(a)].neighbors) {
        if (!tree.nodes[static_cast<std::size_t>(b)].active) continue;
        auto side = skeleton_test::side_taxa(tree, a, b, n);
        if (side[0] != 0) {
          for (char& in : side) in = in != 0 ? 0 : 1;
        }
        out.insert(std::move(side));
      }
    }
    return out;
  };
  int distinguishing = 0;
  for (int trial = 0; trial < 40; ++trial) {
    const int n = 24 + trial % 9;
    const auto labels = numbered_labels(n);
    SolverContext ctx;
    ctx.options = test_options(1);
    ctx.options.nni = false;
    ctx.options.spr_refine = true;
    ctx.options.spr_min_clade = 11;
    ctx.taxa = qraft::TaxonRegistry(labels);
    for (int g = 0; g < 6 + trial % 8; ++g) {
      ctx.genes.push_back(index_gene(
          qraft::parse_newick(random_newick_over(rng, labels, false)), ctx.taxa));
    }
    ctx.workers = configure_threads(ctx.options);
    const SolverTree initial = solver_tree_from(
        qraft::parse_newick(random_newick_over(rng, labels, false)), ctx.taxa);
    const std::vector<int> taxa = iota_taxa(n);
    SolverTree whole = initial;
    ctx.options.spr_rounds = 3;
    refine_spr(ctx, whole, taxa);
    SolverTree stepwise = initial;
    ctx.options.spr_rounds = 1;
    for (int round = 0; round < 3; ++round) {
      if (!refine_spr(ctx, stepwise, taxa)) break;
      if (round + 1 < 3) {
        SolverTree at_floor = stepwise, at_ten = stepwise;
        refine_spr(ctx, at_floor, taxa);
        ctx.options.spr_min_clade = 10;
        refine_spr(ctx, at_ten, taxa);
        ctx.options.spr_min_clade = 11;
        if (splits(at_floor, n) != splits(at_ten, n)) ++distinguishing;
      }
    }
    require(splits(whole, n) == splits(stepwise, n),
            "SPR rounds of one call differ from one-round calls at the same floor");
  }
  require(distinguishing > 0, "SPR floor campaign has no distinguishing trial");

  // Move the intact AB clade from beside EF back beside CD, NNI on.
  SolverContext ctx;
  ctx.options = test_options(1);
  ctx.options.spr_refine = true;
  ctx.options.spr_min_clade = 2;
  ctx.workers = configure_threads(ctx.options);
  const auto good = qraft::parse_newick("(((A,B),(C,D)),((E,F),(G,H)));");
  ctx.taxa = qraft::TaxonRegistry({"A", "B", "C", "D", "E", "F", "G", "H"});
  std::vector<PhysicalLeaves> genes;
  for (int g = 0; g < 10; ++g) {
    ctx.genes.push_back(index_gene(good, ctx.taxa));
    genes.push_back(physical_leaves(good, ctx.taxa));
  }
  SolverTree bad = solver_tree_from(
      qraft::parse_newick("((C,D),(((A,B),(E,F)),(G,H)));"), ctx.taxa);
  const double optimum = enumerated_quartet_score(solver_tree_from(good, ctx.taxa), genes, 8);
  require(enumerated_quartet_score(bad, genes, 8) < optimum, "planted SPR is not worse");
  ctx.options.spr_rounds = 5;
  reset_run_stats();
  set_run_stats_enabled(true);
  require(refine_spr(ctx, bad, iota_taxa(8)), "SPR missed the planted improvement");
  require_close(enumerated_quartet_score(bad, genes, 8), optimum,
                "SPR did not restore the planted optimum");
  // The round after the move finds none and ends the refinement.
  require(run_stats().spr_moves.load() == 1 && run_stats().spr_rounds.load() == 2,
          "SPR did not stop after a round without a move");
  set_run_stats_enabled(false);
  // A star has no edge between two three-arm vertices: nothing to move.
  SolverTree star = base_tree(iota_taxa(8));
  const SolverTree star_before = star;
  require(!refine_spr(ctx, star, iota_taxa(8)), "SPR moved a star");
  for (std::size_t i = 0; i < star.nodes.size(); ++i) {
    require(star.nodes[i].neighbors == star_before.nodes[i].neighbors,
            "SPR changed a star");
  }
  ctx.options = test_options(1);
  require(!ctx.options.spr_refine, "SPR refinement is on by default");
  require(!refine_spr(ctx, star, iota_taxa(8)), "disabled SPR changed a tree");
  std::cout << "SPR exactness: 210 random inputs, 1/4 workers, " << accepted
            << " applied moves checked (" << single_moves
            << " one per call), " << distinguishing
            << " floor-distinguishing rounds; planted optimum restored\n";
}

// The SPR refinement on species trees with polytomies (a starting tree kept
// under --no-divide), NNI off: every applied move is exact against
// enumeration, the polytomies and the number of resolved edges stay, and the
// result does not depend on the worker count.
void test_spr_polytomy_exactness() {
  using namespace qraft::internal;
  std::mt19937 rng(20261025u);
  const auto degrees = [](const SolverTree& tree) {
    std::vector<int> out;
    int nodes = 0;
    for (int u = 0; u < static_cast<int>(tree.nodes.size()); ++u) {
      if (!tree.nodes[static_cast<std::size_t>(u)].active) continue;
      ++nodes;
      if (tree.nodes[static_cast<std::size_t>(u)].taxon >= 0) continue;
      out.push_back(static_cast<int>(active_neighbors(tree, u).size()));
    }
    std::sort(out.begin(), out.end());
    out.push_back(nodes);
    return out;
  };
  int accepted = 0, with_polytomy = 0, polytomy_single_moves = 0;
  for (int trial = 0; trial < 160; ++trial) {
    const int n = 9 + trial % 8;
    const auto labels = numbered_labels(n);
    SolverContext ctx;
    ctx.options = test_options(1);
    ctx.options.nni = false;
    ctx.options.spr_refine = true;
    ctx.options.spr_rounds = 1;
    ctx.options.spr_min_clade = 2;
    ctx.options.spr_reads = trial % 3 == 0 ? 2 + trial % 4 : 128;
    ctx.options.spr_candidates = trial % 5 == 0 ? 1 : 150;
    ctx.taxa = qraft::TaxonRegistry(labels);
    std::vector<PhysicalLeaves> genes;
    for (int g = 0; g < 5 + trial % 12; ++g) {
      const auto gene = qraft::parse_newick(g % 2
          ? random_incomplete_gene(rng, labels, g)
          : random_newick_over(rng, labels, false));
      ctx.genes.push_back(index_gene(gene, ctx.taxa));
      genes.push_back(physical_leaves(gene, ctx.taxa));
    }
    const SolverTree initial = solver_tree_from(
        qraft::parse_newick(random_newick_over(rng, labels, true)), ctx.taxa);
    const auto initial_degrees = degrees(initial);
    const bool polytomous = initial_degrees[initial_degrees.size() - 2] > 3;
    if (polytomous) ++with_polytomy;
    const std::vector<int> taxa = iota_taxa(n);
    SolverTree reference;
    for (int workers : {1, 4}) {
      ctx.options.threads = workers;
      ctx.workers = configure_threads(ctx.options);
      reset_run_stats();
      set_run_stats_enabled(true);
      SolverTree tree = initial;
      double before = enumerated_quartet_score(tree, genes, n);
      const int calls = ctx.options.spr_candidates == 1 ? 10 : 3;
      for (int call = 0; call < calls; ++call) {
        const RegraftCount previous_gain = run_stats().spr_gain;
        const bool moved = refine_spr(ctx, tree, taxa);
        const double after = enumerated_quartet_score(tree, genes, n);
        require(after >= before, "SPR with polytomies decreased the quartet score");
        require(run_stats().spr_gain_exact &&
                    static_cast<RegraftCount>(after - before) ==
                        run_stats().spr_gain - previous_gain,
                "SPR gain with polytomies disagrees with enumeration");
        require(degrees(tree) == initial_degrees,
                "SPR changed the polytomies or the number of vertices");
        if (!moved) break;
        require(after > before, "SPR with polytomies applied a move that does not improve");
        if (workers == 1) {
          ++accepted;
          // With one candidate a call applies one move, checked on its own.
          if (polytomous && ctx.options.spr_candidates == 1) ++polytomy_single_moves;
        }
        before = after;
      }
      if (workers == 1) {
        reference = tree;
      } else {
        for (std::size_t i = 0; i < tree.nodes.size(); ++i) {
          require(reference.nodes[i].neighbors == tree.nodes[i].neighbors &&
                      reference.nodes[i].taxon == tree.nodes[i].taxon,
                  "SPR with polytomies differs between one and four workers");
        }
      }
      set_run_stats_enabled(false);
    }
  }
  require(with_polytomy > 100, "SPR polytomy campaign drew too few polytomies");
  require(accepted > 50, "SPR polytomy campaign applied too few moves");
  require(polytomy_single_moves > 10,
          "SPR polytomy campaign checked too few single moves beside polytomies");
}

// --guide-tree: the starting tree loses exactly its edges that conflict with
// the guide on the shared taxa, keeps every other edge, and gains only guide
// bipartitions. Three fixed cases pin the guided trees, children in order.
std::string guide_test_newick(std::mt19937& rng, std::vector<std::string> labels, double contract) {
  std::vector<std::string> parts = std::move(labels);
  std::vector<std::vector<std::string>> kids(parts.size());  // children of each part, empty for a leaf
  while (parts.size() > 1) {
    std::uniform_int_distribution<std::size_t> pick(0, parts.size() - 1);
    const std::size_t a = pick(rng);
    std::size_t b = pick(rng);
    while (b == a) b = pick(rng);
    std::vector<std::string> joined;
    for (std::size_t i : {a, b}) {
      // A child that is itself a join may be flattened into this node.
      if (!kids[i].empty() && std::uniform_real_distribution<double>(0, 1)(rng) < contract) {
        joined.insert(joined.end(), kids[i].begin(), kids[i].end());
      } else {
        joined.push_back(parts[i]);
      }
    }
    std::string text = "(";
    for (std::size_t i = 0; i < joined.size(); ++i) text += (i ? "," : "") + joined[i];
    text += ")";
    for (std::size_t i : {std::max(a, b), std::min(a, b)}) {
      parts.erase(parts.begin() + static_cast<std::ptrdiff_t>(i));
      kids.erase(kids.begin() + static_cast<std::ptrdiff_t>(i));
    }
    parts.push_back(text);
    kids.push_back(joined);
  }
  return parts[0] + ";";
}

// The nontrivial splits of `tree` restricted to `keep`, each as the sorted
// side without the smallest kept label.
std::set<std::vector<std::string>> guide_test_splits(const qraft::Tree& tree, const std::set<std::string>& keep) {
  std::set<std::vector<std::string>> out;
  const std::string anchor = *keep.begin();
  for (std::size_t u = 0; u < tree.nodes.size(); ++u) {
    for (int v : tree.nodes[u].neighbors) {
      std::vector<char> seen(tree.nodes.size(), 0);
      seen[u] = 1;
      std::vector<int> stack{v};
      std::set<std::string> side;
      while (!stack.empty()) {
        const int x = stack.back();
        stack.pop_back();
        if (seen[static_cast<std::size_t>(x)]) continue;
        seen[static_cast<std::size_t>(x)] = 1;
        const std::string& label = tree.nodes[static_cast<std::size_t>(x)].label;
        if (!label.empty() && keep.count(label)) side.insert(label);
        for (int y : tree.nodes[static_cast<std::size_t>(x)].neighbors) stack.push_back(y);
      }
      if (side.count(anchor)) {
        std::set<std::string> other;
        for (const auto& l : keep) {
          if (!side.count(l)) other.insert(l);
        }
        side = std::move(other);
      }
      if (side.size() >= 2 && keep.size() - side.size() >= 2) out.insert({side.begin(), side.end()});
    }
  }
  return out;
}

bool guide_test_compatible(const std::vector<std::string>& a, const std::vector<std::string>& b,
                           const std::set<std::string>& keep) {
  // Both sides exclude the smallest kept label, so the splits are compatible
  // iff the sides are nested or disjoint.
  const std::set<std::string> x(a.begin(), a.end()), y(b.begin(), b.end());
  bool inter = false, x_in_y = true, y_in_x = true;
  for (const auto& l : x) {
    if (y.count(l)) inter = true; else x_in_y = false;
  }
  for (const auto& l : y) {
    if (!x.count(l)) y_in_x = false;
  }
  (void)keep;
  return !inter || x_in_y || y_in_x;
}

std::string guided(const std::string& start, const std::string& guide, const qraft::TaxonRegistry& taxa) {
  using namespace qraft::internal;
  const StartingTree s = starting_tree_from(qraft::parse_newick(start), taxa);
  return starting_tree_newick(apply_guide_tree(s, qraft::parse_newick(guide), taxa), taxa);
}

void test_guide_tree() {
  {
    const qraft::TaxonRegistry taxa({"a", "b", "c", "d", "e", "f", "g", "h"});
    const qraft::TaxonRegistry six({"a", "b", "c", "d", "e", "f"});
    require(guided("(a,b,c,d,e,f);", "((a,b),(c,d),(e,f));", six) == "(a,b,((c,d),(e,f)));\n",
           "guide: a star is refined by the guide");
    const qraft::TaxonRegistry seven({"a", "b", "c", "d", "e", "f", "g"});
    require(guided("((a,b,c),(d,e),f,g);", "((a,d),(b,c),(e,(f,g)),x);", seven) ==
               "(a,d,((e,(f,g)),(b,c)));\n",
           "guide: conflicting edges go, a label outside the run is ignored");
    require(guided("(((a,b),c),(d,(e,f)),g,h);", "((a,c),(b,(d,e)),(f,g),h);", taxa) ==
               "(a,c,(h,(b,(d,e)),(f,g)));\n",
           "guide: contraction then refinement");
    // A guide that agrees with a binary start changes nothing, and a guide
    // without run taxa is ignored.
    require(guided("((a,b),(c,d),((e,f),(g,h)));", "((a,b),(c,d),((e,f),(g,h)));", taxa) ==
               "((a,b),(c,d),((e,f),(g,h)));\n",
           "guide: an agreeing guide changes nothing");
    require(guided("(a,b,c,d,e,f,g,h);", "((x,y),(z,w));", taxa) == "(a,b,c,d,e,f,g,h);\n",
           "guide: no shared taxa, no change");
  }
  std::mt19937 rng(2026102212);
  int changed = 0;
  for (int trial = 0; trial < 300; ++trial) {
    const int n = 6 + static_cast<int>(rng() % 15);
    std::vector<std::string> labels;
    for (int i = 0; i < n; ++i) labels.push_back((i < 10 ? "t0" : "t") + std::to_string(i));
    std::shuffle(labels.begin(), labels.end(), rng);
    const double cs = std::array<double, 3>{0.0, 0.3, 0.6}[rng() % 3];
    const std::string start = guide_test_newick(rng, labels, cs);
    const int k = 4 + static_cast<int>(rng() % static_cast<unsigned>(n - 3));
    std::vector<std::string> sub(labels.begin(), labels.begin() + k);
    for (int i = 0; i < static_cast<int>(rng() % 3); ++i) sub.push_back("x" + std::to_string(i));
    std::shuffle(sub.begin(), sub.end(), rng);
    const std::string guide = guide_test_newick(rng, sub, std::array<double, 3>{0.0, 0.2, 0.4}[rng() % 3]);
    const qraft::TaxonRegistry taxa(labels);
    const std::string out = guided(start, guide, taxa);
    if (out != start + "\n") ++changed;
    const qraft::Tree s = qraft::parse_newick(start), g = qraft::parse_newick(guide), r = qraft::parse_newick(out);
    require(r.leaf_labels().size() == labels.size(), "guide: every taxon kept once");
    const std::set<std::string> all(labels.begin(), labels.end());
    const std::set<std::string> shared(labels.begin(), labels.begin() + k);
    const auto s_all = guide_test_splits(s, all), r_all = guide_test_splits(r, all);
    const auto g_sh = guide_test_splits(g, shared);
    for (const auto& split : s_all) {
      // The split's side restricted to the shared taxa, as the side without
      // their smallest label.
      std::set<std::string> side;
      for (const auto& l : split) {
        if (shared.count(l)) side.insert(l);
      }
      if (side.count(*shared.begin())) {
        std::set<std::string> other;
        for (const auto& l : shared) {
          if (!side.count(l)) other.insert(l);
        }
        side = std::move(other);
      }
      bool conflicts = false;
      if (side.size() >= 2 && shared.size() - side.size() >= 2) {
        const std::vector<std::string> a(side.begin(), side.end());
        for (const auto& gs : g_sh) conflicts = conflicts || !guide_test_compatible(a, gs, shared);
      }
      require(r_all.count(split) == (conflicts ? 0u : 1u),
             "guide: a starting-tree split survives iff it does not conflict with the guide");
    }
    const auto r_sh = guide_test_splits(r, shared);
    for (const auto& split : r_all) {
      if (s_all.count(split)) continue;
      std::set<std::string> side;
      for (const auto& l : split) {
        if (shared.count(l)) side.insert(l);
      }
      if (side.count(*shared.begin())) {
        std::set<std::string> other;
        for (const auto& l : shared) {
          if (!side.count(l)) other.insert(l);
        }
        side = std::move(other);
      }
      if (side.size() >= 2 && shared.size() - side.size() >= 2) {
        require(g_sh.count({side.begin(), side.end()}) == 1, "guide: an added split is a guide split on the shared taxa");
      }
    }
    for (const auto& split : r_sh) {
      for (const auto& gs : g_sh) require(guide_test_compatible(split, gs, shared), "guide: the result agrees with the guide");
    }
  }
  require(changed > 200, "guide: most random cases change the starting tree");
  std::cout << "guide tree: 300 random cases, " << changed << " changed; 5 fixed cases match\n";
}

// --guide-tree seeding: the MRE built with the guide filter agrees with the guide on
// the shared taxa, so the guide contracts none of its edges, and after the
// refinement the starting tree displays every guide split on those taxa.
void test_guide_seed() {
  using namespace qraft::internal;
  qraft::Options options = test_options(1);
  require(configure_threads(options) == 1, "failed to configure one worker");
  std::mt19937 rng(2026102215);
  int filtered = 0;
  for (int trial = 0; trial < 200; ++trial) {
    const int n = 6 + static_cast<int>(rng() % 15);
    std::vector<std::string> labels;
    for (int i = 0; i < n; ++i) labels.push_back((i < 10 ? "t0" : "t") + std::to_string(i));
    std::shuffle(labels.begin(), labels.end(), rng);
    const qraft::TaxonRegistry taxa(labels);
    // Genes on every taxon or on a random subset of at least four.
    std::vector<GeneIndex> genes;
    const int k = 3 + static_cast<int>(rng() % 8);
    for (int g = 0; g < k; ++g) {
      std::vector<std::string> sub(labels);
      std::shuffle(sub.begin(), sub.end(), rng);
      if (rng() % 2) sub.resize(4 + rng() % static_cast<unsigned>(n - 3));
      genes.push_back(index_gene(qraft::parse_newick(guide_test_newick(rng, sub, 0.2)), taxa));
    }
    const int m = 4 + static_cast<int>(rng() % static_cast<unsigned>(n - 3));
    std::vector<std::string> sub(labels.begin(), labels.begin() + m);
    for (int i = 0; i < static_cast<int>(rng() % 3); ++i) sub.push_back("x" + std::to_string(i));
    std::shuffle(sub.begin(), sub.end(), rng);
    const qraft::Tree guide = qraft::parse_newick(guide_test_newick(rng, sub, std::array<double, 3>{0.0, 0.2, 0.4}[rng() % 3]));
    const GuideFilter filter(guide, taxa);
    const StartingTree plain = build_mre_tree(genes, taxa.size(), options);
    const StartingTree seeded = build_mre_tree(genes, taxa.size(), options, &filter);
    const qraft::Tree p = qraft::parse_newick(starting_tree_newick(plain, taxa));
    const qraft::Tree s = qraft::parse_newick(starting_tree_newick(seeded, taxa));
    const qraft::Tree r = qraft::parse_newick(starting_tree_newick(apply_guide_tree(seeded, guide, taxa), taxa));
    const std::set<std::string> all(labels.begin(), labels.end());
    const std::set<std::string> shared(labels.begin(), labels.begin() + m);
    const auto p_all = guide_test_splits(p, all), s_all = guide_test_splits(s, all), r_all = guide_test_splits(r, all);
    if (s_all != p_all) ++filtered;
    for (const auto& split : s_all) require(r_all.count(split) == 1, "guide seed: the guide contracts no seeded edge");
    const auto g_sh = guide_test_splits(guide, shared);
    for (const auto& split : guide_test_splits(s, shared)) {
      for (const auto& gs : g_sh) require(guide_test_compatible(split, gs, shared), "guide seed: the seeded MRE agrees with the guide");
    }
    const auto r_sh = guide_test_splits(r, shared);
    for (const auto& gs : g_sh) require(r_sh.count(gs) == 1, "guide seed: every guide split is displayed");
  }
  require(filtered > 50, "guide seed: most random cases filter the MRE");
  std::cout << "guide seed: 200 random cases, " << filtered << " filtered\n";
}

void test_nni_gains_match_enumeration() {
  using namespace qraft::internal;
  // Every scored NNI alternative's gain is the change in the tree's quartet
  // score, the number of gene quartets it displays; half the genes miss taxa.
  std::mt19937 rng(20260928u);
  for (int trial = 0; trial < 6; ++trial) {
    const int taxon_count = 7 + trial % 3;
    const std::vector<std::string> labels = numbered_labels(taxon_count);
    SolverContext ctx;
    ctx.options = test_options(2);
    ctx.workers = configure_threads(ctx.options);
    ctx.taxa = qraft::TaxonRegistry(labels);
    std::vector<PhysicalLeaves> genes;
    for (int g = 0; g < 14; ++g) {
      const qraft::Tree tree = qraft::parse_newick(
          g % 2 == 0 ? random_newick_over(rng, labels, false)
                     : random_incomplete_gene(rng, labels, g));
      ctx.genes.push_back(index_gene(tree, ctx.taxa));
      genes.push_back(physical_leaves(tree, ctx.taxa));
    }
    const SolverTree tree = solver_tree_from(
        qraft::parse_newick(random_newick_over(rng, labels, false)), ctx.taxa);
    std::vector<int> taxa(static_cast<std::size_t>(taxon_count));
    std::iota(taxa.begin(), taxa.end(), 0);
    const std::vector<NniMove> moves =
        nni_scored_moves_for_tests(ctx, tree, taxa);
    // A binary tree on n taxa has n-3 internal edges, two alternatives each.
    require(static_cast<int>(moves.size()) == 2 * (taxon_count - 3),
            "NNI did not score every internal-edge alternative");
    const double base = enumerated_quartet_score(tree, genes, taxon_count);
    for (const NniMove& move : moves) {
      SolverTree moved = tree;
      remove_solver_edge(moved, move.edge_a, move.swap_a_neighbor);
      remove_solver_edge(moved, move.edge_b, move.swap_b_neighbor);
      add_solver_edge(moved, move.edge_a, move.swap_b_neighbor);
      add_solver_edge(moved, move.edge_b, move.swap_a_neighbor);
      require_close(move.gain,
                    enumerated_quartet_score(moved, genes, taxon_count) - base,
                    "NNI gain disagrees with the enumerated quartet score");
    }
  }
}

// --support localpp: every internal edge's quartet frequencies equal the
// enumeration over the quartets with one taxon in each of its four arms, for
// genes that miss taxa or have polytomies, summed raw (with the product of
// the arm sizes) or normalized per gene. The posterior and length formulas
// match ASTRAL 5.15.5 (-t 2) on edges of the avian data set.
void test_branch_support_matches_enumeration() {
  using namespace qraft::internal;
  std::mt19937 rng(20261022u);
  const auto near = [](double x, double y, double tolerance) {
    return std::fabs(x - y) <= tolerance * std::max(1.0, std::fabs(y));
  };
  for (int trial = 0; trial < 8; ++trial) {
    const int taxon_count = 7 + trial % 4;
    const std::vector<std::string> labels = numbered_labels(taxon_count);
    SolverContext ctx;
    ctx.options = test_options(1 + trial % 2);
    ctx.workers = configure_threads(ctx.options);
    ctx.taxa = qraft::TaxonRegistry(labels);
    std::vector<PhysicalLeaves> genes;
    for (int g = 0; g < 14; ++g) {
      const qraft::Tree tree = qraft::parse_newick(
          g % 2 == 0 ? random_newick_over(rng, labels, g % 4 == 0)
                     : random_incomplete_gene(rng, labels, g));
      ctx.genes.push_back(index_gene(tree, ctx.taxa));
      genes.push_back(physical_leaves(tree, ctx.taxa));
    }
    const SolverTree tree = solver_tree_from(
        qraft::parse_newick(random_newick_over(rng, labels, false)), ctx.taxa);
    const std::vector<int> taxa = iota_taxa(taxon_count);
    // The taxa on v's side of the edge {u, v}.
    const auto side = [&](int u, int v) {
      std::vector<int> out, stack{v};
      std::vector<char> seen(tree.nodes.size(), 0);
      seen[static_cast<std::size_t>(u)] = 1;
      while (!stack.empty()) {
        const int x = stack.back();
        stack.pop_back();
        if (seen[static_cast<std::size_t>(x)]) continue;
        seen[static_cast<std::size_t>(x)] = 1;
        if (tree.nodes[static_cast<std::size_t>(x)].taxon >= 0) {
          out.push_back(tree.nodes[static_cast<std::size_t>(x)].taxon);
        }
        for (int y : tree.nodes[static_cast<std::size_t>(x)].neighbors) stack.push_back(y);
      }
      return out;
    };
    for (bool per_gene : {false, true}) {
      const std::vector<BranchQuartetFrequencies> branches =
          branch_quartet_frequencies(ctx, tree, taxa, per_gene);
      require(static_cast<int>(branches.size()) == taxon_count - 3,
              "support: not every internal edge has frequencies");
      for (const BranchQuartetFrequencies& f : branches) {
        // A and B at edge_a, C and D at edge_b: the edge's topology is AB|CD.
        std::vector<std::vector<int>> arms;
        for (int end : {f.edge_a, f.edge_b}) {
          const int other = end == f.edge_a ? f.edge_b : f.edge_a;
          for (int nb : tree.nodes[static_cast<std::size_t>(end)].neighbors) {
            if (nb != other) arms.push_back(side(end, nb));
          }
        }
        require(arms.size() == 4, "support: an edge end is not binary");
        double product = 1.0;
        for (const auto& arm : arms) product *= static_cast<double>(arm.size());
        require(near(f.quartets_around, product, 1e-12),
                "support: the quartets around an edge are not the arm product");
        std::array<double, 3> expected{0.0, 0.0, 0.0};
        for (const PhysicalLeaves& gene : genes) {
          const std::vector<int> leaf = leaf_of_taxon(gene, taxon_count);
          std::array<double, 3> counts{0.0, 0.0, 0.0};
          for (int a : arms[0])
            for (int b : arms[1])
              for (int c : arms[2])
                for (int d : arms[3]) {
                  const int r = quartet_resolution(
                      gene, leaf[static_cast<std::size_t>(a)], leaf[static_cast<std::size_t>(b)],
                      leaf[static_cast<std::size_t>(c)], leaf[static_cast<std::size_t>(d)]);
                  if (r >= 0) counts[static_cast<std::size_t>(r)] += 1.0;
                }
          const double total = counts[0] + counts[1] + counts[2];
          for (int i = 0; i < 3; ++i) {
            expected[static_cast<std::size_t>(i)] +=
                per_gene ? (total > 0 ? counts[static_cast<std::size_t>(i)] / total : 0.0)
                         : counts[static_cast<std::size_t>(i)];
          }
        }
        require(near(f.main, expected[0], 1e-9),
                "support: the main topology's frequency disagrees with enumeration");
        const double low = std::min(f.alternative1, f.alternative2);
        const double high = std::max(f.alternative1, f.alternative2);
        require(near(low, std::min(expected[1], expected[2]), 1e-9) &&
                    near(high, std::max(expected[1], expected[2]), 1e-9),
                "support: an alternative's frequency disagrees with enumeration");
      }
    }
  }
  // ASTRAL 5.15.5 -t 2 on the avian data set: f1, f2, f3 -> pp1, pp2, pp3
  // and the length.
  struct Expected {
    std::array<double, 3> f;
    std::array<double, 3> pp;
    double length;
  };
  const Expected rows[] = {
      {{4406.606119728906, 4271.914800348541, 4508.47907992262},
       {0.06367428344763769, 0.01629879717556038, 0.9200269193768018},
       0.001207065759626291},
      {{4750.399249076201, 4701.293869176065, 4676.30688174768},
       {0.563649419044178, 0.2510263489550199, 0.18532423200164516},
       0.004333736460375526},
      {{4893.25106612055, 4871.405182431321, 4681.343751448167},
       {0.6082431383541551, 0.35362147079068557, 0.038135390854625535},
       0.00808803645596889},
      {{12940.053843575597, 417.3152541414592, 252.6309022829439},
       {1.0, 0.0, 0.0},
       2.604479570250765},
  };
  for (const Expected& row : rows) {
    const std::array<double, 3> pp = local_posterior(row.f);
    for (int i = 0; i < 3; ++i) {
      require(std::fabs(pp[static_cast<std::size_t>(i)] - row.pp[static_cast<std::size_t>(i)]) < 1e-9,
              "support: the local posterior disagrees with ASTRAL");
    }
    require(near(coalescent_length(row.f[0], row.f[0] + row.f[1] + row.f[2]), row.length, 1e-9),
            "support: the coalescent length disagrees with ASTRAL");
  }
  std::cout << "branch support: frequencies match enumeration in both modes; "
               "posteriors and lengths match ASTRAL\n";
}

// --support localpp end to end: a fixed starting tree written with a
// degree-two root keeps its topology and gets one label on each of its five
// internal edges, in both label formats and alike at 1 and 4 workers. The
// edge {A,B} has the frequencies counted by hand below. No gene holds both G
// and H, so the edge {G,H} has no information: posterior 1/3 and no length.
// A starting tree with a polytomy leaves the edge into it unlabelled.
void test_branch_support_output() {
  const auto start_path =
      std::filesystem::temp_directory_path() / "qraft-support-start.tre";
  {
    std::ofstream out(start_path);
    out << "((((A,B),C),D),((E,F),(G,H)));\n";
  }
  // Around {A,B} the arms are A, B, C and D..H. Per gene, normalized: genes
  // 1-3 and 5 resolve every such quartet as AB|C.; gene 4 as AC|B.; gene 6
  // has one AC|BD and two of the form A.|BC, so (0, 1/3, 2/3); the polytomy
  // of gene 7 resolves none. So f = (4, 4/3, 2/3) over EN = 6 genes, and the
  // length is ln((6 + 1) / (1.5 (6 - 4 + 1))) = 0.441833.
  const std::vector<std::string> gene_text{
      "((A,B),(C,(D,(E,F))));", "((A,B),(C,(D,(E,F))));",
      "((A,B),(C,(D,(E,F))));", "((A,C),(B,(D,E)));",
      "(((A,B),(G,E)),(C,F));", "((A,(H,F)),((B,D),C));",
      "(A,B,C,(D,E));"};
  std::ostringstream pp_text;
  pp_text << std::setprecision(6)
          << qraft::internal::local_posterior({4.0, 4.0 / 3.0, 2.0 / 3.0})[0];
  const std::string pp = pp_text.str();
  const std::string all = "ABCDEFGH";
  // Each internal edge of the output as its side without A, with its label.
  const auto labelled_edges = [&](const std::string& newick) {
    std::map<std::string, std::string> out;
    std::vector<std::string> open;
    std::size_t i = 0;
    while (i < newick.size()) {
      const char c = newick[i];
      if (c == '(') {
        open.emplace_back();
        ++i;
      } else if (c == ')') {
        std::string side = open.back();
        open.pop_back();
        ++i;
        std::string text;
        if (i < newick.size() && newick[i] == '\'') {
          const std::size_t close = newick.find('\'', i + 1);
          require(close != std::string::npos, "support: unterminated label");
          text = newick.substr(i, close + 1 - i);
          i = close + 1;
        }
        while (i < newick.size() && newick[i] != ',' && newick[i] != ')' &&
               newick[i] != ';') {
          text += newick[i++];
        }
        if (open.empty()) continue;
        open.back() += side;
        if (side.size() < 2 || side.size() > all.size() - 2) continue;
        if (side.find('A') != std::string::npos) {
          std::string other;
          for (char x : all) {
            if (side.find(x) == std::string::npos) other += x;
          }
          side = other;
        }
        std::sort(side.begin(), side.end());
        require(out.emplace(side, text).second, "support: an edge appears twice");
      } else {
        if (std::isupper(static_cast<unsigned char>(c))) open.back() += c;
        ++i;
      }
    }
    return out;
  };
  for (bool detail : {false, true}) {
    std::string canonical;
    for (int workers : {1, 4}) {
      std::vector<qraft::GeneTree> genes;
      for (const std::string& text : gene_text) genes.push_back(qraft::parse_newick(text));
      qraft::Options options;
      options.threads = workers;
      options.progress = qraft::ProgressMode::Off;
      options.starting_tree_path = start_path.string();
      options.divide = false;
      options.nni = false;
      options.support = qraft::SupportMode::LocalPP;
      options.support_detail = detail;
      const std::string text =
          qraft::to_newick(qraft::infer_species_tree(std::move(genes), options));
      if (canonical.empty()) {
        canonical = text;
      } else {
        require(text == canonical, "support: the labelled tree changed with worker count");
      }
    }
    const std::map<std::string, std::string> edges = labelled_edges(canonical);
    require(edges.size() == 5, "support: the output does not have five internal edges");
    for (const char* side : {"CDEFGH", "DEFGH", "EFGH", "EF", "GH"}) {
      const auto it = edges.find(side);
      require(it != edges.end() && !it->second.empty(),
              "support: an internal edge has no label");
    }
    const std::string ab = edges.at("CDEFGH");
    const std::string gh = edges.at("GH");
    if (!detail) {
      require(ab == pp + ":0.441833", "support: the {A,B} label is wrong");
      require(gh == "0.333333", "support: an edge without information is labelled wrongly");
    } else {
      const std::string ab_end = ";CULength=0.441833]':0.441833";
      require(ab.find("'[f1=4;") == 0 && ab.find(";EN=6;q1=0.666667;pp1=" + pp + ";") !=
                  std::string::npos &&
                  (ab.find("f2=1.33333;f3=0.666667;") != std::string::npos ||
                   ab.find("f2=0.666667;f3=1.33333;") != std::string::npos) &&
                  ab.size() > ab_end.size() &&
                  ab.compare(ab.size() - ab_end.size(), ab_end.size(), ab_end) == 0,
              "support: the {A,B} detail label is wrong");
      require(gh.find(";EN=0;") != std::string::npos &&
                  gh.find(";pp1=0.333333;") != std::string::npos && gh.back() == '\'',
              "support: an edge without information has a detail label or length");
    }
  }
  // A polytomy of the starting tree stays: the edge into it has an end that
  // is not binary and gets no label, the edges with binary ends keep theirs.
  {
    std::ofstream out(start_path);
    out << "((A,B),(C,D,E),(F,(G,H)));\n";
  }
  std::vector<qraft::GeneTree> genes;
  for (const std::string& text : gene_text) genes.push_back(qraft::parse_newick(text));
  qraft::Options options;
  options.threads = 1;
  options.progress = qraft::ProgressMode::Off;
  options.starting_tree_path = start_path.string();
  options.divide = false;
  options.nni = false;
  options.support = qraft::SupportMode::LocalPP;
  const std::map<std::string, std::string> edges = labelled_edges(
      qraft::to_newick(qraft::infer_species_tree(std::move(genes), options)));
  require(edges.size() == 4 && edges.count("CDE") == 1 && edges.at("CDE").empty() &&
              !edges.at("CDEFGH").empty() && !edges.at("FGH").empty() &&
              !edges.at("GH").empty(),
          "support: a polytomy's edge is labelled or a binary edge is not");
  std::filesystem::remove(start_path);
  std::cout << "branch support output: every internal edge labelled, both formats, "
               "worker-count invariant, polytomies left unlabelled\n";
}

void test_nni_binary_rooted_starting_tree() {
  using namespace qraft::internal;
  // A starting tree with a binary root has an unlabelled vertex of degree two.
  // NNI scores only edges whose ends have three arms, so the root is
  // suppressed first: every internal edge is scored, and a wrong split across
  // the root edge is repaired.
  std::mt19937 rng(20261019u);
  for (int trial = 0; trial < 6; ++trial) {
    const int taxon_count = 6 + trial;
    const std::vector<std::string> labels = numbered_labels(taxon_count);
    SolverContext ctx;
    ctx.options = test_options(1 + trial % 2);
    ctx.workers = configure_threads(ctx.options);
    ctx.taxa = qraft::TaxonRegistry(labels);
    for (int g = 0; g < 8; ++g) {
      ctx.genes.push_back(index_gene(
          qraft::parse_newick(random_newick_over(rng, labels, false)), ctx.taxa));
    }
    // A binary root: join two random rooted halves.
    std::vector<std::string> left(labels.begin(), labels.begin() + taxon_count / 2);
    std::vector<std::string> right(labels.begin() + taxon_count / 2, labels.end());
    const auto rooted = [&](std::vector<std::string> parts) {
      while (parts.size() > 1) {
        std::shuffle(parts.begin(), parts.end(), rng);
        parts[0] = "(" + parts[0] + "," + parts[1] + ")";
        parts.erase(parts.begin() + 1);
      }
      return parts[0];
    };
    const qraft::Tree binary_rooted =
        qraft::parse_newick("(" + rooted(left) + "," + rooted(right) + ");");
    SolverTree tree = solver_tree_from(binary_rooted, ctx.taxa);
    std::vector<int> taxa(static_cast<std::size_t>(taxon_count));
    std::iota(taxa.begin(), taxa.end(), 0);
    require(static_cast<int>(nni_scored_moves_for_tests(ctx, tree, taxa).size()) <
                2 * (taxon_count - 3),
            "the binary-rooted fixture has every edge scored before suppression");
    suppress_solver_degree_two(tree);
    require(static_cast<int>(nni_scored_moves_for_tests(ctx, tree, taxa).size()) ==
                2 * (taxon_count - 3),
            "NNI did not score every internal edge of a binary-rooted starting tree");
  }

  // Eight genes all show ((A,B),(C,D)) | ((E,F),(G,H)); the binary-rooted
  // starting tree swaps (C,D) and (E,F) across its root edge.
  SolverContext ctx;
  ctx.options = test_options(1);
  ctx.workers = configure_threads(ctx.options);
  const std::vector<std::string> labels{"A", "B", "C", "D", "E", "F", "G", "H"};
  ctx.taxa = qraft::TaxonRegistry(labels);
  const qraft::Tree truth = qraft::parse_newick("(((A,B),(C,D)),((E,F),(G,H)));");
  std::vector<PhysicalLeaves> genes;
  for (int g = 0; g < 8; ++g) {
    ctx.genes.push_back(index_gene(truth, ctx.taxa));
    genes.push_back(physical_leaves(truth, ctx.taxa));
  }
  ctx.starting_tree = starting_tree_from(
      qraft::parse_newick("(((A,B),(E,F)),((C,D),(G,H)));"), ctx.taxa);
  const double optimum = enumerated_quartet_score(solver_tree_from(truth, ctx.taxa), genes, 8);
  const SolverTree start = refine_starting_tree(ctx, false);
  require(enumerated_quartet_score(start, genes, 8) < optimum,
          "the swapped starting tree is not worse than the truth");
  const SolverTree refined = refine_starting_tree(ctx, true);
  require_close(enumerated_quartet_score(refined, genes, 8), optimum,
                "NNI did not repair a split across a binary root");
  std::cout << "NNI on binary-rooted starting trees: 6 random trees fully scored; root-edge swap repaired\n";
}

void test_nni_local_pass_matches_full_pass() {
  using namespace qraft::internal;
  // A local pass scores its edges on genes collapsed to the owners of their
  // neighbourhood; every gain must equal the full pass's exactly, with
  // duplicated and missing taxa, for single edges, random edge sets and all
  // edges.
  std::mt19937 rng(20261001u);
  for (int trial = 0; trial < 8; ++trial) {
    const int taxon_count = 8 + 2 * (trial % 4);
    const std::vector<std::string> labels = numbered_labels(taxon_count);
    SolverContext ctx;
    ctx.options = test_options(2);
    ctx.workers = configure_threads(ctx.options);
    ctx.taxa = qraft::TaxonRegistry(labels);
    for (int g = 0; g < 16; ++g) {
      std::string newick;
      if (g % 4 == 3) {
        std::vector<std::string> gene_labels = labels;
        gene_labels.push_back(labels[rng() % labels.size()]);
        gene_labels.push_back(labels[rng() % labels.size()]);
        newick = random_newick_over(rng, gene_labels, g % 8 == 3);
      } else {
        newick = g % 2 == 0 ? random_newick_over(rng, labels, false)
                            : random_incomplete_gene(rng, labels, g);
      }
      ctx.genes.push_back(index_gene(qraft::parse_newick(newick), ctx.taxa));
    }
    const SolverTree tree = solver_tree_from(
        qraft::parse_newick(random_newick_over(rng, labels, false)), ctx.taxa);
    std::vector<int> taxa(static_cast<std::size_t>(taxon_count));
    std::iota(taxa.begin(), taxa.end(), 0);
    const std::vector<NniMove> full = nni_scored_moves_for_tests(ctx, tree, taxa);
    std::vector<std::pair<int, int>> edges;
    for (const NniMove& move : full) {
      if (edges.empty() || edges.back() != std::pair{move.edge_a, move.edge_b}) {
        edges.emplace_back(move.edge_a, move.edge_b);
      }
    }
    std::vector<std::vector<std::pair<int, int>>> subsets;
    for (const auto& edge : edges) subsets.push_back({edge});
    for (int s = 0; s < 6; ++s) {
      std::vector<std::pair<int, int>> subset = edges;
      std::shuffle(subset.begin(), subset.end(), rng);
      subset.resize(1 + rng() % subset.size());
      subsets.push_back(subset);
    }
    subsets.push_back(edges);
    for (const auto& subset : subsets) {
      const std::vector<NniMove> local =
          nni_local_scored_moves_for_tests(ctx, tree, taxa, subset);
      require(local.size() == 2 * subset.size(),
              "NNI local pass did not score every requested alternative");
      for (const NniMove& move : local) {
        const auto match = std::find_if(full.begin(), full.end(), [&](const NniMove& other) {
          return other.edge_a == move.edge_a && other.edge_b == move.edge_b &&
                 other.swap_a_neighbor == move.swap_a_neighbor &&
                 other.swap_b_neighbor == move.swap_b_neighbor &&
                 other.alternative == move.alternative;
        });
        require(match != full.end() && match->gain == move.gain,
                "NNI local-pass gain differs from the full pass");
      }
    }
  }
}

void test_nni_batch_gains_add() {
  using namespace qraft::internal;
  // One --nni-batch pass applies a set of moves whose gains add: the
  // enumerated score change of applying them all equals their summed gain.
  // The set is also a best one under the rule it implements (no shared
  // vertex; of two central edges one edge apart, both moves keep that edge,
  // each move using either of its two implementations), checked by brute
  // force. Species trees are binary, have polytomies, or are rooted (a
  // vertex of degree two); half the genes miss taxa.
  std::mt19937 rng(20261013u);
  int multi_move_batches = 0;
  for (int trial = 0; trial < 18; ++trial) {
    const int taxon_count = 8 + trial % 3;
    const std::vector<std::string> labels = numbered_labels(taxon_count);
    SolverContext ctx;
    ctx.options = test_options(2);
    ctx.workers = configure_threads(ctx.options);
    ctx.taxa = qraft::TaxonRegistry(labels);
    std::vector<PhysicalLeaves> genes;
    for (int g = 0; g < 14; ++g) {
      const qraft::Tree tree = qraft::parse_newick(
          g % 2 == 0 ? random_newick_over(rng, labels, false)
                     : random_incomplete_gene(rng, labels, g));
      ctx.genes.push_back(index_gene(tree, ctx.taxa));
      genes.push_back(physical_leaves(tree, ctx.taxa));
    }
    std::string newick;
    if (trial % 3 == 2) {
      std::vector<std::string> parts = labels;
      while (parts.size() > 2) {
        std::shuffle(parts.begin(), parts.end(), rng);
        parts[1] = "(" + parts[0] + "," + parts[1] + ")";
        parts.erase(parts.begin());
      }
      newick = "(" + parts[0] + "," + parts[1] + ");";
    } else {
      newick = random_newick_over(rng, labels, trial % 3 == 1);
    }
    const SolverTree tree =
        solver_tree_from(qraft::parse_newick(newick), ctx.taxa);
    std::vector<int> taxa(static_cast<std::size_t>(taxon_count));
    std::iota(taxa.begin(), taxa.end(), 0);
    const auto adjacent = [](const SolverTree& t, int a, int b) {
      const std::vector<int>& nbs = t.nodes[static_cast<std::size_t>(a)].neighbors;
      return std::find(nbs.begin(), nbs.end(), b) != nbs.end();
    };

    const std::vector<NniMove> batch = nni_batch_for_tests(ctx, tree, taxa);
    SolverTree moved = tree;
    double batch_gain = 0.0;
    for (const NniMove& move : batch) {
      require(adjacent(moved, move.edge_a, move.edge_b) &&
                  adjacent(moved, move.edge_a, move.swap_a_neighbor) &&
                  adjacent(moved, move.edge_b, move.swap_b_neighbor),
              "NNI batch move does not apply after the earlier ones");
      remove_solver_edge(moved, move.edge_a, move.swap_a_neighbor);
      remove_solver_edge(moved, move.edge_b, move.swap_b_neighbor);
      add_solver_edge(moved, move.edge_a, move.swap_b_neighbor);
      add_solver_edge(moved, move.edge_b, move.swap_a_neighbor);
      batch_gain += move.gain;
    }
    if (batch.size() >= 2) ++multi_move_batches;
    require_close(enumerated_quartet_score(moved, genes, taxon_count) -
                      enumerated_quartet_score(tree, genes, taxon_count),
                  batch_gain, "NNI batch gains do not add");

    // Brute force over every central edge: no move, or one of its
    // alternatives of positive gain.
    std::vector<std::vector<NniMove>> by_edge;
    for (const NniMove& move : nni_scored_moves_for_tests(ctx, tree, taxa)) {
      if (by_edge.empty() ||
          by_edge.back().front().edge_a != move.edge_a ||
          by_edge.back().front().edge_b != move.edge_b) {
        by_edge.push_back({});
      }
      by_edge.back().push_back(move);
    }
    // The side neighbour of `end` cut by implementation `which` of `move`.
    const auto cut_at = [&](const NniMove& move, int end, int which) {
      const int other = end == move.edge_a ? move.edge_b : move.edge_a;
      const int swapped =
          end == move.edge_a ? move.swap_a_neighbor : move.swap_b_neighbor;
      if (which == 0) return swapped;
      for (int nb : tree.nodes[static_cast<std::size_t>(end)].neighbors) {
        if (nb != other && nb != swapped) return nb;
      }
      return -1;
    };
    const auto allowed = [&](const std::vector<NniMove>& set) {
      for (std::size_t i = 0; i < set.size(); ++i) {
        for (std::size_t j = i + 1; j < set.size(); ++j) {
          for (int x : {set[i].edge_a, set[i].edge_b}) {
            if (x == set[j].edge_a || x == set[j].edge_b) return false;
          }
        }
      }
      for (unsigned mask = 0; mask < (1u << set.size()); ++mask) {
        bool keeps = true;
        for (std::size_t i = 0; i < set.size() && keeps; ++i) {
          for (std::size_t j = 0; j < set.size() && keeps; ++j) {
            if (i == j) continue;
            for (int x : {set[i].edge_a, set[i].edge_b}) {
              for (int y : {set[j].edge_a, set[j].edge_b}) {
                if (adjacent(tree, x, y) &&
                    cut_at(set[i], x, static_cast<int>((mask >> i) & 1u)) == y) {
                  keeps = false;
                }
              }
            }
          }
        }
        if (keeps) return true;
      }
      return false;
    };
    double best = 0.0;
    std::vector<NniMove> chosen;
    const auto search = [&](const auto& self, std::size_t edge) -> void {
      if (edge == by_edge.size()) {
        if (!allowed(chosen)) return;
        double total = 0.0;
        for (const NniMove& move : chosen) total += move.gain;
        best = std::max(best, total);
        return;
      }
      self(self, edge + 1);
      for (const NniMove& move : by_edge[edge]) {
        if (!(move.gain > 1e-9)) continue;
        chosen.push_back(move);
        self(self, edge + 1);
        chosen.pop_back();
      }
    };
    search(search, 0);
    require_close(batch_gain, best,
                  "NNI batch is not a best set of moves whose gains add");
  }
  require(multi_move_batches >= 3,
          "NNI batch fixtures did not select several moves at once");
}

void test_wide_count_exact() {
  using namespace qraft::internal;
  // NNI gains are summed in 128 bits and written as exact decimals: sums
  // past 2^63 and 2^64 stay exact where a double or a long long would not.
  require(wide_decimal(WideCount{}) == "0", "wide count: zero");
  WideCount past_63;
  wide_add(past_63, 9223372036854775808.0);  // 2^63
  for (int i = 0; i < 5; ++i) wide_add(past_63, 1.0);
  require(wide_decimal(past_63) == "9223372036854775813",
          "wide count: unit steps past 2^63");
  WideCount carry;
  wide_add(carry, 18446744073709549568.0);  // 2^64 - 2048
  wide_add(carry, 4096.0);
  require(wide_decimal(carry) == "18446744073709553664",
          "wide count: carry into the high word");
  WideCount repeated;
  for (int i = 0; i < 11; ++i) wide_add(repeated, 9e18);
  require(wide_decimal(repeated) == "99000000000000000000",
          "wide count: eleven times 9e18");
  WideCount large;
  wide_add(large, std::ldexp(1.0, 100));
  wide_add(large, std::ldexp(3.0, 60));
  require(wide_decimal(large) == "1267650600231688166010523746304",
          "wide count: addends past 2^64");
  WideCount rounded;
  wide_add(rounded, 2.6);
  require(wide_decimal(rounded) == "3", "wide count: rounds to the nearest integer");
  const auto throws = [](double value, WideCount start) {
    try {
      wide_add(start, value);
    } catch (const std::overflow_error&) {
      return true;
    }
    return false;
  };
  require(throws(-1.0, WideCount{}), "wide count: accepted a negative addend");
  require(throws(std::nan(""), WideCount{}), "wide count: accepted NaN");
  require(throws(2.0, WideCount{~0ull, ~0ull - 1}), "wide count: passed 2^128");
}

void test_bipartition_test_incomplete_genes() {
  using namespace qraft::internal;
  // On genes that miss taxa the test is J on shared genes: x in A and y in B
  // are compared on the genes holding a, b, x and y, against t(k_xy, R_xy).
  // Checks the index, the anchor choice, the running maximum and minimum, the
  // star of pairs, the per-pair threshold and the decision against an
  // independent implementation on breadth-first distances, on sparse random
  // single-copy genes.
  {
    // A gene holding one taxon is a single node of degree 0; it must not turn
    // the test off for the run.
    const qraft::TaxonRegistry taxa({"A", "B", "C", "D"});
    const auto index = make_bipartition_test_index(
        {index_gene(qraft::parse_newick("((A,B),(C,D));"), taxa),
         index_gene(qraft::parse_newick("A;"), taxa)},
        4);
    require(index->single_copy && !index->complete_single_copy &&
                index->gene_taxa == std::vector<int>({4, 1}),
            "a one-taxon gene disabled the incomplete-gene test");
  }
  {
    // An empty gene holds no taxa; the test still runs on the others. With
    // anchors A and C, B is compared with D on the one gene holding all four.
    SolverContext ctx;
    ctx.options = test_options(1);
    ctx.workers = configure_threads(ctx.options);
    ctx.taxa = qraft::TaxonRegistry({"A", "B", "C", "D"});
    ctx.genes.push_back(index_gene(qraft::parse_newick("((A,B),(C,D));"), ctx.taxa));
    ctx.genes.push_back(index_gene(qraft::GeneTree{}, ctx.taxa));
    ctx.genes.push_back(index_gene(qraft::parse_newick("((A,B),C);"), ctx.taxa));
    ctx.bipartition_test_index = make_bipartition_test_index(ctx.genes, 4);
    require(ctx.bipartition_test_index->single_copy &&
                !ctx.bipartition_test_index->complete_single_copy &&
                ctx.bipartition_test_index->taxon_genes == std::vector<int>({2, 2, 2, 1}) &&
                ctx.bipartition_test_index->gene_taxa == std::vector<int>({4, 0, 3}),
            "an empty gene broke the incomplete-gene test index");
    Bitset split(4);
    split.set(0);
    split.set(1);
    const BipartitionTestResult result = test_bipartition(ctx, {0, 1, 2, 3}, split);
    require(result.available && result.anchor_a == 0 && result.anchor_b == 2 &&
                result.threshold_genes == 1 && result.gap == 1.0,
            "the incomplete-gene test did not run beside an empty gene");
    // A side holding only its anchor has no pair to test and accepts.
    Bitset one(4);
    one.set(3);
    const BipartitionTestResult single = test_bipartition(ctx, {0, 1, 2, 3}, one);
    require(single.available && single.accepted && single.threshold_genes == 0,
            "a one-taxon side was not accepted by the incomplete-gene test");
  }
  {
    // The star is not the all-pairs rule on genes that miss taxa. With
    // anchors a and b, x1 and y1 are the champions and the three pairs tested
    // have gap 1 on two genes each, while the untested pair (x2, y2) attaches
    // in reverse on its two shared genes, with gap -1.
    SolverContext ctx;
    ctx.options = test_options(2);
    ctx.workers = configure_threads(ctx.options);
    ctx.taxa = qraft::TaxonRegistry({"a", "x1", "x2", "b", "y1", "y2"});
    for (const char* tree : {"((a,x2),(b,x1));", "((a,y1),(b,y2));", "((a,x1),(b,y1));",
                             "((a,x1),(b,y2));", "((a,x2),(b,y1));", "((a,y2),(b,x2));"}) {
      for (int copy = 0; copy < 2; ++copy) {
        ctx.genes.push_back(index_gene(qraft::parse_newick(tree), ctx.taxa));
      }
    }
    ctx.bipartition_test_index = make_bipartition_test_index(ctx.genes, 6);
    std::vector<int> active;
    for (const char* label : {"a", "x1", "x2", "b", "y1", "y2"}) {
      active.push_back(ctx.taxa.id_for(label));
    }
    Bitset split(6);
    for (int i = 0; i < 3; ++i) split.set(i);
    const BipartitionTestResult result = test_bipartition(ctx, active, split);
    require(result.available && result.accepted && result.anchor_a == active[0] &&
                result.anchor_b == active[3] && result.gap == 1.0 &&
                result.threshold_genes == 2,
            "the shared-gene test did not test the star of the champions");
  }
  // t(k, R) = C R / sqrt(k) * (sqrt(log(k + e)) - sqrt(log(K0 + e))).
  const auto threshold_of = [](long long k, double range, double c, double k0) {
    const long double e = std::numbers::e_v<long double>;
    return static_cast<double>(
        c * static_cast<long double>(range) / std::sqrt(static_cast<long double>(k)) *
        (std::sqrt(std::log(k + e)) - std::sqrt(std::log(k0 + e))));
  };
  constexpr double kCoefficient = 0.05;
  constexpr double kCrossover = 20.0;
  std::mt19937 rng(20261013u);
  int tested = 0;
  int vacuous = 0;
  int skipped_first = 0;
  int untestable = 0;
  int rejected = 0;
  for (int trial = 0; trial < 80; ++trial) {
    const int n = 7 + trial % 3;
    const bool sparse = trial % 4 == 3;
    const std::vector<std::string> labels = numbered_labels(n);
    SolverContext ctx;
    ctx.options = test_options(1 + trial % 3);
    ctx.workers = configure_threads(ctx.options);
    ctx.options.bipartition_test.crossover = kCrossover;
    ctx.options.bipartition_test.coefficient = kCoefficient;
    ctx.taxa = qraft::TaxonRegistry(labels);
    const int gene_count = sparse ? 3 + trial % 4 : 6 + trial % 25;
    for (int g = 0; g < gene_count; ++g) {
      std::vector<std::string> kept = labels;
      std::shuffle(kept.begin(), kept.end(), rng);
      kept.resize(static_cast<std::size_t>(sparse ? 4 : 4 + rng() % (n - 3)));
      ctx.genes.push_back(index_gene(
          qraft::parse_newick(random_newick_over(rng, kept, g % 4 == 0)), ctx.taxa));
    }
    ctx.bipartition_test_index = make_bipartition_test_index(ctx.genes, n);
    const BipartitionTestIndex& index = *ctx.bipartition_test_index;
    require(index.single_copy, "incomplete single-copy genes have no test index");
    if (index.complete_single_copy) continue;
    // Breadth-first distances in every gene, between taxa it holds.
    std::vector<std::vector<std::vector<int>>> dist(ctx.genes.size());
    for (std::size_t g = 0; g < ctx.genes.size(); ++g) {
      const GeneIndex& gene = ctx.genes[g];
      dist[g].assign(static_cast<std::size_t>(n), std::vector<int>(static_cast<std::size_t>(n), -1));
      for (int a = 0; a < n; ++a) {
        const int start = gene.taxon_to_node[static_cast<std::size_t>(a)];
        if (start < 0) continue;
        std::vector<int> row(gene.nodes.size(), -1);
        std::queue<int> pending;
        row[static_cast<std::size_t>(start)] = 0;
        pending.push(start);
        while (!pending.empty()) {
          const int u = pending.front();
          pending.pop();
          for (int v : gene_neighbors(gene, u)) {
            if (row[static_cast<std::size_t>(v)] >= 0) continue;
            row[static_cast<std::size_t>(v)] = row[static_cast<std::size_t>(u)] + 1;
            pending.push(v);
          }
        }
        for (int b = 0; b < n; ++b) {
          const int node = gene.taxon_to_node[static_cast<std::size_t>(b)];
          if (node >= 0) dist[g][static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] = row[static_cast<std::size_t>(node)];
        }
      }
    }
    std::vector<int> held(static_cast<std::size_t>(n), 0);
    std::vector<int> gene_taxa(ctx.genes.size(), 0);
    for (std::size_t g = 0; g < ctx.genes.size(); ++g) {
      for (int t = 0; t < n; ++t) {
        if (ctx.genes[g].taxon_to_node[static_cast<std::size_t>(t)] < 0) continue;
        ++held[static_cast<std::size_t>(t)];
        ++gene_taxa[g];
      }
    }
    require(index.taxon_genes == held, "incomplete test index counts genes per taxon wrongly");
    require(index.gene_taxa == gene_taxa, "incomplete test index counts taxa per gene wrongly");
    // h_g(z) with anchors a and b, or -1 when gene g lacks a, b or z.
    const auto position = [&](std::size_t g, int a, int b, int z) {
      const auto& d = dist[g];
      const auto at = [&](int x, int y) {
        return d[static_cast<std::size_t>(x)][static_cast<std::size_t>(y)];
      };
      if (at(a, a) < 0 || at(b, b) < 0 || at(z, z) < 0) return -1;
      return (at(a, z) + at(a, b) - at(b, z)) / 2;
    };
    for (int rep = 0; rep < 8; ++rep) {
      std::vector<int> active(static_cast<std::size_t>(n));
      std::iota(active.begin(), active.end(), 0);
      std::shuffle(active.begin(), active.end(), rng);
      active.resize(static_cast<std::size_t>(4 + rng() % (n - 3)));
      const int m = static_cast<int>(active.size());
      Bitset split(m);
      const int left = 1 + static_cast<int>(rng() % static_cast<unsigned>(m - 1));
      for (int i = 0; i < left; ++i) split.set(i);
      // The two taxa of each side held by the most genes, ties by position.
      const auto top_two = [&](bool side) {
        std::vector<int> positions;
        for (int i = 0; i < m; ++i) {
          if (split.test(i) == side) positions.push_back(i);
        }
        std::sort(positions.begin(), positions.end(), [&](int x, int y) {
          const int hx = held[static_cast<std::size_t>(active[static_cast<std::size_t>(x)])];
          const int hy = held[static_cast<std::size_t>(active[static_cast<std::size_t>(y)])];
          return hx != hy ? hx > hy : x < y;
        });
        std::array<int, 2> top{-1, -1};
        for (std::size_t r = 0; r < 2 && r < positions.size(); ++r) {
          top[r] = active[static_cast<std::size_t>(positions[r])];
        }
        return top;
      };
      const std::array<int, 2> top_a = top_two(true);
      const std::array<int, 2> top_b = top_two(false);
      const std::array<std::pair<int, int>, 4> pairs{{
          {top_a[0], top_b[0]}, {top_a[1], top_b[1]}, {top_a[0], top_b[1]}, {top_a[1], top_b[0]}}};
      struct Expected {
        bool usable = false;
        double gap = 0.0;
        double threshold = 0.0;
        int genes = 0;
      };
      // One anchor pair, from the definition: the running maximum x* of A and
      // minimum y* of B, then every pair (x*, y) and (x, y*) sharing a gene.
      const auto expected_for = [&](int a, int b) {
        Expected out;
        std::vector<int> side_a;
        std::vector<int> side_b;
        for (int i = 0; i < m; ++i) {
          const int z = active[static_cast<std::size_t>(i)];
          if (split.test(i) ? z != a : z != b) (split.test(i) ? side_a : side_b).push_back(z);
        }
        if (side_a.empty() || side_b.empty()) {
          out.usable = true;
          out.gap = std::numeric_limits<double>::infinity();
          return out;
        }
        struct Sums {
          long long difference = 0;
          long long genes = 0;
          long long taxa = 0;
        };
        const auto shared = [&](int x, int y) {
          Sums s;
          for (std::size_t g = 0; g < ctx.genes.size(); ++g) {
            const int hx = position(g, a, b, x);
            const int hy = position(g, a, b, y);
            if (hx < 0 || hy < 0) continue;
            s.difference += hy - hx;
            ++s.genes;
            s.taxa += gene_taxa[g];
          }
          return s;
        };
        const auto champion = [&](const std::vector<int>& side, bool maximum) {
          int c = -1;
          for (int z : side) {
            if (c < 0) {
              for (std::size_t g = 0; g < ctx.genes.size() && c < 0; ++g) {
                if (position(g, a, b, z) >= 0) c = z;
              }
              continue;
            }
            const long long difference = shared(c, z).difference;
            if (maximum ? difference > 0 : difference < 0) c = z;
          }
          return c;
        };
        const int xs = champion(side_a, true);
        const int ys = champion(side_b, false);
        if (xs < 0 || ys < 0) return out;
        std::vector<std::pair<int, int>> tested_pairs{{xs, ys}};
        for (int y : side_b) {
          if (y != ys) tested_pairs.emplace_back(xs, y);
        }
        for (int x : side_a) {
          if (x != xs) tested_pairs.emplace_back(x, ys);
        }
        double best = std::numeric_limits<double>::infinity();
        for (const auto& [x, y] : tested_pairs) {
          const Sums s = shared(x, y);
          if (s.genes == 0) {
            if (x == xs && y == ys) return out;
            continue;
          }
          const double gap = static_cast<double>(s.difference) / static_cast<double>(s.genes);
          const double range =
              2.0 * (static_cast<double>(s.taxa) / static_cast<double>(s.genes) - 1.0);
          const double threshold = threshold_of(s.genes, range, kCoefficient, kCrossover);
          if (gap - threshold < best) {
            best = gap - threshold;
            out.gap = gap;
            out.threshold = threshold;
            out.genes = static_cast<int>(s.genes);
          }
        }
        out.usable = true;
        return out;
      };
      // With one anchor pair, the first usable pair decides.
      Expected expect;
      std::pair<int, int> expect_pair{-1, -1};
      bool first_valid = true;
      for (int p = 0; p < 4 && !expect.usable; ++p) {
        const auto anchors = pairs[static_cast<std::size_t>(p)];
        if (anchors.first < 0 || anchors.second < 0 ||
            std::find(pairs.begin(), pairs.begin() + p, anchors) != pairs.begin() + p) {
          continue;
        }
        expect = expected_for(anchors.first, anchors.second);
        if (expect.usable) {
          expect_pair = anchors;
          if (!first_valid) ++skipped_first;
        }
        first_valid = false;
      }
      ctx.options.bipartition_test.anchor_pairs = 1;
      const BipartitionTestResult result = test_bipartition(ctx, active, split);
      require(result.available == expect.usable,
              "incomplete-gene test availability disagrees with the definition");
      if (!expect.usable) {
        ++untestable;
        continue;
      }
      require(result.anchor_a == expect_pair.first && result.anchor_b == expect_pair.second,
              "incomplete-gene test chose the wrong anchors");
      require(result.threshold_genes == expect.genes,
              "incomplete-gene test reported the wrong shared-gene count");
      if (std::isinf(expect.gap)) {
        ++vacuous;
        require(std::isinf(result.gap) && result.gap > 0.0 && result.accepted,
                "a side holding only its anchor was not accepted");
        continue;
      }
      ++tested;
      require_close(result.gap, expect.gap,
                    "incomplete-gene gap disagrees with breadth-first attachment positions");
      require_close(result.threshold, expect.threshold,
                    "incomplete-gene threshold is not t(k_xy, R_xy) of the deciding pair");
      require(result.accepted == (expect.gap > expect.threshold),
              "incomplete-gene decision disagrees with gap and threshold");
      if (!result.accepted) ++rejected;
      // More anchor pairs can only add accepting pairs, and are skipped the same way.
      ctx.options.bipartition_test.anchor_pairs = 4;
      const BipartitionTestResult wide = test_bipartition(ctx, active, split);
      require(wide.available &&
                  wide.gap - wide.threshold >= result.gap - result.threshold &&
                  (!result.accepted || wide.accepted) &&
                  wide.accepted == (wide.gap > wide.threshold),
              "more anchor pairs weakened the incomplete-gene test");
    }
  }
  require(tested > 100 && rejected > 0 && vacuous > 0 && skipped_first > 0 && untestable > 0,
          "incomplete-gene fixtures did not cover tested, rejected, one-taxon, skipped and "
          "untestable cases: tested=" + std::to_string(tested) +
              " rejected=" + std::to_string(rejected) + " vacuous=" + std::to_string(vacuous) +
              " skipped=" + std::to_string(skipped_first) +
              " untestable=" + std::to_string(untestable));

  // A full run on genes that miss taxa uses the test and does not depend on
  // the worker count.
  std::vector<qraft::GeneTree> genes;
  std::mt19937 gene_rng(20261014u);
  const std::vector<std::string> labels = numbered_labels(24);
  for (int g = 0; g < 60; ++g) {
    std::vector<std::string> drop = labels;
    std::shuffle(drop.begin(), drop.end(), gene_rng);
    drop.resize(6);
    std::vector<std::string> keep;
    for (const std::string& label : labels) {
      if (std::find(drop.begin(), drop.end(), label) == drop.end()) keep.push_back(label);
    }
    genes.push_back(qraft::parse_newick(random_newick_over(gene_rng, keep, false)));
  }
  const auto stats_path = std::filesystem::temp_directory_path() / "qraft-incomplete-test.json";
  std::string first;
  for (int workers : {1, 4}) {
    qraft::Options options = test_options(workers);
    options.stats_json = stats_path.string();
    const std::string text = qraft::to_newick(qraft::infer_species_tree(genes, options));
    if (first.empty()) first = text;
    else require(text == first, "incomplete-gene test output changed with workers");
    const std::string stats = read_file(stats_path);
    require(stats.find("\"bipartition_test_available\": true") != std::string::npos &&
                stats.find("\"bipartition_test_complete_genes\": false") != std::string::npos &&
                stats.find("\"bipartition_test_incomplete\": 0,") == std::string::npos &&
                stats.find("\"bipartition_test_shared_genes\": 0,") == std::string::npos,
            "a run on genes that miss taxa did not use the bipartition test");
  }
  std::filesystem::remove(stats_path);
}

void test_bipartition_test_shared_genes_reduce_to_j() {
  using namespace qraft::internal;
  // On complete genes, J on shared genes is the J test: Delta(x, y) is
  // (J(y) - J(x)) / 2 on every gene, every pair has t_k, and x* and y* are the
  // extremes of J. Both paths must give the same decision, separation (to the
  // last bit), threshold, anchors and gene count.
  std::mt19937 rng(20261017u);
  int compared = 0;
  int accepted = 0;
  int rejected = 0;
  for (int trial = 0; trial < 45; ++trial) {
    const int n = 6 + trial % 7;
    const std::vector<std::string> labels = numbered_labels(n);
    SolverContext ctx;
    ctx.options = test_options(1 + trial % 3);
    ctx.workers = configure_threads(ctx.options);
    ctx.options.bipartition_test.crossover = std::array<double, 3>{2.0, 20.0, 10000.0}[trial % 3];
    ctx.options.bipartition_test.anchor_pairs = 1;
    ctx.taxa = qraft::TaxonRegistry(labels);
    const int gene_count = 3 + trial % 30;
    for (int g = 0; g < gene_count; ++g) {
      ctx.genes.push_back(index_gene(
          qraft::parse_newick(random_newick_over(rng, labels, g % 3 == 0)), ctx.taxa));
    }
    ctx.bipartition_test_index = make_bipartition_test_index(ctx.genes, n);
    require(ctx.bipartition_test_index->complete_single_copy,
            "complete random genes have no complete test index");
    for (int rep = 0; rep < 10; ++rep) {
      std::vector<int> active(static_cast<std::size_t>(n));
      std::iota(active.begin(), active.end(), 0);
      std::shuffle(active.begin(), active.end(), rng);
      active.resize(static_cast<std::size_t>(4 + rng() % (n - 3)));
      const int m = static_cast<int>(active.size());
      Bitset split(m);
      const int left = 2 + static_cast<int>(rng() % static_cast<unsigned>(m - 3));
      for (int i = 0; i < left; ++i) split.set(i);
      const BipartitionTestResult plain = test_bipartition(ctx, active, split);
      const BipartitionTestResult shared =
          test_bipartition_on_shared_genes_for_tests(ctx, active, split);
      require(plain.available && shared.available && plain.accepted == shared.accepted &&
                  plain.gap == shared.gap && plain.threshold == shared.threshold &&
                  plain.threshold_genes == shared.threshold_genes &&
                  plain.anchor_a == shared.anchor_a && plain.anchor_b == shared.anchor_b,
              "J on shared genes is not the J test on complete genes");
      ++compared;
      ++(plain.accepted ? accepted : rejected);
    }
  }
  require(compared == 450 && accepted > 0 && rejected > 0,
          "the reduction to J did not see both decisions");
}

void test_bipartition_test_shared_genes_noise_free() {
  using namespace qraft::internal;
  // Genes that are restrictions of one binary species tree to the taxa they
  // hold, with every four taxa sharing a gene and a small positive threshold:
  // every pair across a split attaches at least one edge apart in every gene,
  // and the running maximum and minimum find the witness of a false split, so
  // the test accepts exactly the splits of the species tree.
  std::mt19937 rng(20261018u);
  int true_seen = 0;
  int false_seen = 0;
  for (int trial = 0; trial < 30; ++trial) {
    const int n = 7 + trial % 4;
    const std::vector<std::string> labels = numbered_labels(n);
    // A random rooted binary species tree: leaves 0..n-1, internal nodes after.
    std::vector<std::array<int, 2>> children(static_cast<std::size_t>(n), {-1, -1});
    std::vector<int> roots(static_cast<std::size_t>(n));
    std::iota(roots.begin(), roots.end(), 0);
    while (roots.size() > 1) {
      std::shuffle(roots.begin(), roots.end(), rng);
      children.push_back({roots[0], roots[1]});
      roots.erase(roots.begin(), roots.begin() + 2);
      roots.push_back(static_cast<int>(children.size()) - 1);
    }
    const int root = roots.front();
    const int nodes = static_cast<int>(children.size());
    // The species tree restricted to `kept`, unary vertices suppressed.
    const std::function<std::string(int, const std::vector<char>&)> restricted =
        [&](int v, const std::vector<char>& kept) -> std::string {
      if (v < n) return kept[static_cast<std::size_t>(v)] ? labels[static_cast<std::size_t>(v)] : "";
      const std::string l = restricted(children[static_cast<std::size_t>(v)][0], kept);
      const std::string r = restricted(children[static_cast<std::size_t>(v)][1], kept);
      if (l.empty()) return r;
      if (r.empty()) return l;
      return "(" + l + "," + r + ")";
    };
    // Edge-count distances between the leaves of the species tree.
    std::vector<std::vector<int>> adjacency(static_cast<std::size_t>(nodes));
    for (int v = n; v < nodes; ++v) {
      for (int c : children[static_cast<std::size_t>(v)]) {
        adjacency[static_cast<std::size_t>(v)].push_back(c);
        adjacency[static_cast<std::size_t>(c)].push_back(v);
      }
    }
    std::vector<std::vector<int>> species(static_cast<std::size_t>(n));
    for (int a = 0; a < n; ++a) {
      std::vector<int> row(static_cast<std::size_t>(nodes), -1);
      std::queue<int> pending;
      row[static_cast<std::size_t>(a)] = 0;
      pending.push(a);
      while (!pending.empty()) {
        const int u = pending.front();
        pending.pop();
        for (int v : adjacency[static_cast<std::size_t>(u)]) {
          if (row[static_cast<std::size_t>(v)] >= 0) continue;
          row[static_cast<std::size_t>(v)] = row[static_cast<std::size_t>(u)] + 1;
          pending.push(v);
        }
      }
      species[static_cast<std::size_t>(a)].assign(row.begin(), row.begin() + n);
    }
    const auto d = [&](int x, int y) {
      return species[static_cast<std::size_t>(x)][static_cast<std::size_t>(y)];
    };
    // The species tree displays xx'|yy'.
    const auto displays = [&](int x, int x2, int y, int y2) {
      return d(x, x2) + d(y, y2) < std::min(d(x, y) + d(x2, y2), d(x, y2) + d(x2, y));
    };
    SolverContext ctx;
    ctx.options = test_options(1 + trial % 3);
    ctx.workers = configure_threads(ctx.options);
    ctx.options.bipartition_test.crossover = 0.0;
    ctx.options.bipartition_test.coefficient = 1e-6;
    ctx.taxa = qraft::TaxonRegistry(labels);
    std::vector<std::vector<char>> holds;
    for (int g = 0; g < 40; ++g) {
      std::vector<char> kept(static_cast<std::size_t>(n), 1);
      const int drop = 1 + static_cast<int>(rng() % 2);
      for (int k = 0; k < drop; ++k) kept[static_cast<std::size_t>(rng() % n)] = 0;
      holds.push_back(kept);
      ctx.genes.push_back(index_gene(qraft::parse_newick(restricted(root, kept) + ";"), ctx.taxa));
    }
    ctx.bipartition_test_index = make_bipartition_test_index(ctx.genes, n);
    require(ctx.bipartition_test_index->single_copy &&
                !ctx.bipartition_test_index->complete_single_copy,
            "restricted species trees did not give an incomplete single-copy index");
    for (int rep = 0; rep < 12; ++rep) {
      std::vector<int> active(static_cast<std::size_t>(n));
      std::iota(active.begin(), active.end(), 0);
      std::shuffle(active.begin(), active.end(), rng);
      active.resize(static_cast<std::size_t>(4 + rng() % (n - 3)));
      const int m = static_cast<int>(active.size());
      // Every four active taxa must share a gene.
      bool covered = true;
      for (int i = 0; i < m && covered; ++i) {
        for (int j = i + 1; j < m && covered; ++j) {
          for (int k = j + 1; k < m && covered; ++k) {
            for (int l = k + 1; l < m && covered; ++l) {
              covered = std::any_of(holds.begin(), holds.end(), [&](const std::vector<char>& h) {
                return h[static_cast<std::size_t>(active[static_cast<std::size_t>(i)])] &&
                       h[static_cast<std::size_t>(active[static_cast<std::size_t>(j)])] &&
                       h[static_cast<std::size_t>(active[static_cast<std::size_t>(k)])] &&
                       h[static_cast<std::size_t>(active[static_cast<std::size_t>(l)])];
              });
            }
          }
        }
      }
      if (!covered) continue;
      Bitset split(m);
      std::vector<int> side_a;
      std::vector<int> side_b;
      for (int i = 0; i < m; ++i) {
        if (rng() % 2) split.set(i);
      }
      if (split.count() == 0) split.set(0);
      if (split.count() == m) split.reset(0);
      for (int i = 0; i < m; ++i) {
        (split.test(i) ? side_a : side_b).push_back(active[static_cast<std::size_t>(i)]);
      }
      bool truth = true;
      for (int x : side_a) {
        for (int x2 : side_a) {
          for (int y : side_b) {
            for (int y2 : side_b) {
              if (x < x2 && y < y2 && !displays(x, x2, y, y2)) truth = false;
            }
          }
        }
      }
      for (int pairs : {1, 4}) {
        ctx.options.bipartition_test.anchor_pairs = pairs;
        const BipartitionTestResult result = test_bipartition(ctx, active, split);
        require(result.available && result.accepted == truth,
                std::string("noise-free genes: the shared-gene test ") +
                    (truth ? "rejected a split of" : "accepted a bipartition outside") +
                    " the species tree");
      }
      ++(truth ? true_seen : false_seen);
    }
  }
  require(true_seen > 20 && false_seen > 20,
          "noise-free fixtures did not cover true and false bipartitions: true=" +
              std::to_string(true_seen) + " false=" + std::to_string(false_seen));
}

void test_bipartition_test_threshold_values() {
  using namespace qraft::internal;
  // t_k = C R / sqrt(k) * (sqrt(log(k+e)) - sqrt(log(K0+e))), R = 2(n-1).
  const auto expected = [](int k, int n, double c, double k0) {
    const long double e = std::numbers::e_v<long double>;
    return static_cast<double>(
        c * 2.0L * (n - 1) / std::sqrt(static_cast<long double>(k)) *
        (std::sqrt(std::log(k + e)) - std::sqrt(std::log(k0 + e))));
  };
  for (const auto& [k, n] : std::vector<std::pair<int, int>>{
           {1, 4}, {50, 200}, {1000, 1001}, {9999, 10}, {10001, 37}, {100000, 50}}) {
    require_close(bipartition_test_threshold(k, n, 0.05, 10000.0),
                  expected(k, n, 0.05, 10000.0),
                  "bipartition-test threshold at k=" + std::to_string(k), 0.0, 1e-12);
  }
  require(bipartition_test_threshold(10000, 100, 0.05, 10000.0) == 0.0,
          "bipartition-test threshold must vanish at the crossover");
  require(bipartition_test_threshold(1000, 100, 0.05, 10000.0) < 0.0 &&
              bipartition_test_threshold(20000, 100, 0.05, 10000.0) > 0.0,
          "bipartition-test threshold has the wrong sign around the crossover");
}

// Gene trees under the multispecies coalescent on a fixed species tree.
class CoalescentSimulator {
 public:
  struct Branch {
    int left = -1;
    int right = -1;
    std::string label;   // leaves only
    double length = 0.0; // coalescent units above this node; root: unbounded
  };

  CoalescentSimulator(std::vector<Branch> nodes, int root, std::uint64_t seed)
      : nodes_(std::move(nodes)), root_(root), state_(seed) {}

  std::string gene_tree() {
    std::vector<std::string> lineages = lineages_above(root_);
    return lineages.front() + ";";
  }

 private:
  // xorshift64*: identical on every platform, unlike std::distributions.
  double uniform() {
    state_ ^= state_ >> 12;
    state_ ^= state_ << 25;
    state_ ^= state_ >> 27;
    return static_cast<double>((state_ * 0x2545F4914F6CDD1DULL) >> 11) * 0x1.0p-53;
  }

  std::vector<std::string> lineages_above(int node) {
    const Branch& branch = nodes_[static_cast<std::size_t>(node)];
    std::vector<std::string> lineages;
    if (branch.left < 0) {
      lineages.push_back(branch.label);
    } else {
      lineages = lineages_above(branch.left);
      const std::vector<std::string> right = lineages_above(branch.right);
      lineages.insert(lineages.end(), right.begin(), right.end());
    }
    double time = 0.0;
    while (lineages.size() > 1) {
      const double l = static_cast<double>(lineages.size());
      time += -std::log(1.0 - uniform()) / (l * (l - 1.0) / 2.0);
      if (node != root_ && time > branch.length) break;
      const std::size_t a = static_cast<std::size_t>(uniform() * l);
      std::size_t b = static_cast<std::size_t>(uniform() * (l - 1.0));
      if (b >= a) ++b;
      const std::string joined = "(" + lineages[a] + "," + lineages[b] + ")";
      lineages.erase(lineages.begin() + static_cast<std::ptrdiff_t>(std::max(a, b)));
      lineages.erase(lineages.begin() + static_cast<std::ptrdiff_t>(std::min(a, b)));
      lineages.push_back(joined);
    }
    return lineages;
  }

  std::vector<Branch> nodes_;
  int root_;
  std::uint64_t state_;
};

void test_msc_simulation_recovers_species_tree() {
  // True gene trees from the multispecies coalescent on a nine-taxon species
  // tree with short internal branches: every pipeline variant that carries
  // the consistency theorem must return the species tree.
  using Branch = CoalescentSimulator::Branch;
  // ((((a,b),c),(d,e)),((f,g),(h,i))) with internal lengths in [0.2, 0.6].
  std::vector<Branch> nodes;
  const auto leaf = [&](const std::string& label) {
    nodes.push_back(Branch{-1, -1, label, 0.0});
    return static_cast<int>(nodes.size()) - 1;
  };
  const auto join = [&](int left, int right, double length) {
    nodes.push_back(Branch{left, right, "", length});
    return static_cast<int>(nodes.size()) - 1;
  };
  const int ab = join(leaf("a"), leaf("b"), 0.3);
  const int abc = join(ab, leaf("c"), 0.2);
  const int de = join(leaf("d"), leaf("e"), 0.5);
  const int abcde = join(abc, de, 0.4);
  const int fg = join(leaf("f"), leaf("g"), 0.25);
  const int hi = join(leaf("h"), leaf("i"), 0.6);
  const int fghi = join(fg, hi, 0.35);
  const int root = join(abcde, fghi, 0.0);
  const qraft::Tree truth =
      qraft::parse_newick("((((a,b),c),(d,e)),((f,g),(h,i)));");

  CoalescentSimulator simulator(nodes, root, 0x9e3779b97f4a7c15ULL);
  std::vector<qraft::GeneTree> genes;
  for (int g = 0; g < 400; ++g) genes.push_back(qraft::parse_newick(simulator.gene_tree()));

  qraft::Options defaults = test_options(4);
  qraft::Options no_nni = defaults;
  no_nni.nni = false;
  for (const qraft::Options& options : {defaults, no_nni}) {
    const qraft::SpeciesTree species = qraft::infer_species_tree(genes, options);
    require(qraft::rf_distance(species, truth).first == 0,
            "coalescent simulation did not recover the species tree");
  }
}

void test_trivial_bipartition_completion() {
  const auto four_taxon_starting_tree_path =
      std::filesystem::temp_directory_path() / "qraft-singleton-four-star.tre";
  {
    std::ofstream out(four_taxon_starting_tree_path);
    out << "(A,B,C,D);\n";
  }
  const qraft::Tree expected = qraft::parse_newick("((A,C),(B,D));");
  std::string canonical;
  for (int workers : {1, 4}) {
    std::vector<qraft::GeneTree> genes{
        qraft::parse_newick("((A,C),(B,D));"),
        qraft::parse_newick("(A,B,C,D);"),
        qraft::parse_newick("(A,B,C,D);"),
    };
    qraft::Options options;
    options.threads = workers;
    options.progress = qraft::ProgressMode::Off;
    options.starting_tree_path = four_taxon_starting_tree_path.string();
    options.nni = false;
    const std::filesystem::path stats_path =
        std::filesystem::temp_directory_path() /
        ("qraft-singleton-" + std::to_string(workers) + ".json");
    options.stats_json = stats_path.string();
    const qraft::SpeciesTree inferred =
        qraft::infer_species_tree(std::move(genes), options);
    const auto [rf, nrf] = qraft::rf_distance(inferred, expected);
    (void)nrf;
    require(rf == 0,
            "trivial bipartition and regraft did not recover the supported "
            "four-taxon split");
    const std::string serialized = qraft::to_newick(inferred);
    if (canonical.empty()) {
      canonical = serialized;
    } else {
      require(serialized == canonical,
              "singleton output changed with worker count");
    }
    const qraft::internal::RunStats& stats =
        qraft::internal::run_stats();
    require(stats.mre_accepted_splits.load(std::memory_order_relaxed) == 0,
            "external four-taxon fixture did not keep a star starting tree");
    require(stats.recursion_skewed_splits.load(
                std::memory_order_relaxed) == 1,
            "unresolved four-taxon query did not use one singleton completion");
    require(stats.bipartitions_scored.load(std::memory_order_relaxed) == 0,
            "singleton completion unnecessarily entered K=2 scoring");
    std::filesystem::remove(stats_path);
  }

  std::filesystem::remove(four_taxon_starting_tree_path);

  // Five taxa and a star starting tree force two trivial bipartitions in a
  // row, independently of MRE. The second regraft must not meet an unlabeled
  // degree-two connector left by the first: such a connector is not a vertex
  // with three arms and would give its corridor incomplete scores.
  const auto five_taxon_starting_tree_path =
      std::filesystem::temp_directory_path() / "qraft-singleton-five-star.tre";
  {
    std::ofstream out(five_taxon_starting_tree_path);
    out << "(81,83,84,91,98);\n";
  }
  for (qraft::RegraftScore regraft_mode :
       {qraft::RegraftScore::Raw, qraft::RegraftScore::Normalized}) {
    const std::vector<const char*> newicks{
        "((81,83),(91,98),84);",
        "((81,84),(83,91),98);",
        "((81,84),(83,98),91);",
        "((81,91),(83,84),98);",
        "((81,98),(84,91),83);",
        "((81,98),(84,91),83);",
    };
    std::vector<qraft::GeneTree> genes;
    genes.reserve(newicks.size());
    for (const char* newick : newicks) {
      genes.push_back(qraft::parse_newick(newick));
    }
    qraft::Options options;
    options.threads = 1;
    options.progress = qraft::ProgressMode::Off;
    options.starting_tree_path = five_taxon_starting_tree_path.string();
    options.regraft_score = regraft_mode;
    options.nni = false;
    const std::filesystem::path stats_path =
        std::filesystem::temp_directory_path() /
        (std::string("qraft-singleton-five-taxon-") +
         (regraft_mode == qraft::RegraftScore::Raw ? "raw" : "normalized") +
         ".json");
    options.stats_json = stats_path.string();
    qraft::internal::reset_run_stats();
    const qraft::SpeciesTree inferred =
        qraft::infer_species_tree(std::move(genes), options);
    const qraft::Tree five_taxon_truth =
        qraft::parse_newick("((81,98),(84,91),83);");
    const auto [rf, nrf] = qraft::rf_distance(inferred, five_taxon_truth);
    (void)nrf;
    require(rf == 0,
            "multi-step singleton completion did not maximize the "
            "five-taxon attachment objective");
    const qraft::internal::RunStats& stats =
        qraft::internal::run_stats();
    require(stats.mre_accepted_splits.load(std::memory_order_relaxed) == 0 &&
                stats.recursion_skewed_splits.load(
                    std::memory_order_relaxed) == 2 &&
                stats.recursion_max_depth.load(std::memory_order_relaxed) == 2 &&
                stats.bipartitions_scored.load(
                    std::memory_order_relaxed) == 0 &&
                stats.regraft_edge_ties.load(std::memory_order_relaxed) == 0,
            "five-taxon regression did not exercise two singleton completions");
    std::filesystem::remove(stats_path);
  }

  std::filesystem::remove(five_taxon_starting_tree_path);

  {
    const std::filesystem::path starting_tree_path =
        std::filesystem::temp_directory_path() /
        "qraft-default-external-star-start.tre";
    {
      std::ofstream starting_tree_out(starting_tree_path);
      starting_tree_out << "(T0,T1,T2,T3,T4,T5,T6,T7);\n";
    }
    std::string canonical;
    for (qraft::DivideScore split_mode :
         {qraft::DivideScore::Quartet,
          qraft::DivideScore::Balanced}) {
      const char* mode_name =
          split_mode == qraft::DivideScore::Quartet ? "raw" : "balanced";
      for (int workers : {1, 4}) {
        std::vector<qraft::GeneTree> genes{qraft::parse_newick(
            "(((T0,T1),(T2,T3)),((T4,T5),(T6,T7)));")};
        qraft::Options options;
        options.threads = workers;
        options.progress = qraft::ProgressMode::Off;
        options.starting_tree_path = starting_tree_path.string();
        options.divide_score = split_mode;
        options.nni = false;
        const std::filesystem::path stats_path =
            std::filesystem::temp_directory_path() /
            ("qraft-default-external-star-" + std::string(mode_name) + "-" +
             std::to_string(workers) + ".json");
        options.stats_json = stats_path.string();
        const qraft::SpeciesTree inferred =
            qraft::infer_species_tree(std::move(genes), options);
        require(inferred.leaf_labels().size() == 8,
                "default external-star singleton completion lost a taxon");
        const std::string serialized = qraft::to_newick(inferred);
        if (canonical.empty()) {
          canonical = serialized;
        } else {
          require(serialized == canonical,
                  "default external-star singleton completion changed with "
                  "worker count or retained-edge selector");
        }
        const qraft::internal::RunStats& stats =
            qraft::internal::run_stats();
        require(
            stats.recursion_nonbase_calls.load(std::memory_order_relaxed) ==
                    5 &&
                stats.recursion_skewed_splits.load(
                    std::memory_order_relaxed) == 5 &&
                stats.recursion_balanced_splits.load(
                    std::memory_order_relaxed) == 0,
            "default external star did not realize the exact singleton "
            "recurrence");
        require(stats.recursion_taxon_mass.load(
                    std::memory_order_relaxed) == 38 &&
                    stats.recursion_max_depth.load(std::memory_order_relaxed) ==
                        5,
                "default external-star singleton completion changed skew mass "
                "or depth");
        require(stats.bipartitions_scored.load(
                    std::memory_order_relaxed) == 0,
                "default external-star singleton completion entered K=2 "
                "scoring");
        require(stats.complexity_contract_satisfied.load(
                    std::memory_order_relaxed) == 1,
                "default external-star singleton completion violated the "
                "complexity contract");
        std::filesystem::remove(stats_path);
      }
    }
    std::filesystem::remove(starting_tree_path);
  }

  {
    std::vector<qraft::GeneTree> genes{
        qraft::parse_newick("(T0,T1,T2,T3,T4,T5,T6,T7);")};
    qraft::Options options;
    options.threads = 1;
    options.progress = qraft::ProgressMode::Off;
    options.nni = false;
    const std::filesystem::path stats_path =
        std::filesystem::temp_directory_path() /
        "qraft-singleton-star-chain.json";
    options.stats_json = stats_path.string();
    const qraft::SpeciesTree inferred =
        qraft::infer_species_tree(std::move(genes), options);
    require(inferred.leaf_labels().size() == 8,
            "singleton completion chain lost a taxon");
    const qraft::internal::RunStats& stats =
        qraft::internal::run_stats();
    require(stats.recursion_nonbase_calls.load(std::memory_order_relaxed) == 5 &&
                stats.recursion_skewed_splits.load(
                    std::memory_order_relaxed) == 5 &&
                stats.recursion_balanced_splits.load(std::memory_order_relaxed) ==
                    0,
            "all-star input did not realize the exact singleton recurrence");
    require(stats.recursion_taxon_mass.load(
                std::memory_order_relaxed) == 38 &&
                stats.recursion_max_depth.load(std::memory_order_relaxed) == 5,
            "singleton completion chain changed skew active mass or depth");
    require(stats.bipartitions_scored.load(std::memory_order_relaxed) == 0,
            "all-star singleton chain entered K=2 scoring");
    require(stats.complexity_contract_satisfied.load(
                std::memory_order_relaxed) == 1,
            "singleton completion chain violated the complexity contract");
    std::filesystem::remove(stats_path);
  }

  {
    const std::vector<const char*> newicks{
        "((((A,B),C),D),(((E,F),G),H));",
        "((((A,B),C),D),(((E,F),G),H));",
        "((((A,B),C),D),(((E,F),G),H));",
    };
    std::string canonical_balanced;
    for (int workers : {1, 4}) {
      std::vector<qraft::GeneTree> genes;
      genes.reserve(newicks.size());
      for (const char* newick : newicks) {
        genes.push_back(qraft::parse_newick(newick));
      }
      qraft::Options options;
      options.threads = workers;
      options.progress = qraft::ProgressMode::Off;
      options.divide_score = qraft::DivideScore::Balanced;
      options.nni = false;
      const std::filesystem::path stats_path =
          std::filesystem::temp_directory_path() /
          ("qraft-retained-balanced-" + std::to_string(workers) + ".json");
      options.stats_json = stats_path.string();
      const qraft::SpeciesTree inferred =
          qraft::infer_species_tree(std::move(genes), options);
      const qraft::Tree expected =
          qraft::parse_newick("((((A,B),C),D),(((E,F),G),H));");
      const auto [rf, nrf] = qraft::rf_distance(inferred, expected);
      (void)nrf;
      require(rf == 0,
              "balanced retained-edge selection changed a fully "
              "supported tree");
      const std::string serialized = qraft::to_newick(inferred);
      if (canonical_balanced.empty()) {
        canonical_balanced = serialized;
      } else {
        require(serialized == canonical_balanced,
                "balanced retained-edge output changed with worker "
                "count");
      }
      const qraft::internal::RunStats& stats =
          qraft::internal::run_stats();
      require(stats.recursion_balanced_splits.load(
                  std::memory_order_relaxed) >= 1,
              "balanced retained-edge selection missed the available "
              "four-four starting-tree edge");
      require(stats.bipartitions_scored.load(
                  std::memory_order_relaxed) == 0,
              "balanced retained-edge selection entered K=2 S-V "
              "scoring");
      require(stats.complexity_contract_satisfied.load(
                  std::memory_order_relaxed) == 1,
              "balanced retained-edge selection violated the "
              "complexity contract");
      std::filesystem::remove(stats_path);
    }
  }
}

void test_child_gene_set_partition_bound() {
  using namespace qraft::internal;

  // A parsed gene may have one unary root, and a singleton gene has no
  // internal vertex. These cases check the general I <= L argument behind the
  // 3/2 partition bound, not only the binary case I <= L - 2.
  for (bool unary_root : {false, true}) {
    GeneIndex source;
    source.nodes.resize(unary_root ? 2 : 1);
    const int leaf = unary_root ? 1 : 0;
    source.nodes[static_cast<std::size_t>(leaf)].set_local(0);
    qraft::internal::set_gene_leaf_mass(source, leaf, 1.0);
    if (unary_root) {
      set_gene_adjacency(source, {{1}, {0}});
    } else {
      set_gene_adjacency(source, {{}});
    }
    source.root = 0;
    source.reference_root = leaf;
    source.root_traversal = traversal(source, source.root);

    const GeneIndex selected =
        restrict_gene_by_local(source, std::vector<int>{0});
    const GeneIndex empty = restrict_gene_by_local(source, std::vector<int>{-1});
    const long long input_nodes = static_cast<long long>(source.nodes.size());
    const long long selected_nodes =
        static_cast<long long>(selected.nodes.size());
    const long long empty_nodes = static_cast<long long>(empty.nodes.size());
    require(2 * (selected_nodes + empty_nodes) <= 3 * input_nodes,
            "singleton/unary child outputs exceed the 3/2 partition bound");
    const long long retained_limit =
        input_nodes - (input_nodes + 3) / 4;
    require(selected_nodes <= retained_limit || empty_nodes <= retained_limit,
            "singleton/unary children violate the factor-3/4 gate");
  }

  for (int leaf_count : {4, 5, 8, 17, 32, 65}) {
    for (unsigned seed = 0; seed < 32; ++seed) {
      std::mt19937 rng(seed * 65537u + static_cast<unsigned>(leaf_count));
      GeneIndex root = random_arbitrary_gene_index(leaf_count, rng());
      const int owner_count =
          2 + static_cast<int>(rng() % static_cast<unsigned>(
                  std::max(1, std::min(16, leaf_count) - 1)));
      std::vector<int> leaves = root.taxon_to_node;
      std::shuffle(leaves.begin(), leaves.end(), rng);
      for (int physical = 0; physical < leaf_count; ++physical) {
        GeneIndex::Node& leaf = root.nodes[static_cast<std::size_t>(
            leaves[static_cast<std::size_t>(physical)])];
        leaf.set_taxon(-1);
        leaf.set_local(physical % owner_count);
      }
      root.taxon_to_node.clear();

      std::vector<int> identity_map(static_cast<std::size_t>(owner_count));
      std::iota(identity_map.begin(), identity_map.end(), 0);
      const GeneIndex normalized = restrict_gene_by_local(root, identity_map);

      const std::array<const GeneIndex*, 2> sources{&root, &normalized};
      for (const GeneIndex* source : sources) {
        std::vector<int> left_map(static_cast<std::size_t>(owner_count), -1);
        std::vector<int> right_map(static_cast<std::size_t>(owner_count), -1);
        int left_count = 0;
        int right_count = 0;
        for (int owner = 0; owner < owner_count; ++owner) {
          // Owners zero and one force nonempty sides. The others may go to
          // neither side, like the inactive leaves a heavy reuse keeps.
          const int disposition =
              owner == 0 ? 0 : (owner == 1 ? 1 : static_cast<int>(rng() % 3));
          if (disposition == 0) {
            left_map[static_cast<std::size_t>(owner)] = left_count++;
          } else if (disposition == 1) {
            right_map[static_cast<std::size_t>(owner)] = right_count++;
          }
        }
        const GeneIndex left = restrict_gene_by_local(*source, left_map);
        const GeneIndex right = restrict_gene_by_local(*source, right_map);

        const long long input_nodes =
            static_cast<long long>(source->nodes.size());
        const long long left_nodes = static_cast<long long>(left.nodes.size());
        const long long right_nodes =
            static_cast<long long>(right.nodes.size());
        require(2 * (left_nodes + right_nodes) <= 3 * input_nodes,
                "child physical outputs exceed the 3/2 partition bound");
        const long long retained_limit =
            input_nodes - (input_nodes + 3) / 4;
        require(left_nodes <= retained_limit || right_nodes <= retained_limit,
                "both child physical outputs violate the factor-3/4 gate");
      }
    }
  }
}

void test_mre_complexity_primitives() {
  int previous_height = 0;
  for (int item_count : {257, 4096, 32768}) {
    const qraft::internal::MreComplexityProbeResult result =
        qraft::internal::mre_complexity_probe_for_tests(item_count);
    require(result.radix_matches_oracle,
            "MRE radix deduplication/aggregation diverged from its oracle");
    require(result.dense_split_ids,
            "MRE dense split ids are not contiguous key ranks");
    require(result.occurrence_radix_passes == 8 &&
                result.priority_radix_passes == 8 &&
                result.radix_passes == 16,
            "MRE radix schedules must execute exactly 8+8 stable passes");
    require(result.occurrence_radix_items ==
                8LL * result.emitted_occurrences &&
                result.priority_radix_items ==
                    8LL * result.distinct_splits &&
                result.radix_items == result.occurrence_radix_items +
                                          result.priority_radix_items,
            "MRE radix item touches are not fixed-linear in record counts");
    require(result.deduplicated_gene_splits >= result.distinct_splits &&
                result.distinct_splits == item_count &&
                result.emitted_occurrences >=
                    result.deduplicated_gene_splits,
            "MRE record deduplication/aggregation cardinalities are invalid");
    require(result.avl_invariants_hold,
            "MRE adversarial child-index sequence broke an AVL invariant");
    require(result.range_aggregates_match,
            "MRE AVL range/complement aggregate disagreed with a scan");
    require(result.max_avl_height <= result.avl_height_limit,
            "MRE child-index height exceeded its logarithmic bound");
    require(result.empty_complement_visits <=
                result.empty_complement_visit_limit,
            "MRE empty complement enumeration traversed the large side");
    require(result.flat_topology_matches,
            "MRE flat topology diverged from vector swap-pop semantics");
    require(result.topology_peak_slots <= result.topology_slot_limit,
            "MRE flat topology exceeded its linear slot bound");
    require(result.max_avl_height >= previous_height,
            "MRE adversarial height unexpectedly decreased with input size");
    previous_height = result.max_avl_height;
  }

  qraft::Options options;
  for (int invalid_taxa : {-1, std::numeric_limits<int>::max()}) {
    bool rejected = false;
    try {
      (void)qraft::internal::build_mre_tree({}, invalid_taxa, options);
    } catch (const std::length_error&) {
      rejected = true;
    }
    require(rejected,
            "MRE accepted dimensions outside its integer domain");
  }
}

void test_qraft_smoke() {
  qraft::Options options = test_options();
  const qraft::SpeciesTree inferred =
      qraft::infer_species_tree(small_gene_set(), options);
  require(sorted_labels(inferred) ==
              std::vector<std::string>({"A", "B", "C", "D", "E", "F"}),
          "QRAFT should infer a tree over all input taxa");
}

void test_worker_cap_by_input_size() {
  std::vector<qraft::GeneTree> genes = {
      qraft::parse_newick("((A,B),(C,D));")};
  const int physical_nodes = static_cast<int>(genes.front().nodes.size());
  require(physical_nodes > 1,
          "worker-cap fixture must contain multiple physical nodes");
  qraft::Options options = test_options(physical_nodes + 32);
  (void)qraft::infer_species_tree(std::move(genes), options);
  require(qraft::internal::get_num_workers() == physical_nodes,
          "effective executor workers must not exceed physical input work");
}

void test_worker_count_determinism() {
  require(qraft::internal::gene_partition_count(0) == 0 &&
              qraft::internal::gene_partition_count(1) == 1 &&
              qraft::internal::gene_partition_count(64) == 64 &&
              qraft::internal::gene_partition_count(65) == 65 &&
              qraft::internal::gene_partition_count(257) == 257,
          "score partitions must not impose a fixed gene/core ceiling");
  qraft::Options one = test_options(1);
  qraft::Options four = test_options(4);
  const std::string a = qraft::to_newick(qraft::infer_species_tree(small_gene_set(), one));
  const std::string b = qraft::to_newick(qraft::infer_species_tree(small_gene_set(), four));
  require(a == b, "QRAFT output should be deterministic across thread counts");

  const std::vector<qraft::GeneTree> base = small_gene_set();
  std::vector<qraft::GeneTree> above_old_limit;
  above_old_limit.reserve(65);
  for (int gene = 0; gene < 65; ++gene) {
    above_old_limit.push_back(
        base[static_cast<std::size_t>(gene) % base.size()]);
  }
  const std::string serial_above_old_limit = qraft::to_newick(
      qraft::infer_species_tree(above_old_limit, one));
  const std::string parallel_above_old_limit = qraft::to_newick(
      qraft::infer_species_tree(std::move(above_old_limit), four));
  require(serial_above_old_limit == parallel_above_old_limit,
          "per-gene score reduction should be deterministic above 64 "
          "gene partitions");
}

void test_pipeline_determinism() {
  for (bool divide : {true, false}) {
    std::string reference;
    for (int threads : {1, 2, 8, 16}) {
      qraft::Options options = test_options(threads);
      options.divide = divide;
      const std::string actual = qraft::to_newick(
          qraft::infer_species_tree(small_gene_set(), options));
      if (reference.empty()) {
        reference = actual;
      } else {
        require(actual == reference,
                "pipeline output should be deterministic across workers");
      }
    }
  }
}

void test_ablation_smoke() {
  // QRAFT, the starting tree alone, and the starting tree refined by NNI.
  for (const auto& [divide, nni] : std::vector<std::pair<bool, bool>>{
           {true, true}, {false, false}, {false, true}}) {
    qraft::Options options = test_options();
    options.divide = divide;
    options.nni = nni;
    const qraft::SpeciesTree inferred =
        qraft::infer_species_tree(small_gene_set(), options);
    require(sorted_labels(inferred) ==
                std::vector<std::string>({"A", "B", "C", "D", "E", "F"}),
            "ablation variant should retain all taxa");
  }
}

void test_rf_identity() {
  const qraft::Tree tree = qraft::parse_newick("((A,B),(C,D));");
  const auto [rf, nrf] = qraft::rf_distance(tree, tree);
  require(rf == 0, "RF of a tree with itself should be zero");
  require_close(nrf, 0.0, "NRF of a tree with itself should be zero");
}

void test_fixture_smoke() {
  const std::filesystem::path root = QRAFT_TEST_TREE_DIR;
  const std::filesystem::path genes_path = root / "astral2_200taxa_1000genes.tre";
  qraft::Options options = test_options(
      std::max(1u, std::min(4u, std::thread::hardware_concurrency())));
  const qraft::SpeciesTree inferred =
      qraft::infer_species_tree(parse_gene_tree_file(genes_path), options);
  require(inferred.leaf_labels().size() == 201,
          "1000-gene fixture should infer all 201 taxa");
}

void test_hdt_structural_bounds() {
  using namespace qraft::internal;
  int max_leaves = 512;
  if (const char* requested = std::getenv("QRAFT_TEST_STRUCTURAL_MAX")) {
    max_leaves = std::max(32, std::atoi(requested));
  }
  int fixed_live_caterpillar_units = -1;

  for (int nominal = 32; nominal <= max_leaves; nominal *= 2) {
    GeneIndex bit_reversed = balanced_gene_index(nominal);
    bit_reverse_gene_taxa(bit_reversed);
    std::vector<std::pair<std::string, GeneIndex>> fixtures;
    fixtures.emplace_back("balanced", balanced_gene_index(nominal));
    fixtures.emplace_back("caterpillar", caterpillar_gene_index(nominal));
    fixtures.emplace_back("bit-reversed", std::move(bit_reversed));
    fixtures.emplace_back("star", star_gene_index(nominal));
    fixtures.emplace_back(
        "alternating-corridor",
        alternating_corridor_gene_index(std::max(3, nominal / 2)));
    fixtures.emplace_back(
        "random-arbitrary-degree",
        random_linear_gene_index(nominal, 0x5eed0000u +
                                              static_cast<unsigned>(nominal)));
    fixtures.emplace_back("broom", broom_gene_index(nominal / 2));

    for (auto& [name, gene] : fixtures) {
      const int leaves = static_cast<int>(gene.taxon_to_node.size());
      const OwnerProjection owners = identity_owner_map(leaves);
      HdtColorPolicy colors;
      colors.current_by_owner.resize(static_cast<std::size_t>(leaves));
      colors.background_by_owner.assign(static_cast<std::size_t>(leaves), 0);
      for (int owner = 0; owner < leaves; ++owner) {
        colors.current_by_owner[static_cast<std::size_t>(owner)] =
            static_cast<std::uint8_t>(owner & 1);
      }
      const HdtState base = build_hdt_state(
          build_base_component_graph(gene, owners, 2, colors));
      for (int wanted : {1, 2, 4, 8, 16}) {
        if (wanted > leaves) continue;
        std::vector<int> selected;
        selected.reserve(static_cast<std::size_t>(wanted));
        for (int i = 0; i < wanted; ++i) {
          selected.push_back(static_cast<int>(
              static_cast<long long>(i) * leaves / wanted));
        }
        const HdtContractResult contracted = hdt_contract(
            base, hdt_owner_selector(leaves, selected),
            HdtContractMode::LightCopy);
        const int live = contracted.stats.selected_live_atoms;
        const int r = contracted.state.components.component_units();
        const int a = contracted.state.hdt.stats.initial_atoms;
        const int h = contracted.state.hdt.stats.reachable_nodes;
        const std::string label = name + " n=" + std::to_string(leaves) +
                                  " live=" + std::to_string(live);
        require(r <= 4 * std::max(1, live) + 3,
                label + " violated R <= 4*max(1,live)+3");
        require(a <= 2 * r, label + " violated A <= 2R");
        require(h == 2 * a - 1, label + " violated H == 2A-1");
        require(h <= 16 * std::max(1, live) + 16,
                label + " violated H live-set bound");
        require(contracted.state.hdt.stats.max_recompute_degree == 2,
                label + " rebuilt HDT is not binary");
        require(contracted.state.hdt.stats.local_balance_violations == 0,
                label + " rebuilt HDT is not locally balanced");
        require(contracted.state.hdt.stats.node_count_ok &&
                    contracted.state.hdt.stats.atom_bound_ok &&
                    contracted.stats.size_bound_ok,
                label + " structural bound flags failed");
        const long long work =
            contracted.stats.input_units + contracted.stats.fragment_hdt_nodes +
            contracted.stats.cc_ops + contracted.stats.gg_ops +
            contracted.stats.ig_ops + contracted.stats.c_to_g_ops;
        require(work <= 256LL * contracted.stats.input_units + 256,
                label + " Contract operation count is not linear");
      }

      if (name == "caterpillar" && nominal >= 64) {
        const std::vector<int> endpoints{0, leaves - 1};
        const HdtContractResult fixed = hdt_contract(
            base, hdt_owner_selector(leaves, endpoints),
            HdtContractMode::LightCopy);
        const int units = fixed.state.components.component_units();
        if (fixed_live_caterpillar_units < 0) {
          fixed_live_caterpillar_units = units;
        } else {
          require(units == fixed_live_caterpillar_units,
                  "fixed-live caterpillar Contract output grew with source n");
        }
      }
    }

    if (nominal > max_leaves / 2) break;
  }
}

void test_hdt_heavy_threshold() {
  using namespace qraft::internal;
  // The rooted three-leaf query has a two-leaf heavy child, so the Contract
  // bound is 4*2+3 = 11 and a heavy Contract runs only above 4*11 = 44 input
  // component units. A star gene has one center plus N leaves.
  struct Boundary {
    int leaves;
    long long expected_heavy_calls;
    const char* position;
  };
  for (const Boundary boundary : {
           Boundary{42, 0, "immediately below"},
           Boundary{43, 0, "at"},
           Boundary{44, 1, "immediately above"}}) {
    const GeneIndex gene = star_gene_index(boundary.leaves);
    const OwnerProjection owners = identity_owner_map(boundary.leaves);
    HdtColorPolicy colors;
    colors.current_by_owner.assign(static_cast<std::size_t>(boundary.leaves),
                                   1);
    colors.background_by_owner.assign(
        static_cast<std::size_t>(boundary.leaves), 0);
    HdtState state = build_hdt_state(
        build_base_component_graph(gene, owners, 2, colors));
    require(state.components.component_units() == boundary.leaves + 1,
            "heavy-threshold fixture has an unexpected input R");

    const ActiveTaxa active = active_set({0, 1, 2}, boundary.leaves);
    const QueryTree query = induce_query_tree(
        balanced_starting_tree(3), active);
    require(query.root >= 0 &&
                query.heavy_child[static_cast<std::size_t>(query.root)] >= 0,
            "heavy-threshold query did not expose a heavy child");
    const int heavy = query.heavy_child[static_cast<std::size_t>(query.root)];
    require(query.subtree_leaf_count[static_cast<std::size_t>(heavy)] == 2,
            "three-leaf threshold query heavy child should contain two atoms");
    std::vector<unsigned char> score_node(query.nodes.size(), 0);
    score_node[static_cast<std::size_t>(heavy)] = 1;
    HdtExtractWorkspace workspace;
    HdtTraversalStats stats;
    int visits = 0;
    count_bipartitions(
        std::move(state), query, query.leaf_active_order, score_node, workspace,
        stats,
        [&](int node, HdtState& current) {
          require(node == heavy, "threshold traversal visited an unmarked node");
          const std::array<double, 3> counts = hdt_root_counts(current);
          require_finite(counts[0], "threshold traversal support channel");
          require_finite(counts[1], "threshold traversal violation channel");
          ++visits;
        });
    require(visits == 1, "threshold traversal did not visit its heavy leaf");
    require(stats.heavy_calls == boundary.expected_heavy_calls,
            std::string("heavy gate fired ") + boundary.position +
                " its strict threshold");
    require(stats.contract_heavy_calls == boundary.expected_heavy_calls &&
                stats.heavy_shrink_violations == 0 &&
                stats.complexity_contract_satisfied,
            std::string("heavy threshold contract failed ") +
                boundary.position + " boundary");
  }
}

void test_hdt_units() {
  using namespace qraft::internal;
  struct Fixture {
    const char* name;
    GeneIndex gene;
    int leaves;
  };
  std::vector<Fixture> fixtures;
  fixtures.push_back({"star", star_gene_index(64), 64});
  fixtures.push_back({"caterpillar", caterpillar_gene_index(96), 96});
  fixtures.push_back({"broom", broom_gene_index(48), 96});
  fixtures.push_back({"random", random_arbitrary_gene_index(80, 73), 80});

  for (int palette : {2, 4}) {
    for (const Fixture& fixture : fixtures) {
      OwnerProjection owners = identity_owner_map(fixture.leaves);
      HdtColorPolicy colors;
      colors.current_by_owner.resize(static_cast<std::size_t>(fixture.leaves));
      colors.background_by_owner.resize(
          static_cast<std::size_t>(fixture.leaves));
      for (int owner = 0; owner < fixture.leaves; ++owner) {
        colors.current_by_owner[static_cast<std::size_t>(owner)] =
            static_cast<std::uint8_t>(owner % palette);
        colors.background_by_owner[static_cast<std::size_t>(owner)] =
            static_cast<std::uint8_t>(palette == 4 && owner % 7 == 0 ? 3 : 0);
      }
      HdtComponentGraph graph = build_base_component_graph(
          fixture.gene, owners, palette, colors);
      std::string error;
      require(validate_component_graph(graph, &error),
              std::string(fixture.name) + " component graph invalid: " + error);
      require(graph.component_units() ==
                  static_cast<int>(fixture.gene.nodes.size()),
              std::string(fixture.name) + " component unit count mismatch");
      require(static_cast<int>(graph.edges.size()) + 1 ==
                  graph.component_units(),
              std::string(fixture.name) + " component graph not a tree");
      require(static_cast<int>(std::count_if(
                  graph.units.begin(), graph.units.end(),
                  [](const HdtComponentUnit& unit) {
                    return unit.kind == HdtComponentKind::LiveLeaf;
                  })) == fixture.leaves,
              std::string(fixture.name) + " owner range lost leaves");
      require(graph.root_unit == fixture.gene.root,
              std::string(fixture.name) + " root key mismatch");

      const HdtComponentGraph second = build_base_component_graph(
          fixture.gene, owners, palette, colors);
      require(second.edges.size() == graph.edges.size(),
              std::string(fixture.name) + " nondeterministic edge count");
      for (int edge = 0; edge < static_cast<int>(graph.edges.size()); ++edge) {
        const HdtIncidenceEdge& a =
            graph.edges[static_cast<std::size_t>(edge)];
        const HdtIncidenceEdge& b =
            second.edges[static_cast<std::size_t>(edge)];
        require(a.a == b.a && a.b == b.b && a.a_port == b.a_port &&
                    a.b_port == b.b_port,
                std::string(fixture.name) + " unstable source keys");
      }
    }
  }
}

void test_hdt_construction() {
  using namespace qraft::internal;
  for (int palette : {2, 4}) {
    HdtComponentGraph empty;
    empty.palette = palette;
    const HdtState state = build_hdt_state(std::move(empty));
    std::string error;
    require(state.hdt.root == -1 && state.hdt.nodes.empty() &&
                state.hdt.stats.node_count_ok &&
                state.hdt.stats.atom_bound_ok &&
                validate_hdt(state.components, state.hdt, &error),
            "empty HDT did not satisfy its structural base case: " +
                error);
  }

  // Extracting one arm of a degree-three root must rebuild a valid binary HDT,
  // not leave a high-degree recompute node.
  for (int palette : {2, 4}) {
    const GeneIndex gene = three_leaf_star_gene_index();
    const OwnerProjection owners = identity_owner_map(3);
    HdtColorPolicy colors;
    colors.current_by_owner =
        palette == 2 ? std::vector<std::uint8_t>{1, 0, 1}
                     : std::vector<std::uint8_t>{1, 2, 3};
    colors.background_by_owner.assign(3, 0);
    const HdtState base = build_hdt_state(
        build_base_component_graph(gene, owners, palette, colors));
    const HdtExtractContractResult extracted = hdt_extract_contract(
        base, hdt_owner_selector(3, {0}));
    std::string error;
    require(extracted.contract.size_bound_ok &&
                extracted.state.hdt.stats.max_recompute_degree == 2 &&
                extracted.state.hdt.stats.local_balance_violations == 0 &&
                validate_hdt(extracted.state.components, extracted.state.hdt,
                             &error),
            "resolved-root ExtractContract did not rebuild a valid binary "
            "HDT: " +
                error);
    require_quartet_counts(hdt_root_counts(extracted.state),
                           {0.0, 0.0, 0.0}, palette,
                           "resolved-root three-leaf base case");
  }

  struct Fixture {
    const char* name;
    GeneIndex gene;
    int leaves;
  };
  std::vector<Fixture> fixtures;
  fixtures.push_back({"star", star_gene_index(18), 18});
  fixtures.push_back({"caterpillar", caterpillar_gene_index(24), 24});
  fixtures.push_back({"broom", broom_gene_index(14), 28});
  fixtures.push_back({"random", random_arbitrary_gene_index(24, 91), 24});

  auto require_counts = [](const std::array<double, 3>& actual,
                           const std::array<double, 3>& expected,
                           int channels,
                           const std::string& label) {
    for (int channel = 0; channel < channels; ++channel) {
      const double a = actual[static_cast<std::size_t>(channel)];
      const double e = expected[static_cast<std::size_t>(channel)];
      if (!finite_close(a, e, 1e-8, 1e-8)) {
        fail(label + " channel=" + std::to_string(channel) +
             " actual=" + std::to_string(a) +
             " expected=" + std::to_string(e));
      }
    }
  };

  for (const Fixture& fixture : fixtures) {
    const OwnerProjection owners = identity_owner_map(fixture.leaves);
    for (int palette : {2, 4}) {
      HdtColorPolicy colors;
      colors.current_by_owner.resize(static_cast<std::size_t>(fixture.leaves));
      colors.background_by_owner.resize(
          static_cast<std::size_t>(fixture.leaves), 0);
      for (int owner = 0; owner < fixture.leaves; ++owner) {
        colors.current_by_owner[static_cast<std::size_t>(owner)] =
            static_cast<std::uint8_t>(owner % palette);
        if (palette == 4 && owner % 5 == 0) {
          colors.background_by_owner[static_cast<std::size_t>(owner)] = 3;
        }
      }
      HdtState state = build_hdt_state(
          build_base_component_graph(fixture.gene, owners, palette, colors));
      require(state.hdt.stats.node_count_ok,
              std::string(fixture.name) + " H != 2A-1");
      require(state.hdt.stats.atom_bound_ok,
              std::string(fixture.name) + " A > 2R");
      require(state.hdt.stats.max_recompute_degree == 2,
              std::string(fixture.name) + " HDT is not binary");
      require(state.hdt.stats.local_balance_violations == 0,
              std::string(fixture.name) + " HDT not locally balanced");

      require_counts(
          hdt_root_counts(state),
          independent_gene_quartet_counts(
              fixture.gene, owners, palette, colors.current_by_owner),
          palette == 2 ? 2 : 3,
          std::string(fixture.name) + " fixed-vector base");

      std::vector<int> changed_owners;
      for (int owner = 1; owner < fixture.leaves; owner += 3) {
        changed_owners.push_back(owner);
      }
      hdt_recolor_owners(state, changed_owners,
                               static_cast<std::uint8_t>(palette - 1));
      HdtColorPolicy recolored = colors;
      for (int owner : changed_owners) {
        recolored.current_by_owner[static_cast<std::size_t>(owner)] =
            static_cast<std::uint8_t>(palette - 1);
      }
      const HdtState rebuilt = build_hdt_state(
          build_base_component_graph(fixture.gene, owners, palette, recolored));
      require_counts(hdt_root_counts(state), hdt_root_counts(rebuilt),
                     palette == 2 ? 2 : 3,
                     std::string(fixture.name) + " recolor");
    }
  }

  // A child keeps a restricted gene set only after a fixed-factor shrink, so a
  // scoring gene can contain taxa outside the current owner map. Those taxa
  // are absent from the objective, not background-colored leaves.
  const GeneIndex extra_taxa_gene = random_arbitrary_gene_index(24, 177);
  OwnerProjection subset_owners;
  subset_owners.taxon_to_owner.assign(24, -1);
  subset_owners.live_owner.assign(12, 1);
  subset_owners.owner_count = 12;
  for (int owner = 0; owner < subset_owners.owner_count; ++owner) {
    const int taxon = 2 * owner;
    subset_owners.taxon_to_owner[static_cast<std::size_t>(taxon)] = owner;
    subset_owners.live_taxa.push_back(taxon);
  }
  for (int palette : {2, 4}) {
    HdtColorPolicy colors;
    colors.default_current = static_cast<std::uint8_t>(palette - 1);
    colors.default_background = static_cast<std::uint8_t>(palette - 1);
    colors.current_by_owner.resize(12);
    colors.background_by_owner.assign(12, 0);
    for (int owner = 0; owner < 12; ++owner) {
      colors.current_by_owner[static_cast<std::size_t>(owner)] =
          static_cast<std::uint8_t>(owner % palette);
    }
    HdtComponentGraph base = build_base_component_graph(
        extra_taxa_gene, subset_owners, palette, colors);
    const int base_units = base.component_units();
    const HdtState uncontracted = build_hdt_state(base);
    HdtInitialStateResult initial =
        hdt_build_initial_state(std::move(base));
    require(initial.contracted,
            "extra-taxa initial graph did not take compact Contract path");
    require(initial.contract.size_bound_ok &&
                initial.state.components.component_units() <= 4 * 12 + 3 &&
                initial.state.components.component_units() <= base_units,
            "extra-taxa initial Contract violated its size bound");
    require_counts(hdt_root_counts(initial.state),
                   hdt_root_counts(uncontracted),
                   palette == 2 ? 2 : 3,
                   "extra-taxa initial Contract parity");
    HdtInitialStateResult induced = hdt_build_initial_state(
        build_induced_base_component_graph(
            extra_taxa_gene, subset_owners, palette, colors));
    require(!induced.contracted,
            "direct source induction retained an absent source leaf");
    require(induced.state.components.component_units() < base_units,
            "direct source induction did not reduce its source topology");
    require_counts(hdt_root_counts(induced.state),
                   hdt_root_counts(uncontracted),
                   palette == 2 ? 2 : 3,
                   "extra-taxa direct-induction parity");
    const HdtState& state = initial.state;
    require_counts(
        hdt_root_counts(state),
        independent_gene_quartet_counts(
            extra_taxa_gene, subset_owners, palette,
            colors.current_by_owner),
        palette == 2 ? 2 : 3, "extra-taxa fixed-vector base");
  }

  // A scoring gene can be rooted at a taxon absent from the current owner map.
  // Its component is then a zero-mass C1Opaque root, not a structural Contract
  // anchor.
  GeneIndex opaque_root_gene = extra_taxa_gene;
  opaque_root_gene.root = opaque_root_gene.taxon_to_node[23];
  opaque_root_gene.root_traversal =
      traversal(opaque_root_gene, opaque_root_gene.root);
  for (int palette : {2, 4}) {
    HdtColorPolicy colors;
    colors.default_current = static_cast<std::uint8_t>(palette - 1);
    colors.default_background = static_cast<std::uint8_t>(palette - 1);
    colors.current_by_owner.resize(12);
    colors.background_by_owner.assign(12, 0);
    for (int owner = 0; owner < 12; ++owner) {
      colors.current_by_owner[static_cast<std::size_t>(owner)] =
          static_cast<std::uint8_t>(owner % palette);
    }

    HdtComponentGraph base = build_base_component_graph(
        opaque_root_gene, subset_owners, palette, colors);
    const HdtComponentUnit& source_root =
        base.units[static_cast<std::size_t>(base.root_unit)];
    require(source_root.kind == HdtComponentKind::C1Opaque &&
                source_root.weight == 0.0 && source_root.incidence_count == 1,
            "opaque-root fixture did not construct a zero-mass opaque root");
    HdtInitialStateResult initial =
        hdt_build_initial_state(std::move(base));
    require(initial.contracted && initial.contract.size_bound_ok,
            "opaque-root initial Contract failed its size bound");
    const HdtComponentKind contracted_root =
        initial.state.components.units[static_cast<std::size_t>(
            initial.state.components.root_unit)]
            .kind;
    require(contracted_root == HdtComponentKind::IAnchor ||
                contracted_root == HdtComponentKind::LiveLeaf,
            "opaque-root Contract did not emit a structural root");
    require_quartet_counts(
        hdt_root_counts(initial.state),
        independent_gene_quartet_counts(extra_taxa_gene, subset_owners,
                                         palette, colors.current_by_owner),
        palette, "opaque-root initial Contract parity");
    HdtInitialStateResult induced = hdt_build_initial_state(
        build_induced_base_component_graph(
            opaque_root_gene, subset_owners, palette, colors));
    require(!induced.contracted,
            "opaque-root direct induction retained its absent root");
    require_quartet_counts(
        hdt_root_counts(induced.state), hdt_root_counts(initial.state),
        palette, "opaque-root direct-induction parity");

    std::vector<int> changed{1, 4, 7, 10};
    hdt_recolor_owners(initial.state, changed,
                             static_cast<std::uint8_t>(palette - 1));
    HdtColorPolicy recolored = colors;
    for (int owner : changed) {
      recolored.current_by_owner[static_cast<std::size_t>(owner)] =
          static_cast<std::uint8_t>(palette - 1);
    }
    require_quartet_counts(
        hdt_root_counts(initial.state),
        independent_gene_quartet_counts(extra_taxa_gene, subset_owners,
                                         palette,
                                         recolored.current_by_owner),
        palette, "opaque-root recolor parity");

    const std::vector<int> selected{0, 2, 4, 6, 8, 10};
    const HdtLiveSelector selector =
        hdt_owner_selector(12, selected);
    const HdtExtractContractResult extracted =
        hdt_extract_contract(initial.state, selector);
    for (int owner = 0; owner < 12; ++owner) {
      if (!selector.contains_owner(owner)) {
        recolored.current_by_owner[static_cast<std::size_t>(owner)] =
            recolored.background_by_owner[static_cast<std::size_t>(owner)];
      }
    }
    require(extracted.contract.size_bound_ok,
            "opaque-root ExtractContract exceeded its size bound");
    require_quartet_counts(
        hdt_root_counts(extracted.state),
        independent_gene_quartet_counts(extra_taxa_gene, subset_owners,
                                         palette,
                                         recolored.current_by_owner),
        palette, "opaque-root ExtractContract parity");
  }

  // A two-leaf gene has no I anchor. Contract must fall through from an
  // absent opaque root to the selected live leaf without inventing a topology
  // contribution. An HDT of the two leaves themselves would join two C1
  // values at their top ports, which no composition does, so building one
  // fails closed; the scorers skip genes with fewer than four live leaves.
  GeneIndex two_leaf_gene;
  two_leaf_gene.nodes.resize(2);
  two_leaf_gene.nodes[0].set_taxon(0);
  two_leaf_gene.nodes[1].set_taxon(1);
  set_gene_adjacency(two_leaf_gene, {{1}, {0}});
  two_leaf_gene.root = 1;
  two_leaf_gene.reference_root = 0;
  two_leaf_gene.taxon_to_node = {0, 1};
  two_leaf_gene.root_traversal = traversal(two_leaf_gene, two_leaf_gene.root);
  OwnerProjection one_owner;
  one_owner.taxon_to_owner = {0, -1};
  one_owner.live_owner = {1};
  one_owner.live_taxa = {0};
  one_owner.owner_count = 1;
  for (int palette : {2, 4}) {
    HdtColorPolicy colors;
    colors.current_by_owner = {static_cast<std::uint8_t>(palette - 1)};
    colors.background_by_owner = {0};
    HdtComponentGraph base = build_base_component_graph(
        two_leaf_gene, one_owner, palette, colors);
    bool rejected = false;
    try {
      static_cast<void>(build_hdt_state(base));
    } catch (const std::logic_error&) {
      rejected = true;
    }
    require(rejected, "an HDT of two joined leaves was built");
    const HdtInitialStateResult initial =
        hdt_build_initial_state(std::move(base));
    require(initial.contracted && initial.contract.size_bound_ok,
            "two-leaf opaque-root Contract failed its size bound");
    require(initial.state.components.units[static_cast<std::size_t>(
                initial.state.components.root_unit)]
                .kind == HdtComponentKind::LiveLeaf,
            "two-leaf opaque-root Contract did not select its live root");
    require_counts(hdt_root_counts(initial.state), {0.0, 0.0, 0.0},
                   palette == 2 ? 2 : 3,
                   "two-leaf opaque-root Contract counted a quartet");
  }

  const OwnerProjection two_owners = identity_owner_map(2);
  for (int palette : {2, 4}) {
    HdtColorPolicy colors;
    colors.current_by_owner = {0, static_cast<std::uint8_t>(palette - 1)};
    colors.background_by_owner = {0, 0};
    bool rejected = false;
    try {
      static_cast<void>(build_hdt_state(build_base_component_graph(
          two_leaf_gene, two_owners, palette, colors)));
    } catch (const std::logic_error&) {
      rejected = true;
    }
    require(rejected, "an HDT of two joined live leaves was built");
  }
}

void test_hdt_contract() {
  using namespace qraft::internal;
  struct Fixture {
    const char* name;
    GeneIndex gene;
    int leaves;
  };
  std::vector<Fixture> fixtures;
  fixtures.push_back({"star", star_gene_index(64), 64});
  fixtures.push_back({"caterpillar", caterpillar_gene_index(96), 96});
  fixtures.push_back({"broom", broom_gene_index(48), 96});
  fixtures.push_back({"random", random_arbitrary_gene_index(80, 119), 80});

  auto spread = [](int leaves, int count) {
    std::vector<int> owners;
    for (int i = 0; i < count; ++i) {
      owners.push_back(static_cast<int>(
          static_cast<long long>(i) * leaves / std::max(1, count)));
    }
    std::sort(owners.begin(), owners.end());
    owners.erase(std::unique(owners.begin(), owners.end()), owners.end());
    return owners;
  };
  auto require_counts = [](const std::array<double, 3>& actual,
                           const std::array<double, 3>& expected,
                           int channels,
                           const std::string& label) {
    for (int channel = 0; channel < channels; ++channel) {
      const double a = actual[static_cast<std::size_t>(channel)];
      const double e = expected[static_cast<std::size_t>(channel)];
      if (!finite_close(a, e, 1e-8, 1e-8)) {
        fail(label + " channel=" + std::to_string(channel) +
             " actual=" + std::to_string(a) +
             " expected=" + std::to_string(e));
      }
    }
  };

  for (const Fixture& fixture : fixtures) {
    const OwnerProjection owners = identity_owner_map(fixture.leaves);
    for (int palette : {2, 4}) {
      HdtColorPolicy base_colors;
      base_colors.current_by_owner.resize(
          static_cast<std::size_t>(fixture.leaves));
      base_colors.background_by_owner.resize(
          static_cast<std::size_t>(fixture.leaves));
      for (int owner = 0; owner < fixture.leaves; ++owner) {
        base_colors.current_by_owner[static_cast<std::size_t>(owner)] =
            static_cast<std::uint8_t>((owner * 3 + 1) % palette);
        base_colors.background_by_owner[static_cast<std::size_t>(owner)] =
            static_cast<std::uint8_t>(palette == 4 && owner % 9 == 0 ? 3 : 0);
      }
      const HdtState base = build_hdt_state(
          build_base_component_graph(fixture.gene, owners, palette,
                                     base_colors));
      for (int wanted : {0, 1, 2, 4, 8}) {
        const std::vector<int> selected = spread(fixture.leaves, wanted);
        const HdtLiveSelector selector =
            hdt_owner_selector(fixture.leaves, selected);
        const HdtContractResult contracted = hdt_contract(
            base, selector, HdtContractMode::LightCopy);
        require(contracted.stats.size_bound_ok,
                std::string(fixture.name) + " Contract size bound failed");
        require(contracted.stats.output_units <=
                    4 * std::max(1, contracted.stats.selected_live_atoms) + 3,
                std::string(fixture.name) + " Contract R bound failed");
        require(contracted.state.hdt.stats.node_count_ok &&
                    contracted.state.hdt.stats.atom_bound_ok &&
                    contracted.state.hdt.stats.max_recompute_degree <= 2,
                std::string(fixture.name) + " Contract rebuild invalid");
        require(
            contracted.stats.scratch_peak_bytes ==
                contracted.stats.persistent_scratch_bytes +
                    contracted.stats.ephemeral_scratch_peak_bytes &&
                contracted.stats.persistent_scratch_bytes > 0 &&
                contracted.stats.ephemeral_scratch_peak_bytes > 0,
            std::string(fixture.name) +
                " Contract scratch high-water composition is incomplete");

        HdtColorPolicy expected_colors = base_colors;
        for (int owner = 0; owner < fixture.leaves; ++owner) {
          if (selector[static_cast<std::size_t>(owner)] == 0) {
            expected_colors.current_by_owner[static_cast<std::size_t>(owner)] =
                expected_colors.background_by_owner[static_cast<std::size_t>(owner)];
          }
        }
        const HdtState expected = build_hdt_state(
            build_base_component_graph(fixture.gene, owners, palette,
                                       expected_colors));
        require_counts(hdt_root_counts(contracted.state),
                       hdt_root_counts(expected), palette == 2 ? 2 : 3,
                       std::string(fixture.name) + " Contract parity");

        if (selected.size() >= 2u) {
          std::vector<int> nested;
          for (std::size_t i = 0; i < selected.size(); i += 2u) {
            nested.push_back(selected[i]);
          }
          const HdtLiveSelector nested_selector =
              hdt_owner_selector(fixture.leaves, nested);
          const HdtContractResult second = hdt_contract(
              contracted.state, nested_selector,
              HdtContractMode::HeavyReplace);
          HdtColorPolicy nested_colors = base_colors;
          for (int owner = 0; owner < fixture.leaves; ++owner) {
            if (nested_selector[static_cast<std::size_t>(owner)] == 0) {
              nested_colors.current_by_owner[static_cast<std::size_t>(owner)] =
                  nested_colors.background_by_owner[static_cast<std::size_t>(owner)];
            }
          }
          const HdtState nested_expected = build_hdt_state(
              build_base_component_graph(fixture.gene, owners, palette,
                                         nested_colors));
          require_counts(hdt_root_counts(second.state),
                         hdt_root_counts(nested_expected),
                         palette == 2 ? 2 : 3,
                         std::string(fixture.name) + " nested Contract parity");
        }
      }
    }
  }
}

void test_hdt_extract() {
  using namespace qraft::internal;
  struct Fixture {
    const char* name;
    GeneIndex gene;
    int leaves;
  };
  std::vector<Fixture> fixtures;
  fixtures.push_back({"star", star_gene_index(72), 72});
  fixtures.push_back({"caterpillar", caterpillar_gene_index(128), 128});
  fixtures.push_back({"broom", broom_gene_index(64), 128});
  fixtures.push_back({"random", random_arbitrary_gene_index(96, 211), 96});

  auto spread = [](int leaves, int count, int offset) {
    std::vector<int> owners;
    for (int i = 0; i < count; ++i) {
      owners.push_back((offset + static_cast<int>(
          static_cast<long long>(i) * leaves / std::max(1, count))) % leaves);
    }
    std::sort(owners.begin(), owners.end());
    owners.erase(std::unique(owners.begin(), owners.end()), owners.end());
    return owners;
  };
  auto require_counts = [](const std::array<double, 3>& actual,
                           const std::array<double, 3>& expected,
                           int channels,
                           const std::string& label) {
    for (int channel = 0; channel < channels; ++channel) {
      const double a = actual[static_cast<std::size_t>(channel)];
      const double e = expected[static_cast<std::size_t>(channel)];
      if (!finite_close(a, e, 1e-8, 1e-8)) {
        fail(label + " channel=" + std::to_string(channel) +
             " actual=" + std::to_string(a) +
             " expected=" + std::to_string(e));
      }
    }
  };

  for (const Fixture& fixture : fixtures) {
    const OwnerProjection owners = identity_owner_map(fixture.leaves);
    for (int palette : {2, 4}) {
      HdtColorPolicy colors;
      colors.current_by_owner.resize(static_cast<std::size_t>(fixture.leaves));
      colors.background_by_owner.resize(
          static_cast<std::size_t>(fixture.leaves));
      for (int owner = 0; owner < fixture.leaves; ++owner) {
        colors.current_by_owner[static_cast<std::size_t>(owner)] =
            static_cast<std::uint8_t>((owner + 1) % palette);
        colors.background_by_owner[static_cast<std::size_t>(owner)] =
            static_cast<std::uint8_t>(palette == 4 && owner % 11 == 0 ? 3 : 0);
      }
      const HdtState base = build_hdt_state(
          build_base_component_graph(fixture.gene, owners, palette, colors));
      for (int wanted : {1, 2, 4, 8}) {
        const std::vector<int> selected = spread(fixture.leaves, wanted, 1);
        const HdtLiveSelector selector =
            hdt_owner_selector(fixture.leaves, selected);
        const std::array<double, 3> parent_before =
            hdt_root_counts(base);
        const HdtExtractContractResult extracted =
            hdt_extract_contract(base, selector);
        require_counts(hdt_root_counts(base), parent_before,
                       palette == 2 ? 2 : 3,
                       std::string(fixture.name) + " Extract mutated parent");
        require(extracted.extract.marked_hdt_nodes <=
                    extracted.extract.source_hdt_nodes,
                std::string(fixture.name) + " invalid marked-node count");
        require(extracted.contract.size_bound_ok &&
                    extracted.state.components.component_units() <=
                        4 * std::max(1, extracted.extract.selected_live_atoms) + 3,
                std::string(fixture.name) + " ExtractContract final bound");

        HdtColorPolicy expected_colors = colors;
        for (int owner = 0; owner < fixture.leaves; ++owner) {
          if (selector[static_cast<std::size_t>(owner)] == 0) {
            expected_colors.current_by_owner[static_cast<std::size_t>(owner)] =
                expected_colors.background_by_owner[static_cast<std::size_t>(owner)];
          }
        }
        const HdtState expected = build_hdt_state(
            build_base_component_graph(fixture.gene, owners, palette,
                                       expected_colors));
        require_counts(hdt_root_counts(extracted.state),
                       hdt_root_counts(expected), palette == 2 ? 2 : 3,
                       std::string(fixture.name) + " ExtractContract parity");

        const std::vector<int> parent_selected =
            spread(fixture.leaves, std::min(16, fixture.leaves), 0);
        const HdtContractResult restricted_parent = hdt_contract(
            base, hdt_owner_selector(fixture.leaves, parent_selected),
            HdtContractMode::HeavyReplace);
        std::vector<int> nested;
        for (std::size_t i = 0; i < parent_selected.size(); i += 2u) {
          nested.push_back(parent_selected[i]);
        }
        const HdtLiveSelector nested_selector =
            hdt_owner_selector(fixture.leaves, nested);
        const HdtExtractContractResult nested_result =
            hdt_extract_contract(restricted_parent.state, nested_selector);
        HdtColorPolicy nested_colors = colors;
        for (int owner = 0; owner < fixture.leaves; ++owner) {
          if (nested_selector[static_cast<std::size_t>(owner)] == 0) {
            nested_colors.current_by_owner[static_cast<std::size_t>(owner)] =
                nested_colors.background_by_owner[static_cast<std::size_t>(owner)];
          }
        }
        const HdtState nested_expected = build_hdt_state(
            build_base_component_graph(fixture.gene, owners, palette,
                                       nested_colors));
        require_counts(hdt_root_counts(nested_result.state),
                       hdt_root_counts(nested_expected),
                       palette == 2 ? 2 : 3,
                       std::string(fixture.name) + " nested Extract parity");
      }
    }
  }
}

void test_count_traversal_oracle() {
  using namespace qraft::internal;
  constexpr int kLeaves = 12;
  const OwnerProjection owners = identity_owner_map(kLeaves);
  const ActiveTaxa all_active = [] {
    std::vector<int> taxa(kLeaves);
    std::iota(taxa.begin(), taxa.end(), 0);
    return active_set(std::move(taxa), kLeaves);
  }();
  const QueryTree active_query = induce_query_tree(
      balanced_starting_tree(kLeaves), all_active);

  std::vector<GeneIndex> genes;
  genes.push_back(star_gene_index(kLeaves));
  genes.push_back(caterpillar_gene_index(kLeaves));
  genes.push_back(broom_gene_index(kLeaves / 2));
  genes.push_back(random_arbitrary_gene_index(kLeaves, 917));

  auto close = [](double actual, double expected, const std::string& label) {
    if (!finite_close(actual, expected, 1e-8, 1e-8)) {
      fail(label + " actual=" + std::to_string(actual) +
           " expected=" + std::to_string(expected));
    }
  };

  for (int fixture = 0; fixture < static_cast<int>(genes.size()); ++fixture) {
    const GeneIndex& gene = genes[static_cast<std::size_t>(fixture)];
    HdtColorPolicy colors;
    colors.default_current = 0;
    colors.default_background = 0;
    colors.current_by_owner.assign(kLeaves, 1);
    colors.background_by_owner.assign(kLeaves, 0);
    HdtState state = build_hdt_state(
        build_base_component_graph(gene, owners, 2, colors));
    HdtExtractWorkspace workspace;
    HdtTraversalStats stats;
    int visits = 0;
    count_bipartitions(
        std::move(state), active_query, active_query.leaf_active_order, {},
        workspace, stats,
        [&](int query_node, HdtState& current) {
          std::vector<std::uint8_t> owner_colors(kLeaves, 0);
          const int begin = active_query.subtree_leaf_first[
              static_cast<std::size_t>(query_node)];
          const int end = begin + active_query.subtree_leaf_count[
                                      static_cast<std::size_t>(query_node)];
          for (int position = begin; position < end; ++position) {
            owner_colors[static_cast<std::size_t>(
                active_query.leaf_active_order[
                    static_cast<std::size_t>(position)])] = 1;
          }
          const std::array<double, 3> expected =
              independent_gene_quartet_counts(gene, owners, 2, owner_colors);
          const std::array<double, 3> actual = hdt_root_counts(current);
          close(actual[0], expected[0], "active support fixture=" +
                                                std::to_string(fixture));
          close(actual[1], expected[1], "active violation fixture=" +
                                                std::to_string(fixture));
          ++visits;
        });
    require(visits > 0 && stats.complexity_contract_satisfied,
            "active oracle traversal did not complete");
  }

  std::vector<int> dynamic_taxa(10);
  std::iota(dynamic_taxa.begin(), dynamic_taxa.end(), 0);
  const ActiveTaxa dynamic = active_set(dynamic_taxa, kLeaves);
  const QueryTree query = induce_query_tree(
      balanced_starting_tree(kLeaves), dynamic);
  std::vector<int> active_to_owner(10);
  std::iota(active_to_owner.begin(), active_to_owner.end(), 0);
  const std::vector<int> owner_order =
      query_leaf_owner_order(query, active_to_owner);
  const std::vector<int> context{10, 11};

  for (int fixture = 0; fixture < static_cast<int>(genes.size()); ++fixture) {
    const GeneIndex& gene = genes[static_cast<std::size_t>(fixture)];
    HdtColorPolicy colors;
    colors.default_current = 0;
    colors.default_background = 0;
    colors.current_by_owner.assign(kLeaves, 0);
    colors.background_by_owner.assign(kLeaves, 0);
    for (int owner : owner_order) {
      colors.current_by_owner[static_cast<std::size_t>(owner)] = 1;
    }
    for (int owner : context) {
      colors.current_by_owner[static_cast<std::size_t>(owner)] = 3;
      colors.background_by_owner[static_cast<std::size_t>(owner)] = 3;
    }
    HdtState state = build_hdt_state(
        build_base_component_graph(gene, owners, 4, colors));
    auto oracle_counts = [&](const std::vector<int>& owner_colors) {
      std::vector<std::uint8_t> compact(owner_colors.begin(),
                                        owner_colors.end());
      return independent_gene_quartet_counts(gene, owners, 4, compact);
    };

    HdtExtractWorkspace workspace;
    HdtTraversalStats stats;
    int visits = 0;
    count_vertices(
        std::move(state), query, owner_order, workspace, stats,
        [&](int vertex, const std::array<int, 3>& branch_nodes,
            const std::array<std::uint8_t, 3>& branch_colors,
            HdtState& current) {
          std::vector<int> owner_colors(kLeaves, 0);
          owner_colors[10] = 3;
          owner_colors[11] = 3;
          std::vector<int> branch_of_position(owner_order.size(), -1);
          const int vertex_begin = query.subtree_leaf_first[
              static_cast<std::size_t>(vertex)];
          const int vertex_end = vertex_begin + query.subtree_leaf_count[
                                                  static_cast<std::size_t>(vertex)];
          for (int branch = 0; branch < 3; ++branch) {
            const int branch_node =
                branch_nodes[static_cast<std::size_t>(branch)];
            if (branch_node == query.parent[static_cast<std::size_t>(vertex)]) {
              for (int position = 0;
                   position < static_cast<int>(owner_order.size()); ++position) {
                if (position < vertex_begin || position >= vertex_end) {
                  branch_of_position[static_cast<std::size_t>(position)] =
                      branch;
                }
              }
              continue;
            }
            const int begin = query.subtree_leaf_first[
                static_cast<std::size_t>(branch_node)];
            const int end = begin + query.subtree_leaf_count[
                                        static_cast<std::size_t>(branch_node)];
            for (int position = begin; position < end; ++position) {
              branch_of_position[static_cast<std::size_t>(position)] = branch;
            }
          }
          for (int position = 0;
               position < static_cast<int>(owner_order.size()); ++position) {
            const int branch =
                branch_of_position[static_cast<std::size_t>(position)];
            require(branch >= 0, "vertex branch partition missed a query leaf");
            owner_colors[static_cast<std::size_t>(
                owner_order[static_cast<std::size_t>(position)])] =
                branch_colors[static_cast<std::size_t>(branch)];
          }
          const std::array<double, 3> expected = oracle_counts(owner_colors);
          const std::array<double, 3> actual = hdt_root_counts(current);
          std::string regraft_case =
              " fixture=" + std::to_string(fixture) +
              " vertex=" + std::to_string(vertex);
          regraft_case += " index=" +
                        std::to_string(current.hdt.owner_position_begin) +
                        ":" +
                        std::to_string(current.hdt.owner_offsets.size()) +
                        " expected=";
          for (int position = 0;
               position < static_cast<int>(owner_order.size()); ++position) {
            regraft_case += std::to_string(owner_colors[static_cast<std::size_t>(
                owner_order[static_cast<std::size_t>(position)])]);
          }
          regraft_case += " actual=";
          for (int position = 0;
               position < static_cast<int>(owner_order.size()); ++position) {
            const int owner = owner_order[static_cast<std::size_t>(position)];
            const int slot = current.hdt.owner_index_by_position
                                 ? position -
                                       current.hdt.owner_position_begin
                                 : owner;
            int actual_color = -1;
            if (slot >= 0 && slot + 1 <
                                 static_cast<int>(
                                     current.hdt.owner_offsets.size())) {
              const int begin = current.hdt.owner_offsets[
                  static_cast<std::size_t>(slot)];
              const int end = current.hdt.owner_offsets[
                  static_cast<std::size_t>(slot) + 1u];
              if (begin < end) {
                const HdtNode& atom_node = current.hdt.nodes[
                    static_cast<std::size_t>(current.hdt.owner_atom_nodes[
                        static_cast<std::size_t>(begin)])];
                actual_color = current.hdt.atoms[static_cast<std::size_t>(
                    atom_node.atom_descriptor)].current_color;
              }
            }
            regraft_case += std::to_string(actual_color);
          }
          close(actual[0], expected[0], "regraft q01" + regraft_case);
          close(actual[1], expected[1], "regraft q02" + regraft_case);
          close(actual[2], expected[2], "regraft q03" + regraft_case);

          for (int branch = 0; branch < 3; ++branch) {
            const int branch_node =
                branch_nodes[static_cast<std::size_t>(branch)];
            if (branch_node < 0 ||
                branch_node == query.parent[static_cast<std::size_t>(vertex)]) {
              continue;
            }
            const int begin = query.subtree_leaf_first[
                static_cast<std::size_t>(branch_node)];
            const int end = begin + query.subtree_leaf_count[
                                        static_cast<std::size_t>(branch_node)];
            hdt_recolor_owner_range(current, owner_order, begin, end, 3);
            std::vector<int> nni_colors = owner_colors;
            for (int position = begin; position < end; ++position) {
              nni_colors[static_cast<std::size_t>(
                  owner_order[static_cast<std::size_t>(position)])] = 3;
            }
            const std::array<double, 3> nni_expected =
                oracle_counts(nni_colors);
            const std::array<double, 3> nni_actual =
                hdt_root_counts(current);
            close(nni_actual[0], nni_expected[0],
                  "NNI q01 fixture=" + std::to_string(fixture));
            close(nni_actual[1], nni_expected[1],
                  "NNI q02 fixture=" + std::to_string(fixture));
            close(nni_actual[2], nni_expected[2],
                  "NNI q03 fixture=" + std::to_string(fixture));
            hdt_recolor_owner_range(
                current, owner_order, begin, end,
                branch_colors[static_cast<std::size_t>(branch)]);
            break;
          }
          ++visits;
        });
    require(visits > 0 && stats.complexity_contract_satisfied,
            "regraft/NNI oracle traversal did not complete");
  }
}

void test_count_traversal_scaling() {
  using namespace qraft::internal;
  for (int leaves : {32, 64, 128, 256}) {
    const OwnerProjection owners = identity_owner_map(leaves);
    const GeneIndex gene = random_arbitrary_gene_index(
        leaves, static_cast<unsigned>(3100 + leaves));
    HdtColorPolicy colors;
    colors.default_current = 0;
    colors.default_background = 0;
    colors.current_by_owner.assign(static_cast<std::size_t>(leaves), 1);
    colors.background_by_owner.assign(static_cast<std::size_t>(leaves), 0);
    const HdtState base = build_hdt_state(
        build_base_component_graph(gene, owners, 4, colors));
    const std::vector<int> selected{leaves / 3, (2 * leaves) / 3};
    const HdtLiveSelector selector =
        hdt_owner_selector(leaves, selected);
    HdtExtractWorkspace workspace;
    const HdtExtractContractResult extracted =
        hdt_extract_contract(base, selector, workspace);
    const double path_envelope =
        static_cast<double>(selected.size()) *
        (1.0 + std::log2(static_cast<double>(
                   base.hdt.stats.initial_atoms) / selected.size()));
    require_finite(path_envelope, "fixed-live Extract path envelope");
    require(extracted.extract.marked_hdt_nodes <=
                    static_cast<int>(64.0 * path_envelope + 16.0) &&
                extracted.extract.intermediate_units <=
                    static_cast<int>(128.0 * path_envelope + 32.0),
            "fixed-live Extract path union exceeded s(1+log(N/s)) envelope");
    require(extracted.state.components.component_units() <= 11 &&
                extracted.contract.size_bound_ok,
            "fixed-live ExtractContract output grew with source size");

    const HdtContractResult contracted = hdt_contract(
        base, selector, HdtContractMode::LightCopy);
    const long long contract_work =
        contracted.stats.input_units + contracted.stats.fragment_hdt_nodes +
        contracted.stats.cc_ops + contracted.stats.gg_ops +
        contracted.stats.ig_ops + contracted.stats.c_to_g_ops;
    require(contract_work <= 64LL * base.components.component_units() + 64,
            "Contract operation count is not linear in input units");
  }

  std::vector<double> normalized;
  for (int leaves : {16, 32, 64, 128}) {
    const OwnerProjection owners = identity_owner_map(leaves);
    const GeneIndex gene = random_arbitrary_gene_index(
        leaves, static_cast<unsigned>(700 + leaves));
    std::vector<int> taxa(static_cast<std::size_t>(leaves));
    std::iota(taxa.begin(), taxa.end(), 0);
    const ActiveTaxa active = active_set(std::move(taxa), leaves);
    const QueryTree query = induce_query_tree(
        balanced_starting_tree(leaves), active);
    HdtColorPolicy colors;
    colors.default_current = 0;
    colors.default_background = 0;
    colors.current_by_owner.assign(static_cast<std::size_t>(leaves), 1);
    colors.background_by_owner.assign(static_cast<std::size_t>(leaves), 0);
    HdtState state = build_hdt_state(
        build_base_component_graph(gene, owners, 2, colors));
    HdtExtractWorkspace workspace;
    HdtTraversalStats stats;
    count_bipartitions(
        std::move(state), query, query.leaf_active_order, {}, workspace, stats,
        [](int, HdtState& current) {
          (void)hdt_root_counts(current);
        });
    const long long structural_work =
        stats.extract_marked_hdt_nodes + stats.extract_intermediate_units +
        stats.contract_input_units + stats.hdt_rebuild_nodes +
        stats.recolor_marked_nodes;
    const double envelope = static_cast<double>(leaves) *
                            (1.0 + std::log2(static_cast<double>(leaves)));
    require_finite(envelope, "balanced scaling envelope");
    const double normalized_work =
        static_cast<double>(structural_work) / envelope;
    require_finite(normalized_work, "balanced normalized work");
    normalized.push_back(normalized_work);
    require(stats.complexity_contract_satisfied &&
                stats.contract_bound_violations == 0 &&
                stats.heavy_shrink_violations == 0,
            "scaling structural invariant failed");
    require(structural_work <= static_cast<long long>(5000.0 * envelope),
            "traversal exceeded n log n work envelope");
    require(stats.pending_state_units_peak <= 32LL * leaves + 64 &&
                stats.pending_state_hdt_nodes_peak <= 128LL * leaves + 256,
            "traversal retained superlinear logical state");
    require(stats.extract_scratch_peak_bytes <=
                static_cast<std::size_t>(1024LL * leaves + 4096),
            "Extract workspace exceeded linear envelope");
  }
  require_finite(normalized.front(), "first normalized work");
  require_finite(normalized.back(), "last normalized work");
  require(normalized.back() <= 4.0 * normalized.front() + 1.0,
          "normalized work grows faster than n log n");

  constexpr int kSkewLeaves = 128;
  const OwnerProjection owners = identity_owner_map(kSkewLeaves);
  const GeneIndex gene = caterpillar_gene_index(kSkewLeaves);
  std::vector<int> taxa(kSkewLeaves);
  std::iota(taxa.begin(), taxa.end(), 0);
  const QueryTree query = induce_query_tree(
      caterpillar_starting_tree(kSkewLeaves),
      active_set(std::move(taxa), kSkewLeaves));
  HdtColorPolicy colors;
  colors.default_current = 0;
  colors.default_background = 0;
  colors.current_by_owner.assign(kSkewLeaves, 1);
  colors.background_by_owner.assign(kSkewLeaves, 0);
  HdtExtractWorkspace workspace;
  HdtTraversalStats stats;
  count_bipartitions(
      build_hdt_state(
          build_base_component_graph(gene, owners, 2, colors)),
      query, query.leaf_active_order, {}, workspace, stats,
      [](int, HdtState& current) {
        (void)hdt_root_counts(current);
      });
  const long long skew_work =
      stats.extract_marked_hdt_nodes + stats.extract_intermediate_units +
      stats.contract_input_units + stats.hdt_rebuild_nodes +
      stats.recolor_marked_nodes;
  const double skew_envelope = kSkewLeaves *
                               (1.0 + std::log2(kSkewLeaves));
  require_finite(skew_envelope, "skew scaling envelope");
  require(stats.complexity_contract_satisfied &&
              skew_work <= static_cast<long long>(5000.0 * skew_envelope),
          "skew traversal exceeded n log n work envelope");
  require(stats.extract_calls == stats.contract_calls &&
              stats.contract_calls ==
                  stats.contract_light_calls + stats.contract_heavy_calls &&
              stats.heavy_calls == stats.contract_heavy_calls &&
              stats.heavy_calls > 0,
          "counters do not distinguish fused light/heavy "
          "ExtractContract calls");
  require(stats.extract_source_hdt_nodes >=
                  stats.extract_marked_hdt_nodes &&
              stats.dp_transition_evals > 0 && stats.dp_scalar_fmas > 0 &&
              stats.hdt_max_recompute_degree == 2 &&
              stats.contract_bound_violations == 0 &&
              stats.heavy_shrink_violations == 0,
          "proof-counter deltas are inconsistent");
}

// ---------------------------------------------------------------------------
// Counting-scheme parity on random rooted trees: folding the generated counter
// transitions (leaf, C-to-G, GG, IG) bottom-up must reproduce brute-force
// weighted quartet counts under both a current and a background coloring (for
// K=4 the background mixes colors 0 and 3).
namespace counting_parity {

struct BTree {
  std::vector<std::vector<int>> ch;
  std::vector<int> color;   // -1 internal
  std::vector<double> weight;
  int root = 0;
  int add(int c, double w) {
    ch.push_back({}); color.push_back(c); weight.push_back(w);
    return static_cast<int>(ch.size()) - 1;
  }
};

BTree random_tree(std::mt19937& rng, int nlv, int K) {
  BTree t;
  int root = t.add(-1, 0.0);
  t.root = root;
  std::vector<int> internal{root};
  std::uniform_int_distribution<int> col(0, K - 1), w(1, 4), coin(0, 1);
  for (int i = 0; i < nlv; ++i) {
    int par = internal[std::uniform_int_distribution<std::size_t>(0, internal.size() - 1)(rng)];
    int lf = t.add(col(rng), 1.0 + w(rng));
    t.ch[par].push_back(lf);
    if (coin(rng) == 0) { int nv = t.add(-1, 0.0); t.ch[par].push_back(nv); internal.push_back(nv); }
  }
  for (int v = 0; v < static_cast<int>(t.ch.size()); ++v)
    if (v != root && t.color[v] < 0 && t.ch[v].empty()) {
      int lf = t.add(col(rng), 1.0); t.ch[v].push_back(lf);
    }
  return t;
}

int channel_idx(int u, int v, int x, int y, int K) {
  int a = std::min(u, v), b = std::max(u, v), c = std::min(x, y), d = std::max(x, y);
  if (a > c || (a == c && b > d)) { std::swap(a, c); std::swap(b, d); }
  if (K == 2) {
    if (a == 0 && b == 0 && c == 1 && d == 1) return 0;
    if (a == 0 && b == 1 && c == 0 && d == 1) return 1;
    return -1;
  }
  if (a == 0 && b == 1 && c == 2 && d == 3) return 0;
  if (a == 0 && b == 2 && c == 1 && d == 3) return 1;
  if (a == 0 && b == 3 && c == 1 && d == 2) return 2;
  return -1;
}

std::array<double, 3> brute(const BTree& t, const std::vector<int>& color, int K) {
  const int n = static_cast<int>(t.ch.size());
  std::vector<std::vector<int>> adj(n);
  for (int v = 0; v < n; ++v)
    for (int c : t.ch[v]) { adj[v].push_back(c); adj[c].push_back(v); }
  std::vector<int> leaves;
  for (int v = 0; v < n; ++v) if (t.color[v] >= 0) leaves.push_back(v);
  auto pathset = [&](int a, int b) {
    std::vector<int> par(n, -2); std::queue<int> q; q.push(a); par[a] = -1;
    while (!q.empty()) { int x = q.front(); q.pop(); if (x == b) break;
      for (int y : adj[x]) if (par[y] == -2) { par[y] = x; q.push(y); } }
    std::set<int> s; int x = b; while (x != -1) { s.insert(x); x = par[x]; } return s;
  };
  std::array<double, 3> out{0, 0, 0};
  const int L = static_cast<int>(leaves.size());
  for (int i = 0; i < L; ++i) for (int j = i + 1; j < L; ++j)
    for (int k = j + 1; k < L; ++k) for (int m = k + 1; m < L; ++m) {
      int q4[4] = {leaves[i], leaves[j], leaves[k], leaves[m]};
      double w = 1; for (int z = 0; z < 4; ++z) w *= t.weight[q4[z]];
      int pr[3][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
      for (int p = 0; p < 3; ++p) {
        int A = q4[pr[p][0]], B = q4[pr[p][1]], C = q4[pr[p][2]], D = q4[pr[p][3]];
        auto s1 = pathset(A, B), s2 = pathset(C, D);
        bool disj = true; for (int z : s1) if (s2.count(z)) { disj = false; break; }
        if (disj) { int ci = channel_idx(color[A], color[B], color[C], color[D], K);
          if (ci >= 0) out[ci] += w; break; }
      }
    }
  return out;
}

qraft::CounterVec full_star_g(const BTree& t, int K, int v,
                             const std::vector<int>& color) {
  qraft::CounterVec g;
  bool first = true;
  for (int c : t.ch[v]) {
    qraft::CounterVec c1 =
        t.color[c] >= 0
            ? qraft::counters_leaf(K, color[c], t.weight[c])
            : qraft::counters_ig_tf(K, full_star_g(t, K, c, color));
    qraft::CounterVec arm = qraft::counters_c_to_g(K, c1);
    g = first ? arm : qraft::counters_gg(K, g, arm);
    first = false;
  }
  return g;
}

std::array<double, 3> full_dp(const BTree& t, int K,
                              const std::vector<int>& color) {
  return qraft::counters_readout(
      K, qraft::counters_ig_ff(K, full_star_g(t, K, t.root, color)));
}

void require_full_vec_close(const qraft::CounterVec& a,
                            const qraft::CounterVec& b,
                            const std::string& label) {
  require(a.interface == b.interface && a.orientation == b.orientation &&
              a.dimension == b.dimension,
          label + " type mismatch");
  for (int i = 0; i < a.dimension; ++i) {
    const double av = a.values[static_cast<std::size_t>(i)];
    const double bv = b.values[static_cast<std::size_t>(i)];
    if (!finite_close(av, bv, 1e-10, 1e-10)) {
      fail(label + " coordinate mismatch");
    }
  }
}

}  // namespace bdp_parity

void test_counters_packed_color_masks() {
  using qraft::CounterInterface;
  using qraft::CounterOrientation;
  using qraft::CounterVec;

  // The constant packed-layout table must agree with the validating
  // generated-table functions for every palette, packed sort and mask.
  for (int palette : {2, 4}) {
    const std::pair<CounterInterface, CounterOrientation> sorts[] = {
        {CounterInterface::Closed, CounterOrientation::None},
        {CounterInterface::C1, CounterOrientation::U},
        {CounterInterface::C1, CounterOrientation::D},
        {CounterInterface::C2, CounterOrientation::None},
        {CounterInterface::G, CounterOrientation::None}};
    for (const auto& [interface, orientation] : sorts) {
      for (int raw_mask = 0; raw_mask < (1 << palette); ++raw_mask) {
        const auto mask = static_cast<std::uint8_t>(raw_mask);
        const qraft::PackedCounterLayout& layout = qraft::counters_packed_layout(
            palette, interface, orientation, mask);
        const std::uint64_t bits = qraft::counters_packed_coordinate_bits(
            palette, interface, orientation, mask);
        require(layout.bits == bits, "packed layout bits mismatch");
        require(layout.width == qraft::counters_packed_coordinate_count(
                                    palette, interface, orientation, mask),
                "packed layout width mismatch");
        int rank = 0;
        for (int coordinate = 0;
             coordinate < static_cast<int>(qraft::kMaxCounters);
             ++coordinate) {
          const bool stored = ((bits >> coordinate) & 1u) != 0;
          require(layout.rank[static_cast<std::size_t>(coordinate)] ==
                      (stored ? static_cast<std::uint8_t>(rank)
                              : qraft::kCounterNotPacked),
                  "packed layout rank mismatch");
          if (stored) ++rank;
        }
      }
    }
  }

  const auto make_value = [](int palette, CounterInterface interface,
                             CounterOrientation orientation,
                             std::uint8_t mask, int salt) {
    CounterVec value;
    value.interface = interface;
    value.orientation = orientation;
    value.dimension = static_cast<std::uint8_t>(
        qraft::counters_dimension(palette, interface));
    value.values.fill(0.0);
    const std::uint64_t packed_bits = qraft::counters_packed_coordinate_bits(
        palette, interface, orientation, mask);
    for (int coordinate = 0; coordinate < value.dimension; ++coordinate) {
      if ((packed_bits >> coordinate) & 1u) {
        require(qraft::counters_coordinate_active(
                    palette, interface, coordinate, mask),
                "packed coordinate outside the active colour mask");
        value.values[static_cast<std::size_t>(coordinate)] =
            static_cast<double>((coordinate + 3) * (salt + 5)) / 17.0;
      }
    }
    return value;
  };
  const auto pack = [](int palette, const CounterVec& value,
                       std::uint8_t mask) {
    std::array<double, qraft::kMaxCounters> packed{};
    const std::uint64_t packed_bits = qraft::counters_packed_coordinate_bits(
        palette, value.interface, value.orientation, mask);
    int rank = 0;
    for (int coordinate = 0; coordinate < value.dimension; ++coordinate) {
      if ((packed_bits >> coordinate) & 1u) {
        packed[static_cast<std::size_t>(rank++)] =
            value.values[static_cast<std::size_t>(coordinate)];
      } else if (value.values[static_cast<std::size_t>(coordinate)] != 0.0) {
        fail("unpacked coordinate carries a nonzero value");
      }
    }
    require(rank == qraft::counters_packed_coordinate_count(
                        palette, value.interface, value.orientation, mask),
            "packed coordinate count mismatch");
    return packed;
  };
  const auto require_packed = [](int palette, const CounterVec& expected,
                                 std::uint8_t mask, const double* actual,
                                 const std::string& leg) {
    int rank = 0;
    const std::uint64_t packed_bits = qraft::counters_packed_coordinate_bits(
        palette, expected.interface, expected.orientation, mask);
    for (int coordinate = 0; coordinate < expected.dimension; ++coordinate) {
      if (!((packed_bits >> coordinate) & 1u)) {
        if (expected.values[static_cast<std::size_t>(coordinate)] != 0.0) {
          fail(leg + " dense value outside the packed block mask=" +
               std::to_string(mask) + " coordinate=" +
               std::to_string(coordinate));
        }
        continue;
      }
      const double expected_value =
          expected.values[static_cast<std::size_t>(coordinate)];
      const double actual_value = actual[rank++];
      if (!finite_close(actual_value, expected_value, 1e-12, 1e-12)) {
        fail(leg + " mask=" + std::to_string(mask) +
             " coordinate=" + std::to_string(coordinate) +
             " rank=" + std::to_string(rank - 1) +
             " actual=" + std::to_string(actual_value) +
             " expected=" + std::to_string(expected_value));
      }
    }
  };
  const auto fresh_output = [] {
    std::array<double, qraft::kMaxCounters> output;
    output.fill(std::numeric_limits<double>::quiet_NaN());
    return output;
  };

  for (int palette : {2, 4}) {
    const int mask_count = 1 << palette;
    for (int raw_mask = 0; raw_mask < mask_count; ++raw_mask) {
      const auto mask = static_cast<std::uint8_t>(raw_mask);
      const CounterVec c1_u = make_value(
          palette, CounterInterface::C1, CounterOrientation::U, mask, 1);
      const CounterVec c2 = make_value(
          palette, CounterInterface::C2, CounterOrientation::None, mask, 2);
      const CounterVec g = make_value(
          palette, CounterInterface::G, CounterOrientation::None, mask, 3);
      const auto packed_c1_u = pack(palette, c1_u, mask);
      const auto packed_c2 = pack(palette, c2, mask);
      const auto packed_g = pack(palette, g, mask);

      auto output = fresh_output();
      qraft::counters_c_to_g_color_mask_into(
          palette, mask, packed_c1_u.data(), output.data());
      require_packed(palette, qraft::counters_c_to_g(palette, c1_u), mask,
                     output.data(), "packed C-to-G");

      output = fresh_output();
      qraft::counters_close_color_mask_into(
          palette, mask, packed_c1_u.data(), output.data());
      const CounterVec closed = qraft::counters_close(palette, c1_u);
      require_packed(palette, closed, mask, output.data(), "packed close");
      const auto packed_closed = pack(palette, closed, mask);
      const auto packed_counts = qraft::counters_color_mask_readout(
          palette, mask, packed_closed.data());
      const auto dense_counts = qraft::counters_readout(palette, closed);
      for (int channel = 0; channel < 3; ++channel) {
        if (!finite_close(
                packed_counts[static_cast<std::size_t>(channel)],
                dense_counts[static_cast<std::size_t>(channel)], 1e-12,
                1e-12)) {
          fail("packed root readout mismatch palette=" +
               std::to_string(palette) + " mask=" + std::to_string(mask) +
               " channel=" + std::to_string(channel) + " actual=" +
               std::to_string(
                   packed_counts[static_cast<std::size_t>(channel)]) +
               " expected=" + std::to_string(
                   dense_counts[static_cast<std::size_t>(channel)]));
        }
      }

      output = fresh_output();
      qraft::counters_reverse_c1_color_mask_into(
          palette, CounterOrientation::U, mask, packed_c1_u.data(),
          output.data());
      require_packed(palette, qraft::counters_reverse(palette, c1_u), mask,
                     output.data(), "packed C1 reversal U->D");
      {
        const CounterVec c1_d = make_value(
            palette, CounterInterface::C1, CounterOrientation::D, mask, 4);
        const auto packed_c1_d = pack(palette, c1_d, mask);
        output = fresh_output();
        qraft::counters_reverse_c1_color_mask_into(
            palette, CounterOrientation::D, mask, packed_c1_d.data(),
            output.data());
        require_packed(palette, qraft::counters_reverse(palette, c1_d), mask,
                       output.data(), "packed C1 reversal D->U");
      }

      output = fresh_output();
      qraft::counters_reverse_c2_color_mask_into(
          palette, mask, packed_c2.data(), output.data());
      require_packed(palette, qraft::counters_reverse(palette, c2), mask,
                     output.data(), "packed C2 reversal");

      output = fresh_output();
      qraft::counters_ig_ff_color_mask_into(
          palette, mask, packed_g.data(), output.data());
      require_packed(palette, qraft::counters_ig_ff(palette, g), mask,
                     output.data(), "packed IG-FF");

      output = fresh_output();
      qraft::counters_ig_tf_color_mask_into(
          palette, mask, packed_g.data(), output.data());
      require_packed(palette, qraft::counters_ig_tf(palette, g), mask,
                     output.data(), "packed IG-TF");

      output = fresh_output();
      qraft::counters_ig_tt_color_mask_into(
          palette, mask, packed_g.data(), output.data());
      require_packed(palette, qraft::counters_ig_tt(palette, g), mask,
                     output.data(), "packed IG-TT");

      for (int color = 0; color < palette; ++color) {
        const auto leaf_mask = static_cast<std::uint8_t>(1u << color);
        output = fresh_output();
        qraft::counters_leaf_color_mask_into(
            palette, color, 2.75, leaf_mask, output.data());
        require_packed(palette, qraft::counters_leaf(palette, color, 2.75),
                       leaf_mask, output.data(), "packed leaf");
      }

      for (int raw_right_mask = 0; raw_right_mask < mask_count;
           ++raw_right_mask) {
        const auto right_mask =
            static_cast<std::uint8_t>(raw_right_mask);
        const auto union_mask =
            static_cast<std::uint8_t>(mask | right_mask);
        const CounterVec right_g = make_value(
            palette, CounterInterface::G, CounterOrientation::None,
            right_mask, 7);
        const auto packed_right_g = pack(palette, right_g, right_mask);
        output = fresh_output();
        qraft::counters_gg_color_mask_into(
            palette, mask, right_mask, packed_g.data(),
            packed_right_g.data(), output.data());
        require_packed(palette, qraft::counters_gg(palette, g, right_g),
                       union_mask, output.data(),
                       "packed GG palette=" + std::to_string(palette) +
                           " left=" + std::to_string(mask) + " right=" +
                           std::to_string(right_mask));

        for (bool higher_c1 : {false, true}) {
          for (bool lower_c1 : {false, true}) {
            const CounterVec higher = make_value(
                palette,
                higher_c1 ? CounterInterface::C1 : CounterInterface::C2,
                higher_c1 ? CounterOrientation::D
                          : CounterOrientation::None,
                mask, 11);
            const CounterVec lower = make_value(
                palette,
                lower_c1 ? CounterInterface::C1 : CounterInterface::C2,
                lower_c1 ? CounterOrientation::U
                         : CounterOrientation::None,
                right_mask, 13);
            const auto packed_higher = pack(palette, higher, mask);
            const auto packed_lower = pack(palette, lower, right_mask);
            output = fresh_output();
            qraft::counters_cc_color_mask_into(
                palette, higher.interface, lower.interface, mask, right_mask,
                packed_higher.data(), packed_lower.data(), output.data());
            require_packed(
                palette, qraft::counters_cc(palette, higher, lower), union_mask,
                output.data(),
                "packed CC palette=" + std::to_string(palette) +
                    " left=" + std::to_string(mask) + " right=" +
                    std::to_string(right_mask) + " interfaces=" +
                    std::to_string(higher_c1) + std::to_string(lower_c1));
          }
        }
      }
    }
  }
}

void test_counters_match_enumeration() {
  using namespace counting_parity;
  std::mt19937 rng(20260712u);
  for (int K : {2, 4}) {
    const int bg_hi = (K == 4) ? 3 : 1;  // background colors: {0, bg_hi}
    std::uniform_int_distribution<int> bgcol(0, 1);
    for (int trial = 0; trial < 400; ++trial) {
      BTree t = random_tree(rng, 4 + static_cast<int>(rng() % 5), K);
      // current coloring
      std::vector<int> cur(t.ch.size());
      for (int v = 0; v < static_cast<int>(t.ch.size()); ++v) cur[v] = t.color[v];
      // background coloring: each leaf gets 0 or bg_hi
      std::vector<int> bg(t.ch.size(), -1);
      for (int v = 0; v < static_cast<int>(t.ch.size()); ++v)
        if (t.color[v] >= 0) bg[v] = bgcol(rng) ? bg_hi : 0;
      for (const std::vector<int>* col : {&cur, &bg}) {
        std::array<double, 3> b = brute(t, *col, K);
        std::array<double, 3> fd = full_dp(t, K, *col);
        for (int c = 0; c < (K == 2 ? 2 : 3); ++c)
          if (!finite_close(b[c], fd[c], 1e-6, 1e-6))
            fail("counting-scheme parity mismatch (K=" + std::to_string(K) +
                 ", ch=" + std::to_string(c) + ", bf=" + std::to_string(b[c]) +
                 ", counters=" + std::to_string(fd[c]) + ")");
      }
    }

    const qraft::CounterVec leaf0 = qraft::counters_leaf(K, 0, 2.0);
    const qraft::CounterPair absent =
        qraft::counters_leaf_dual(K, K - 1, 0, 0.0);
    require(absent.current_equals_background,
            "zero-mass leaf was not recognized as frozen");
    require_full_vec_close(absent.current, absent.background,
                           "zero-mass leaf dual parity");
    const qraft::CounterVec arm0 = qraft::counters_c_to_g(K, leaf0);
    require_full_vec_close(
        leaf0, qraft::counters_reverse(
                   K, qraft::counters_reverse(K, leaf0)),
        "C1 reverse involution");
    const qraft::CounterVec identity = qraft::counters_identity_c2(K);
    require_full_vec_close(
        identity, qraft::counters_reverse(K, identity),
        "C2 identity reversal");

    auto center_piece = [&](int color, double weight) {
      return qraft::counters_ig_tt(
          K, qraft::counters_c_to_g(
                 K, qraft::counters_leaf(K, color, weight)));
    };
    const qraft::CounterVec a = center_piece(0, 2.0);
    const qraft::CounterVec b = center_piece(1, 3.0);
    const qraft::CounterVec c = center_piece(K == 4 ? 2 : 0, 5.0);
    const qraft::CounterVec c1u = qraft::counters_ig_tf(K, arm0);
    const qraft::CounterVec c1d = qraft::counters_reverse(K, c1u);
    require_full_vec_close(
        qraft::counters_cc(K, qraft::counters_cc(K, a, b), c),
        qraft::counters_cc(K, a, qraft::counters_cc(K, b, c)),
        "C2 corridor associativity");
    require_full_vec_close(
        qraft::counters_reverse(K, qraft::counters_cc(K, a, b)),
        qraft::counters_cc(K, qraft::counters_reverse(K, b),
                          qraft::counters_reverse(K, a)),
        "C2 corridor opposite-orientation equivalence");

    const qraft::CounterPair dual =
        qraft::counters_leaf_dual(K, 1, K == 4 ? 3 : 0, 7.0);
    require_full_vec_close(dual.current, qraft::counters_leaf(K, 1, 7.0),
                           "dual current leaf");
    require_full_vec_close(
        dual.background, qraft::counters_leaf(K, K == 4 ? 3 : 0, 7.0),
        "dual background leaf");
  }
}

struct ProductionSmallTreeCase {
  std::uint64_t case_index = 0;
  int leaf_count = 0;
  int palette = 0;
  std::string weight_kind;
  std::uint64_t selector_mask = 0;
  std::uint64_t nested_mask = 0;
  int recolor_owner = -1;
  int recolor_color = 0;
  int node_count = 0;
  int root = -1;
  std::string topology;
  std::vector<std::uint8_t> colors;
  std::vector<std::uint8_t> background;
  std::vector<std::pair<long long, long long>> rational_weights;
  std::vector<std::pair<int, int>> edges;
  // current, background, selected, selected+recolor, nested.
  std::array<std::array<double, 3>, 5> expected{};
};

void require_protocol_token(std::istream& input,
                            const std::string& expected) {
  std::string actual;
  require(static_cast<bool>(input >> actual),
          "production-small protocol ended before " + expected);
  require(actual == expected,
          "production-small protocol expected " + expected + " but read " +
              actual);
}

std::vector<int> owners_from_mask(int leaf_count, std::uint64_t mask) {
  std::vector<int> owners;
  for (int owner = 0; owner < leaf_count; ++owner) {
    if ((mask & (std::uint64_t{1} << owner)) != 0) owners.push_back(owner);
  }
  return owners;
}

std::string small_tree_case_context(const ProductionSmallTreeCase& item,
                                    std::uint64_t seed,
                                    int shard_index,
                                    int shard_count) {
  auto bytes = [](const std::vector<std::uint8_t>& values) {
    std::ostringstream out;
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i) {
      if (i != 0) out << ',';
      out << static_cast<int>(values[i]);
    }
    return out.str() + ']';
  };
  std::ostringstream weights;
  weights << '[';
  for (std::size_t i = 0; i < item.rational_weights.size(); ++i) {
    if (i != 0) weights << ',';
    weights << item.rational_weights[i].first << '/'
            << item.rational_weights[i].second;
  }
  weights << ']';

  std::ostringstream out;
  out << "seed=" << seed << " shard=" << shard_index << '/' << shard_count
      << " case=" << item.case_index << " n=" << item.leaf_count
      << " K=" << item.palette << " topology=" << item.topology
      << " coloring=" << bytes(item.colors)
      << " background=" << bytes(item.background)
      << " weights=" << weights.str() << " weight_kind=" << item.weight_kind
      << " selector=0x" << std::hex << item.selector_mask
      << " nested=0x" << item.nested_mask << std::dec
      << " recolor=" << item.recolor_owner << "->" << item.recolor_color;
  return out.str();
}

qraft::internal::GeneIndex production_small_gene(
    const ProductionSmallTreeCase& item) {
  using qraft::internal::GeneIndex;
  require(item.node_count >= item.leaf_count,
          "protocol node count is smaller than its leaf count");
  require(item.root >= item.leaf_count && item.root < item.node_count,
          "protocol root is not an internal node");
  require(static_cast<int>(item.rational_weights.size()) == item.leaf_count,
          "protocol weight count mismatch");

  GeneIndex gene;
  gene.nodes.resize(static_cast<std::size_t>(item.node_count));
  gene.taxon_to_node.resize(static_cast<std::size_t>(item.leaf_count));
  for (int leaf = 0; leaf < item.leaf_count; ++leaf) {
    const auto [numerator, denominator] =
        item.rational_weights[static_cast<std::size_t>(leaf)];
    require(denominator > 0, "protocol weight denominator is not positive");
    const double weight = static_cast<double>(numerator) /
                          static_cast<double>(denominator);
    require(std::isfinite(weight) && weight >= 0.0,
            "protocol injected a non-finite or negative weight");
    gene.nodes[static_cast<std::size_t>(leaf)].set_taxon(leaf);
    qraft::internal::set_gene_leaf_mass(gene, leaf, weight);
    gene.taxon_to_node[static_cast<std::size_t>(leaf)] = leaf;
  }
  std::vector<std::vector<int>> adjacency(gene.nodes.size());
  for (const auto& [a, b] : item.edges) {
    require(a >= 0 && a < item.node_count && b >= 0 && b < item.node_count &&
                a != b,
            "protocol edge endpoint is invalid");
    adjacency[static_cast<std::size_t>(a)].push_back(b);
    adjacency[static_cast<std::size_t>(b)].push_back(a);
  }
  qraft::internal::set_gene_adjacency(gene, adjacency);
  require(static_cast<int>(item.edges.size()) == item.node_count - 1,
          "protocol topology is not a tree by edge count");
  gene.root = item.root;
  gene.reference_root = 0;
  gene.root_traversal = qraft::internal::traversal(gene, gene.root);
  require(static_cast<int>(gene.root_traversal.preorder.size()) ==
              item.node_count,
          "protocol topology is disconnected");
  return gene;
}

counting_parity::BTree production_small_boundary_tree(
    const qraft::internal::GeneIndex& gene,
    const std::vector<std::uint8_t>& colors) {
  counting_parity::BTree tree;
  tree.root = gene.root;
  tree.ch.resize(gene.nodes.size());
  tree.color.assign(gene.nodes.size(), -1);
  tree.weight.assign(gene.nodes.size(), 0.0);
  auto orient = [&](auto&& self, int node, int parent) -> void {
    const auto& source = gene.nodes[static_cast<std::size_t>(node)];
    if (source.taxon() >= 0) {
      require(source.taxon() < static_cast<int>(colors.size()),
              "counter fold taxon is outside the coloring");
      tree.color[static_cast<std::size_t>(node)] =
          colors[static_cast<std::size_t>(source.taxon())];
      tree.weight[static_cast<std::size_t>(node)] = qraft::internal::gene_leaf_mass(gene, node);
      return;
    }
    for (int child : qraft::internal::gene_neighbors(gene, node)) {
      if (child == parent) continue;
      tree.ch[static_cast<std::size_t>(node)].push_back(child);
      self(self, child, node);
    }
  };
  orient(orient, gene.root, -1);
  return tree;
}

void require_production_small_counts(
    const std::array<double, 3>& actual,
    const std::array<double, 3>& expected,
    int palette,
    const std::string& leg) {
  const int channels = palette == 2 ? 2 : 3;
  for (int channel = 0; channel < channels; ++channel) {
    require_close(actual[static_cast<std::size_t>(channel)],
                  expected[static_cast<std::size_t>(channel)],
                  leg + " channel=" + std::to_string(channel), 1e-9, 1e-9);
  }
}

void validate_production_small_state(
    const qraft::internal::HdtState& state,
    const std::string& leg) {
  for (const auto& unit : state.components.units) {
    require_finite(unit.weight, leg + " component weight");
  }
  for (double value : state.components.current_payload) {
    require_finite(value, leg + " component current payload");
  }
  for (double value : state.components.background_payload) {
    require_finite(value, leg + " component background payload");
  }
  for (double value : state.hdt.current_payload.values) {
    require_finite(value, leg + " HDT current payload");
  }
  for (double value : state.hdt.background_payload.values) {
    require_finite(value, leg + " HDT background payload");
  }
  std::string error;
  require(qraft::internal::validate_component_graph(state.components, &error),
          leg + " component graph invalid: " + error);
  error.clear();
  require(qraft::internal::validate_hdt(state.components, state.hdt,
                                                 &error),
          leg + " HDT invalid: " + error);
}

struct ProductionSmallTreeCoverage {
  std::uint64_t cases = 0;
  std::uint64_t quotient_folds = 0;
  std::uint64_t base_builds = 0;
  std::uint64_t contract_calls = 0;
  std::uint64_t extract_contract_calls = 0;
  std::uint64_t recolor_calls = 0;
  std::uint64_t oracle_comparisons = 0;
  long long gg_ops = 0;
  long long ig_ops = 0;
  long long cc_ops = 0;
};

void run_production_small_tree_case(const ProductionSmallTreeCase& item,
                                    ProductionSmallTreeCoverage& coverage,
                                    std::string& operation_trace) {
  using namespace qraft::internal;
  require(item.palette == 2 || item.palette == 4,
          "protocol palette must be K=2 or K=4");
  require(item.leaf_count >= 4 && item.leaf_count <= 8,
          "protocol leaf count must be in [4,8]");
  require(item.recolor_owner >= 0 && item.recolor_owner < item.leaf_count &&
              item.recolor_color >= 0 && item.recolor_color < item.palette,
          "protocol recolor is out of range");
  require(static_cast<int>(item.colors.size()) == item.leaf_count &&
              static_cast<int>(item.background.size()) == item.leaf_count,
          "protocol coloring size mismatch");
  const std::uint64_t valid_mask =
      (std::uint64_t{1} << item.leaf_count) - 1;
  require((item.selector_mask & ~valid_mask) == 0 &&
              (item.nested_mask & ~item.selector_mask) == 0,
          "protocol selector masks are invalid");
  for (std::uint8_t color : item.colors) {
    require(color < item.palette, "protocol current color is out of range");
  }
  for (std::uint8_t color : item.background) {
    require(color < item.palette, "protocol background color is out of range");
  }
  for (const auto& expected : item.expected) {
    for (double value : expected) {
      require_finite(value, "protocol independent-oracle value");
    }
  }

  GeneIndex gene = production_small_gene(item);
  const OwnerProjection owners = identity_owner_map(item.leaf_count);
  const counting_parity::BTree boundary =
      production_small_boundary_tree(gene, item.colors);
  std::vector<int> quotient_current(boundary.color.begin(),
                                    boundary.color.end());
  std::vector<int> quotient_background(gene.nodes.size(), -1);
  for (int node = 0; node < static_cast<int>(gene.nodes.size()); ++node) {
    const int taxon = gene.nodes[static_cast<std::size_t>(node)].taxon();
    if (taxon >= 0) {
      quotient_background[static_cast<std::size_t>(node)] =
          item.background[static_cast<std::size_t>(taxon)];
    }
  }

  operation_trace = "quotient(current)";
  require_production_small_counts(
      counting_parity::full_dp(boundary, item.palette, quotient_current),
      item.expected[0], item.palette, operation_trace);
  operation_trace += ">quotient(background)";
  require_production_small_counts(
      counting_parity::full_dp(boundary, item.palette, quotient_background),
      item.expected[1], item.palette, operation_trace);
  coverage.quotient_folds += 2;
  coverage.oracle_comparisons += 2;

  HdtColorPolicy policy;
  policy.current_by_owner = item.colors;
  policy.background_by_owner = item.background;
  operation_trace += ">base>merge";
  HdtComponentGraph graph =
      build_base_component_graph(gene, owners, item.palette, policy);
  std::string graph_error;
  require(validate_component_graph(graph, &graph_error),
          "base component graph invalid: " + graph_error);
  HdtState base = build_hdt_state(std::move(graph));
  validate_production_small_state(base, "base");
  require(base.hdt.stats.gg_ops + base.hdt.stats.ig_ops +
                  base.hdt.stats.cc_ops >
              0,
          "base build did not execute a merge");
  require_production_small_counts(hdt_root_counts(base), item.expected[0],
                                  item.palette, "base current");
  const qraft::CounterVec base_background =
      hdt_node_background(base.hdt, base.hdt.root);
  require_production_small_counts(
      qraft::counters_readout(item.palette, base_background), item.expected[1],
      item.palette, "base background");
  ++coverage.base_builds;
  coverage.oracle_comparisons += 2;
  coverage.gg_ops += base.hdt.stats.gg_ops;
  coverage.ig_ops += base.hdt.stats.ig_ops;
  coverage.cc_ops += base.hdt.stats.cc_ops;

  const HdtLiveSelector selector = hdt_owner_selector(
      item.leaf_count, owners_from_mask(item.leaf_count, item.selector_mask));
  operation_trace += ">Contract(light)";
  HdtContractResult contracted = hdt_contract(
      base, selector, HdtContractMode::LightCopy);
  validate_production_small_state(contracted.state, "Contract(light)");
  require_production_small_counts(hdt_root_counts(contracted.state),
                                  item.expected[2], item.palette,
                                  "Contract(light)");
  ++coverage.contract_calls;
  ++coverage.oracle_comparisons;
  coverage.gg_ops += contracted.stats.gg_ops;
  coverage.ig_ops += contracted.stats.ig_ops;
  coverage.cc_ops += contracted.stats.cc_ops;

  HdtExtractWorkspace workspace;
  operation_trace += ">ExtractContract";
  HdtExtractContractResult extracted =
      hdt_extract_contract(base, selector, workspace);
  validate_production_small_state(extracted.state, "ExtractContract");
  require_production_small_counts(hdt_root_counts(extracted.state),
                                  item.expected[2], item.palette,
                                  "ExtractContract");
  ++coverage.extract_contract_calls;
  ++coverage.oracle_comparisons;
  coverage.gg_ops += extracted.contract.gg_ops;
  coverage.ig_ops += extracted.contract.ig_ops;
  coverage.cc_ops += extracted.contract.cc_ops;

  operation_trace += ">recolor(Contract)";
  const int contract_changed = hdt_recolor_owner(
      contracted.state, item.recolor_owner,
      static_cast<std::uint8_t>(item.recolor_color));
  operation_trace += ">recolor(ExtractContract)";
  const int extract_changed = hdt_recolor_owner(
      extracted.state, item.recolor_owner,
      static_cast<std::uint8_t>(item.recolor_color));
  const int expected_changed =
      item.recolor_owner >= 0 &&
              (item.selector_mask &
               (std::uint64_t{1} << item.recolor_owner)) != 0
          ? 1
          : 0;
  require(contract_changed == expected_changed &&
              extract_changed == expected_changed,
          "contextual recolor changed an unexpected number of live leaves");
  require_production_small_counts(hdt_root_counts(contracted.state),
                                  item.expected[3], item.palette,
                                  "Contract contextual recolor");
  require_production_small_counts(hdt_root_counts(extracted.state),
                                  item.expected[3], item.palette,
                                  "ExtractContract contextual recolor");
  coverage.recolor_calls += 2;
  coverage.oracle_comparisons += 2;

  const HdtLiveSelector nested_selector = hdt_owner_selector(
      item.leaf_count, owners_from_mask(item.leaf_count, item.nested_mask));
  operation_trace += ">ExtractContract(after Contract+recolor)";
  HdtExtractContractResult nested_extract = hdt_extract_contract(
      contracted.state, nested_selector, workspace);
  validate_production_small_state(
      nested_extract.state, "ExtractContract after Contract+recolor");
  require_production_small_counts(
      hdt_root_counts(nested_extract.state), item.expected[4],
      item.palette, "ExtractContract after Contract+recolor");
  ++coverage.extract_contract_calls;
  ++coverage.oracle_comparisons;
  coverage.gg_ops += nested_extract.contract.gg_ops;
  coverage.ig_ops += nested_extract.contract.ig_ops;
  coverage.cc_ops += nested_extract.contract.cc_ops;

  operation_trace += ">Contract(heavy after ExtractContract+recolor)";
  HdtContractResult nested_contract = hdt_contract(
      extracted.state, nested_selector, HdtContractMode::HeavyReplace);
  validate_production_small_state(
      nested_contract.state, "Contract after ExtractContract+recolor");
  require_production_small_counts(
      hdt_root_counts(nested_contract.state), item.expected[4],
      item.palette, "Contract after ExtractContract+recolor");
  ++coverage.contract_calls;
  ++coverage.oracle_comparisons;
  coverage.gg_ops += nested_contract.stats.gg_ops;
  coverage.ig_ops += nested_contract.stats.ig_ops;
  coverage.cc_ops += nested_contract.stats.cc_ops;
  ++coverage.cases;
}

void test_production_small_trees_protocol() {
  require_protocol_token(std::cin, "QRAFT_PRODUCTION_SMALL_TREES_V1");
  std::uint64_t seed = 0;
  int shard_index = -1;
  int shard_count = 0;
  int max_leaves = 0;
  require(static_cast<bool>(std::cin >> seed >> shard_index >> shard_count >>
                            max_leaves),
          "invalid production-small protocol header");
  require(max_leaves >= 4 && max_leaves <= 8 && shard_count > 0 &&
              shard_index >= 0 && shard_index < shard_count,
          "invalid production-small protocol configuration");

  ProductionSmallTreeCoverage coverage;
  std::uint64_t declared_global_cases = 0;
  std::uint64_t declared_emitted_cases = 0;
  for (;;) {
    std::string token;
    require(static_cast<bool>(std::cin >> token),
            "production-small protocol ended without DONE");
    if (token == "DONE") {
      require(static_cast<bool>(std::cin >> declared_global_cases >>
                                declared_emitted_cases),
              "invalid production-small DONE record");
      break;
    }
    require(token == "CASE",
            "production-small protocol expected CASE but read " + token);

    ProductionSmallTreeCase item;
    require(static_cast<bool>(
                std::cin >> item.case_index >> item.leaf_count >> item.palette >>
                item.weight_kind >> item.selector_mask >> item.nested_mask >>
                item.recolor_owner >> item.recolor_color >> item.node_count >>
                item.root >> item.topology),
            "invalid production-small CASE record");
    require(item.case_index % static_cast<std::uint64_t>(shard_count) ==
                static_cast<std::uint64_t>(shard_index),
            "protocol emitted a case into the wrong shard");

    require_protocol_token(std::cin, "COLORS");
    item.colors.resize(static_cast<std::size_t>(item.leaf_count));
    for (std::uint8_t& color : item.colors) {
      int value = 0;
      require(static_cast<bool>(std::cin >> value),
              "invalid protocol current coloring");
      color = static_cast<std::uint8_t>(value);
    }
    require_protocol_token(std::cin, "BACKGROUND");
    item.background.resize(static_cast<std::size_t>(item.leaf_count));
    for (std::uint8_t& color : item.background) {
      int value = 0;
      require(static_cast<bool>(std::cin >> value),
              "invalid protocol background coloring");
      color = static_cast<std::uint8_t>(value);
    }
    require_protocol_token(std::cin, "WEIGHTS");
    item.rational_weights.resize(static_cast<std::size_t>(item.leaf_count));
    for (auto& weight : item.rational_weights) {
      require(static_cast<bool>(std::cin >> weight.first >> weight.second),
              "invalid protocol rational weight");
    }
    require_protocol_token(std::cin, "EDGES");
    int edge_count = 0;
    require(static_cast<bool>(std::cin >> edge_count) && edge_count >= 0,
            "invalid protocol edge count");
    item.edges.resize(static_cast<std::size_t>(edge_count));
    for (auto& edge : item.edges) {
      require(static_cast<bool>(std::cin >> edge.first >> edge.second),
              "invalid protocol edge");
    }
    require_protocol_token(std::cin, "EXPECTED");
    for (auto& expected : item.expected) {
      for (double& value : expected) {
        require(static_cast<bool>(std::cin >> value),
                "invalid protocol expected count");
      }
    }
    require_protocol_token(std::cin, "END");

    std::string operation_trace = "parse";
    try {
      run_production_small_tree_case(item, coverage, operation_trace);
    } catch (const std::exception& ex) {
      fail(small_tree_case_context(item, seed, shard_index, shard_count) +
           " operation_trace=" + operation_trace + " cause=" + ex.what());
    }
  }

  require(coverage.cases == declared_emitted_cases,
          "production-small protocol emitted-count mismatch");
  require(coverage.cases > 0, "production-small shard contained no cases");
  std::cout << "production-small exhaustive shard passed: seed=" << seed
            << " shard=" << shard_index << '/' << shard_count
            << " max_leaves=" << max_leaves << " checked=" << coverage.cases
            << '/' << declared_global_cases
            << " quotient_folds=" << coverage.quotient_folds
            << " base_builds=" << coverage.base_builds
            << " Contract=" << coverage.contract_calls
            << " ExtractContract=" << coverage.extract_contract_calls
            << " recolors=" << coverage.recolor_calls
            << " oracle_comparisons=" << coverage.oracle_comparisons
            << " merges=(GG:" << coverage.gg_ops << ",IG:" << coverage.ig_ops
            << ",CC:" << coverage.cc_ops << ")\n";
}

using TestFn = void (*)();

const std::map<std::string, TestFn>& tests() {
  static const std::map<std::string, TestFn> kTests = {
      {"finite_score_guards",
       test_finite_score_guards},
      {"tie_band_selection_permutations",
       test_tie_band_selection_permutations},
      {"bipartition_near_ties",
       test_bipartition_near_ties},
      {"newick_parser", test_newick_parser},
      {"newick_structure", test_newick_structure},
      {"newick_deep_trees", test_newick_deep_trees},
      {"gene_structure_validation",
       test_gene_structure_validation},
      {"checked_arithmetic_guards",
       test_checked_arithmetic_guards},
      {"hdt_fail_closed_guards",
       test_hdt_fail_closed_guards},
      {"taxon_registry_scaling",
       test_taxon_registry_scaling},
      {"mre_duplicate_and_incomplete_genes",
       test_mre_duplicate_and_incomplete_genes},
      {"mre_exact_cluster_ids", test_mre_exact_cluster_ids},
      {"mre_exact_cluster_ids_parallel", test_mre_exact_cluster_ids_parallel},
      {"mre_exact_split_entries", test_mre_exact_split_entries},
      {"mre_update_stream_stop", test_mre_update_stream_stop},
      {"incomplete_gene_trees", test_incomplete_gene_trees},
      {"bipartition_test", test_bipartition_test},
      {"bipartition_test_distance_index", test_bipartition_test_distance_index},
      {"bipartition_test_incomplete_genes", test_bipartition_test_incomplete_genes},
      {"bipartition_test_shared_genes_reduce_to_j",
       test_bipartition_test_shared_genes_reduce_to_j},
      {"bipartition_test_shared_genes_noise_free",
       test_bipartition_test_shared_genes_noise_free},
      {"nni_split_support", test_nni_split_support},
      {"nni_exact_split_support", test_nni_exact_split_support},
      {"bipartition_scores_match_enumeration",
       test_bipartition_scores_match_enumeration},
      {"regraft_scores_match_enumeration",
       test_regraft_scores_match_enumeration},
      {"skeleton_gains_exactness", test_skeleton_gains_exactness},
      {"spr_exactness", test_spr_exactness},
      {"spr_polytomy_exactness", test_spr_polytomy_exactness},
      {"spr_gain_exactness", test_spr_gain_exactness},
      {"guide_tree", test_guide_tree},
      {"guide_seed", test_guide_seed},
      {"nni_gains_match_enumeration",
       test_nni_gains_match_enumeration},
      {"branch_support_matches_enumeration",
       test_branch_support_matches_enumeration},
      {"branch_support_output", test_branch_support_output},
      {"nni_binary_rooted_starting_tree",
       test_nni_binary_rooted_starting_tree},
      {"nni_local_pass_matches_full_pass",
       test_nni_local_pass_matches_full_pass},
      {"nni_batch_gains_add", test_nni_batch_gains_add},
      {"wide_count_exact", test_wide_count_exact},
      {"bipartition_test_threshold_values",
       test_bipartition_test_threshold_values},
      {"msc_simulation_recovers_species_tree",
       test_msc_simulation_recovers_species_tree},
      {"trivial_bipartition_completion", test_trivial_bipartition_completion},
      {"child_gene_set_partition_bound",
       test_child_gene_set_partition_bound},
      {"mre_complexity_primitives",
       test_mre_complexity_primitives},
      {"counters_packed_color_masks",
       test_counters_packed_color_masks},
      {"counters_match_enumeration", test_counters_match_enumeration},
      {"hdt_weighted_fixtures",
       test_hdt_weighted_fixtures},
      {"hdt_duplicate_owner_contexts",
       test_hdt_duplicate_owner_contexts},
      {"regraft_near_ties",
       test_regraft_near_ties},
      {"hdt_units",
       test_hdt_units},
      {"hdt_construction", test_hdt_construction},
      {"hdt_contract", test_hdt_contract},
      {"hdt_extract",
       test_hdt_extract},
      {"count_traversal_oracle",
       test_count_traversal_oracle},
      {"count_traversal_scaling",
       test_count_traversal_scaling},
      {"qraft_smoke", test_qraft_smoke},
      {"worker_cap_by_input_size",
       test_worker_cap_by_input_size},
      {"worker_count_determinism", test_worker_count_determinism},
      {"pipeline_determinism",
       test_pipeline_determinism},
      {"ablation_smoke", test_ablation_smoke},
      {"hdt_structural_bounds", test_hdt_structural_bounds},
      {"hdt_heavy_threshold",
       test_hdt_heavy_threshold},
      {"production_small_trees_protocol",
       test_production_small_trees_protocol},
      {"rf_identity", test_rf_identity},
      {"fixture_smoke", test_fixture_smoke},
  };
  return kTests;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: " << argv[0] << " <test-name>\n";
    return 2;
  }
  const auto it = tests().find(argv[1]);
  if (it == tests().end()) {
    std::cerr << "unknown test: " << argv[1] << "\n";
    return 2;
  }
  try {
    it->second();
  } catch (const std::exception& ex) {
    std::cerr << argv[1] << " failed: " << ex.what() << "\n";
    return 1;
  }
  return 0;
}
