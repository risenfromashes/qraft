// Branch support of the output tree (--support localpp): the local posterior
// probability of each internal edge (Sayyari and Mirarab 2016) and its length
// in coalescent units, from the gene-tree quartet frequencies f1, f2, f3 of
// the edge's own topology and its two alternatives. One K=4 counting pass
// over the final tree gives every edge's frequencies (branch_quartet_
// frequencies), as a full NNI scoring pass does.
//
// With n = f1 + f2 + f3 and the Yule prior rate lambda = 0.5 of ASTRAL and
// ASTER, topology i has posterior weight
//   w_i = (1 - I_{1/3}(f_i + 1, n - f_i + 2 lambda))
//         * 2^{f_i} * Gamma(f_i + 1) * Gamma(n - f_i + 2 lambda),
// with I the regularized incomplete beta function, and the local posterior
// probability of the tree's topology is w_1 / (w_1 + w_2 + w_3). The length
// is the posterior mode of t under the same prior, -ln(3/2 - 3/2 f1 /
// (n + 2 lambda)), when f1 > (n + 2 lambda) / 3, and 0 otherwise.

#include "internal.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace qraft::internal {
namespace {

constexpr double kYuleLambda = 0.5;

// The continued fraction of I_x(a, b) (modified Lentz), which converges fast
// for x < (a + 1) / (a + b + 2).
double beta_continued_fraction(double a, double b, double x) {
  constexpr int kMaxIterations = 1000000;
  constexpr double kEpsilon = 1e-16;
  constexpr double kTiny = 1e-300;
  const double qab = a + b, qap = a + 1.0, qam = a - 1.0;
  double c = 1.0;
  double d = 1.0 - qab * x / qap;
  if (std::fabs(d) < kTiny) d = kTiny;
  d = 1.0 / d;
  double h = d;
  for (int m = 1; m <= kMaxIterations; ++m) {
    const double m2 = 2.0 * m;
    double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d;
    if (std::fabs(d) < kTiny) d = kTiny;
    c = 1.0 + aa / c;
    if (std::fabs(c) < kTiny) c = kTiny;
    d = 1.0 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1.0 + aa * d;
    if (std::fabs(d) < kTiny) d = kTiny;
    c = 1.0 + aa / c;
    if (std::fabs(c) < kTiny) c = kTiny;
    d = 1.0 / d;
    const double step = d * c;
    h *= step;
    if (std::fabs(step - 1.0) < kEpsilon) return h;
  }
  throw std::runtime_error("incomplete beta continued fraction did not converge");
}

}  // namespace

double log_regularized_incomplete_beta(double a, double b, double x) {
  if (x <= 0.0) return -std::numeric_limits<double>::infinity();
  if (x >= 1.0) return 0.0;
  // Gamma(a + b) / (Gamma(a) Gamma(b)) x^a (1 - x)^b, shared by I_x(a, b)
  // and its complement I_{1-x}(b, a).
  const double log_front = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
                           a * std::log(x) + b * std::log1p(-x);
  if (x < (a + 1.0) / (a + b + 2.0)) {
    return log_front + std::log(beta_continued_fraction(a, b, x) / a);
  }
  return std::log1p(-std::exp(log_front) * beta_continued_fraction(b, a, 1.0 - x) / b);
}

std::array<double, 3> local_posterior(const std::array<double, 3>& f) {
  const double n = f[0] + f[1] + f[2];
  if (!(n > 0.0)) return {1.0 / 3, 1.0 / 3, 1.0 / 3};
  std::array<double, 3> log_weight{};
  for (int i = 0; i < 3; ++i) {
    const double a = f[static_cast<std::size_t>(i)] + 1.0;
    const double b = n - f[static_cast<std::size_t>(i)] + 2.0 * kYuleLambda;
    // 1 - I_{1/3}(a, b) = I_{2/3}(b, a), computed directly so that it keeps
    // its precision when it is tiny.
    log_weight[static_cast<std::size_t>(i)] =
        log_regularized_incomplete_beta(b, a, 2.0 / 3.0) +
        f[static_cast<std::size_t>(i)] * std::log(2.0) + std::lgamma(a) + std::lgamma(b);
  }
  const double top = *std::max_element(log_weight.begin(), log_weight.end());
  std::array<double, 3> pp{};
  double sum = 0.0;
  for (int i = 0; i < 3; ++i) {
    pp[static_cast<std::size_t>(i)] = std::exp(log_weight[static_cast<std::size_t>(i)] - top);
    sum += pp[static_cast<std::size_t>(i)];
  }
  for (double& p : pp) p /= sum;
  return pp;
}

double coalescent_length(double f1, double n) {
  if (!(n > 0.0)) return std::numeric_limits<double>::quiet_NaN();
  if (!(3.0 * f1 > n)) return 0.0;
  // -ln(3/2 (n - f1 + 2 lambda) / (n + 2 lambda)) as a difference of logs,
  // which keeps the length finite when f1 = n is too large for n + 2 lambda
  // to differ from n.
  return std::max(0.0, std::log(n + 2.0 * kYuleLambda) -
                           std::log(n - f1 + 2.0 * kYuleLambda) - std::log(1.5));
}

void annotate_branch_support(const SolverContext& ctx,
                             const SolverTree& tree,
                             const std::vector<int>& species_node,
                             SpeciesTree& species) {
  const auto started = std::chrono::steady_clock::now();
  std::vector<int> taxa(static_cast<std::size_t>(context_taxa(ctx).size()));
  std::iota(taxa.begin(), taxa.end(), 0);
  const bool per_gene = ctx.options.support_frequencies == SupportFrequencies::Gene;
  const std::vector<BranchQuartetFrequencies> branches =
      branch_quartet_frequencies(ctx, tree, taxa, per_gene);
  for (const BranchQuartetFrequencies& branch : branches) {
    std::array<double, 3> f{branch.main, branch.alternative1, branch.alternative2};
    if (!per_gene) {
      // The average over the quartets around the edge.
      if (!(branch.quartets_around > 0.0)) continue;
      for (double& x : f) x /= branch.quartets_around;
    }
    const int a = species_node[static_cast<std::size_t>(branch.edge_a)];
    const int b = species_node[static_cast<std::size_t>(branch.edge_b)];
    if (a < 0 || b < 0) continue;
    const double n = f[0] + f[1] + f[2];
    const std::array<double, 3> pp = local_posterior(f);
    const double length = coalescent_length(f[0], n);
    std::ostringstream text;
    text << std::setprecision(6);
    if (ctx.options.support_detail) {
      text << "'[f1=" << f[0] << ";f2=" << f[1] << ";f3=" << f[2] << ";EN=" << n
           << ";q1=" << (n > 0 ? f[0] / n : 0.0) << ";pp1=" << pp[0] << ";pp2=" << pp[1]
           << ";pp3=" << pp[2] << ";CULength=" << length << "]'";
    } else {
      text << pp[0];
    }
    if (std::isfinite(length)) text << ':' << length;
    species.edge_labels.push_back(
        SpeciesTree::EdgeLabel{std::min(a, b), std::max(a, b), text.str()});
  }
  std::sort(species.edge_labels.begin(), species.edge_labels.end(),
            [](const SpeciesTree::EdgeLabel& x, const SpeciesTree::EdgeLabel& y) {
              return std::pair<int, int>{x.a, x.b} < std::pair<int, int>{y.a, y.b};
            });
  stats_add(run_stats().support_edges, static_cast<long long>(species.edge_labels.size()));
  stats_add_time(run_stats().support_ns,
                 std::chrono::duration_cast<std::chrono::nanoseconds>(
                     std::chrono::steady_clock::now() - started));
}

}  // namespace qraft::internal
