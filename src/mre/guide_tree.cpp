// Guide tree (--guide-tree): edits the starting tree with a guide tree on
// some or all of the taxa. The run builds the MRE without the clusters that
// conflict with the guide (GuideFilter), so the first step below finds
// nothing there; it contracts every starting-tree edge whose bipartition,
// restricted to the taxa both trees hold, conflicts with a guide bipartition.
// Then each polytomy is refined by the guide clusters that keep its arms
// whole and are bipartitions of the whole guide on the shared taxa. The
// result keeps every other starting-tree bipartition and adds only guide
// bipartitions; the divide then runs from it as from any starting tree.
//
// Every step applies its groups and keeps children in a fixed order, so the
// guided tree depends only on the starting tree and the guide.

#include "internal.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace qraft::internal {
namespace {

// A rooted tree as children lists. A leaf carries its taxon, or -1 when its
// label is not a taxon of the run; internal nodes carry -1.
struct Rooted {
  std::vector<std::vector<int>> kids;
  std::vector<int> taxon;
  int root = 0;
};

std::vector<int> preorder(const Rooted& t, int root) {
  std::vector<int> out, stack{root};
  while (!stack.empty()) {
    const int v = stack.back();
    stack.pop_back();
    out.push_back(v);
    const auto& k = t.kids[static_cast<std::size_t>(v)];
    for (auto it = k.rbegin(); it != k.rend(); ++it) stack.push_back(*it);
  }
  return out;
}

// Lowest common ancestors by an Euler tour and a sparse table, ancestors at a
// given depth by binary lifting, and the number of marked leaves below each
// node.
class Lca {
 public:
  explicit Lca(const Rooted& t) {
    const std::size_t n = t.kids.size();
    tin.assign(n, 0);
    depth.assign(n, 0);
    parent.assign(n, -1);
    std::vector<int> euler;
    std::vector<std::pair<int, std::size_t>> stack{{t.root, 0}};
    while (!stack.empty()) {
      const auto [v, k] = stack.back();
      stack.pop_back();
      if (k == 0) tin[static_cast<std::size_t>(v)] = static_cast<int>(euler.size());
      euler.push_back(v);
      const auto& kids = t.kids[static_cast<std::size_t>(v)];
      if (k < kids.size()) {
        stack.push_back({v, k + 1});
        const int c = kids[k];
        depth[static_cast<std::size_t>(c)] = depth[static_cast<std::size_t>(v)] + 1;
        parent[static_cast<std::size_t>(c)] = v;
        stack.push_back({c, 0});
      }
    }
    up.push_back(std::vector<int>(n));
    for (std::size_t v = 0; v < n; ++v) up[0][v] = parent[v] >= 0 ? parent[v] : static_cast<int>(v);
    while ((std::size_t{1} << up.size()) < n) {
      const auto& u = up.back();
      std::vector<int> next(n);
      for (std::size_t v = 0; v < n; ++v) next[v] = u[static_cast<std::size_t>(u[v])];
      up.push_back(std::move(next));
    }
    sparse.push_back(euler);
    for (std::size_t j = 1; (std::size_t{1} << j) <= euler.size(); ++j) {
      const auto& p = sparse.back();
      const std::size_t h = std::size_t{1} << (j - 1);
      std::vector<int> next(euler.size() - (std::size_t{1} << j) + 1);
      for (std::size_t i = 0; i < next.size(); ++i) next[i] = shallower(p[i], p[i + h]);
      sparse.push_back(std::move(next));
    }
  }

  int operator()(int a, int b) const {
    int l = tin[static_cast<std::size_t>(a)], r = tin[static_cast<std::size_t>(b)];
    if (l > r) std::swap(l, r);
    int j = 0;
    while ((2 << j) <= r - l + 1) ++j;
    return shallower(sparse[static_cast<std::size_t>(j)][static_cast<std::size_t>(l)],
                     sparse[static_cast<std::size_t>(j)][static_cast<std::size_t>(r - (1 << j) + 1)]);
  }

  int ancestor_at(int v, int target_depth) const {
    int k = depth[static_cast<std::size_t>(v)] - target_depth;
    for (std::size_t j = 0; k != 0; k >>= 1, ++j) {
      if (k & 1) v = up[j][static_cast<std::size_t>(v)];
    }
    return v;
  }

  void count_below(const Rooted& t, const std::vector<char>& marked) {
    below.assign(t.kids.size(), 0);
    const std::vector<int> order = preorder(t, t.root);
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
      const auto v = static_cast<std::size_t>(*it);
      int c = marked[v] ? 1 : 0;
      for (int u : t.kids[v]) c += below[static_cast<std::size_t>(u)];
      below[v] = c;
    }
    total = below[static_cast<std::size_t>(t.root)];
  }

  std::vector<int> tin, depth, parent, below;
  int total = 0;

 private:
  int shallower(int a, int b) const {
    return depth[static_cast<std::size_t>(a)] <= depth[static_cast<std::size_t>(b)] ? a : b;
  }
  std::vector<std::vector<int>> up, sparse;
};

Rooted rooted_from(const StartingTree& s) {
  Rooted t;
  t.root = s.root;
  t.kids.resize(s.nodes.size());
  t.taxon.resize(s.nodes.size());
  for (std::size_t i = 0; i < s.nodes.size(); ++i) {
    t.kids[i] = s.nodes[i].children;
    t.taxon[i] = s.nodes[i].taxon;
  }
  return t;
}

// The guide as parsed, rooted at its Newick root; labels outside the run get
// taxon -1.
Rooted rooted_from(const Tree& tree, const TaxonRegistry& taxa) {
  std::vector<std::pair<std::string, int>> by_label;
  by_label.reserve(static_cast<std::size_t>(taxa.size()));
  for (int id = 0; id < taxa.size(); ++id) by_label.push_back({taxa.label_for(id), id});
  std::sort(by_label.begin(), by_label.end());
  Rooted t;
  t.root = tree.root >= 0 ? tree.root : 0;
  t.kids.resize(tree.nodes.size());
  t.taxon.assign(tree.nodes.size(), -1);
  std::vector<int> parent(tree.nodes.size(), -2);
  std::vector<int> stack{t.root};
  parent[static_cast<std::size_t>(t.root)] = -1;
  while (!stack.empty()) {
    const int v = stack.back();
    stack.pop_back();
    const auto& node = tree.nodes[static_cast<std::size_t>(v)];
    for (int u : node.neighbors) {
      if (u == parent[static_cast<std::size_t>(v)]) continue;
      if (parent[static_cast<std::size_t>(u)] != -2) throw std::invalid_argument("guide tree is not acyclic");
      parent[static_cast<std::size_t>(u)] = v;
      t.kids[static_cast<std::size_t>(v)].push_back(u);
      stack.push_back(u);
    }
    if (t.kids[static_cast<std::size_t>(v)].empty()) {
      const auto it = std::lower_bound(by_label.begin(), by_label.end(),
                                       std::pair<std::string, int>{node.label, -1});
      if (it != by_label.end() && it->first == node.label) t.taxon[static_cast<std::size_t>(v)] = it->second;
    }
  }
  return t;
}

std::string escape_label(const std::string& label) {
  bool quote = label.empty();
  for (char ch : label) {
    if (std::isspace(static_cast<unsigned char>(ch)) || ch == '(' || ch == ')' || ch == ',' || ch == ':' ||
        ch == ';' || ch == '\'' || ch == '[' || ch == ']') {
      quote = true;
      break;
    }
  }
  if (!quote) return label;
  std::string out = "'";
  for (char ch : label) out += ch == '\'' ? std::string("''") : std::string(1, ch);
  return out + "'";
}

std::string newick(const Rooted& t, const TaxonRegistry& taxa) {
  std::string out;
  std::vector<std::pair<int, int>> stack{{t.root, -1}};  // node, or -1 with a literal in the second slot
  while (!stack.empty()) {
    const auto [v, lit] = stack.back();
    stack.pop_back();
    if (v < 0) {
      out += static_cast<char>(lit);
      continue;
    }
    const auto& k = t.kids[static_cast<std::size_t>(v)];
    if (k.empty()) {
      out += escape_label(taxa.label_for(t.taxon[static_cast<std::size_t>(v)]));
      continue;
    }
    out += '(';
    stack.push_back({-1, ')'});
    for (std::size_t j = k.size(); j-- > 0;) {
      stack.push_back({k[j], -1});
      if (j > 0) stack.push_back({-1, ','});
    }
  }
  return out + ";\n";
}

// Whether a bipartition conflicts with the guide on the shared taxa. With
// the guide rooted at the neighbour of a shared taxon b0 and A the restricted
// side without b0, the bipartition is compatible iff every child subtree of
// the guide node LCA(A) holds only A taxa or no A taxa among the shared ones.
class ConflictTest {
 public:
  ConflictTest(const Rooted& guide, const std::vector<int>& gleaf, const std::vector<int>& shared)
      : gleaf_(gleaf), shared_(shared), b0_(gleaf[static_cast<std::size_t>(shared.front())]),
        g_(rerooted(guide, b0_)), lca_(g_), in_a_(g_.kids.size(), 0) {
    std::vector<char> marked(g_.kids.size(), 0);
    for (int t : shared_) marked[static_cast<std::size_t>(gleaf_[static_cast<std::size_t>(t)])] = 1;
    lca_.count_below(g_, marked);
  }

  int b0() const { return b0_; }

  // `side`: the guide leaves of the shared taxa on one side; `has_b0`
  // whether b0 is among them.
  bool conflicts(const std::vector<int>& side, bool has_b0) const {
    std::vector<int> a;
    if (has_b0) {
      for (int x : side) in_a_[static_cast<std::size_t>(x)] = 1;
      for (int t : shared_) {
        const int x = gleaf_[static_cast<std::size_t>(t)];
        if (!in_a_[static_cast<std::size_t>(x)]) a.push_back(x);
      }
      for (int x : side) in_a_[static_cast<std::size_t>(x)] = 0;
    } else {
      a = side;
    }
    if (a.size() < 2 || shared_.size() - a.size() < 2) return false;
    std::sort(a.begin(), a.end(), [&](int x, int y) {
      return lca_.tin[static_cast<std::size_t>(x)] < lca_.tin[static_cast<std::size_t>(y)];
    });
    int top = a.front();
    for (std::size_t i = 1; i < a.size(); ++i) top = lca_(top, a[i]);
    std::map<int, int> per_child;
    const int dtop = lca_.depth[static_cast<std::size_t>(top)];
    for (int x : a) ++per_child[lca_.ancestor_at(x, dtop + 1)];
    for (const auto& [c, k] : per_child) {
      if (k != lca_.below[static_cast<std::size_t>(c)]) return true;
    }
    return false;
  }

 private:
  static Rooted rerooted(const Rooted& guide, int b0) {
    std::vector<std::vector<int>> adj(guide.kids.size());
    for (std::size_t v = 0; v < guide.kids.size(); ++v) {
      for (int c : guide.kids[v]) {
        adj[v].push_back(c);
        adj[static_cast<std::size_t>(c)].push_back(static_cast<int>(v));
      }
    }
    if (adj[static_cast<std::size_t>(b0)].size() != 1) throw std::invalid_argument("guide tree leaf has several neighbours");
    Rooted g;
    g.root = adj[static_cast<std::size_t>(b0)][0];
    g.kids.resize(guide.kids.size());
    g.taxon = guide.taxon;
    std::vector<char> seen(guide.kids.size(), 0);
    std::vector<int> stack{g.root};
    seen[static_cast<std::size_t>(g.root)] = 1;
    while (!stack.empty()) {
      const int v = stack.back();
      stack.pop_back();
      for (int u : adj[static_cast<std::size_t>(v)]) {
        if (seen[static_cast<std::size_t>(u)]) continue;
        seen[static_cast<std::size_t>(u)] = 1;
        g.kids[static_cast<std::size_t>(v)].push_back(u);
        stack.push_back(u);
      }
    }
    return g;
  }

  std::vector<int> gleaf_, shared_;
  int b0_;
  Rooted g_;
  Lca lca_;
  mutable std::vector<char> in_a_;
};

// Contracts every edge of `s` whose bipartition,
// restricted to the shared taxa, conflicts with the guide.
long long collapse_conflicts(Rooted& s, const Rooted& guide, const std::vector<int>& gleaf,
                             const std::vector<int>& shared) {
  const ConflictTest test(guide, gleaf, shared);
  const int b0 = test.b0();
  const std::vector<int> order = preorder(s, s.root);
  std::vector<std::vector<int>> leaves_of(s.kids.size());  // shared taxa below, as guide leaves
  std::vector<char> has_b0(s.kids.size(), 0), conflict(s.kids.size(), 0);
  long long contracted = 0;
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    const auto v = static_cast<std::size_t>(*it);
    if (s.kids[v].empty()) {
      const int t = s.taxon[v];
      if (t >= 0 && gleaf[static_cast<std::size_t>(t)] >= 0) {
        leaves_of[v].push_back(gleaf[static_cast<std::size_t>(t)]);
        has_b0[v] = gleaf[static_cast<std::size_t>(t)] == b0;
      }
      continue;
    }
    std::vector<int> acc;
    for (int c : s.kids[v]) {
      auto& b = leaves_of[static_cast<std::size_t>(c)];
      if (b.size() > acc.size()) std::swap(acc, b);
      acc.insert(acc.end(), b.begin(), b.end());
      std::vector<int>().swap(b);
      has_b0[v] = has_b0[v] || has_b0[static_cast<std::size_t>(c)];
    }
    leaves_of[v] = std::move(acc);
    if (static_cast<int>(v) == s.root) continue;
    if (test.conflicts(leaves_of[v], has_b0[v])) {
      conflict[v] = 1;
      ++contracted;
    }
  }
  // Contract: the children of a contracted node take its place, in order.
  for (int v : order) {
    if (conflict[static_cast<std::size_t>(v)]) continue;
    auto& k = s.kids[static_cast<std::size_t>(v)];
    bool any = false;
    for (int c : k) any = any || conflict[static_cast<std::size_t>(c)];
    if (!any) continue;
    std::vector<int> out, stack(k.rbegin(), k.rend());
    while (!stack.empty()) {
      const int c = stack.back();
      stack.pop_back();
      if (!conflict[static_cast<std::size_t>(c)]) {
        out.push_back(c);
        continue;
      }
      const auto& ck = s.kids[static_cast<std::size_t>(c)];
      for (auto it = ck.rbegin(); it != ck.rend(); ++it) stack.push_back(*it);
    }
    k = std::move(out);
  }
  return contracted;
}

// The groups of arms that the guide,
// restricted to arm_of's taxa and rooted at r, puts in one cluster that
// splits no arm and holds no shared taxon from outside the polytomy.
std::set<std::vector<int>> polytomy_groups(const Lca& lca, const std::map<int, int>& arm_of,
                                           const std::map<int, int>& tot, int r,
                                           const std::vector<int>& groupable) {
  std::vector<int> l;
  for (const auto& [x, a] : arm_of) l.push_back(x);
  const auto by_tin = [&](int x, int y) { return lca.tin[static_cast<std::size_t>(x)] < lca.tin[static_cast<std::size_t>(y)]; };
  std::sort(l.begin(), l.end(), by_tin);
  std::vector<int> nodes(l);
  for (std::size_t i = 0; i + 1 < l.size(); ++i) nodes.push_back(lca(l[i], l[i + 1]));
  std::sort(nodes.begin(), nodes.end());
  nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
  std::sort(nodes.begin(), nodes.end(), by_tin);
  std::map<int, int> vpar;
  std::vector<int> stack;
  for (int x : nodes) {
    while (!stack.empty() && lca(stack.back(), x) != stack.back()) stack.pop_back();
    vpar[x] = stack.empty() ? -1 : stack.back();
    stack.push_back(x);
  }
  const std::set<int> gset(groupable.begin(), groupable.end());
  std::map<int, std::map<int, int>> cnt;
  std::map<int, char> has_r;
  std::set<std::vector<int>> found;
  for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
    const int x = *it;
    std::map<int, int> c;
    if (const auto f = cnt.find(x); f != cnt.end()) {
      c = std::move(f->second);
      cnt.erase(f);
    }
    if (const auto f = arm_of.find(x); f != arm_of.end()) ++c[f->second];
    bool hr = x == r;
    if (const auto f = has_r.find(x); f != has_r.end()) {
      hr = hr || f->second;
      has_r.erase(f);
    }
    const int px = vpar[x];
    if (px < 0) continue;
    int n_below = 0;
    for (const auto& [a, k] : c) n_below += k;
    bool whole;
    if (hr) {  // the side without r lies outside the topmost guide node on the path from px to x
      const int top = lca.ancestor_at(x, lca.depth[static_cast<std::size_t>(px)] + 1);
      whole = lca.total - lca.below[static_cast<std::size_t>(top)] == static_cast<int>(l.size()) - n_below;
    } else {
      whole = lca.below[static_cast<std::size_t>(x)] == n_below;
    }
    bool arms_whole = true;
    for (const auto& [a, k] : c) arms_whole = arms_whole && k == tot.at(a);
    if (whole && arms_whole) {
      std::vector<int> g;
      for (int a : gset) {
        const bool present = c.count(a) != 0;
        if (present != hr) g.push_back(a);
      }
      if (g.size() >= 2 && g.size() < groupable.size()) found.insert(std::move(g));
    }
    auto& pc = cnt[px];
    if (pc.size() < c.size()) std::swap(pc, c);
    for (const auto& [a, k] : c) pc[a] += k;
    if (hr) has_r[px] = 1;
  }
  return found;
}

// Refines each polytomy (the root with at least four arms, any
// other node with at least three) by the guide clusters that keep its arms
// whole and are bipartitions of the whole guide on the shared taxa.
long long refine_polytomies(Rooted& s, const Rooted& guide, const std::vector<int>& gleaf) {
  Lca lca(guide);
  const std::vector<int> order = preorder(s, s.root);
  std::vector<int> lo(s.kids.size(), 0), hi(s.kids.size(), 0);
  int pos = 0;
  for (int v : order) {
    lo[static_cast<std::size_t>(v)] = pos;
    if (s.kids[static_cast<std::size_t>(v)].empty()) ++pos;
  }
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    const auto v = static_cast<std::size_t>(*it);
    if (s.kids[v].empty()) {
      hi[v] = lo[v] + 1;
    } else {
      int h = 0;
      for (int c : s.kids[v]) h = std::max(h, hi[static_cast<std::size_t>(c)]);
      hi[v] = h;
    }
  }
  std::vector<int> bb_pos, bb_leaf;  // shared leaves in preorder: position, guide leaf
  std::vector<char> marked(guide.kids.size(), 0);
  int i = 0;
  for (int v : order) {
    if (!s.kids[static_cast<std::size_t>(v)].empty()) continue;
    const int t = s.taxon[static_cast<std::size_t>(v)];
    if (t >= 0 && gleaf[static_cast<std::size_t>(t)] >= 0) {
      bb_pos.push_back(i);
      bb_leaf.push_back(gleaf[static_cast<std::size_t>(t)]);
      marked[static_cast<std::size_t>(gleaf[static_cast<std::size_t>(t)])] = 1;
    }
    ++i;
  }
  lca.count_below(guide, marked);
  const auto span = [&](int v) {
    return std::pair<int, int>{
        static_cast<int>(std::lower_bound(bb_pos.begin(), bb_pos.end(), lo[static_cast<std::size_t>(v)]) - bb_pos.begin()),
        static_cast<int>(std::lower_bound(bb_pos.begin(), bb_pos.end(), hi[static_cast<std::size_t>(v)]) - bb_pos.begin())};
  };
  std::vector<std::pair<int, std::set<std::vector<int>>>> groups_at;
  for (int v : order) {
    const int d = static_cast<int>(s.kids[static_cast<std::size_t>(v)].size());
    if (d == 0 || d < (v == s.root ? 4 : 3)) continue;
    std::map<int, int> arm_of, tot;
    std::vector<int> arms_bb;
    std::vector<std::pair<int, int>> first_leaf;  // arm -> its first guide leaf, in preorder
    for (int a = 0; a < d; ++a) {
      const auto [a0, b0] = span(s.kids[static_cast<std::size_t>(v)][static_cast<std::size_t>(a)]);
      if (b0 > a0) {
        arms_bb.push_back(a);
        tot[a] = b0 - a0;
        first_leaf.push_back({a, bb_leaf[static_cast<std::size_t>(a0)]});
        for (int j = a0; j < b0; ++j) arm_of[bb_leaf[static_cast<std::size_t>(j)]] = a;
      }
    }
    if (arms_bb.size() < 3) continue;
    const auto [a0, b0] = span(v);
    int r, ref;
    if (v != s.root && (a0 > 0 || b0 < static_cast<int>(bb_leaf.size()))) {
      r = a0 > 0 ? bb_leaf[static_cast<std::size_t>(a0 - 1)] : bb_leaf[static_cast<std::size_t>(b0)];
      arm_of[r] = -1;
      tot[-1] = 1;
      ref = -1;
    } else {
      ref = arms_bb.front();
      for (int a : arms_bb) {
        if (tot[a] > tot[ref]) ref = a;
      }
      r = -1;
      for (const auto& [a, x] : first_leaf) {
        if (a == ref) r = x;
      }
    }
    std::vector<int> groupable;
    for (int a : arms_bb) {
      if (a != ref) groupable.push_back(a);
    }
    auto found = polytomy_groups(lca, arm_of, tot, r, groupable);
    if (!found.empty()) groups_at.push_back({v, std::move(found)});
  }
  long long added = 0;
  for (auto& [v, found] : groups_at) {  // larger groups first; arms move down as groups nest
    const std::vector<int> arm_nodes = s.kids[static_cast<std::size_t>(v)];
    std::vector<int> cur(arm_nodes.size(), v);
    std::vector<std::vector<int>> groups(found.begin(), found.end());
    std::stable_sort(groups.begin(), groups.end(),
                     [](const std::vector<int>& x, const std::vector<int>& y) { return x.size() > y.size(); });
    for (const auto& g : groups) {
      const int p = cur[static_cast<std::size_t>(g.front())];
      for (int a : g) {
        if (cur[static_cast<std::size_t>(a)] != p) throw std::logic_error("guide groups of one polytomy are not nested");
      }
      const int u = static_cast<int>(s.kids.size());
      s.kids.emplace_back();
      s.taxon.push_back(-1);
      for (int a : g) {
        auto& pk = s.kids[static_cast<std::size_t>(p)];
        pk.erase(std::find(pk.begin(), pk.end(), arm_nodes[static_cast<std::size_t>(a)]));
        s.kids[static_cast<std::size_t>(u)].push_back(arm_nodes[static_cast<std::size_t>(a)]);
        cur[static_cast<std::size_t>(a)] = u;
      }
      s.kids[static_cast<std::size_t>(p)].push_back(u);
      ++added;
    }
  }
  return added;
}

}  // namespace

struct GuideFilter::Impl {
  std::vector<int> gleaf;
  std::optional<ConflictTest> test;
};

GuideFilter::GuideFilter(const Tree& guide_tree, const TaxonRegistry& taxa) : impl_(std::make_unique<Impl>()) {
  const Rooted guide = rooted_from(guide_tree, taxa);
  impl_->gleaf.assign(static_cast<std::size_t>(taxa.size()), -1);
  for (std::size_t v = 0; v < guide.kids.size(); ++v) {
    if (guide.taxon[v] >= 0) impl_->gleaf[static_cast<std::size_t>(guide.taxon[v])] = static_cast<int>(v);
  }
  std::vector<int> shared;
  for (int t = 0; t < taxa.size(); ++t) {
    if (impl_->gleaf[static_cast<std::size_t>(t)] >= 0) shared.push_back(t);
  }
  std::sort(shared.begin(), shared.end(), [&](int x, int y) { return taxa.label_for(x) < taxa.label_for(y); });
  if (shared.size() >= 4) impl_->test.emplace(guide, impl_->gleaf, shared);
}

GuideFilter::~GuideFilter() = default;

bool GuideFilter::conflicts(const std::vector<int>& cluster_taxa) const {
  if (!impl_->test) return false;
  std::vector<int> side;
  bool has_b0 = false;
  for (int t : cluster_taxa) {
    const int x = impl_->gleaf[static_cast<std::size_t>(t)];
    if (x < 0) continue;
    side.push_back(x);
    has_b0 = has_b0 || x == impl_->test->b0();
  }
  return impl_->test->conflicts(side, has_b0);
}

std::string starting_tree_newick(const StartingTree& tree, const TaxonRegistry& taxa) {
  return newick(rooted_from(tree), taxa);
}

StartingTree apply_guide_tree(const StartingTree& tree, const Tree& guide_tree, const TaxonRegistry& taxa) {
  const auto started = std::chrono::steady_clock::now();
  Rooted s = rooted_from(tree);
  const Rooted guide = rooted_from(guide_tree, taxa);
  std::vector<int> gleaf(static_cast<std::size_t>(taxa.size()), -1);
  for (std::size_t v = 0; v < guide.kids.size(); ++v) {
    const int t = guide.taxon[v];
    if (t < 0) continue;
    if (gleaf[static_cast<std::size_t>(t)] >= 0) {
      throw std::invalid_argument("duplicate taxon in guide tree: " + taxa.label_for(t));
    }
    gleaf[static_cast<std::size_t>(t)] = static_cast<int>(v);
  }
  std::vector<int> shared;  // ordered by label: b0 is the first
  for (int t = 0; t < taxa.size(); ++t) {
    if (gleaf[static_cast<std::size_t>(t)] >= 0) shared.push_back(t);
  }
  std::sort(shared.begin(), shared.end(),
            [&](int x, int y) { return taxa.label_for(x) < taxa.label_for(y); });
  long long contracted = 0, added = 0;
  if (shared.size() >= 4) {
    contracted = collapse_conflicts(s, guide, gleaf, shared);
    added = refine_polytomies(s, guide, gleaf);
  }
  stats_add(run_stats().guide_contracted_edges, contracted);
  stats_add(run_stats().guide_added_groups, added);
  stats_add(run_stats().guide_shared_taxa, static_cast<long long>(shared.size()));
  // Through the starting-tree reader, so the run sees the tree exactly as if
  // it were given with --starting-tree.
  StartingTree out = starting_tree_from(parse_newick(newick(s, taxa)), taxa);
  stats_add_time(run_stats().guide_ns,
                 std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - started));
  return out;
}

}  // namespace qraft::internal
