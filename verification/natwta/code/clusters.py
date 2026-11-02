"""
Clusters with boundary vertices (top-tree style), the HDT signature of
Brodal et al. (SODA'13) in cluster form, fragments, the fragment automaton,
and the exact computation of minimal counter dimensions.

Marked trees / forests are stored as (kinds, adj):
    kinds : list of node kinds; a kind is one of
            'T'  (top boundary vertex),  'B' (bottom boundary vertex),
            ('L', c) (leaf of colour c), 'I' (unmarked internal vertex)
    adj   : list of sets (undirected adjacency)
Everything is tiny (<= 12 vertices), so we do not optimise.
"""
from fractions import Fraction
import itertools, random

D = 4          # pattern size (quartets)

# ----------------------------------------------------------------------
# basic graph helpers
# ----------------------------------------------------------------------
def copy_graph(g):
    kinds, adj = g
    return (list(kinds), [set(a) for a in adj])

def is_leaf_kind(k):
    return isinstance(k, tuple) and k[0] == 'L'

def marked(k):
    return k == 'T' or k == 'B' or k == 'R' or is_leaf_kind(k)

def glue(g1, g2, pairs, newkind):
    """Disjoint union of g1, g2, then identify vertices u (in g1) with v (in g2)
    for each (u, v) in pairs; the identified vertex gets kind newkind[(u,v)]
    (or keeps u's kind if not given).  Returns (graph, map1, map2)."""
    k1, a1 = g1
    k2, a2 = g2
    n1 = len(k1)
    kinds = list(k1) + list(k2)
    adj = [set(a) for a in a1] + [set(x + n1 for x in a) for a in a2]
    # identification: redirect v -> u
    ident = {}
    for (u, v) in pairs:
        ident[v + n1] = u
        kinds[u] = newkind.get((u, v), kinds[u])
    def rep(x):
        return ident.get(x, x)
    m = len(kinds)
    for v in list(ident.keys()):
        u = ident[v]
        for w in adj[v]:
            adj[w].discard(v)
            adj[w].add(u)
            adj[u].add(w)
        adj[v] = set()
    # compact
    keep = [i for i in range(m) if i not in ident]
    newidx = {old: i for i, old in enumerate(keep)}
    nk = [kinds[i] for i in keep]
    na = [set(newidx[rep(w)] for w in adj[i]) for i in keep]
    for i in range(len(na)):
        na[i].discard(i)
    map1 = {i: newidx[i] for i in range(n1)}
    map2 = {j: newidx[rep(j + n1)] for j in range(len(k2))}
    return (nk, na), map1, map2

def relabel_kind(g, v, kind):
    g = copy_graph(g)
    g[0][v] = kind
    return g

def reduce_marked(g):
    """Restriction to the marked vertices: prune unmarked vertices of degree <= 1
    repeatedly, then suppress unmarked vertices of degree 2.  Works on forests
    (each component is treated the same way; components without marks vanish)."""
    kinds, adj = copy_graph(g)
    n = len(kinds)
    alive = [True] * n
    changed = True
    while changed:
        changed = False
        for v in range(n):
            if alive[v] and not marked(kinds[v]) and len(adj[v]) <= 1:
                for w in adj[v]:
                    adj[w].discard(v)
                adj[v] = set()
                alive[v] = False
                changed = True
    for v in range(n):
        if alive[v] and not marked(kinds[v]) and len(adj[v]) == 2:
            a, b = tuple(adj[v])
            adj[a].discard(v); adj[b].discard(v)
            adj[a].add(b); adj[b].add(a)
            adj[v] = set()
            alive[v] = False
    keep = [v for v in range(n) if alive[v]]
    idx = {v: i for i, v in enumerate(keep)}
    nk = [kinds[v] for v in keep]
    na = [set(idx[w] for w in adj[v]) for v in keep]
    return (nk, na)

def components(g):
    kinds, adj = g
    seen = [False] * len(kinds)
    comps = []
    for s in range(len(kinds)):
        if seen[s]:
            continue
        stack = [s]; seen[s] = True; comp = []
        while stack:
            v = stack.pop(); comp.append(v)
            for w in adj[v]:
                if not seen[w]:
                    seen[w] = True; stack.append(w)
        comps.append(sorted(comp))
    return comps

def kind_key(k):
    if k == 'T': return (0,)
    if k == 'B': return (1,)
    if k == 'I': return (2,)
    if k == 'R': return (4,)
    return (3, k[1])

def rooted_code(g, v, parent=None):
    kinds, adj = g
    return (kind_key(kinds[v]), tuple(sorted(rooted_code(g, w, v) for w in adj[v] if w != parent)))

def canon(g):
    """Canonical code of a marked forest.  Each component is rooted at its
    'T' vertex if present, else at its 'B' vertex, else (closed pattern) the
    minimum over all rootings.  The forest code is the sorted tuple of
    component codes.  The empty forest has code ()."""
    kinds, adj = g
    codes = []
    for comp in components(g):
        roots = [v for v in comp if kinds[v] == 'R']
        tops = [v for v in comp if kinds[v] == 'T']
        bots = [v for v in comp if kinds[v] == 'B']
        if roots:
            codes.append(rooted_code(g, roots[0]))
        elif tops:
            codes.append(rooted_code(g, tops[0]))
        elif bots:
            codes.append(rooted_code(g, bots[0]))
        else:
            codes.append(min(rooted_code(g, v) for v in comp))
    return tuple(sorted(codes))

def graph_from_code(code):
    """Inverse of canon (up to isomorphism)."""
    kinds, adj = [], []
    def build(c, parent):
        kk, children = c
        if kk[0] == 0: kind = 'T'
        elif kk[0] == 1: kind = 'B'
        elif kk[0] == 2: kind = 'I'
        elif kk[0] == 4: kind = 'R'
        else: kind = ('L', kk[1])
        v = len(kinds); kinds.append(kind); adj.append(set())
        if parent is not None:
            adj[v].add(parent); adj[parent].add(v)
        for ch in children:
            build(ch, v)
        return v
    for c in code:
        build(c, None)
    return (kinds, adj)

def degree(code):
    """number of leaves in a fragment code"""
    def cnt(c):
        kk, children = c
        return (1 if kk[0] == 3 else 0) + sum(cnt(ch) for ch in children)
    return sum(cnt(c) for c in code)

def colours_of(code):
    def col(c):
        kk, children = c
        s = set([kk[1]]) if kk[0] == 3 else set()
        for ch in children:
            s |= col(ch)
        return s
    s = set()
    for c in code:
        s |= col(c)
    return s

def pretty(code):
    """human readable code: T(...) with arms; leaves as colours; I as '(...)'"""
    def pr(c):
        kk, children = c
        inner = ",".join(pr(ch) for ch in children)
        if kk[0] == 0: return "T(" + inner + ")" if inner else "T"
        if kk[0] == 1: return "B(" + inner + ")" if inner else "B"
        if kk[0] == 2: return "(" + inner + ")"
        if kk[0] == 4: return "R(" + inner + ")" if inner else "R"
        return str(kk[1])
    return " | ".join(pr(c) for c in code) if code else "empty"

# ----------------------------------------------------------------------
# concrete clusters (for brute-force validation)
# ----------------------------------------------------------------------
class Cluster:
    """A concrete tree with boundary vertices; leaves carry weights."""
    def __init__(self, kinds, adj, weights):
        self.g = (list(kinds), [set(a) for a in adj])
        self.w = dict(weights)      # leaf vertex -> weight
    @property
    def sort(self):
        kinds = self.g[0]
        hasT = 'T' in kinds; hasB = 'B' in kinds
        if hasT and hasB: return 'C2'
        if hasT:
            return 'C1' if len(self.g[1][kinds.index('T')]) == 1 else 'G'
        return 'T'
    def find(self, kind):
        return self.g[0].index(kind)
    def leaves(self):
        return [v for v, k in enumerate(self.g[0]) if is_leaf_kind(k)]

def leaf_atom(c, w):
    return Cluster(['T', ('L', c)], [{1}, {0}], {1: w})

def edge_atom():
    return Cluster(['T', 'B'], [{1}, {0}], {})

def merge_clusters(A, B, pairs_kinds, relabel):
    """generic: pairs_kinds = list of (kindA, kindB, newkind); relabel = dict oldkind->newkind
    applied to A's and B's remaining boundary vertices (after gluing)."""
    pairs = []; newkind = {}
    for ka, kb, nk in pairs_kinds:
        u = A.find(ka); v = B.find(kb)
        pairs.append((u, v)); newkind[(u, v)] = nk
    g, m1, m2 = glue(A.g, B.g, pairs, newkind)
    kinds, adj = g
    for (side, old, new) in relabel:
        src = A if side == 'A' else B
        m = m1 if side == 'A' else m2
        v = src.find(old)
        if (side == 'A' and any(u == v for (u, _) in pairs)) or (side == 'B' and any(vv == v for (_, vv) in pairs)):
            continue
        kinds[m[v]] = new
    weights = {}
    for v, w in A.w.items(): weights[m1[v]] = w
    for v, w in B.w.items(): weights[m2[v]] = w
    return Cluster(kinds, adj, weights)

def op_sib(A, B):    # P x P -> P : identify tops, top stays
    return merge_clusters(A, B, [('T', 'T', 'T')], [])
def op_hang(A, B):   # L x P -> L : A.bot = B.top, stays boundary (bottom)
    return merge_clusters(A, B, [('B', 'T', 'B')], [])
def op_cat(A, B):    # L x L -> L : A.bot = B.top closes; result top = A.top, bottom = B.bot
    return merge_clusters(A, B, [('B', 'T', 'I')], [])
def op_cap(A, B):    # L x P -> P : A.bot = B.top closes
    return merge_clusters(A, B, [('B', 'T', 'I')], [])
def op_close(A):     # P -> C : top becomes internal
    C = Cluster(*A.g, A.w)
    C.g[0][A.find('T')] = 'I'
    return C
def op_up_open(G):   # P -> L : attach the parent edge of the top p; p becomes the bottom (one pending edge)
    return op_hang(edge_atom(), G)
def op_up_closed(G): # P -> P : attach the parent edge of the top p; p closes
    return op_cap(edge_atom(), G)

# The HDT signature of Brodal et al. in cluster form.  'hang' (general top-tree
# merge keeping the shared vertex open) is NOT part of it; it is kept only for the
# comparison remark (general top trees).
def op_toG(A):       # C1 -> G : retyping (identity on clusters)
    return Cluster(*A.g, A.w)

# Sorts: T (closed tree), G (point cluster, any number of arms at the top = SODA'13 G),
# C1 (point cluster with exactly one arm at the top = downward-closed C, incl. leaves),
# C2 (path cluster; top has one arm, bottom has exactly one pending edge = downward-open C).
OPS = {  # name: (arity, input sorts, output sort, function)
    'toG':       (1, ('C1',), 'G', op_toG),
    'sib':       (2, ('G', 'G'), 'G', op_sib),
    'up_open':   (1, ('G',), 'C2', op_up_open),
    'up_closed': (1, ('G',), 'C1', op_up_closed),
    'cat':       (2, ('C2', 'C2'), 'C2', op_cat),
    'cap':       (2, ('C2', 'C1'), 'C1', op_cap),
    'close':     (1, ('G',), 'T', op_close),
}
GENERAL_TOPTREE = False   # if True, completions allow arbitrary material at the bottom (hang)

def fragment(A, X):
    """fragment of selection X (set of leaf vertices) in cluster A: restriction
    to X together with the boundary vertices; returns a canonical code."""
    kinds = list(A.g[0])
    for v in A.leaves():
        if v not in X:
            kinds[v] = 'I'          # unselected leaf becomes unmarked (pruned)
    return canon(reduce_marked((kinds, A.g[1])))

def fragment_vector(A, d=D):
    """h(A): dict code -> total weight of selections of size <= d"""
    h = {}
    lv = A.leaves()
    for r in range(0, d + 1):
        for X in itertools.combinations(lv, r):
            w = Fraction(1)
            for x in X:
                w *= A.w[x]
            s = fragment(A, set(X))
            h[s] = h.get(s, 0) + w
    return {s: v for s, v in h.items() if v != 0}

# ----------------------------------------------------------------------
# the fragment automaton (transitions on shapes, computed by glue + reduce)
# ----------------------------------------------------------------------
def shape_graph(code):
    return graph_from_code(code)

def find_kind(g, kind):
    return g[0].index(kind)

def ghat(op, codes):
    """transition on shapes: returns the output shape code (or None if degree > D)"""
    if sum(degree(c) for c in codes) > D:
        return None
    if op == 'toG':
        return codes[0]
    if op == 'close':
        g = graph_from_code(codes[0])
        g = relabel_kind(g, find_kind(g, 'T'), 'I')
        return canon(reduce_marked(g))
    if op in ('up_open', 'up_closed'):
        gG = graph_from_code(codes[0])
        gE = (['T', 'B'], [{1}, {0}])
        u = 1; v = find_kind(gG, 'T')
        g, _, _ = glue(gE, gG, [(u, v)], {(u, v): 'B' if op == 'up_open' else 'I'})
        return canon(reduce_marked(g))
    gA = graph_from_code(codes[0]); gB = graph_from_code(codes[1])
    if op == 'sib':
        u = find_kind(gA, 'T'); v = find_kind(gB, 'T')
        g, _, _ = glue(gA, gB, [(u, v)], {(u, v): 'T'})
    elif op == 'hang':
        u = find_kind(gA, 'B'); v = find_kind(gB, 'T')
        g, _, _ = glue(gA, gB, [(u, v)], {(u, v): 'B'})
    elif op in ('cat', 'cap', 'capG'):
        u = find_kind(gA, 'B'); v = find_kind(gB, 'T')
        g, _, _ = glue(gA, gB, [(u, v)], {(u, v): 'I'})
    else:
        raise ValueError(op)
    return canon(reduce_marked(g))

def apply_op(op, states):
    """multilinear map on state dicts"""
    out = {}
    if len(states) == 1:
        for s, x in states[0].items():
            t = ghat(op, (s,))
            if t is not None:
                out[t] = out.get(t, 0) + x
    else:
        for s1, x1 in states[0].items():
            for s2, x2 in states[1].items():
                t = ghat(op, (s1, s2))
                if t is not None:
                    out[t] = out.get(t, 0) + x1 * x2
    return {s: v for s, v in out.items() if v != 0}

def atom_state(kind, c=None, w=None):
    if kind == 'leaf':
        empty = canon((['T'], [set()]))
        one = canon((['T', ('L', c)], [{1}, {0}]))
        st = {empty: Fraction(1)}
        if w != 0:
            st[one] = Fraction(w)
        return st
    if kind == 'edge':
        return {canon((['T', 'B'], [{1}, {0}])): Fraction(1)}

# ----------------------------------------------------------------------
# random terms, brute-force validation
# ----------------------------------------------------------------------
def random_term(sort, K, rng, depth=0, weights=None):
    """returns (term, cluster, state) for a random term of the given sort.
    weights: callable giving leaf weight."""
    if weights is None:
        weights = lambda: Fraction(rng.randint(1, 3))
    stop = rng.random() < 0.35 + 0.1 * depth
    if sort == 'C1':
        if stop or depth > 6:
            c = rng.randrange(K); w = weights()
            return ('leaf', c, w), leaf_atom(c, w), atom_state('leaf', c, w)
        if rng.random() < 0.5:
            t1, A, s1 = random_term('G', K, rng, depth + 1, weights)
            return ('up_closed', t1), op_up_closed(A), apply_op('up_closed', (s1,))
        t1, A, s1 = random_term('C2', K, rng, depth + 1, weights)
        t2, B, s2 = random_term('C1', K, rng, depth + 1, weights)
        return ('cap', t1, t2), op_cap(A, B), apply_op('cap', (s1, s2))
    if sort == 'G':
        if stop or depth > 6:
            t1, A, s1 = random_term('C1', K, rng, depth + 1, weights)
            return ('toG', t1), op_toG(A), apply_op('toG', (s1,))
        t1, A, s1 = random_term('G', K, rng, depth + 1, weights)
        t2, B, s2 = random_term('G', K, rng, depth + 1, weights)
        return ('sib', t1, t2), op_sib(A, B), apply_op('sib', (s1, s2))
    if sort == 'C2':
        if stop or depth > 6:
            t1, A, s1 = random_term('G', K, rng, depth + 1, weights)
            return ('up_open', t1), op_up_open(A), apply_op('up_open', (s1,))
        t1, A, s1 = random_term('C2', K, rng, depth + 1, weights)
        t2, B, s2 = random_term('C2', K, rng, depth + 1, weights)
        return ('cat', t1, t2), op_cat(A, B), apply_op('cat', (s1, s2))
    if sort == 'T':
        t1, A, s1 = random_term('G', K, rng, depth, weights)
        return ('close', t1), op_close(A), apply_op('close', (s1,))
    raise ValueError(sort)

def validate(K=3, trials=200, seed=1):
    """the fragment automaton computes the fragment vector of every term value,
    and the cluster has the declared sort"""
    rng = random.Random(seed)
    for _ in range(trials):
        sort = rng.choice(['C1', 'G', 'C2', 'T'])
        t, A, s = random_term(sort, K, rng)
        assert A.sort == sort or (sort == 'G' and A.sort == 'C1'), (sort, A.sort, t)
        h = fragment_vector(A)
        assert h == s, (t, h, s)
    return True

# ----------------------------------------------------------------------
# enumeration of fragments and completions (abstract generators)
# ----------------------------------------------------------------------
def arms(r, K):
    """all rooted 'arm' trees with exactly r coloured leaves, as codes rooted at
    the arm root; an arm root is a leaf (r=1) or an unmarked vertex with >= 2 arms."""
    if r == 1:
        return [((3, c), ()) for c in range(K)]
    out = set()
    for ms in multisets_of_arms(r, K, minparts=2, maxpart=r - 1):
        out.add(((2,), tuple(sorted(ms))))
    return sorted(out)

_arm_cache = {}
def arms_cached(r, K):
    key = (r, K)
    if key not in _arm_cache:
        _arm_cache[key] = arms(r, K)
    return _arm_cache[key]

def multisets_of_arms(r, K, minparts, maxpart=None):
    """multisets of arms (codes) with total leaf count r and at least minparts arms;
    every arm has at most maxpart leaves"""
    if maxpart is None:
        maxpart = r
    results = set()
    def rec(remaining, minpart_size, chosen):
        if remaining == 0:
            if len(chosen) >= minparts:
                results.add(tuple(sorted(chosen)))
            return
        for size in range(minpart_size, min(remaining, maxpart) + 1):
            for a in arms_cached(size, K):
                rec(remaining - size, size, chosen + [a])
    rec(r, 1, [])
    return results

def marked_trees(mark, r, K, mindeg=0):
    """codes of trees with one marked vertex `mark` (0='T', 1='B') of degree >= mindeg
    carrying a multiset of arms with r leaves in total."""
    out = set()
    for ms in multisets_of_arms(r, K, minparts=mindeg):
        out.add(((mark,), tuple(sorted(ms))))
    if r == 0 and mindeg == 0:
        out.add(((mark,), ()))
    return sorted(out)

def fragments(sort, r, K):
    """fragment shapes of the given sort with exactly r leaves.
    G : one marked vertex (top) of any degree carrying arms;
    C1: as G but the top has exactly one arm (degree 0 only for r = 0);
    C2: top of degree 1, path top - u_1 - ... - u_k - bottom, each u_i with >= 1 arm,
        bottom with >= 0 arms;
    T : closed patterns."""
    if sort == 'G':
        return [(t,) for t in marked_trees(0, r, K)]
    if sort == 'C1':
        return [(t,) for t in marked_trees(0, r, K) if len(t[1]) <= 1]
    if sort == 'T':
        out = set()
        for t in marked_trees(0, r, K):
            g = graph_from_code((t,))
            g = relabel_kind(g, find_kind(g, 'T'), 'I')
            out.add(canon(reduce_marked(g)))
        return sorted(out)
    if sort == 'C2':
        out = set()
        def seqs(remaining):
            res = []
            def rec2(rem, interior):
                for bms in (multisets_of_arms(rem, K, minparts=0) if rem > 0 else [()]):
                    res.append((list(interior), bms))
                for use in range(1, rem + 1):
                    for ims in multisets_of_arms(use, K, minparts=1):
                        rec2(rem - use, interior + [ims])
            rec2(remaining, [])
            return res
        for interior, bms in seqs(r):
            node = ((1,), tuple(sorted(bms)))
            for ims in reversed(interior):
                node = ((2,), tuple(sorted(list(ims) + [node])))
            out.add((((0,), (node,)),))
        return sorted(out)
    raise ValueError(sort)

def rooted_completion_trees(r, K, hole_kind=0):
    """Completions above a boundary vertex v (kind hole_kind: 0='T') when the output
    sort is a rooted G (root R marked): rooted trees at R containing v.  Either R = v
    (then a single vertex of kind R carrying v's other arms), or R is above v: a path
    R - u_1 - ... - u_k - v whose interior vertices carry >= 1 arm, R carrying >= 0
    further arms, v carrying >= 0 arms; r leaves in total.  Codes are rooted at R."""
    out = set()
    for t in marked_trees(4, r, K):             # case R = v
        out.add(t)
    for rv in range(r + 1):                     # case R above v; rv leaves in v's arms
        for vms in (multisets_of_arms(rv, K, minparts=0) if rv > 0 else [()]):
            vnode = ((hole_kind,), tuple(sorted(vms)))
            def chains(rem, node):
                res = []
                for rms in (multisets_of_arms(rem, K, minparts=0) if rem > 0 else [()]):
                    res.append(((4,), tuple(sorted(list(rms) + [node]))))
                for use in range(1, rem + 1):
                    for ims in multisets_of_arms(use, K, minparts=1):
                        res += chains(rem - use, ((2,), tuple(sorted(list(ims) + [node]))))
                return res
            for code in chains(r - rv, vnode):
                out.add(code)
    return sorted(out)

def completions(sort, r, K, out='T'):
    """completions for a hole of the given sort with exactly r context leaves, for
    a task whose output sort is `out`:
      out = 'T': the tree is closed at the end (unrooted patterns);
      out = 'G': the final component is the root point cluster (rooted patterns, root R).
    G, C1: material at the top (any degree);
    C2   : material at the top and, at the bottom, a single arm or nothing
           (the bottom has exactly one pending edge) -- with GENERAL_TOPTREE any tree."""
    if sort == 'T':
        return [()] if (r == 0 and out == 'T') else []
    if out == 'T':
        tops = lambda rr: marked_trees(0, rr, K)
    else:
        tops = lambda rr: rooted_completion_trees(rr, K)
    if sort in ('G', 'C1'):
        return [(t,) for t in tops(r)]
    if sort == 'C2':
        res = set()
        for r1 in range(r + 1):
            for t1 in tops(r1):
                for t2 in marked_trees(1, r - r1, K):
                    if not GENERAL_TOPTREE and len(t2[1]) > 1:
                        continue
                    res.add(canon(graph_from_code((t1, t2))))
        return sorted(res)
    raise ValueError(sort)

def complete(frag, comp):
    """glue fragment and completion at the boundary vertices (top with top-or-root,
    bottom with bottom), unmark them, keep the root mark R if present, and reduce.
    Returns the code of the resulting (rooted or unrooted) pattern."""
    gF = graph_from_code(frag); gC = graph_from_code(comp)
    pairs = []; newkind = {}
    kindsC = gC[0]
    if 'T' in gF[0]:
        u = gF[0].index('T')
        if 'T' in kindsC:
            v = kindsC.index('T'); pairs.append((u, v)); newkind[(u, v)] = 'I'
        else:                       # completion has R = v
            v = kindsC.index('R'); pairs.append((u, v)); newkind[(u, v)] = 'R'
    if 'B' in gF[0]:
        u = gF[0].index('B'); v = kindsC.index('B')
        pairs.append((u, v)); newkind[(u, v)] = 'I'
    g, _, _ = glue(gF, gC, pairs, newkind)
    return canon(reduce_marked(g))

# ----------------------------------------------------------------------
# closed patterns and channels
# ----------------------------------------------------------------------
def split_pattern(a, b, c, d):
    """closed pattern ab|cd (code)"""
    kinds = ['I', 'I', ('L', a), ('L', b), ('L', c), ('L', d)]
    adj = [{1, 2, 3}, {0, 4, 5}, {0}, {0}, {1}, {1}]
    return canon((kinds, adj))

def star_pattern(a, b, c, d):
    kinds = ['I', ('L', a), ('L', b), ('L', c), ('L', d)]
    adj = [{1, 2, 3, 4}, {0}, {0}, {0}, {0}]
    return canon((kinds, adj))

def channel_from_patterns(pats):
    ch = {}
    for p in pats:
        ch[p] = 1
    return ch

def qraft_channels(K):
    if K == 2:
        return {'00|11': channel_from_patterns([split_pattern(0, 0, 1, 1)]),
                '01|01': channel_from_patterns([split_pattern(0, 1, 0, 1)])}
    if K == 4:
        return {'01|23': channel_from_patterns([split_pattern(0, 1, 2, 3)]),
                '02|13': channel_from_patterns([split_pattern(0, 2, 1, 3)]),
                '03|12': channel_from_patterns([split_pattern(0, 3, 1, 2)])}

def soda_channels(x, which):
    """SODA'13 families with colours 0 (outside) and 1..x (children of v).
    A: shared resolved (ii|jj, ii|jk, ii|0j); B: disagreeing resolved (ij|ij, ij|ik, i0|ij);
    E: unresolved in both (star ijkl, star 0ijk); D: star in T1, resolved in T2 (ij|kl, 0i|jk...)."""
    nz = range(1, x + 1)
    pats = set()
    if which == 'A':
        for i in nz:
            for j in nz:
                if j != i: pats.add(split_pattern(i, i, j, j)); pats.add(split_pattern(i, i, 0, j))
                for k in nz:
                    if len({i, j, k}) == 3: pats.add(split_pattern(i, i, j, k))
    elif which == 'B':
        for i in nz:
            for j in nz:
                if j != i: pats.add(split_pattern(i, j, i, j)); pats.add(split_pattern(i, 0, i, j))
                for k in nz:
                    if len({i, j, k}) == 3: pats.add(split_pattern(i, j, i, k))
    elif which == 'E':
        for i, j, k, l in itertools.combinations(nz, 4):
            pats.add(star_pattern(i, j, k, l))
        for i, j, k in itertools.combinations(nz, 3):
            pats.add(star_pattern(0, i, j, k))
    elif which == 'C':   # resolved in T1 (anchored at v), unresolved in T2
        for i in nz:
            for j in nz:
                if j != i: pats.add(star_pattern(i, i, j, j)); pats.add(star_pattern(i, i, 0, j))
                for k in nz:
                    if len({i, j, k}) == 3: pats.add(star_pattern(i, i, j, k))
    elif which == 'D':   # unresolved in T1 (4 distinct directions), resolved in T2
        for i, j, k, l in itertools.combinations(nz, 4):
            for (a, b, c, d) in ((i, j, k, l), (i, k, j, l), (i, l, j, k)):
                pats.add(split_pattern(a, b, c, d))
        for i, j, k in itertools.combinations(nz, 3):
            for (a, b, c, d) in ((0, i, j, k), (0, j, i, k), (0, k, i, j)):
                pats.add(split_pattern(a, b, c, d))
    else:
        raise ValueError(which)
    return channel_from_patterns(pats)

# ----------------------------------------------------------------------
# exact linear algebra
# ----------------------------------------------------------------------
def rank_mod(rows, p):
    """rank over GF(p) of a list of integer rows (lists)"""
    M = [[x % p for x in r] for r in rows]
    if not M:
        return 0
    ncols = len(M[0]); rank = 0; r = 0
    for c in range(ncols):
        piv = None
        for i in range(r, len(M)):
            if M[i][c]:
                piv = i; break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        r += 1
        if r == len(M):
            break
    return r

def rref_mod(rows, p):
    """reduced row echelon basis over GF(p)"""
    M = [[x % p for x in r] for r in rows]
    if not M:
        return []
    ncols = len(M[0]); r = 0
    for c in range(ncols):
        piv = None
        for i in range(r, len(M)):
            if M[i][c]:
                piv = i; break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        r += 1
        if r == len(M):
            break
    return M[:r]

def rank_exact(rows):
    """rank over Q (Fraction Gaussian elimination)"""
    M = [[Fraction(x) for x in r] for r in rows]
    if not M:
        return 0
    ncols = len(M[0]); r = 0
    for c in range(ncols):
        piv = None
        for i in range(r, len(M)):
            if M[i][c] != 0:
                piv = i; break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = 1 / M[r][c]
        M[r] = [x * inv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
        if r == len(M):
            break
    return r

# ----------------------------------------------------------------------
# minimal dimensions via the completion matrix
# ----------------------------------------------------------------------
P1 = (1 << 31) - 1
P2 = 999983

def completion_matrix(sort, m, K, channels, out='T'):
    """rows: (completion kappa of degree D-m, channel j); columns: fragments of degree m.
    entry = channel_j(pattern of fragment completed by kappa)."""
    frags = fragments(sort, m, K)
    comps = completions(sort, D - m, K, out)
    rows = []
    rowlabels = []
    for kappa in comps:
        pats = [complete(s, kappa) for s in frags]
        for name, ch in channels.items():
            rows.append([ch.get(p, 0) for p in pats])
            rowlabels.append((kappa, name))
    return frags, rows, rowlabels

def minimal_dims(sort, K, channels, exact=True, p=P1, out='T'):
    """dict m -> (#fragments, rank), over Q by default.

    exact=False computes rank over GF(p), only a lower bound for rank over Q.
    Matching answers at several primes do not certify the rational rank.
    """
    res = {}
    for m in range(D + 1):
        frags, rows, _ = completion_matrix(sort, m, K, channels, out)
        if exact:
            rk = rank_exact(rows)
        else:
            rk = rank_mod(rows, p)
        res[m] = (len(frags), rk)
    return res

def packed_width(sort, K, channels, present, p=P1, out='T', exact=True):
    """dimension of the observable space restricted to fragments whose colours
    lie in `present` (a set of colours), summed over degrees.

    Contexts still use the ambient palette of size K. This need not equal a
    fresh minimum computed with K=len(present).
    """
    total = 0
    for m in range(D + 1):
        frags, rows, _ = completion_matrix(sort, m, K, channels, out)
        cols = [i for i, s in enumerate(frags) if colours_of(s) <= set(present)]
        sub = [[r[i] for i in cols] for r in rows]
        total += (rank_exact(sub) if exact else rank_mod(sub, p)) if cols else 0
    return total

# ----------------------------------------------------------------------
# Brodal et al.'s per-run readout: rooted quartet forms
# ----------------------------------------------------------------------
from collections import Counter

def rooted_patterns(K):
    """all rooted 4-leaf patterns: fragments of the root G of degree 4, with the top
    re-labelled as the root mark R"""
    out = []
    for (t,) in fragments('G', 4, K):
        out.append((((4,), t[1]),))
    return out

def rooted_quartet_info(code):
    """for a rooted 4-leaf pattern (root R): returns (form, split) where form is one of
    'alpha' (root has two cherry children), 'beta' (cherry + two leaves), 'gamma'
    (leaf + ((xx),y)), 'star4' (four leaves), 'star31' (leaf + 3-star) and split is
    a frozenset of two colour tuples or None."""
    (kk, children), = code
    node = (kk, children)
    while len(node[1]) == 1:          # root above the LCA of the four leaves
        node = node[1][0]
    kids = node[1]
    def leaves(c):
        k, ch = c
        return [k[1]] if k[0] == 3 else sum((leaves(x) for x in ch), [])
    sizes = sorted(((len(leaves(c)), c) for c in kids), key=lambda t: t[0])
    sz = [a for a, _ in sizes]
    if sz == [2, 2]:
        a, b = leaves(sizes[0][1]), leaves(sizes[1][1])
        return 'alpha', (tuple(sorted(a)), tuple(sorted(b)))
    if sz == [1, 1, 2]:
        a = [leaves(sizes[0][1])[0], leaves(sizes[1][1])[0]]; b = leaves(sizes[2][1])
        return 'beta', (tuple(sorted(a)), tuple(sorted(b)))
    if sz == [1, 3]:
        single = leaves(sizes[0][1])[0]
        sub = sizes[1][1]
        subkids = sub[1]
        if len(subkids) == 2:            # ((xx), y)
            s2 = sorted(((len(leaves(c)), c) for c in subkids), key=lambda t: t[0])
            y = leaves(s2[0][1])[0]; xx = leaves(s2[1][1])
            return 'gamma', (tuple(sorted(xx)), tuple(sorted([single, y])))
        return 'star31', None
    if sz == [1, 1, 1, 1]:
        return 'star4', None
    raise ValueError(code)

ORDER = {'alpha': 0, 'beta': 1, 'gamma': 2}

def t1_type(colours):
    """T1 form of a resolved quartet anchored at v, from its colour multiset (0 = outside):
    alpha {i,i,j,j}; beta {i,i,j,k}; gamma {i,i,0,j}; None otherwise."""
    c = Counter(colours)
    nz = [k for k in c if k != 0]
    if c.get(0, 0) == 0:
        if sorted(c.values()) == [2, 2]:
            return 'alpha'
        if sorted(c.values()) == [1, 1, 2]:
            return 'beta'
        return None
    if c[0] == 1 and sorted(c[k] for k in nz) == [1, 2]:
        return 'gamma'
    return None

def brodal_run_channel(x, which, cells='lower'):
    """per-run readout on rooted patterns for family which in {'A','B'}:
    counts quartets whose (T1 form, rooted T2 form) lies in the lower triangle
    (T2 form <= T1 form); the upper triangle is obtained by swapping the trees.
    cells='all' gives the full family (all nine cells)."""
    K = x + 1
    ch = {}
    for code in rooted_patterns(K):
        form, split = rooted_quartet_info(code)
        if split is None:
            continue
        cols = [c for side in split for c in side]
        t1 = t1_type(cols)
        if t1 is None:
            continue
        # 'same_side': a doubled nonzero colour sits on one side of the split (A);
        # otherwise the two equally coloured leaves are separated (B)
        same_side = any(side[0] == side[1] and side[0] != 0 for side in split)
        if which == 'A' and not same_side:
            continue
        if which == 'B' and same_side:
            continue
        if cells == 'lower' and ORDER[form] > ORDER[t1]:
            continue
        ch[code] = 1
    return ch
