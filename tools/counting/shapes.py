#!/usr/bin/env python3
"""Canonical component shapes and the exact reference counting.

A shape is a coloured tree fragment with at most four selected leaves and the
boundary vertices of its component sort (top port U, bottom port D, super
root Z); compositions glue shapes at a boundary vertex and restrict. This
module enumerates the canonical shapes of every sort, composes them, and
checks three facts that the counter tables rest on:

  1. Exactness: counting through shape compositions equals brute-force
     weighted-quartet enumeration on binary and polytomous trees, with
     integer and fractional weights, along two independent decompositions.
  2. Finiteness and completeness: each sort has a constant number of shapes,
     independent of tree size and degree, and composition closes exactly on
     an independent enumeration of every legal reduced shape.
  3. Minimality and reproducibility: an exact pullback closure gives the
     minimal linear counters, and a canonical RREF fixes their coordinates
     independently of Python's hash order.

Run:  python3 tools/counting/shapes.py
"""
from __future__ import annotations
import heapq, itertools, random, sys
from collections import defaultdict, deque
from fractions import Fraction

sys.setrecursionlimit(100000)

# ---------------------------------------------------------------------------
# Graph = colored tree fragment.  nodes: id -> ('leaf',color)|('U',)|('D',)|
# ('Z',)|('I',).  adj: id -> set(id).  Ports U/D/Z are unique per component.
# ---------------------------------------------------------------------------

LEAF, U, D, Z, I = 'leaf', 'U', 'D', 'Z', 'I'
PORTS = (U, D, Z)


class Frag:
    __slots__ = ('kind', 'adj')

    def __init__(self):
        self.kind = {}          # id -> tuple token
        self.adj = defaultdict(set)

    def copy(self):
        f = Frag()
        f.kind = dict(self.kind)
        f.adj = defaultdict(set, {k: set(v) for k, v in self.adj.items()})
        return f

    def add(self, nid, token):
        self.kind[nid] = token
        self.adj.setdefault(nid, set())

    def link(self, a, b):
        self.adj[a].add(b)
        self.adj[b].add(a)

    def unlink(self, a, b):
        self.adj[a].discard(b)
        self.adj[b].discard(a)

    def remove(self, nid):
        for nb in list(self.adj[nid]):
            self.adj[nb].discard(nid)
        del self.adj[nid]
        del self.kind[nid]


def _relabel(frag, base):
    """Return a copy with ids shifted by `base` (to make two frags disjoint)."""
    f = Frag()
    for nid, tok in frag.kind.items():
        f.kind[nid + base] = tok
    for nid, nbs in frag.adj.items():
        f.adj[nid + base] = set(n + base for n in nbs)
    return f


def is_marked(tok):
    return tok[0] in (LEAF, U, D, Z)


# ---------------------------------------------------------------------------
# reduce_and_code: minimal subtree over marked (leaves+ports), suppress
# unlabeled degree-2 I nodes, canonicalize (AHU).  Returns (code, canon_frag).
# ---------------------------------------------------------------------------

def reduce_and_code(frag):
    frag = frag.copy()
    # 1/2/3: prune non-marked degree-1 nodes (leaf-strip to the Steiner subtree).
    changed = True
    while changed:
        changed = False
        for nid in list(frag.kind):
            if is_marked(frag.kind[nid]):
                continue
            if len(frag.adj[nid]) <= 1:
                frag.remove(nid)
                changed = True
    # 4: suppress unlabeled ('I') degree-2 nodes (Z is protected).
    changed = True
    while changed:
        changed = False
        for nid in list(frag.kind):
            if frag.kind[nid][0] != I:
                continue
            nbs = list(frag.adj[nid])
            if len(nbs) == 2:
                a, b = nbs
                frag.remove(nid)
                if a != b:
                    frag.link(a, b)
                changed = True
            elif len(nbs) < 2:
                # degree-0/1 unlabeled leftover -> drop (leaf-less stub).
                frag.remove(nid)
                changed = True
    if not frag.kind:
        return ((('EMPTY',),), frag)
    code, canon = _canonical(frag)
    return (code, canon)


def _rooted_code(frag, v, parent):
    kids = sorted(_rooted_code(frag, c, v) for c in frag.adj[v] if c != parent)
    return (frag.kind[v],) + tuple(kids)


def _tree_centers(frag):
    # iterative leaf removal; last 1 or 2 nodes are the center(s).
    deg = {n: len(frag.adj[n]) for n in frag.kind}
    remaining = set(frag.kind)
    if len(remaining) <= 2:
        return list(remaining)
    leaves = deque(n for n in remaining if deg[n] <= 1)
    cnt = len(remaining)
    while cnt > 2:
        for _ in range(len(leaves)):
            leaf = leaves.popleft()
            remaining.discard(leaf)
            cnt -= 1
            for nb in frag.adj[leaf]:
                if nb in remaining:
                    deg[nb] -= 1
                    if deg[nb] == 1:
                        leaves.append(nb)
    return list(remaining)


def _canonical(frag):
    # choose root: U > D > Z > tree center(s) (min rooted code).
    root_candidates = []
    for pk in (U, D, Z):
        found = [n for n, t in frag.kind.items() if t[0] == pk]
        if found:
            root_candidates = found
            break
    if not root_candidates:
        root_candidates = _tree_centers(frag)
    best = None
    best_root = None
    for r in root_candidates:
        code = _rooted_code(frag, r, None)
        if best is None or code < best:
            best = code
            best_root = r
    # build a fresh canonical frag by DFS preorder following sorted children.
    canon = Frag()
    counter = itertools.count()
    def build(v, parent, pnew):
        nid = next(counter)
        canon.add(nid, frag.kind[v])
        if pnew is not None:
            canon.link(pnew, nid)
        kids = sorted(frag.adj[v] - ({parent} if parent is not None else set()),
                      key=lambda c: _rooted_code(frag, c, v))
        for c in kids:
            build(c, v, nid)
    build(best_root, None, None)
    return best, canon


def sel_count(code):
    # code = (token_tuple, *child_codes); token kind = code[0][0].
    kind = code[0][0]
    n = 1 if kind == LEAF else 0
    for ch in code[1:]:
        n += sel_count(ch)
    return n


def interface_of(code):
    # which ports present in the code tokens.
    toks = set()
    def walk(c):
        k = c[0][0]
        if k in (U, D, Z):
            toks.add(k)
        for ch in c[1:]:
            walk(ch)
    walk(code)
    if Z in toks:
        return 'G'
    if U in toks and D in toks:
        return 'C2'
    if U in toks or D in toks:
        return 'C1'
    return 'Closed'


def state_sort_key(code):
    """Stable state order used by registries, exact bases, and emitted tables."""
    return (sel_count(code), code)


# ---------------------------------------------------------------------------
# Atoms and operations.  Each returns (code, canon_frag).
# ---------------------------------------------------------------------------

def atom_leaf(color):
    """LEAF_TO_C: two states (absent -> U only; present -> U--leaf)."""
    a = Frag(); a.add(0, (U,))
    absent = reduce_and_code(a)
    p = Frag(); p.add(0, (U,)); p.add(1, (LEAF, color)); p.link(0, 1)
    present = reduce_and_code(p)
    return absent, present


def op_c_to_g(comp):
    """Relabel the (unique) U port as center Z."""
    _, frag = comp
    frag = frag.copy()
    for nid, tok in list(frag.kind.items()):
        if tok[0] == U:
            frag.kind[nid] = (Z,)
    return reduce_and_code(frag)


def _merge_ports(fa, fb, kinda, kindb):
    """Disjoint-union fa,fb; return combined frag + the port-node ids of each."""
    base = (max(fa.kind) + 1) if fa.kind else 0
    fb2 = _relabel(fb, base + 1)
    comb = Frag()
    for f in (fa, fb2):
        for nid, tok in f.kind.items():
            comb.add(nid, tok)
        for nid, nbs in f.adj.items():
            for n in nbs:
                comb.adj[nid].add(n)
    pa = [n for n, t in fa.kind.items() if t[0] == kinda]
    pb = [n + base + 1 for n, t in fb.kind.items() if t[0] == kindb]
    return comb, pa, pb


def op_gg_to_g(compa, compb):
    """Identify the two Z centers into one vertex (NOT an edge)."""
    _, fa = compa
    _, fb = compb
    comb, za, zb = _merge_ports(fa, fb, Z, Z)
    za, zb = za[0], zb[0]
    # merge zb into za.
    for nb in list(comb.adj[zb]):
        comb.adj[za].add(nb)
        comb.adj[nb].discard(zb)
        comb.adj[nb].add(za)
    comb.adj[za].discard(za)
    comb.remove(zb)
    comb.kind[za] = (Z,)
    return reduce_and_code(comb)


def op_ig_to_c(comp_g, mask_u, mask_d):
    """Identify I-template source vertex with G.Z, drop Z label, attach ports."""
    _, fg = comp_g
    frag = fg.copy()
    zid = [n for n, t in frag.kind.items() if t[0] == Z][0]
    frag.kind[zid] = (I,)
    nxt = (max(frag.kind) + 1)
    if mask_u:
        frag.add(nxt, (U,)); frag.link(zid, nxt); nxt += 1
    if mask_d:
        frag.add(nxt, (D,)); frag.link(zid, nxt); nxt += 1
    return reduce_and_code(frag)


def op_close_root(comp):
    """CLOSE_ROOT: delete the U port -> Closed."""
    _, frag = comp
    frag = frag.copy()
    for nid, tok in list(frag.kind.items()):
        if tok[0] == U:
            frag.remove(nid)
    return reduce_and_code(frag)


def op_reverse_c(comp):
    """REVERSE_C: swap U<->D."""
    _, frag = comp
    frag = frag.copy()
    for nid, tok in list(frag.kind.items()):
        if tok[0] == U:
            frag.kind[nid] = (D,)
        elif tok[0] == D:
            frag.kind[nid] = (U,)
    return reduce_and_code(frag)


def op_cc_to_c(higher, lower):
    """EDGE_JOIN(higher.D, lower.U): consume the shared edge, connect neighbors.
    Output U=higher.U, D=lower.D."""
    _, fh = higher
    _, fl = lower
    comb, dports, uports = _merge_ports(fh, fl, D, U)
    dnode, unode = dports[0], uports[0]
    na = next(iter(comb.adj[dnode]), None)
    nb = next(iter(comb.adj[unode]), None)
    comb.remove(dnode)
    comb.remove(unode)
    if na is not None and nb is not None:
        comb.link(na, nb)
    return reduce_and_code(comb)


# ---------------------------------------------------------------------------
# Sigma = dict code->Fraction weight + global registry code->canon_frag.
# ---------------------------------------------------------------------------

REG = {}  # code -> canon_frag


def reg(comp):
    code, frag = comp
    if code not in REG:
        REG[code] = frag
    return code


def sig_unary(sig, op, *args):
    out = defaultdict(Fraction)
    for code, w in sig.items():
        oc = reg(op((code, REG[code]), *args))
        out[oc] += w
    return out


def sig_binary(sigA, sigB, op):
    out = defaultdict(Fraction)
    for ca, wa in sigA.items():
        sa = sel_count(ca)
        for cb, wb in sigB.items():
            if sa + sel_count(cb) > 4:
                continue
            oc = reg(op((ca, REG[ca]), (cb, REG[cb])))
            out[oc] += wa * wb
    return out


def leaf_sigma(color, weight):
    absent, present = atom_leaf(color)
    reg(absent); reg(present)
    return {absent[0]: Fraction(1), present[0]: Fraction(weight)}


# ---------------------------------------------------------------------------
# Tree model + two DP evaluators + brute force.
# ---------------------------------------------------------------------------

class Tree:
    """Rooted tree.  node -> children list; leaves carry (color,weight)."""
    def __init__(self):
        self.children = defaultdict(list)
        self.leafinfo = {}   # node -> (color, weight)
        self.root = None

    def leaves(self):
        return list(self.leafinfo)


def eval_center_fold(tree):
    """Root Closed Sigma via star folds at every internal node."""
    def down(v, parent):
        if v in tree.leafinfo:
            c, w = tree.leafinfo[v]
            return leaf_sigma(c, w)
        kids = [c for c in tree.children[v]]
        G = None
        for c in kids:
            childC1 = down(c, v)
            childG = sig_unary(childC1, op_c_to_g)
            G = childG if G is None else sig_binary(G, childG, op_gg_to_g)
        if parent is None:
            return sig_unary(G, op_ig_to_c, False, False)
        return sig_unary(G, op_ig_to_c, True, False)
    return down(tree.root, None)


def eval_cc_spine(tree):
    """Alternative decomposition for a CATERPILLAR: fold the spine with CC_TO_C.
    Only valid when the tree is a caterpillar rooted so the spine is a path.
    Each spine node i has one pendant leaf (or subtree) and one child spine node.
    We build C2 corridor pieces and CC-join them, then close.  Returns Closed."""
    # Build the spine order from root down its single internal child chain.
    spine = []
    v = tree.root
    while True:
        spine.append(v)
        inter = [c for c in tree.children[v] if c not in tree.leafinfo]
        if len(inter) != 1:
            break
        v = inter[0]
    # each spine node -> its pendant leaves (as G at that center).
    # Represent spine node i as a C2 corridor piece: I center with its pendant
    # leaves folded as G, plus U (to prev) and D (to next).  End nodes differ.
    # Build per-node C2 (or C1 at ends) then CC fold; close remaining U.
    def center_c2(node, want_u, want_d):
        pend = [c for c in tree.children[node] if c in tree.leafinfo]
        G = None
        for c in pend:
            cc, ww = tree.leafinfo[c]
            g = sig_unary(leaf_sigma(cc, ww), op_c_to_g)
            G = g if G is None else sig_binary(G, g, op_gg_to_g)
        if G is None:  # no pendant: empty center G (identity center)
            e = Frag(); e.add(0, (Z,))
            G = {reg(reduce_and_code(e)): Fraction(1)}
        return sig_unary(G, op_ig_to_c, want_u, want_d)
    n = len(spine)
    # first node: only D (down); interior: U+D; last node: also its terminal
    # child chain leaves handled below.
    pieces = []
    for idx, node in enumerate(spine):
        want_u = idx != 0
        # last spine node: its children are all leaves -> fold as C1 (only U).
        is_last = (idx == n - 1)
        want_d = not is_last
        pieces.append(center_c2(node, want_u, want_d))
    # CC-fold: piece0 (U? no; D yes) is 'higher', join .D to next.U
    acc = pieces[0]
    for nxt in pieces[1:]:
        acc = sig_binary(acc, nxt, op_cc_to_c)
    # acc is C1 (leftover U at spine[0] end? spine[0] had want_u=False so it's D-... )
    # after all joins, remaining ports: spine[0].U(none) -> actually spine0 has no U,
    # last has no D.  So the corridor is Closed already if both ends portless, else C1.
    iface = None
    for code in acc:
        iface = interface_of(code)
        break
    if iface == 'C1':
        acc = sig_unary(acc, op_close_root)
    return acc


def readout(sig, K):
    """Sum over Closed states with 4 selected leaves -> channels."""
    ch = defaultdict(Fraction)
    for code, w in sig.items():
        if interface_of(code) != 'Closed':
            continue
        if sel_count(code) != 4:
            continue
        split = _four_leaf_split(code)
        if split is None:
            continue  # star / unresolved
        s1, s2 = split
        key = _channel(s1, s2, K)
        if key:
            ch[key] += w
    return ch


def _collect_leaf_colors(code):
    cols = []
    def walk(c):
        if c[0][0] == LEAF:
            cols.append(c[0][1])
        for ch in c[1:]:
            walk(ch)
    walk(code)
    return cols


def _four_leaf_split(code):
    """For a 4-leaf Closed reduced tree, return the unique 2|2 color multisets,
    or None if star (no internal 2|2 edge).  Rebuild frag and test edges."""
    frag = REG[code]
    # find internal edge whose removal gives 2 leaves | 2 leaves.
    leafnodes = [n for n, t in frag.kind.items() if t[0] == LEAF]
    if len(leafnodes) != 4:
        return None
    edges = set()
    for a in frag.adj:
        for b in frag.adj[a]:
            if a < b:
                edges.add((a, b))
    for (a, b) in edges:
        # split components by removing edge (a,b)
        comp = _side(frag, a, b)
        side_leaves = [n for n in comp if frag.kind[n][0] == LEAF]
        if len(side_leaves) == 2:
            other = [n for n in leafnodes if n not in side_leaves]
            s1 = sorted(frag.kind[n][1] for n in side_leaves)
            s2 = sorted(frag.kind[n][1] for n in other)
            return (tuple(s1), tuple(s2))
    return None


def _side(frag, a, b):
    """Nodes reachable from a without crossing edge (a,b)."""
    seen = {a}
    stack = [a]
    while stack:
        x = stack.pop()
        for y in frag.adj[x]:
            if (x, y) == (a, b) or (x, y) == (b, a):
                continue
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return seen


def _channel(s1, s2, K):
    m = tuple(sorted([s1, s2]))
    if K == 2:
        if m == ((0, 0), (1, 1)):
            return 'Q00_11'
        if m == ((0, 1), (0, 1)):
            return 'Q01_01'
        return None
    # K==4 four-distinct channels (plus the two-color ones for completeness)
    if m == ((0, 1), (2, 3)):
        return 'Q01_23'
    if m == ((0, 2), (1, 3)):
        return 'Q02_13'
    if m == ((0, 3), (1, 2)):
        return 'Q03_12'
    if m == ((0, 0), (1, 1)):
        return 'Q00_11'
    if m == ((0, 1), (0, 1)):
        return 'Q01_01'
    return None


# ---- brute force -----------------------------------------------------------

def brute_force(tree, K):
    leaves = tree.leaves()
    # build undirected adjacency for path tests
    adj = defaultdict(set)
    def walk(v):
        for c in tree.children[v]:
            adj[v].add(c); adj[c].add(v)
            walk(c)
    walk(tree.root)
    def path(a, b):
        prev = {a: None}
        dq = deque([a])
        while dq:
            x = dq.popleft()
            if x == b:
                break
            for y in adj[x]:
                if y not in prev:
                    prev[y] = x
                    dq.append(y)
        nodes = set()
        x = b
        while x is not None:
            nodes.add(x)
            x = prev[x]
        return nodes
    ch = defaultdict(Fraction)
    for quad in itertools.combinations(leaves, 4):
        # test the three pairings; resolved xy|zw iff path(x,y) & path(z,w) disjoint
        cols = {l: tree.leafinfo[l][0] for l in quad}
        w = Fraction(1)
        for l in quad:
            w *= Fraction(tree.leafinfo[l][1])
        for (x, y) in itertools.combinations(quad, 2):
            zw = [l for l in quad if l not in (x, y)]
            z, ww = zw
            if not (path(x, y) & path(z, ww)):
                s1 = tuple(sorted((cols[x], cols[y])))
                s2 = tuple(sorted((cols[z], cols[ww])))
                key = _channel(s1, s2, K)
                if key:
                    ch[key] += w
                break  # a 4-leaf tree has at most one 2|2 split
    return ch


# ---------------------------------------------------------------------------
# random tree generators
# ---------------------------------------------------------------------------

def random_tree(n_leaves, K, rng, allow_poly=True, frac=False):
    t = Tree()
    ctr = itertools.count()
    root = next(ctr)
    t.root = root
    t.children[root] = []
    internal = [root]
    for _ in range(n_leaves):
        par = rng.choice(internal)
        leaf = next(ctr)
        col = rng.randrange(K)
        wt = Fraction(rng.randint(1, 4)) if not frac else Fraction(rng.randint(1, 7), rng.randint(1, 5))
        t.leafinfo[leaf] = (col, wt)
        t.children[par].append(leaf)
        # sometimes insert an internal node to grow depth / create resolution
        if rng.random() < 0.5:
            nv = next(ctr)
            t.children[nv] = []
            t.children[par].append(nv)
            internal.append(nv)
    # ensure every internal node has >=1 child; drop empties by attaching a leaf
    for v in list(t.children):
        if v != root and not t.children[v] and v not in t.leafinfo:
            leaf = next(ctr)
            t.leafinfo[leaf] = (rng.randrange(K), Fraction(1))
            t.children[v].append(leaf)
    if not allow_poly:
        t = _binarize(t, rng)
    return t


def caterpillar(n_leaves, K, rng, frac=False):
    """Root spine chain; each spine node has one pendant leaf, last has two."""
    t = Tree()
    ctr = itertools.count()
    prev = None
    spine = []
    for i in range(n_leaves - 1):
        v = next(ctr)
        t.children[v] = []
        if prev is None:
            t.root = v
        else:
            t.children[prev].append(v)
        # pendant leaf
        leaf = next(ctr)
        col = rng.randrange(K)
        wt = Fraction(rng.randint(1, 4)) if not frac else Fraction(rng.randint(1, 7), rng.randint(1, 5))
        t.leafinfo[leaf] = (col, wt)
        t.children[v].append(leaf)
        spine.append(v)
        prev = v
    # final leaf on the last spine node
    leaf = next(ctr)
    col = rng.randrange(K)
    wt = Fraction(rng.randint(1, 4)) if not frac else Fraction(rng.randint(1, 7), rng.randint(1, 5))
    t.leafinfo[leaf] = (col, wt)
    t.children[prev].append(leaf)
    return t


def _binarize(t, rng):
    nt = Tree()
    ctr = itertools.count(10**6)
    def rec(v):
        if v in t.leafinfo:
            nid = next(ctr)
            nt.leafinfo[nid] = t.leafinfo[v]
            return nid
        kids = [rec(c) for c in t.children[v]]
        while len(kids) > 2:
            a = kids.pop(); b = kids.pop()
            nv = next(ctr)
            nt.children[nv] = [a, b]
            kids.append(nv)
        nv = next(ctr)
        nt.children[nv] = kids
        return nv
    nt.root = rec(t.root)
    return nt


# ---------------------------------------------------------------------------
# TEST 1: exactness (DP == brute force), two decompositions
# ---------------------------------------------------------------------------

def test_exactness():
    print("=== TEST 1: DP root readout == brute-force quartets ===")
    rng = random.Random(20260712)
    fails = 0
    trials = 0
    for K in (2, 4):
        for frac in (False, True):
            for kind in ('poly', 'binary', 'caterpillar'):
                for _ in range(40):
                    nlv = rng.randint(4, 8)
                    if kind == 'caterpillar':
                        t = caterpillar(nlv, K, rng, frac=frac)
                    else:
                        t = random_tree(nlv, K, rng,
                                        allow_poly=(kind == 'poly'), frac=frac)
                    bf = brute_force(t, K)
                    cf = readout(eval_center_fold(t), K)
                    trials += 1
                    if dict(bf) != dict(cf):
                        fails += 1
                        if fails <= 5:
                            print(f"  MISMATCH center-fold K={K} frac={frac} {kind}: bf={dict(bf)} dp={dict(cf)}")
                    # CC-spine decomposition only for caterpillars
                    if kind == 'caterpillar':
                        sp = readout(eval_cc_spine(t), K)
                        if dict(sp) != dict(bf):
                            fails += 1
                            if fails <= 8:
                                print(f"  MISMATCH cc-spine K={K} frac={frac}: bf={dict(bf)} dp={dict(sp)}")
    print(f"  trials={trials} fails={fails} -> {'PASS' if fails == 0 else 'FAIL'}")
    return fails == 0




# ---------------------------------------------------------------------------
# PART 2: center-fold state space (Closed/C1/G) + observable quotient rank.
# The scoring combine uses exactly these interfaces (C_TO_G, GG_TO_G, IG);
# the C2 corridor is an Extract/Contract rebuild optimization, finite by the
# same <=4-leaf/<=2-port bound but not the scoring hot-node payload.
# ---------------------------------------------------------------------------

def enumerate_cf(K):
    states = set()

    def add(comp):
        c = reg(comp); states.add(c); return c

    a = Frag(); a.add(0, (U,)); add(reduce_and_code(a))
    for c in range(K):
        p = Frag(); p.add(0, (U,)); p.add(1, (LEAF, c)); p.link(0, 1)
        add(reduce_and_code(p))

    U_C_TO_G, U_IG_FF, U_IG_TF, U_CLOSE = {}, {}, {}, {}
    B_GG = {}
    frontier = set(states)
    while frontier:
        new = set()
        for c in sorted(states, key=state_sort_key):
            iface = interface_of(c)
            if iface == 'C1':
                if c not in U_C_TO_G:
                    oc = reg(op_c_to_g((c, REG[c]))); U_C_TO_G[c] = oc
                    if oc not in states: new.add(oc)
                if c not in U_CLOSE:
                    oc = reg(op_close_root((c, REG[c]))); U_CLOSE[c] = oc
                    if oc not in states: new.add(oc)
            elif iface == 'G':
                if c not in U_IG_FF:
                    oc = reg(op_ig_to_c((c, REG[c]), False, False)); U_IG_FF[c] = oc
                    if oc not in states: new.add(oc)
                if c not in U_IG_TF:
                    oc = reg(op_ig_to_c((c, REG[c]), True, False)); U_IG_TF[c] = oc
                    if oc not in states: new.add(oc)
        # GG over G pairs, gated by sel<=4 (index by sel to avoid O(|G|^2))
        by_sel = defaultdict(list)
        for c in sorted(states, key=state_sort_key):
            if interface_of(c) == 'G':
                by_sel[sel_count(c)].append(c)
        for sx in range(5):
            for sy in range(5 - sx):
                for x in sorted(by_sel[sx], key=state_sort_key):
                    for y in sorted(by_sel[sy], key=state_sort_key):
                        if (x, y) in B_GG:
                            continue
                        oc = reg(op_gg_to_g((x, REG[x]), (y, REG[y]))); B_GG[(x, y)] = oc
                        if oc not in states: new.add(oc)
        states |= new; frontier = new

    by_iface = defaultdict(set)
    for c in states:
        by_iface[interface_of(c)].add(c)
    unary = {'C_TO_G': (U_C_TO_G, 'C1', 'G'), 'IG_FF': (U_IG_FF, 'G', 'Closed'),
             'IG_TF': (U_IG_TF, 'G', 'C1'), 'CLOSE': (U_CLOSE, 'C1', 'Closed')}
    return by_iface, unary, B_GG


class Basis:
    """Incremental exact-rational row-echelon basis; rank = len(rows)."""
    def __init__(self):
        self.rows = []
        self.pivots = {}

    def add(self, row):
        row = {c: Fraction(v) for c, v in row.items() if v != 0}
        while row:
            piv = [c for c in row if c in self.pivots]
            if piv:
                pc = min(piv)
                pr = self.rows[self.pivots[pc]]
                f = row[pc]
                for c, v in pr.items():
                    nv = row.get(c, Fraction(0)) - f * v
                    if nv == 0:
                        row.pop(c, None)
                    else:
                        row[c] = nv
            else:
                lead = min(row); pv = row[lead]
                self.pivots[lead] = len(self.rows)
                self.rows.append({c: v / pv for c, v in row.items()})
                return True
        return False

    def rank(self):
        return len(self.rows)

    def in_span(self, row):
        """True iff row is in the span of the current basis (non-mutating)."""
        row = {c: Fraction(v) for c, v in row.items() if v != 0}
        while row:
            piv = [c for c in row if c in self.pivots]
            if not piv:
                return False
            pc = min(piv)
            pr = self.rows[self.pivots[pc]]
            f = row[pc]
            for c, v in pr.items():
                nv = row.get(c, Fraction(0)) - f * v
                if nv == 0:
                    row.pop(c, None)
                else:
                    row[c] = nv
        return True


def _row_sort_key(row):
    return tuple((c, v.numerator, v.denominator)
                 for c, v in sorted(row.items()))


def canonical_rref(rows):
    """Return the unique exact RREF of a row space in canonical column order."""
    work = [{c: Fraction(v) for c, v in row.items() if v != 0}
            for row in rows]
    work = [row for row in work if row]
    work.sort(key=_row_sort_key)
    columns = sorted({c for row in work for c in row})
    rank = 0
    for col in columns:
        candidates = [i for i in range(rank, len(work))
                      if work[i].get(col, Fraction(0)) != 0]
        if not candidates:
            continue
        pivot = min(candidates, key=lambda i: _row_sort_key(work[i]))
        work[rank], work[pivot] = work[pivot], work[rank]
        pv = work[rank][col]
        work[rank] = {c: v / pv for c, v in work[rank].items()}
        prow = work[rank]
        for i in range(len(work)):
            if i == rank:
                continue
            factor = work[i].get(col, Fraction(0))
            if factor == 0:
                continue
            row = dict(work[i])
            for c, v in prow.items():
                nv = row.get(c, Fraction(0)) - factor * v
                if nv == 0:
                    row.pop(c, None)
                else:
                    row[c] = nv
            work[i] = row
        rank += 1
        if rank == len(work):
            break
    result = [row for row in work if row]
    result.sort(key=lambda row: min(row))
    return result


def canonicalize_basis(basis):
    result = Basis()
    result.rows = canonical_rref(basis.rows)
    result.pivots = {min(row): i for i, row in enumerate(result.rows)}
    return result


def readout_rows(K, closed, channels=None):
    if channels is None:
        # The channels QRAFT reads: the divide (K=2) reads the two 2-colour
        # channels, regrafting and NNI (K=4) the three quartets with four
        # distinct colours. Pass channels=... to override.
        channels = ['Q00_11', 'Q01_01'] if K == 2 else \
            ['Q01_23', 'Q02_13', 'Q03_12']
    rows = {ch: defaultdict(Fraction) for ch in channels}
    for code in sorted(closed, key=state_sort_key):
        if sel_count(code) != 4:
            continue
        sp = _four_leaf_split(code)
        if sp is None:
            continue
        ch = _channel(sp[0], sp[1], K)
        if ch in rows:
            rows[ch][code] += 1
    return [dict(rows[ch]) for ch in channels if rows[ch]]


def quotient(K, ifaces, closed_states, unary_tbls, binary_tbls):
    """§3.16 observable-quotient rank via an EXACT worklist closure (each
    functional pulled back once through every op's inverse; runs to true
    fixpoint, no iteration cap).  Returns dim(W_T) per interface.

    CRITICAL: unary and binary tables can be MIXED-interface (REVERSE has domain
    C1 u C2; CC_TO_C's higher operand is C2 or D-only-C1).  A pulled-back row
    must be a functional on a SINGLE interface, so every table is split into
    single-interface buckets by operand interface before inversion; otherwise a
    row would mix coordinates from two different spaces and corrupt the rank."""
    Bq = {t: Basis() for t in ifaces}
    # unary buckets: (in_iface) -> {out_code: [x, ...]}
    un_buckets = []
    for tbl in unary_tbls:
        by_in = defaultdict(lambda: defaultdict(list))
        for x, ox in sorted(tbl.items(), key=lambda item: state_sort_key(item[0])):
            by_in[interface_of(x)][ox].append(x)
        for _in_if, inv in by_in.items():
            un_buckets.append(inv)
    # binary buckets: (hif, lif) -> {out_code: [(a, b), ...]}
    bin_buckets = []
    for tbl in binary_tbls:
        by_pair = defaultdict(lambda: defaultdict(list))
        for (a, b), oc in sorted(
                tbl.items(),
                key=lambda item: (state_sort_key(item[0][0]),
                                  state_sort_key(item[0][1]))):
            by_pair[(interface_of(a), interface_of(b))][oc].append((a, b))
        for _pair, inv in by_pair.items():
            bin_buckets.append(inv)
    work = deque()

    def push(row):
        row = {c: v for c, v in row.items() if v != 0}
        if not row:
            return
        ks = iter(row)
        t = interface_of(next(ks))
        # invariant: a pullback row lives in one interface (bucketing guarantees it)
        assert all(interface_of(k) == t for k in ks), "mixed-interface pullback row"
        if Bq[t].add(dict(row)):
            work.append(row)

    for r in readout_rows(K, closed_states):
        push(r)
    while work:
        f = work.popleft()
        for inv in un_buckets:                  # unary: row[x] += f[op(x)]
            row = defaultdict(Fraction)
            for oc, val in f.items():
                for x in inv.get(oc, ()):
                    row[x] += val
            push(row)
        for inv in bin_buckets:                 # binary: fix one operand (single-iface bucket)
            left = defaultdict(lambda: defaultdict(Fraction))
            right = defaultdict(lambda: defaultdict(Fraction))
            for oc, val in f.items():
                for (a, b) in inv.get(oc, ()):
                    left[b][a] += val
                    right[a][b] += val
            for _, row in left.items():
                push(row)
            for _, row in right.items():
                push(row)
    Bq = {t: canonicalize_basis(Bq[t]) for t in ifaces}
    return {t: Bq[t].rank() for t in ifaces}, Bq


def _has(code, port):
    def walk(c):
        if c[0][0] == port:
            return True
        return any(walk(ch) for ch in c[1:])
    return walk(code)


def enumerate_full(K):
    """Full grammar incl. the two-port C2 corridor (IG_TT, CC_TO_C, REVERSE)."""
    states = set()
    def add(comp):
        c = reg(comp); states.add(c); return c
    a = Frag(); a.add(0, (U,)); add(reduce_and_code(a))
    for c in range(K):
        p = Frag(); p.add(0, (U,)); p.add(1, (LEAF, c)); p.link(0, 1)
        add(reduce_and_code(p))

    UC_TO_G, UIG_FF, UIG_TF, UIG_TT, UCLOSE, UREV = {}, {}, {}, {}, {}, {}
    B_GG, B_CC = {}, {}
    frontier = set(states)
    while frontier:
        new = set()
        for c in sorted(states, key=state_sort_key):
            iface = interface_of(c)
            if iface == 'C1':
                # C_TO_G and CLOSE_ROOT act on the UP port; on a D-only C1 they are
                # no-op identities (spurious self-loops) -- guard to U-bearing C1s.
                ops = [(UREV, lambda c: op_reverse_c((c, REG[c])))]
                if _has(c, U):
                    ops = [(UC_TO_G, lambda c: op_c_to_g((c, REG[c]))),
                           (UCLOSE, lambda c: op_close_root((c, REG[c])))] + ops
                for tbl, fn in ops:
                    if c not in tbl:
                        oc = reg(fn(c)); tbl[c] = oc
                        if oc not in states: new.add(oc)
            elif iface == 'C2':
                if c not in UREV:
                    oc = reg(op_reverse_c((c, REG[c]))); UREV[c] = oc
                    if oc not in states: new.add(oc)
            elif iface == 'G':
                for tbl, fn in ((UIG_FF, lambda c: op_ig_to_c((c, REG[c]), False, False)),
                                (UIG_TF, lambda c: op_ig_to_c((c, REG[c]), True, False)),
                                (UIG_TT, lambda c: op_ig_to_c((c, REG[c]), True, True))):
                    if c not in tbl:
                        oc = reg(fn(c)); tbl[c] = oc
                        if oc not in states: new.add(oc)
        # GG (G x G), sel-indexed
        gby = defaultdict(list)
        for c in sorted(states, key=state_sort_key):
            if interface_of(c) == 'G':
                gby[sel_count(c)].append(c)
        for sx in range(5):
            for sy in range(5 - sx):
                for x in sorted(gby[sx], key=state_sort_key):
                    for y in sorted(gby[sy], key=state_sort_key):
                        if (x, y) in B_GG: continue
                        oc = reg(op_gg_to_g((x, REG[x]), (y, REG[y]))); B_GG[(x, y)] = oc
                        if oc not in states: new.add(oc)
        # CC (higher has D) x (lower has U), sel-indexed
        hby, lby = defaultdict(list), defaultdict(list)
        for c in sorted(states, key=state_sort_key):
            if _has(c, D): hby[sel_count(c)].append(c)
            if _has(c, U): lby[sel_count(c)].append(c)
        for sx in range(5):
            for sy in range(5 - sx):
                for h in sorted(hby[sx], key=state_sort_key):
                    for l in sorted(lby[sy], key=state_sort_key):
                        if (h, l) in B_CC: continue
                        oc = reg(op_cc_to_c((h, REG[h]), (l, REG[l]))); B_CC[(h, l)] = oc
                        if oc not in states: new.add(oc)
        states |= new; frontier = new

    by_iface = defaultdict(set)
    for c in states:
        by_iface[interface_of(c)].add(c)
    unary = {'C_TO_G': (UC_TO_G,), 'IG_FF': (UIG_FF,), 'IG_TF': (UIG_TF,),
             'IG_TT': (UIG_TT,), 'CLOSE': (UCLOSE,), 'REV': (UREV,)}
    return by_iface, unary, B_GG, B_CC


# ---------------------------------------------------------------------------
# Independent completeness oracle: bounded Pruefer enumeration.
#
# Selected leaves and ordinary ports are terminals, hence never occur in a
# Pruefer sequence. Only unlabeled internal vertices (and protected Z for G)
# occur in the sequence. The worst ordinary case is four leaves + U/D + four
# internal vertices: the degree>=3 filter leaves only 2520 length-eight
# sequences, so this is exhaustive without a large search.
# ---------------------------------------------------------------------------

def _tree_from_prufer(sequence, node_count):
    degree = [1] * node_count
    for value in sequence:
        degree[value] += 1
    leaves = [i for i, value in enumerate(degree) if value == 1]
    heapq.heapify(leaves)
    edges = []
    for value in sequence:
        leaf = heapq.heappop(leaves)
        edges.append((leaf, value))
        degree[leaf] -= 1
        degree[value] -= 1
        if degree[value] == 1:
            heapq.heappush(leaves, value)
    if node_count >= 2:
        a = heapq.heappop(leaves)
        b = heapq.heappop(leaves)
        edges.append((a, b))
    return edges


def _prufer_topologies(interface, leaf_count, c1_port=U):
    if interface == 'Closed':
        port_tokens = []
    elif interface == 'C1':
        port_tokens = [(c1_port,)]
    elif interface == 'C2':
        port_tokens = [(U,), (D,)]
    elif interface == 'G':
        port_tokens = [(Z,)]
    else:
        raise ValueError(f"unknown interface {interface}")

    terminal_count = leaf_count + len(port_tokens)
    if terminal_count == 0:
        return {reduce_and_code(Frag())[0]: Frag()}

    max_internal = (leaf_count - 1) if interface == 'G' else (terminal_count - 2)
    max_internal = max(0, max_internal)
    topologies = {}
    for internal_count in range(max_internal + 1):
        node_count = terminal_count + internal_count
        if node_count == 1:
            f = Frag()
            token = (LEAF, 0) if leaf_count else port_tokens[0]
            f.add(0, token)
            code, canon = reduce_and_code(f)
            topologies[code] = canon
            continue

        first_port = leaf_count
        first_internal = terminal_count
        internal_ids = list(range(first_internal, node_count))
        if interface == 'G':
            z_id = first_port
            alphabet = [z_id] + internal_ids
        else:
            alphabet = internal_ids
        sequence_length = node_count - 2
        for sequence in itertools.product(alphabet, repeat=sequence_length):
            if any(sequence.count(i) < 2 for i in internal_ids):
                continue
            f = Frag()
            for i in range(leaf_count):
                f.add(i, (LEAF, i))
            for offset, token in enumerate(port_tokens):
                f.add(first_port + offset, token)
            for i in internal_ids:
                f.add(i, (I,))
            for a, b in _tree_from_prufer(sequence, node_count):
                f.link(a, b)
            code, canon = reduce_and_code(f)
            topologies[code] = canon
    return topologies


def enumerate_prufer_shapes(K):
    """Enumerate all legal reduced shapes independently of the HDT grammar."""
    result = {t: set() for t in ('Closed', 'C1', 'C2', 'G')}
    for interface in result:
        for leaf_count in range(5):
            orientations = (U, D) if interface == 'C1' else (U,)
            for orientation in orientations:
                for frag in _prufer_topologies(
                        interface, leaf_count, orientation).values():
                    for colors in itertools.product(range(K), repeat=leaf_count):
                        colored = frag.copy()
                        for nid, token in list(colored.kind.items()):
                            if token[0] == LEAF:
                                colored.kind[nid] = (LEAF, colors[token[1]])
                        code, _ = reduce_and_code(colored)
                        result[interface].add(code)
    return result


def verify_state_completeness(K, grammar_states):
    independent = enumerate_prufer_shapes(K)
    for interface in ('Closed', 'C1', 'C2', 'G'):
        missing = independent[interface] - grammar_states[interface]
        extra = grammar_states[interface] - independent[interface]
        if missing or extra:
            raise AssertionError(
                f"state completeness failed K={K} {interface}: "
                f"missing={len(missing)} extra={len(extra)}")
    return {t: len(independent[t]) for t in ('Closed', 'C1', 'C2', 'G')}


if __name__ == '__main__':
    ok = test_exactness()
    if ok:
        for K in (2, 4):
            print(f"\n=== PART 2: K={K} center-fold (scoring interfaces C1/G/Closed) ===")
            by_iface, unary, B_GG = enumerate_cf(K)
            counts = {t: len(by_iface[t]) for t in ('Closed', 'C1', 'G')}
            print(f"  reachable canonical states: {counts}  (GG table {len(B_GG)})")
            q, _ = quotient(K, ('Closed', 'C1', 'G'), by_iface['Closed'],
                            [tbl for (tbl, _i, _o) in unary.values()], [B_GG])
            print(f"  observable quotient rank (min O(1) scoring payload): {q}")
        for K in (2, 4):
            print(f"\n=== PART 3: K={K} FULL grammar (incl. two-port C2 corridor) ===")
            REG.clear()
            by, un, gg, cc = enumerate_full(K)
            counts = verify_state_completeness(K, by)
            print(f"  complete canonical states: {counts}  (GG {len(gg)}, CC {len(cc)})")
            qf, _ = quotient(K, ('Closed', 'C1', 'C2', 'G'), by['Closed'],
                             [tbl for (tbl,) in un.values()], [gg, cc])
            print(f"  observable quotient rank (min O(1) node payload): {qf}")
    sys.exit(0 if ok else 1)
