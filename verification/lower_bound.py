"""Lower bound on one HDT (Lemma 4.1, 'one HDT', of the note).

For a node x of an HDT of T, of sort X at some stage, with inside leaves L(x) and outside leaves
L(T) - L(x), the matrix N_x has
   rows     (S, c)       S subset of L(x),       |S| <= 4,  c a colouring of S,
   columns  (S', c', a)  S' subset of outside,   |S|+|S'| = 4, c' a colouring of S', a a channel,
   entry    1 iff the quartet T|(S u S') is resolved and of type a under c u c'.
Every decoration scheme that is correct on this HDT for all colourings and all 0/1 weights has at
least rank N_x counters of sort X.  We compute the rank of a submatrix (one leaf set per 'shape'
T|(S u boundary vertices), all colourings), which is a lower bound for rank N_x; exact arithmetic
modulo p = 2^31-1 (a lower bound for the rank over Q)."""
import itertools, sys
from kmw_closure import RowSpace

def channels(k):
    if k == 2: return [((0, 0), (1, 1)), ((0, 1), (0, 1))]
    if k == 4: return [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
    raise ValueError(k)

def canon_shape(T, S, marks):
    """canonical code of the smallest subtree of T containing S and the marks, degree-two
    vertices suppressed (marks and leaves kept); marks: dict name -> vertex"""
    lab = {v: 'l' for v in S}
    for name, v in marks.items(): lab[v] = name
    M = list(lab)
    if not M: return ''
    V = set(M)
    for a, b in itertools.combinations(M, 2): V.add(T.lca(a, b))
    V = list(V)
    adj = {v: set() for v in V}
    for u in V:
        best = None
        for w in V:
            if w != u and T.lca(w, u) == w and (best is None or T.depth[w] > T.depth[best]): best = w
        if best is not None: adj[u].add(best); adj[best].add(u)
    changed = True
    while changed:                     # suppress unlabelled vertices of degree two
        changed = False
        for v in list(adj):
            if v not in lab and len(adj[v]) == 2:
                a, b = adj[v]; adj[a].discard(v); adj[b].discard(v); adj[a].add(b); adj[b].add(a)
                del adj[v]; changed = True
    def code(v, p):
        return lab.get(v, 'i') + '(' + ','.join(sorted(code(w, v) for w in adj[v] if w != p)) + ')'
    roots = [marks[n] for n in ('T', 'B', 'v') if n in marks]
    if roots: return code(roots[0], None)
    return min(code(v, None) for v in adj)

def reps(T, leaves, marks, maxd):
    """one leaf set per shape, by number of leaves"""
    out = {d: {} for d in range(maxd + 1)}
    for d in range(min(maxd, len(leaves)) + 1):
        for S in itertools.combinations(leaves, d):
            key = canon_shape(T, S, marks)
            out[d].setdefault(key, S)
    return {d: list(v.values()) for d, v in out.items()}

def split(T, q):
    """the resolved split of four leaves as a pair of index pairs, or None (star)"""
    a, b, c, d = q
    s1 = T.dist(a, b) + T.dist(c, d); s2 = T.dist(a, c) + T.dist(b, d); s3 = T.dist(a, d) + T.dist(b, c)
    if s1 < s2 and s1 < s3: return ((0, 1), (2, 3))
    if s2 < s1 and s2 < s3: return ((0, 2), (1, 3))
    if s3 < s1 and s3 < s2: return ((0, 3), (1, 2))
    return None

def typ(cols, sp):
    A = tuple(sorted((cols[sp[0][0]], cols[sp[0][1]]))); B = tuple(sorted((cols[sp[1][0]], cols[sp[1][1]])))
    return tuple(sorted((A, B)))

def rank_at(T, node, stage, k, verbose=False):
    """rank of the submatrix of N_x for the HDT node at the given stage (index into node.sorts)"""
    sort, marks, _ = node.sorts[stage]
    inside = [l for l in T.leaves if l in node.nodes]
    outside = [l for l in T.leaves if l not in node.nodes]
    R_in = reps(T, inside, marks, 4)
    R_out = reps(T, outside, marks, 4)
    chans = [tuple(sorted(a)) for a in channels(k)]
    prof = []
    for d in range(5):
        rows_S = R_in.get(d, []); cols_S = R_out.get(4 - d, [])
        if not rows_S or not cols_S:
            prof.append(0); continue
        col_index = {}
        for j, Sp in enumerate(cols_S):
            for cp in itertools.product(range(k), repeat=4 - d):
                for ai in range(len(chans)):
                    col_index[(j, cp, ai)] = len(col_index)
        sp = [[split(T, tuple(S) + tuple(Sp)) for Sp in cols_S] for S in rows_S]
        rs = RowSpace(len(col_index))
        for i, S in enumerate(rows_S):
            for c in itertools.product(range(k), repeat=d):
                row = {}
                for j, Sp in enumerate(cols_S):
                    s = sp[i][j]
                    if s is None: continue
                    for cp in itertools.product(range(k), repeat=4 - d):
                        t = typ(c + cp, s)
                        for ai, a in enumerate(chans):
                            if t == a: row[col_index[(j, cp, ai)]] = 1
                if row: rs.insert(row)
        prof.append(rs.dim())
    if verbose:
        print(f'  sort {sort}: |L(x)|={len(inside)}, outside {len(outside)}, '
              f'inside shapes {[len(R_in.get(d, [])) for d in range(5)]}, '
              f'outside shapes {[len(R_out.get(d, [])) for d in range(5)]}, rank profile {prof} = {sum(prof)}')
    return prof

TARGET = {2: {'L': 3, 'I': 1, 'Irho': 1, 'C1': 12, 'C2': 23, 'C3': 12, 'G': 15, 'R': 2},
          4: {'L': 5, 'I': 1, 'Irho': 1, 'C1': 26, 'C2': 56, 'C3': 26, 'G': 32, 'R': 3}}
