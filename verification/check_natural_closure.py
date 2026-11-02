"""Upper bound of Theorem 4.2 (minimal number of counters): the natural counters of Sec. 5, all sorts,
k = 2 and 4, checked directly on the census automaton:

  (a) the counters of each sort are linearly independent vectors over the shapes of that sort;
  (b) each counter has a witness shape on which it alone is non-zero (value 1);
  (c) the gamma^a lie in the span of the counters of sort R;
  (d) closure: for every symbol sigma: X1 X2 -> X, every counter b of sort X and every shape pi2
      of sort X2, the vector pi1 -> b(sigma(pi1, pi2)) lies in the span of the counters of sort
      X1; likewise in the second argument; for the transformations L -> C1 -> G the restriction
      of b lies in the span.

Shapes and the gluing sigma(pi1, pi2) are those of kmw_closure_full.py (sorts L, I, Irho, C1, C2,
C3, G, R; Irho is the I component at the root, which has no top vertex).  Membership is tested
mod 2^31-1 and exactly over Q (a vector is in the span iff it equals the combination given by its
values at the witnesses)."""
import sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE / 'natwta'), str(HERE / 'natwta' / 'code')]
import clusters as C
from natural_wta import NaturalWTA, unrooted
from kmw_closure import RowSpace, P
import kmw_closure_full as KF

def in_span(rs, row):
    row = {c: v % P for c, v in row.items() if v % P}
    for c in sorted(row):
        if row.get(c) and c in rs.piv:
            f = row[c]
            for cc, vv in rs.piv[c].items():
                row[cc] = (row.get(cc, 0) - f * vv) % P
            row = {cc: vv for cc, vv in row.items() if vv}
    return not row

def exact_member(nat_s, wit_s, row):
    """exact test over Q (integers): row lies in the span of the counters iff it equals the
    combination whose coefficients are its values at the witnesses (each witness is 1 on its
    counter and 0 on all others)"""
    comb = {}
    for (lab, vec), w in zip(nat_s, wit_s):
        c = row.get(w, 0)
        if c:
            for i, v in vec.items():
                comb[i] = comb.get(i, 0) + c * v
    comb = {i: v for i, v in comb.items() if v}
    row = {i: v for i, v in row.items() if v}
    return comb == row

def natural_counters(K, sh):
    """sort -> list of (label, {shape index: value})"""
    M = NaturalWTA(K); ch = M.channels
    pairs = [(0, 0), (0, 1), (1, 1)] if K == 2 else [(i, j) for i in range(K) for j in range(i + 1, K)]
    out = {}
    def from_features(sort, fsort, conv=lambda c: c):
        labels = M.labels[fsort]
        vecs = [dict() for _ in labels]
        for i, s in enumerate(sh[sort]):
            for lab, v in M.features(fsort, conv(s)).items():
                vecs[lab][i] = v
        return list(zip(labels, vecs))
    out['C1'] = from_features('C1', 'C1'); out['C2'] = from_features('C2', 'C2'); out['G'] = from_features('G', 'G')
    out['R'] = from_features('R', 'T')
    # L: unit and leaf counters (the C1 counters restricted to the L shapes)
    out['L'] = [(lab, {i: v for i, v in vec.items()}) for lab, vec in from_features('L', 'C1') if vec]
    out['I'] = [(('unit',), {0: 1})]; out['Irho'] = [(('unit',), {0: 1})]
    # C3: 1, n_i, n_[ij], n_down[xbar(yy')], q_a   (shapes of C3 = G shapes with b as super root)
    labels = [('unit',)] + [('leaf', c) for c in range(K)] + [('pair', i, j) for i, j in pairs] \
             + [('triple', 'B', c, a) for c in range(K) for a in ch] + [('quartet', a) for a in ch]
    idx = {l: n for n, l in enumerate(labels)}
    vecs = [dict() for _ in labels]
    for n, s in enumerate(sh['C3']):
        g = C.graph_from_code(s)
        gT = C.canon(C.relabel_kind(g, C.find_kind(g, 'B'), 'T'))        # back to a G code
        d = C.degree(s); f = []
        if d == 0: f = [('unit',)]
        elif d == 1: f = [('leaf', next(iter(C.colours_of(s))))]
        elif d == 2:
            kinds, adj = C.graph_from_code(s)
            cols = tuple(sorted(k[1] for k in kinds if C.is_leaf_kind(k)))
            if ('pair',) + cols in idx: f = [('pair',) + cols]
        elif d == 3:
            for c in range(K):
                closed = unrooted(gT, 'T', c)
                f += [('triple', 'B', c, a) for a, vals in ch.items() if vals.get(closed, 0)]
        else:
            closed = unrooted(gT)
            f = [('quartet', a) for a, vals in ch.items() if vals.get(closed, 0)]
        for l in f: vecs[idx[l]][n] = 1
    out['C3'] = list(zip(labels, vecs))
    return out, ch

def run(K):
    t0 = time.time()
    sh = {s: KF.shapes(s, K) for s in KF.SORTS}
    idx = {s: {c: i for i, c in enumerate(sh[s])} for s in KF.SORTS}
    deg = {s: [C.degree(c) for c in sh[s]] for s in KF.SORTS}
    bydeg = {s: {} for s in KF.SORTS}
    for s in KF.SORTS:
        for i, d in enumerate(deg[s]): bydeg[s].setdefault(d, []).append(i)
    nat, ch = natural_counters(K, sh)
    span = {}; wit = {}
    for s in KF.SORTS:
        rs = RowSpace(len(sh[s]))
        for lab, vec in nat[s]:
            assert vec, ('empty counter', s, lab)
            assert rs.insert(vec), ('dependent counter', s, lab)                       # (a)
        span[s] = rs
        wit[s] = []
        for lab, vec in nat[s]:                                                          # (b)
            ws = [i for i, v in vec.items() if v == 1 and all(i not in v2 for l2, v2 in nat[s] if l2 != lab)]
            assert ws, ('no witness', s, lab)
            wit[s].append(ws[0])
    for a, vals in ch.items():                                                           # (c)
        g = {idx['R'][p]: 1 for p in vals if p in idx['R']}
        assert in_span(span['R'], g) and exact_member(nat['R'], wit['R'], g), ('gamma not in span', a)
    # (d) closure
    cache = {}
    def glue(sym, a, b):
        key = (sym, a, b)
        if key not in cache: cache[key] = KF.ghat(sym, (a, b))
        return cache[key]
    nchk = 0
    for sym, (ins, out) in KF.SIG.items():
        if len(ins) == 1:
            src = ins[0]
            for lab, vec in nat[out]:
                row = {i: vec.get(idx[out].get(sh[src][i], -1), 0) for i in range(len(sh[src]))}
                row = {i: v for i, v in row.items() if v}
                assert in_span(span[src], row) and exact_member(nat[src], wit[src], row), ('unary closure fails', sym, lab); nchk += 1
            continue
        s1, s2 = ins
        for lab, vec in nat[out]:
            d = deg[out][next(iter(vec))]                   # every counter lives on one leaf count
            assert all(deg[out][i] == d for i in vec)
            for (sa, sb, first) in ((s1, s2, True), (s2, s1, False)):
                for db in range(d + 1):
                    for jb in bydeg[sb].get(db, []):
                        row = {}
                        for ja in bydeg[sa].get(d - db, []):
                            r = glue(sym, sh[sa][ja], sh[sb][jb]) if first else glue(sym, sh[sb][jb], sh[sa][ja])
                            if r is None: continue
                            v = vec.get(idx[out].get(r, -1), 0)
                            if v: row[ja] = v
                        if row:
                            assert in_span(span[sa], row) and exact_member(nat[sa], wit[sa], row), ('closure fails', sym, lab, first); nchk += 1
    dims = {s: span[s].dim() for s in KF.SORTS}
    print(f'k={K}: counters per sort {dims}; independent, witnesses, gamma in span, '
          f'closed under all {len(KF.SIG)} symbols ({nchk} non-zero pull-backs checked, exactly over Q via witnesses): ok  ({time.time()-t0:.1f}s)')
    return dims

if __name__ == '__main__':
    for K in (2, 4):
        run(K)
