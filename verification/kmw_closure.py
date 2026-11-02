"""Run the minimisation procedure of Kiefer-Marusic-Worrell (LMCS 2017, Sec. 3) on the
shape automaton A_0 for the QRAFT tasks, literally:

  forward space  F_X  = closure of the leaf vectors under (v1 (x) v2) mu(sigma)   [Prop. 3.1(a)]
  backward space B_X  = closure of the final vectors gamma^a under mu(c) for one-step
                        contexts c = sigma(., t) / sigma(t, .) with t in a spanning set S of F [Prop. 3.2(a)]
  m_X = rank(F_X B_X)                                                           [Lemma 3.3]

Exact arithmetic modulo the prime p = 2^31 - 1 (0/1 matrices; ranks also verified over Q
elsewhere).  Degree bookkeeping (a shape with > 4 leaves is dropped) keeps the closure finite.
"""
import sys, time, itertools
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / 'natwta' / 'code'))
import clusters as C

P = (1 << 31) - 1
SORTS = ('T', 'C1', 'G', 'C2')
READOUT = sys.argv[1] if len(sys.argv) > 1 else 'C1'   # 'C1' (Brodal: final component is a closed C) or 'T'

def shapes(sort, K):
    out = []
    for m in range(5):
        out += C.fragments(sort, m, K)
    return out

class RowSpace:
    """incrementally maintained reduced basis (rows over GF(p)) of a subspace of F^n"""
    def __init__(self, n):
        self.n = n; self.piv = {}          # pivot column -> reduced row (dict col->val)
    def insert(self, row):                 # row: dict col -> val ; returns True if new
        row = {c: v % P for c, v in row.items() if v % P}
        for c in sorted(row):
            if c in self.piv and row.get(c):
                f = row[c]; prow = self.piv[c]
                for cc, vv in prow.items():
                    row[cc] = (row.get(cc, 0) - f * vv) % P
                row = {cc: vv for cc, vv in row.items() if vv}
        if not row:
            return False
        c0 = min(row); inv = pow(row[c0], P - 2, P)
        row = {cc: (vv * inv) % P for cc, vv in row.items()}
        for c, prow in self.piv.items():
            if prow.get(c0):
                f = prow[c0]
                new = dict(prow)
                for cc, vv in row.items():
                    new[cc] = (new.get(cc, 0) - f * vv) % P
                self.piv[c] = {cc: vv for cc, vv in new.items() if vv}
        self.piv[c0] = row
        return True
    def dim(self):
        return len(self.piv)

def forward(K):
    """Prop. 3.1(a): closure of the leaf vectors mu(L_{i,w}) = e_empty + w e_i (w in {0,1})
    under the transitions.  Since e_empty, e_i in F, the closure on basis vectors is the set of
    reachable shapes; F_X is full iff every enumerated shape is reachable."""
    reach = {s: set() for s in SORTS}
    for c in range(K):
        for w in (0, 1):
            reach['C1'] |= set(C.atom_state('leaf', c, w).keys())
    done = set(); changed = True
    while changed:
        changed = False
        for op, (ar, ins, out, _) in C.OPS.items():
            bydeg = [{} for _ in ins]
            for k, s in enumerate(ins):
                for c in reach[s]:
                    bydeg[k].setdefault(C.degree(c), []).append(c)
            if ar == 1:
                combos = [(c,) for d, cs in bydeg[0].items() for c in cs]
            else:
                combos = [(c1, c2) for d1, cs1 in bydeg[0].items() for d2, cs2 in bydeg[1].items()
                          if d1 + d2 <= 4 for c1 in cs1 for c2 in cs2]
            for codes in combos:
                if (op, codes) in done: continue
                done.add((op, codes))
                r = C.ghat(op, codes)
                if r is not None and r not in reach[out]:
                    reach[out].add(r); changed = True
    dims = {}
    for s in SORTS:
        sh = set(shapes(s, K))
        assert reach[s] == sh, (s, len(reach[s]), len(sh))
        dims[s] = len(sh)
    return dims

def backward(K, verbose=True, return_basis=False):
    """Prop. 3.2(a): B is the smallest space containing the gamma^a and closed under
    v -> mu(c) v for one-step contexts c whose other subtree ranges over a spanning set of F
    (here: the unit vectors e_t of all shapes)."""
    ch = C.qraft_channels(K)
    sh = {s: shapes(s, K) for s in SORTS}
    idx = {s: {c: i for i, c in enumerate(sh[s])} for s in SORTS}
    bydeg = {s: {} for s in SORTS}
    for s in SORTS:
        for c in sh[s]:
            bydeg[s].setdefault(C.degree(c), []).append(c)
    B = {s: RowSpace(len(sh[s])) for s in SORTS}
    if READOUT == 'T':
        for name, d in ch.items():
            B['T'].insert({idx['T'][pat]: 1 for pat in d if pat in idx['T']})
    else:   # read-out at the downward closed C (Brodal's final component): q_a on its 4-leaf shapes
        for name, d in ch.items():
            row = {}
            for i, t in enumerate(sh['C1']):
                if C.degree(t) == 4 and d.get(C.ghat('close', (t,)), 0):
                    row[i] = 1
            B['C1'].insert(row)
    cache = {}
    def ghat(op, codes):
        if (op, codes) not in cache:
            cache[(op, codes)] = C.ghat(op, codes)
        return cache[(op, codes)]
    rounds = 0; changed = True
    while changed:
        changed = False; rounds += 1
        for op, (ar, ins, out, _) in C.OPS.items():
            if READOUT != 'T' and op == 'close':
                continue                      # no separate whole-tree type: the answer is read at C1
            outrows = list(B[out].piv.values())
            if not outrows: continue
            if ar == 1:
                src = ins[0]
                for phi in outrows:
                    row = {}
                    for i, t in enumerate(sh[src]):
                        r = ghat(op, (t,))
                        if r is not None and r in idx[out] and phi.get(idx[out][r]):
                            row[i] = phi[idx[out][r]]
                    if row and B[src].insert(row): changed = True
            else:
                s1, s2 = ins
                for phi in outrows:
                    # pull back along position 1, other subtree = e_{t2}
                    for d2, ts2 in bydeg[s2].items():
                        for t2 in ts2:
                            row = {}
                            for d1 in range(0, 5 - d2):
                                for t1 in bydeg[s1].get(d1, []):
                                    r = ghat(op, (t1, t2))
                                    if r is not None and r in idx[out] and phi.get(idx[out][r]):
                                        row[idx[s1][t1]] = phi[idx[out][r]]
                            if row and B[s1].insert(row): changed = True
                    for d1, ts1 in bydeg[s1].items():
                        for t1 in ts1:
                            row = {}
                            for d2 in range(0, 5 - d1):
                                for t2 in bydeg[s2].get(d2, []):
                                    r = ghat(op, (t1, t2))
                                    if r is not None and r in idx[out] and phi.get(idx[out][r]):
                                        row[idx[s2][t2]] = phi[idx[out][r]]
                            if row and B[s2].insert(row): changed = True
        if verbose:
            print('  round', rounds, {s: B[s].dim() for s in SORTS}, flush=True)
    dims = {s: B[s].dim() for s in SORTS}
    prof = {s: [0]*5 for s in SORTS}
    for s in SORTS:
        for col in B[s].piv:
            prof[s][C.degree(sh[s][col])] += 1
    if return_basis:
        return dims, prof, B
    return dims, prof

if __name__ == '__main__':
    for K in (2, 4):
        t0 = time.time()
        fdims = forward(K)
        print(f'k={K}: forward space F_X is the full shape space, dims', fdims, f'({time.time()-t0:.1f}s)')
        t0 = time.time()
        bdims, prof = backward(K)
        print(f'k={K}: backward space dims', bdims, 'degree profiles', prof, f'({time.time()-t0:.1f}s)')
        print(f'k={K}: minimal dimensions m_X = rank(F_X B_X) = dim B_X:', bdims)
