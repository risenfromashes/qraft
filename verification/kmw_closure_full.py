"""KMW17 Sec. 3 on the census automaton for the full sorted alphabet with the sorts
    L, I, I_rho, C1 (closed C), C2 (open C), C3 (contains the root, open below), G, R (whole tree).
Read-out at R.  Shapes are stored as canonical codes of the QRAFT package (clusters.py):
marks 'T' (top boundary vertex), 'B' (bottom boundary vertex), leaves ('L', c).
A C3 shape has a 'B' and no 'T'; an R shape has no marks; the I_rho shape is the single vertex 'B'.

  forward: reachability of every shape from the leaf atoms (proof of Theorem 4.2, 0/1 weights)
  backward: closure of the gamma^a under the one-step operators M^sigma_{.,pi'}, M^sigma_{pi',.}
            (the pull-backs in the proof of Theorem 4.2), exact arithmetic mod 2^31-1
Prints dim B_X and the leaf-count profile for every sort.  ~1 s (k=2), ~30 s (k=4)."""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / 'natwta' / 'code'))
import clusters as C
from kmw_closure import RowSpace, P

SORTS = ['L', 'I', 'Irho', 'C1', 'C2', 'C3', 'G', 'R']
# symbol: (argument sorts, result sort)
SIG = {'LC': (('L',), 'C1'), 'CG': (('C1',), 'G'), 'GG': (('G', 'G'), 'G'),
       'IG1': (('I', 'G'), 'C1'), 'IG2': (('I', 'G'), 'C2'), 'IG3': (('Irho', 'G'), 'C3'), 'IGR': (('Irho', 'G'), 'R'),
       'CC1': (('C2', 'C1'), 'C1'), 'CC2': (('C2', 'C2'), 'C2'), 'CC3': (('C3', 'C2'), 'C3'), 'CCR': (('C3', 'C1'), 'R')}

def to_c3(code):
    """G shape -> C3 shape: the super root becomes the bottom boundary vertex"""
    g = C.graph_from_code(code)
    return C.canon(C.relabel_kind(g, C.find_kind(g, 'T'), 'B'))

def shapes(sort, K):
    if sort == 'L':    return [s for m in (0, 1) for s in C.fragments('C1', m, K)]
    if sort == 'I':    return [C.canon((['T', 'B'], [{1}, {0}]))]
    if sort == 'Irho': return [C.canon((['B'], [set()]))]
    if sort == 'C3':   return [to_c3(s) for m in range(5) for s in C.fragments('G', m, K)]
    if sort == 'R':    return [s for m in range(5) for s in C.fragments('T', m, K)]
    return [s for m in range(5) for s in C.fragments(sort, m, K)]

def ghat(sym, codes):
    """the transition on shapes: glue at the boundary vertices, re-mark, reduce; None if > 4 leaves"""
    if sum(C.degree(c) for c in codes) > 4: return None
    if sym in ('LC', 'CG'): return codes[0]
    if sym in ('IG1', 'IG2', 'IG3', 'IGR'):
        gI = C.graph_from_code(codes[0]); gG = C.graph_from_code(codes[1])
        u = C.find_kind(gI, 'B'); v = C.find_kind(gG, 'T')
        new = {'IG1': 'I', 'IG2': 'B', 'IG3': 'B', 'IGR': 'I'}[sym]
        g, _, _ = C.glue(gI, gG, [(u, v)], {(u, v): new})
        return C.canon(C.reduce_marked(g))
    gA = C.graph_from_code(codes[0]); gB = C.graph_from_code(codes[1])
    if sym == 'GG':
        u = C.find_kind(gA, 'T'); v = C.find_kind(gB, 'T')
        g, _, _ = C.glue(gA, gB, [(u, v)], {(u, v): 'T'})
    else:   # CC1, CC2, CC3, CCR: bottom of the upper part = top of the lower part, becomes internal
        u = C.find_kind(gA, 'B'); v = C.find_kind(gB, 'T')
        g, _, _ = C.glue(gA, gB, [(u, v)], {(u, v): 'I'})
    return C.canon(C.reduce_marked(g))

def run(K, verbose=True):
    sh = {s: shapes(s, K) for s in SORTS}
    idx = {s: {c: i for i, c in enumerate(sh[s])} for s in SORTS}
    bydeg = {s: {} for s in SORTS}
    for s in SORTS:
        for c in sh[s]: bydeg[s].setdefault(C.degree(c), []).append(c)
    # ---- forward: reachability of every shape (unit vectors of the 0/1-weighted leaf atoms)
    reach = {s: set() for s in SORTS}
    for c in range(K):
        for w in (0, 1): reach['L'] |= set(C.atom_state('leaf', c, w).keys())
    reach['I'] = set(sh['I']); reach['Irho'] = set(sh['Irho'])
    done = set(); changed = True
    while changed:
        changed = False
        for sym, (ins, out) in SIG.items():
            bd = [{} for _ in ins]
            for k, a in enumerate(ins):
                for c in reach[a]: bd[k].setdefault(C.degree(c), []).append(c)
            if len(ins) == 1:
                combos = [(c,) for cs in bd[0].values() for c in cs]
            else:
                combos = [(a, b) for d1, cs1 in bd[0].items() for d2, cs2 in bd[1].items()
                          if d1 + d2 <= 4 for a in cs1 for b in cs2]
            for codes in combos:
                if (sym, codes) in done: continue
                done.add((sym, codes))
                r = ghat(sym, codes)
                if r is not None and r not in reach[out]:
                    reach[out].add(r); changed = True
    for s in SORTS:
        assert reach[s] == set(sh[s]), (s, len(reach[s]), len(sh[s]))
    if verbose: print(f'k={K}: every shape is reachable; |Sh_X| =', {s: len(sh[s]) for s in SORTS})
    # ---- backward closure seeded at R
    ch = C.qraft_channels(K)
    B = {s: RowSpace(len(sh[s])) for s in SORTS}
    for name, d in ch.items():
        B['R'].insert({idx['R'][pat]: 1 for pat in d if pat in idx['R']})
    cache = {}
    def gh(sym, codes):
        if (sym, codes) not in cache: cache[(sym, codes)] = ghat(sym, codes)
        return cache[(sym, codes)]
    rounds = 0; changed = True
    while changed:
        changed = False; rounds += 1
        for sym, (ins, out) in SIG.items():
            outrows = list(B[out].piv.values())
            if not outrows: continue
            if len(ins) == 1:
                src = ins[0]
                for phi in outrows:
                    row = {}
                    for i, t in enumerate(sh[src]):
                        r = gh(sym, (t,))
                        if r is not None and phi.get(idx[out].get(r, -1)):
                            row[i] = phi[idx[out][r]]
                    if row and B[src].insert(row): changed = True
            else:
                s1, s2 = ins
                for phi in outrows:
                    for d2, ts2 in bydeg[s2].items():
                        for t2 in ts2:
                            row = {}
                            for d1 in range(0, 5 - d2):
                                for t1 in bydeg[s1].get(d1, []):
                                    r = gh(sym, (t1, t2))
                                    if r is not None and phi.get(idx[out].get(r, -1)):
                                        row[idx[s1][t1]] = phi[idx[out][r]]
                            if row and B[s1].insert(row): changed = True
                    for d1, ts1 in bydeg[s1].items():
                        for t1 in ts1:
                            row = {}
                            for d2 in range(0, 5 - d1):
                                for t2 in bydeg[s2].get(d2, []):
                                    r = gh(sym, (t1, t2))
                                    if r is not None and phi.get(idx[out].get(r, -1)):
                                        row[idx[s2][t2]] = phi[idx[out][r]]
                            if row and B[s2].insert(row): changed = True
        if verbose: print('  pass', rounds, {s: B[s].dim() for s in SORTS}, flush=True)
    dims = {s: B[s].dim() for s in SORTS}
    prof = {s: [0] * 5 for s in SORTS}
    for s in SORTS:
        for col in B[s].piv: prof[s][C.degree(sh[s][col])] += 1
    return sh, idx, B, dims, prof

if __name__ == '__main__':
    for K in (2, 4):
        t0 = time.time()
        sh, idx, B, dims, prof = run(K)
        print(f'k={K}: dim B_X = m_X:', dims)
        print(f'k={K}: leaf-count profiles:', {s: prof[s] for s in SORTS}, f'({time.time()-t0:.1f}s)')
