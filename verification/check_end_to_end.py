"""End-to-end check of Theorem 4.1(a): the natural scheme, with the update rules
mu~(sigma) = (F~ (x) F~) mu(sigma) B of the quotient lemma (F~ = unit rows at the witness shapes),
run on the HDTs that Brodal's construction (brodal_hdt.py) builds for random trees, with random
colours and weights, gives the weighted quartet counts computed by brute force."""
import sys, random, itertools, time
from fractions import Fraction
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE / 'natwta'), str(HERE / 'natwta' / 'code')]
import clusters as C
import kmw_closure_full as KF
from check_natural_closure import natural_counters
from brodal_hdt import Tree, build_hdt, check_invariants
from lower_bound import split, typ, channels

def scheme(K):
    sh = {s: KF.shapes(s, K) for s in KF.SORTS}
    nat, ch = natural_counters(K, sh)
    idx = {s: {c: i for i, c in enumerate(sh[s])} for s in KF.SORTS}
    wit = {}
    for s in KF.SORTS:
        wit[s] = []
        for lab, vec in nat[s]:
            (i,) = [i for i, v in vec.items() if v == 1 and all(i not in v2 for l2, v2 in nat[s] if l2 != lab)][:1]
            wit[s].append(sh[s][i])
    def val(sort, code, l):
        if code is None: return 0
        return nat[sort][l][1].get(idx[sort].get(code, -1), 0)
    rules = {}
    for sym, (ins, out) in KF.SIG.items():
        tab = {}
        for combo in itertools.product(*[range(len(wit[s])) for s in ins]):
            codes = tuple(wit[s][i] for s, i in zip(ins, combo))
            r = KF.ghat(sym, codes)
            for l in range(len(nat[out])):
                v = val(out, r, l)
                if v: tab[combo + (l,)] = v
        rules[sym] = tab
    # leaf initial vectors in the L counters: census e_empty + w e_c, projected
    def leaf(c, w):
        e0 = [s for s in sh['L'] if C.degree(s) == 0][0]
        ec = [s for s in sh['L'] if C.degree(s) == 1 and next(iter(C.colours_of(s))) == c][0]
        return [val('L', e0, l) + w * val('L', ec, l) for l in range(len(nat['L']))]
    # read-out gamma~ = F~_R gamma: gamma evaluated at the witness shapes of R
    gam = {}
    for ai, (a, vals) in enumerate(ch.items()):
        gam[ai] = [vals.get(wit['R'][l], 0) for l in range(len(nat['R']))]
    dims = {s: len(nat[s]) for s in KF.SORTS}
    return rules, leaf, gam, dims

def apply(rule, vecs, dim_out):
    out = [Fraction(0)] * dim_out
    for key, v in rule.items():
        *ins, l = key
        p = Fraction(v)
        for vec, i in zip(vecs, ins):
            p *= vec[i]
            if not p: break
        out[l] += p
    return out

def run_hdt(T, hdt, root, col, wt, rules, leaf, dims):
    state = {}
    def last_sort(h): return hdt[h].sorts[-1][0]
    order = sorted(range(len(hdt)), key=lambda h: len(hdt[h].nodes))
    for h in order:
        node = hdt[h]
        if not node.kids:
            (v,) = node.nodes
            if T.ch[v]:
                cur = ('Irho' if v == 0 else 'I', [Fraction(1)])
            else:
                cur = ('L', leaf(col[v], wt[v]))
        else:
            a, b = node.kids
            sa, va = state[a]; sb, vb = state[b]
            res = node.sorts[0][0]
            if sa in ('I', 'Irho'):
                sym = {'C1': 'IG1', 'C2': 'IG2', 'C3': 'IG3', 'R': 'IGR'}[res]
            elif sa == 'G' and sb == 'G':
                sym = 'GG'
            else:
                sym = {'C1': 'CC1', 'C2': 'CC2', 'C3': 'CC3', 'R': 'CCR'}[res]
            assert KF.SIG[sym][0] == (sa, sb), (sym, sa, sb)
            cur = (res, apply(rules[sym], [va, vb], dims[res]))
        # transformations recorded after the first sort of this node
        for s, m, r in node.sorts[1:] if not node.kids else node.sorts[1:]:
            sym = {('L', 'C1'): 'LC', ('C1', 'G'): 'CG'}[(cur[0], s)]
            cur = (s, apply(rules[sym], [cur[1]], dims[s]))
        state[h] = cur
    return state[root]

def brute(T, col, wt, K):
    chans = [tuple(sorted(a)) for a in channels(K)]
    res = [Fraction(0)] * len(chans)
    for q in itertools.combinations(T.leaves, 4):
        s = split(T, q)
        if s is None: continue
        t = typ(tuple(col[l] for l in q), s)
        for ai, a in enumerate(chans):
            if t == a:
                w = Fraction(1)
                for l in q: w *= wt[l]
                res[ai] += w
    return res

def random_tree(rng, n_internal, maxdeg):
    ch = {0: []}; nxt = 1; internal = [0]
    while len(internal) < n_internal:
        v = rng.choice(internal)
        if len(ch[v]) >= maxdeg: continue
        ch[v].append(nxt); ch[nxt] = []; internal.append(nxt); nxt += 1
    for v in list(internal):
        while len(ch[v]) < 2 or (len(ch[v]) < maxdeg and rng.random() < 0.35):
            ch[v].append(nxt); ch[nxt] = []; nxt += 1
    return {v: c for v, c in ch.items() if c}

if __name__ == '__main__':
    for K in (2, 4):
        t0 = time.time()
        rules, leaf, gam, dims = scheme(K)
        rng = random.Random(7 + K)
        n = 0; sorts_seen = set(); choice_trees = 0
        for trial in range(60):
            T = Tree(random_tree(rng, rng.randint(1, 14), rng.choice([2, 3, 4, 5])))
            hdt, root, rounds, choices = build_hdt(T)
            check_invariants(T, hdt, root)
            choice_trees += choices > 0
            for h in hdt: sorts_seen |= {s for s, m, r in h.sorts}
            col = {l: rng.randrange(K) for l in T.leaves}
            wt = {l: rng.choice([Fraction(0), Fraction(1), Fraction(1), Fraction(2), Fraction(-1), Fraction(2, 3)]) for l in T.leaves}
            sort, vec = run_hdt(T, hdt, root, col, wt, rules, leaf, dims)
            assert sort == 'R'
            got = [sum(v * g for v, g in zip(vec, gam[ai])) for ai in range(len(gam))]
            assert got == brute(T, col, wt, K), (K, trial, got, brute(T, col, wt, K))
            n += 1
        print(f'k={K}: natural scheme = brute force on the HDTs of {n} random trees '
              f'({choice_trees} with pairing choices), sorts met {sorted(sorts_seen)}: ok ({time.time()-t0:.1f}s)')
