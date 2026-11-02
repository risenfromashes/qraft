"""Rules for the sorts C_3 (contains the root, downward open) and R (the whole tree),
k = 2 or 4, in Brodal's notation.  Shapes of C_3 are the G shapes with the bottom vertex b
in the role of the super root (the mark 'T' of a G-code plays the role of b here).
Counters of C_3: 1, n_i, n_[ij], n_{down[xbar(yy')]}, q_a;  of R: q_a.

Symbols:  IG_3 : I G -> C_3   (shape: identity, b = the super root)
          CC_3 : C_3 C_2 -> C_3   (glue b with the top of C_2; new b = its bottom)
          CC_R : C_3 C_1 -> R     (glue b with the top of C_1; closed pattern)
"""
import sys, json
from collections import defaultdict
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE / 'natwta'), str(HERE / 'natwta' / 'code')]
import clusters as C
from natural_wta import NaturalWTA, unrooted
from gen_brodal_rules import sym, triple_pattern   # reuse the Brodal-notation printer

K = int(sys.argv[1]) if len(sys.argv) > 1 else 2
OUT = sys.argv[2] if len(sys.argv) > 2 else f'brodal_rules_c3_k{K}.tex'
model = NaturalWTA(K)
ch = model.channels
pairs = [(0, 0), (0, 1), (1, 1)] if K == 2 else [(i, j) for i in range(K) for j in range(i + 1, K)]

# ---- coordinates of C_rho
labels = [('unit',)] + [('leaf', c) for c in range(K)] + [('pair', i, j, 'all') for i, j in pairs] \
         + [('triple', 'B', c, a) for c in range(K) for a in ch] + [('quartet', a) for a in ch]
index = {l: i for i, l in enumerate(labels)}
RHO = 'C'

def features_rho(code):
    if code is None: return {}
    d = C.degree(code)
    if d == 0: return {index[('unit',)]: 1}
    if d == 1: return {index[('leaf', next(iter(C.colours_of(code))))]: 1}
    if d == 2:
        kinds, adj = C.graph_from_code(code)
        cols = tuple(sorted(k[1] for k in kinds if C.is_leaf_kind(k)))
        return {index[('pair', cols[0], cols[1], 'all')]: 1} if ('pair', cols[0], cols[1], 'all') in index else {}
    if d == 3:
        f = {}
        for c in range(K):
            closed = unrooted(code, 'T', c)            # test leaf at b (the mark 'T' of the G-code)
            for a, vals in ch.items():
                if vals.get(closed, 0): f[index[('triple', 'B', c, a)]] = 1
        return f
    closed = unrooted(code)
    return {index[('quartet', a)]: 1 for a, vals in ch.items() if vals.get(closed, 0)}

shapes = [s for m in range(5) for s in C.fragments('G', m, K)]
wit = [None] * len(labels)
for s in shapes:
    f = features_rho(s)
    if len(f) == 1:
        (i, v), = f.items()
        if v == 1 and wit[i] is None: wit[i] = s
assert all(w is not None for w in wit), [labels[i] for i, w in enumerate(wit) if w is None]

def glue_rho(code_rho, code_c, result_open):
    gA = C.graph_from_code(code_rho); gB = C.graph_from_code(code_c)
    u = C.find_kind(gA, 'T'); v = C.find_kind(gB, 'T')
    g, _, m2 = C.glue(gA, gB, [(u, v)], {(u, v): 'I'})
    kinds, adj = g
    if result_open:                       # the bottom of the lower part becomes the new b, i.e. mark 'T'
        b = kinds.index('B'); kinds[b] = 'T'
    return C.canon(C.reduce_marked((kinds, adj)))

def rules(name_in2, sort2, result, out):
    """binary CC->C with C_rho first; result 'rho' or 'R'"""
    entries = defaultdict(list)
    for i, wi in enumerate(wit):
        for j, wj in enumerate(model.witnesses[sort2]):
            if C.degree(wi) + C.degree(wj) > 4: continue
            r = glue_rho(wi, wj, result == 'rho')
            if r is None: continue
            if result == 'rho':
                for t in features_rho(r): entries[t].append((i, j))
            else:
                closed = unrooted(r)
                for a, vals in ch.items():
                    if vals.get(closed, 0): entries[('quartet', a)].append((i, j))
    lines = []
    if result == 'rho':
        targets = [(t, labels[t]) for t in range(len(labels)) if labels[t][0] != 'unit']
    else:
        targets = [(('quartet', a), ('quartet', a)) for a in ch]
    for t, lab in targets:
        lhs = sym(lab, RHO if result == 'rho' else 'R')
        terms = []
        for (i, j) in entries.get(t, []):
            a = sym(labels[i], "C'"); b = sym(model.labels[sort2][j], "C''")
            terms.append(b if a == '1' else a if b == '1' else a + '\\,' + b)
        rhs = ' \\\\\n &\\quad + '.join(' + '.join(terms[k:k+3]) for k in range(0, len(terms), 3)) if terms else '0'
        lines.append(f'{lhs} &= {rhs}')
    return '\\begin{align*}\n' + ' \\\\\n'.join(lines) + '\n\\end{align*}'

out = []
# IG->C, I is the root with one child outside: G -> C_rho, shape identity
lines = []
for t in range(1, len(labels)):
    f = [j for j, wj in enumerate(model.witnesses['G']) if features_rho(wj).get(t)]
    terms = [sym(model.labels['G'][j], 'G') for j in f]
    lhs = sym(labels[t], RHO)
    lines.append(lhs + ' &= ' + (' + '.join(terms) if terms else '0'))
out.append(r'\subsubsection*{$\mathrm{IG}\to\mathrm{C}$ at the root, result $\mathrm{C}_3$ (one child outside)}')
out.append('\\begin{align*}\n' + ' \\\\\n'.join(lines) + '\n\\end{align*}')
out.append(r'\subsubsection*{$\mathrm{CC}\to\mathrm{C}$, upper part $\mathrm{C}_3$, lower part downward open: result $\mathrm{C}_3$}')
out.append(rules('C_2', 'C2', 'rho', out))
out.append(r'\subsubsection*{$\mathrm{CC}\to\mathrm{C}$, upper part $\mathrm{C}_3$, lower part downward closed: result $\mathrm{R}$}')
out.append(rules('C_2', 'C1', 'R', out))
open(OUT, 'w').write('\n'.join(out) + '\n')
print('wrote', OUT, 'C_3 dimension', len(labels))

# ---- brute-force check of the C_rho rules on random terms:
#      build a C_rho term as a chain: root-open IG->C on a random G, then cat with random C2's,
#      compare the rule-based counters with the counters computed directly from the fragment vector.
import random
from fractions import Fraction
rng = random.Random(7)
def project(raw):
    v = [Fraction(0)] * len(labels)
    for code, w in raw.items():
        for i, x in features_rho(code).items(): v[i] += w * x
    return v
def apply_rules(entries_fn, vec1, vec2, labels2):
    pass
checks = 0
for _ in range(40):
    _, clusterG, rawG = C.random_term('G', K, rng)
    raw_rho = dict(rawG)                                   # IG->C root-open: shapes unchanged
    v_rho = project(raw_rho)
    for step in range(rng.randint(1, 3)):
        _, clusterC, rawC = C.random_term('C2', K, rng)
        # exact new state by gluing shapes
        new_raw = {}
        for s1, x1 in raw_rho.items():
            for s2, x2 in rawC.items():
                if C.degree(s1) + C.degree(s2) > 4: continue
                r = glue_rho(s1, s2, True)
                if r is not None: new_raw[r] = new_raw.get(r, 0) + x1 * x2
        # rule-based state: tensor from witnesses
        v2 = model.project('C2', rawC)
        v_new = [Fraction(0)] * len(labels)
        for i, wi in enumerate(wit):
            for j, wj in enumerate(model.witnesses['C2']):
                if C.degree(wi) + C.degree(wj) > 4: continue
                r = glue_rho(wi, wj, True)
                for t in features_rho(r): v_new[t] += v_rho[i] * v2[j]
        assert v_new == project(new_raw), 'rule mismatch'
        raw_rho, v_rho = new_raw, v_new; checks += 1
    # close with a C1
    _, clusterC1, rawC1 = C.random_term('C1', K, rng)
    v1 = model.project('C1', rawC1)
    q_rule = {a: Fraction(0) for a in ch}
    for i, wi in enumerate(wit):
        for j, wj in enumerate(model.witnesses['C1']):
            if C.degree(wi) + C.degree(wj) > 4: continue
            r = glue_rho(wi, wj, False)
            closed = unrooted(r)
            for a, vals in ch.items():
                if vals.get(closed, 0): q_rule[a] += v_rho[i] * v1[j]
    q_exact = {a: Fraction(0) for a in ch}
    for s1, x1 in raw_rho.items():
        for s2, x2 in rawC1.items():
            if C.degree(s1) + C.degree(s2) > 4: continue
            closed = unrooted(glue_rho(s1, s2, False))
            for a, vals in ch.items():
                if vals.get(closed, 0): q_exact[a] += x1 * x2
    assert q_rule == q_exact, 'root rule mismatch'
    checks += 1
print('C_3 / R rules verified on', checks, 'random compositions')
