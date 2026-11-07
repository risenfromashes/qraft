#!/usr/bin/env python3
"""Adversarial verification of shapes.py.

A passing exactness test does NOT by itself prove the canonical form is a true
graph-isomorphism invariant, nor that the observable-quotient rank is the exact
dimension of the observable subspace. This script attacks both, plus the
brute-force oracle itself, with independent methods:

  T1  Extended exactness battery: adversarial trees (stars, all-one-color,
      deep chains, near-degenerate), up to 9 leaves, integer+fractional, both
      decompositions, thousands of trials.
  T2  INDEPENDENT oracle: reclassify every quartet by the distance four-point
      condition (Buneman) and require it to agree with the path-disjointness
      brute force AND with the DP -- three independent methods must coincide.
  T3  Canonical-form invariance: relabel each reduced fragment by a random node
      permutation and require an identical code (no false split); confirm reduce
      is idempotent (canonical is a fixed point).
  T4  DECISIVE quotient check -- CLOSURE: for every op and every W_out basis
      functional, its pullback to the operand interface(s) already lies in
      span(W_in) (Basis.in_span, non-mutating). If any pullback escaped the span,
      the observable subspace would not be closed and the reported rank would be
      too small. Closed + generated from the readout => W_T is exactly the
      observable subspace => the rank is the exact payload, not an
      underestimate.
  T5  Exhaustive algebra laws: every reachable degree-bounded G triple and
      every legal C corridor triple is parenthesization-invariant; GG is
      permutation-invariant and every CC transition agrees in both corridor
      orientations.

Run:  python3 tools/counting/check_counting_scheme.py
"""
from __future__ import annotations
import itertools, random, sys
from collections import defaultdict, deque
from fractions import Fraction

import shapes as G
from shapes import (Frag, U, D, Z, I, LEAF, reduce_and_code, interface_of,
                             sel_count, _channel, _four_leaf_split, readout,
                             eval_center_fold, eval_cc_spine, brute_force,
                             random_tree, caterpillar, enumerate_full, quotient,
                             _canonical, REG, state_sort_key,
                             verify_state_completeness)

FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        if len(FAIL) <= 20:
            print("  FAIL:", msg)


# ---------------------------------------------------------------------------
# T1: extended exactness battery incl. adversarial shapes
# ---------------------------------------------------------------------------

def star_tree(nlv, K, rng, frac):
    t = G.Tree(); ctr = itertools.count(); r = next(ctr); t.root = r; t.children[r] = []
    for _ in range(nlv):
        lf = next(ctr)
        w = Fraction(rng.randint(1, 5)) if not frac else Fraction(rng.randint(1, 9), rng.randint(1, 4))
        t.leafinfo[lf] = (rng.randrange(K), w); t.children[r].append(lf)
    return t


def deep_chain(nlv, K, rng, frac):
    # every internal node has one leaf and one internal child, plus 2 leaves at tip
    return caterpillar(nlv, K, rng, frac=frac)


def onecolor_tree(nlv, K, rng, frac):
    t = random_tree(nlv, K, rng, allow_poly=True, frac=frac)
    for lf in t.leafinfo:
        c, w = t.leafinfo[lf]; t.leafinfo[lf] = (0, w)  # collapse to one color
    return t


def test_t1_exactness():
    print("=== T1: extended exactness (DP == brute force) ===")
    rng = random.Random(777)
    trials = 0
    corridor_trials = 0
    for K in (2, 4):
        for frac in (False, True):
            gens = [(lambda n: random_tree(n, K, rng, True, frac), False),
                    (lambda n: random_tree(n, K, rng, False, frac), False),
                    (lambda n: caterpillar(n, K, rng, frac=frac), True),
                    (lambda n: star_tree(n, K, rng, frac), False),
                    (lambda n: onecolor_tree(n, K, rng, frac), False)]
            for gen, check_corridor in gens:
                for _ in range(120):
                    nlv = rng.randint(4, 9)
                    t = gen(nlv)
                    bf = dict(brute_force(t, K))
                    cf = dict(readout(eval_center_fold(t), K))
                    trials += 1
                    check(bf == cf, f"center-fold K={K} frac={frac} n={nlv}: {bf} vs {cf}")
                    if check_corridor:
                        cc = dict(readout(eval_cc_spine(t), K))
                        corridor_trials += 1
                        check(bf == cc,
                              f"CC-corridor K={K} frac={frac} n={nlv}: "
                              f"{bf} vs {cc}")
    print(f"  center-fold trials={trials} CC-corridor trials={corridor_trials}")


# ---------------------------------------------------------------------------
# T2: independent four-point (distance) oracle
# ---------------------------------------------------------------------------

def _adj_of(tree):
    adj = defaultdict(set)
    def walk(v):
        for c in tree.children[v]:
            adj[v].add(c); adj[c].add(v); walk(c)
    walk(tree.root)
    return adj


def _dists(adj, src, targets):
    d = {src: 0}; dq = deque([src]); need = set(targets)
    while dq and need:
        x = dq.popleft()
        for y in adj[x]:
            if y not in d:
                d[y] = d[x] + 1; dq.append(y)
                need.discard(y)
    return d


def brute_force_fourpoint(tree, K):
    """Independent oracle: resolution by the tree four-point condition on
    unit-edge distances (Buneman). ab|cd resolved iff d(a,b)+d(c,d) is the
    strict minimum of the three pair-sums; all equal => star (unresolved)."""
    adj = _adj_of(tree)
    leaves = list(tree.leafinfo)
    ch = defaultdict(Fraction)
    for quad in itertools.combinations(leaves, 4):
        dm = {}
        for a in quad:
            da = _dists(adj, a, quad)
            for b in quad:
                dm[(a, b)] = da[b]
        a, b, c, d = quad
        pairings = [((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))]
        sums = [(dm[p1] + dm[p2], (p1, p2)) for (p1, p2) in pairings]
        sums.sort()
        w = Fraction(1)
        for l in quad:
            w *= Fraction(tree.leafinfo[l][1])
        if sums[0][0] < sums[1][0]:               # strict min => resolved
            (p1, p2) = sums[0][1]
            cols = tree.leafinfo
            s1 = tuple(sorted((cols[p1[0]][0], cols[p1[1]][0])))
            s2 = tuple(sorted((cols[p2[0]][0], cols[p2[1]][0])))
            key = _channel(s1, s2, K)
            if key:
                ch[key] += w
        # else all three equal within the group => unresolved star (no channel)
    return ch


def test_t2_independent_oracle():
    print("=== T2: independent four-point oracle == path-disjoint brute force == DP ===")
    rng = random.Random(4242)
    n = 0
    for K in (2, 4):
        for frac in (False, True):
            for _ in range(300):
                nlv = rng.randint(4, 9)
                kind = rng.choice(['poly', 'binary', 'cat'])
                if kind == 'cat':
                    t = caterpillar(nlv, K, rng, frac=frac)
                else:
                    t = random_tree(nlv, K, rng, allow_poly=(kind == 'poly'), frac=frac)
                a = dict(brute_force(t, K))          # path-disjoint
                b = dict(brute_force_fourpoint(t, K))  # distance four-point
                c = dict(readout(eval_center_fold(t), K))  # DP
                n += 1
                check(a == b, f"four-point != path-disjoint K={K}: {a} vs {b}")
                check(a == c, f"DP != brute force K={K}: {a} vs {c}")
    print(f"  trials={n}")


# ---------------------------------------------------------------------------
# T3: canonical-form invariance under relabeling + idempotence
# ---------------------------------------------------------------------------

def _relabel_perm(frag, rng):
    ids = list(frag.kind)
    perm = ids[:]; rng.shuffle(perm)
    m = {old: new for old, new in zip(ids, perm)}
    f = Frag()
    for nid, tok in frag.kind.items():
        f.kind[m[nid]] = tok
    for nid, nbs in frag.adj.items():
        f.adj[m[nid]] = set(m[x] for x in nbs)
    return f


def test_t3_canonical():
    print("=== T3: canonical code invariant under relabeling + idempotent ===")
    rng = random.Random(9)
    n = 0
    states = 0
    for K in (2, 4):
        REG.clear()
        by, _un, _gg, _cc = enumerate_full(K)
        for interface in ('Closed', 'C1', 'C2', 'G'):
            for code in sorted(by[interface], key=state_sort_key):
                states += 1
                frag = REG[code]
                # Idempotence: recanonicalizing the canonical fragment reproduces
                # the code. Relabel every state, not a prefix sample.
                code2, _ = reduce_and_code(frag)
                check(code2 == code, f"not idempotent: {code} -> {code2}")
                for _ in range(3):
                    rf = _relabel_perm(frag, rng)
                    rc, _ = reduce_and_code(rf)
                    check(rc == code, f"relabel changed code: {code} -> {rc}")
                    n += 1
    print(f"  states={states} relabelings checked={n}")


# ---------------------------------------------------------------------------
# T4: quotient congruence (decisive: proves the rank is the true minimum)
# ---------------------------------------------------------------------------

def _pull_unary_row(inv, f):
    row = defaultdict(Fraction)
    for oc, val in f.items():
        for x in inv.get(oc, ()):
            row[x] += val
    return row


def test_t4_closure(K):
    """DECISIVE: verify the observable subspace is CLOSED under every op's
    pullback -- for each op and each W_out basis functional, its pullback to the
    operand interface(s) already lies in span(W_in).  Closed + generated from the
    readout  =>  W_T is exactly the observable subspace  =>  rank is the true
    payload (not an underestimate).  Uses the same interface bucketing as
    the generator, and Basis.in_span (non-mutating)."""
    print(f"=== T4: observable-subspace CLOSURE under every op (K={K}) ===")
    REG.clear()
    by, un, gg, cc = enumerate_full(K)
    counts = verify_state_completeness(K, by)
    ifaces = ('Closed', 'C1', 'C2', 'G')
    ranks, Bq = quotient(K, ifaces, by['Closed'],
                         [tbl for (tbl,) in un.values()], [gg, cc])

    # sanity: readout functionals lie in W_Closed
    for r in G.readout_rows(K, by['Closed']):
        check(Bq['Closed'].in_span(r), "readout functional not in W_Closed")

    checked = 0
    # unary buckets (in-iface) : {out_code -> [x]}
    for tbl in [t for (t,) in un.values()]:
        by_in = defaultdict(lambda: defaultdict(list))
        for x, ox in tbl.items():
            by_in[interface_of(x)][ox].append(x)
        for in_if, inv in by_in.items():
            out_ifs = {interface_of(oc) for oc in inv}
            for out_if in out_ifs:
                for f in Bq[out_if].rows:
                    row = _pull_unary_row(inv, f)
                    row = {c: v for c, v in row.items() if v != 0}
                    if row:
                        check(Bq[in_if].in_span(row),
                              f"unary pullback NOT in span: in={in_if} out={out_if} "
                              f"(=> rank {ranks[in_if]} is too small)")
                        checked += 1
    # binary buckets ((hif,lif)) : {out_code -> [(a,b)]}
    for tbl in (gg, cc):
        by_pair = defaultdict(lambda: defaultdict(list))
        for (a, b), oc in tbl.items():
            by_pair[(interface_of(a), interface_of(b))][oc].append((a, b))
        for (hif, lif), inv in by_pair.items():
            out_ifs = {interface_of(oc) for oc in inv}
            for out_if in out_ifs:
                for f in Bq[out_if].rows:
                    left = defaultdict(lambda: defaultdict(Fraction))
                    right = defaultdict(lambda: defaultdict(Fraction))
                    for oc, val in f.items():
                        for (a, b) in inv.get(oc, ()):
                            left[b][a] += val
                            right[a][b] += val
                    for _, row in left.items():
                        row = {c: v for c, v in row.items() if v != 0}
                        if row:
                            check(Bq[hif].in_span(row),
                                  f"binary(fix-right) pullback NOT in span: h={hif} "
                                  f"(=> rank {ranks[hif]} too small)")
                            checked += 1
                    for _, row in right.items():
                        row = {c: v for c, v in row.items() if v != 0}
                        if row:
                            check(Bq[lif].in_span(row),
                                  f"binary(fix-left) pullback NOT in span: l={lif} "
                                  f"(=> rank {ranks[lif]} too small)")
                            checked += 1
    print(f"  states={counts} ranks={ranks}  "
          f"({checked} pullbacks all in span => closed)")


# ---------------------------------------------------------------------------
# T5: exhaustive parenthesization, permutation, and corridor orientation
# ---------------------------------------------------------------------------

def test_t5_parenthesization(K):
    print(f"=== T5: exhaustive parenthesization/orientation laws (K={K}) ===")
    REG.clear()
    by, unary, gg, cc = enumerate_full(K)
    reverse = unary['REV'][0]

    gg_checks = 0
    for (left, right), output in gg.items():
        check(gg.get((right, left)) == output,
              f"GG permutation mismatch K={K}: {left}, {right}")
        gg_checks += 1

    g_by_degree = defaultdict(list)
    for state in by['G']:
        g_by_degree[sel_count(state)].append(state)
    for left_degree in range(5):
        for middle_degree in range(5 - left_degree):
            for right_degree in range(5 - left_degree - middle_degree):
                for left in g_by_degree[left_degree]:
                    for middle in g_by_degree[middle_degree]:
                        left_fold = gg.get((left, middle))
                        for right in g_by_degree[right_degree]:
                            right_fold = gg.get((middle, right))
                            lhs = gg.get((left_fold, right))
                            rhs = gg.get((left, right_fold))
                            check(lhs is not None and lhs == rhs,
                                  f"GG parenthesization mismatch K={K}")
                            gg_checks += 1

    cc_checks = 0
    for (higher, lower), output in cc.items():
        reversed_pair = cc.get((reverse[lower], reverse[higher]))
        expected = reverse.get(output, output)
        check(reversed_pair == expected,
              f"CC orientation mismatch K={K}: {higher}, {lower}")
        cc_checks += 1

    identities = [state for state in by['C2'] if sel_count(state) == 0]
    check(len(identities) == 1, f"K={K} expected one degree-zero C2 identity")
    identity = identities[0]
    for (higher, lower), output in cc.items():
        if lower == identity:
            check(output == higher, f"right C2 identity failed K={K}")
            cc_checks += 1
        if higher == identity:
            check(output == lower, f"left C2 identity failed K={K}")
            cc_checks += 1

    before = defaultdict(list)
    after = defaultdict(list)
    for (higher, lower), output in cc.items():
        before[lower].append((higher, output))
        after[higher].append((lower, output))
    for middle in before.keys() & after.keys():
        # Triples with the unique C2 identity are already covered exactly by
        # both identity laws above; avoiding their Cartesian product saves
        # ~19 million redundant checks for K=4.
        if middle == identity:
            continue
        for left, left_fold in before[middle]:
            for right, right_fold in after[middle]:
                lhs = cc.get((left_fold, right))
                rhs = cc.get((left, right_fold))
                check(lhs == rhs,
                      f"CC parenthesization mismatch K={K}: "
                      f"{left}, {middle}, {right}")
                cc_checks += 1
    print(f"  GG laws={gg_checks} CC laws={cc_checks}")


if __name__ == '__main__':
    test_t1_exactness()
    test_t2_independent_oracle()
    test_t3_canonical()
    test_t4_closure(2)
    test_t4_closure(4)
    test_t5_parenthesization(2)
    test_t5_parenthesization(4)
    print()
    if FAIL:
        print(f"ADVERSARIAL VERIFICATION FAILED: {len(FAIL)} issue(s)")
        sys.exit(1)
    print("ALL ADVERSARIAL CHECKS PASSED")
    sys.exit(0)
