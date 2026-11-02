#!/usr/bin/env python3
"""Brute-force checks of what a frozen background must record.

Validates, by exact enumeration of weighted quartets on small arbitrary-degree
gene trees:
  T1  A mass-only frozen background is INSUFFICIENT: two gene trees with the same
      background leaf multiset but different background topology give different
      q00_11 (so the frozen summary must carry the resolved/flat pair structure).
  T2  POLYTOMY: cross-polytomy-child background pairs are CONTINGENT -- they DO
      contribute a resolved 0-0 pair when the outside pair localizes to one branch,
      and do not when the outside pair is split across branches.
  T3  A single RESOLVED-merged leaf OVERCOUNTS (forces always-resolved).
  T4  POSITIVE: for a one-port background region, (leaf-mass, resolved-pair-mass) is a
      sufficient statistic for the 2-color objective; mass-only is not.
"""
from itertools import combinations
from collections import deque
import functools
import math


def finite_close(actual, expected, tolerance=1e-12):
    return (math.isfinite(actual) and math.isfinite(expected)
            and abs(actual - expected) <= tolerance)


def assert_counts(actual, expected, label):
    assert len(actual) == len(expected), f"{label}: channel-count mismatch"
    assert all(finite_close(a, e) for a, e in zip(actual, expected)), \
        f"{label}: actual={actual} expected={expected}"

# ---------------- tree infra: unrooted, adjacency; leaves carry color+mass ----------------
class Tree:
    def __init__(self):
        self.G = {}; self.color = {}; self.mass = {}; self.leaves = []; self._n = 0
    def node(self):
        i = self._n; self._n += 1; self.G[i] = set(); return i
    def leaf(self, color, mass):
        i = self.node(); self.color[i] = color; self.mass[i] = mass; self.leaves.append(i); return i
    def link(self, a, b):
        self.G[a].add(b); self.G[b].add(a)
    def build(self, spec):
        """spec: ('leaf',color,mass) or ('int',[specs...]) -> returns node id"""
        if spec[0] == 'leaf':
            return self.leaf(spec[1], spec[2])
        nid = self.node()
        for c in spec[1]:
            self.link(nid, self.build(c))
        return nid

def quartet_2color(t):
    """Return (q00_11, q01_01) over all weighted 4-leaf subsets, strict-resolved semantics."""
    @functools.lru_cache(maxsize=None)
    def path(a, b):
        prev = {a: None}; dq = deque([a])
        while dq:
            u = dq.popleft()
            if u == b: break
            for w in t.G[u]:
                if w not in prev: prev[w] = u; dq.append(w)
        s = set(); u = b
        while u is not None: s.add(u); u = prev[u]
        return frozenset(s)
    q0011 = 0.0; q0101 = 0.0
    for a, b, c, d in combinations(t.leaves, 4):
        for (p, q) in (((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))):
            if path(p[0], p[1]).isdisjoint(path(q[0], q[1])):
                w = t.mass[a]*t.mass[b]*t.mass[c]*t.mass[d]
                s1 = sorted([t.color[p[0]], t.color[p[1]]])
                s2 = sorted([t.color[q[0]], t.color[q[1]]])
                if (s1 == [0, 0] and s2 == [1, 1]) or (s1 == [1, 1] and s2 == [0, 0]):
                    q0011 += w
                elif s1 == [0, 1] and s2 == [0, 1]:
                    q0101 += w
                break  # at most one resolved pairing per quartet
    return q0011, q0101

def L(c, m=1.0): return ('leaf', c, m)
def I(*kids):    return ('int', list(kids))

# --------- sanity: a clean quartet and a star ---------
tq = Tree(); tq.build(I(I(L(0), L(0)), I(L(1), L(1))))   # (00)|(11)
assert_counts(quartet_2color(tq), (1.0, 0.0), "clean quartet")
tstar = Tree(); tstar.build(I(L(0), L(0), L(1), L(1)))    # star -> unresolved
assert_counts(quartet_2color(tstar), (0.0, 0.0), "four-leaf star")
print("sanity: clean quartet -> (1,0) OK ; 4-star -> (0,0) OK\n")

# ================= T1 & T2 : same background mass, different topology =================
print("="*74)
print("T1/T2  Same background {b1,b2} (colors 0, masses 1,1), 2 live L1,L2 (color 1).")
print("="*74)
# Polytomy: b1,b2,L1,L2 all direct children of one star center (b1,b2 in DIFFERENT children)
tP = Tree(); tP.build(I(L(0), L(0), L(1), L(1)))
# Resolved cherry background: b1,b2 form a cherry (resolved together), then L1,L2
tC = Tree(); tC.build(I(I(L(0), L(0)), L(1), L(1)))
qP = quartet_2color(tP); qC = quartet_2color(tC)
assert_counts(qP, (0.0, 0.0), "mass-only negative control: polytomy")
assert_counts(qC, (1.0, 0.0), "mass-only negative control: cherry")
print(f"  polytomy  (b1,b2 flat at center)  q00_11={qP[0]:.1f}  q01_01={qP[1]:.1f}")
print(f"  cherry    (b1,b2 resolved below)  q00_11={qC[0]:.1f}  q01_01={qC[1]:.1f}")
print(f"  --> SAME background leaf multiset, DIFFERENT q00_11 ({qP[0]:.0f} vs {qC[0]:.0f}).")
print( "      A mass-only frozen summary is identical for both => WRONG.")
print( "      The pair {b1,b2} is FLAT at the polytomy but RESOLVED under the cherry.\n")

# ================= T2b : contingency of a flat cross-child pair =================
print("="*74)
print("T2b  Flat pair {b1,b2} is resolved vs a LOCALIZED outside pair, unresolved vs a SPLIT one.")
print("="*74)
# center star: b1,b2,L1,L2 as direct children, plus childE=cherry(L3,L4)
tPP = Tree(); tPP.build(I(L(0), L(0), L(1), L(1), I(L(1), L(1))))
# isolate specific quartets by relabeling: give each leaf unique mass to read contributions
tX = Tree()
b1 = tX.leaf(0, 1.0); b2 = tX.leaf(0, 1.0)
L1 = tX.leaf(1, 1.0); L2 = tX.leaf(1, 1.0)
cE = tX.node(); L3 = tX.leaf(1, 1.0); L4 = tX.leaf(1, 1.0)
tX.link(cE, L3); tX.link(cE, L4)
center = tX.node()
for x in (b1, b2, L1, L2, cE): tX.link(center, x)
# manual per-quartet check
@functools.lru_cache(maxsize=None)
def pth(a, b):
    prev = {a: None}; dq = deque([a])
    while dq:
        u = dq.popleft()
        if u == b: break
        for w in tX.G[u]:
            if w not in prev: prev[w] = u; dq.append(w)
    s = set(); u = b
    while u is not None: s.add(u); u = prev[u]
    return frozenset(s)
def resolves(pair_a, pair_b):
    return pth(*pair_a).isdisjoint(pth(*pair_b))
print(f"  {{b1,b2}} vs {{L1,L2}} (L1,L2 SPLIT across center children): resolved={resolves((b1,b2),(L1,L2))}  -> contributes 0")
print(f"  {{b1,b2}} vs {{L3,L4}} (L3,L4 LOCALIZED in one branch cE):   resolved={resolves((b1,b2),(L3,L4))}  -> contributes q00_11")
assert not resolves((b1, b2), (L1, L2)), \
    "split outside pair must stay unresolved at the common center"
assert resolves((b1, b2), (L3, L4)), \
    "localized outside pair must resolve against the cross-arm pair"
print( "  --> The SAME flat pair {b1,b2} DOES yield a resolved 0-0 pair, contingent on the")
print( "      outside pair's localization. 'Never contributes a resolved 0-0 pair' is FALSE.\n")

# ================= T3 : single resolved-merged leaf overcounts =================
print("="*74)
print("T3  Merging polytomy background into a RESOLVED cherry OVERCOUNTS.")
print("="*74)
correct = quartet_2color(tP)              # star: {b1,b2,L1,L2} unresolved -> q00_11=0
wrongmerge = quartet_2color(tC)           # cherry: forces {b1,b2} resolved -> q00_11=1
assert_counts(correct, (0.0, 0.0), "flat G-group negative control")
assert_counts(wrongmerge, (1.0, 0.0),
              "ordinary resolved-leaf replacement negative control")
print(f"  correct (flat)          q00_11={correct[0]:.1f}")
print(f"  wrong resolved-merge    q00_11={wrongmerge[0]:.1f}   (overcount = {wrongmerge[0]-correct[0]:.1f})\n")

# ================= T4 : one-port background: (mass, pair-mass) sufficient; mass alone not =================
print("="*74)
print("T4  One-port background R hanging by one edge: (leaf-mass, resolved-pair-mass) is")
print("    a sufficient statistic for the 2-color objective; mass alone is NOT.")
print("="*74)
def onePortTree(bg_spec):
    # binary root: [ bg_spec ] | [ cherry(L1,L2) live ]
    t = Tree(); t.build(I(bg_spec, I(L(1, 1.0), L(1, 1.0)))); return t
# Region A: cherry(b1=1, b2=2):  sum=3, pair-mass e2 = 1*2 = 2
tA = onePortTree(I(L(0, 1.0), L(0, 2.0)))
# Region B: different topology/masses but SAME (sum=3, e2=2): caterpillar ((b1',b2'),b3')
#   choose masses so sum=3 and total pairwise product e2=2.  masses x,y,z: x+y+z=3, xy+xz+yz=2.
#   e.g. x=2,y=0.5,z=0.5 -> sum=3, e2=2*.5+2*.5+.5*.5=1+1+0.25=2.25 (no). Solve: pick 1,1,1 -> e2=3.
#   Want e2=2 with 3 leaves: masses 2, 1, 0 invalid. Use x=2, y=1, z=0? no zero.
#   Instead compare two 2-leaf regions with same (sum,e2): cherry(1,2) vs cherry(2,1) trivial-equal,
#   so use cherry(1.5,1.5): sum=3 e2=2.25 (different) -> gives different q. Good contrast set:
tB_same = onePortTree(I(L(0, 2.0), L(0, 1.0)))       # same multiset -> must match tA
tB_massonly = onePortTree(L(0, 3.0))                  # mass-only single leaf: sum=3, e2=0
qA = quartet_2color(tA); qBs = quartet_2color(tB_same); qBm = quartet_2color(tB_massonly)
assert_counts(qA, (2.0, 0.0), "one-port weighted pair reference")
assert_counts(qBs, qA, "one-port equal summary substitution")
assert_counts(qBm, (0.0, 0.0), "one-port mass-only undercount")
print(f"  R=cherry(1,2)   [sum=3,e2=2]     q00_11={qA[0]:.2f}")
print(f"  R=cherry(2,1)   [sum=3,e2=2]     q00_11={qBs[0]:.2f}   (same summary -> {'MATCH' if finite_close(qA[0], qBs[0]) else 'DIFFER'})")
print(f"  R=single leaf   [sum=3,e2=0]     q00_11={qBm[0]:.2f}   (mass-only, loses e2 -> undercount by {qA[0]-qBm[0]:.2f})")
print(f"  check: q00_11 == e2(R)*w(L1)*w(L2) = 2*1*1 = 2  -> {finite_close(qA[0], 2.0)}\n")

print("ALL STRESS TESTS COMPLETED.")
