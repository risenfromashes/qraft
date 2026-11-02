"""Lower bound of the note on Brodal's own HDTs.

The tree T* (Fig. 5 of the note; node 0 is the root rho) has a unique HDT: the construction of
[SODA'13, Sec. 4] never meets an I component with three or more G children, so the pairing in (i)
is forced.  For the designated nodes x of this HDT the rank of the matrix N_x of Lemma 4.1
(lower_bound.py) equals the minimal number of counters m_X, with the same leaf-count profile, for
k = 2 and k = 4."""
import sys, time
from brodal_hdt import Tree, build_hdt, check_invariants
from lower_bound import rank_at, TARGET

# T*: rho = 0 has the leaves 10, 11 and the child u = 1;  u has the children 3 (cherry 12,13),
# 7 (leaf 19 and cherry 9 = {22,23}) and the spine 2 - 4 - 5 - 6:  2 has the cherry 8 = {20,21},
# 4 the leaf 14, 5 the leaves 15, 16, and 6 is the cherry {17,18}.
TSTAR = {0: [1, 10, 11], 1: [2, 3, 7], 2: [4, 8], 3: [12, 13], 4: [5, 14], 5: [6, 15, 16],
         6: [17, 18], 7: [9, 19], 8: [20, 21], 9: [22, 23]}
# designated components: node sets (the HDT node is found by its node set) and expected boundary
DESIGNATED = {   # label: (sort in brodal_hdt.py, node set, expected boundary)
    'L':    ('L',  {14}, {'T': 4}),
    'I':    ('I',  {1}, {'v': 1}),          # an internal node other than the root
    'Irho': ('I',  {0}, {'v': 0}),          # the root
    'G':    ('G',  {3, 12, 13, 7, 9, 19, 22, 23}, {'T': 1}),
    'C2':   ('C2', {2, 8, 20, 21, 4, 14, 5, 15, 16}, {'T': 1, 'B': 5}),
    'C1':   ('C1', {2, 8, 20, 21, 4, 14, 5, 15, 16, 6, 17, 18}, {'T': 1}),
    'C3':   ('C3', {0, 10, 11, 1, 3, 12, 13, 7, 9, 19, 22, 23}, {'B': 1}),
    'R':    ('R',  set(range(24)), {}),
}

def main():
    T = Tree(TSTAR)
    hdt, root, rounds, choices = build_hdt(T)
    check_invariants(T, hdt, root)
    assert choices == 0, 'the HDT of T* is not unique'
    print(f'T*: {len(T.leaves)} leaves, {len(T.nodes)} nodes; HDT: {len(hdt)} nodes, {rounds} rounds, '
          f'no pairing choice (unique HDT)')
    for k in (2, 4):
        t0 = time.time()
        for label, (sort, nodes, marks) in DESIGNATED.items():
            (h,) = [h for h in hdt if h.nodes == frozenset(nodes)]
            (st,) = [i for i, (s, m, r) in enumerate(h.sorts) if s == sort]
            assert h.sorts[st][1] == marks, (sort, h.sorts[st])
            prof = rank_at(T, h, st, k)
            assert sum(prof) == TARGET[k][label], (k, label, prof)
            print(f'  k={k} {label:4s}: rank N_x = {sum(prof):2d}, profile {prof}  = m_{label}')
        print(f'  ({time.time()-t0:.1f}s)')

if __name__ == '__main__':
    main()
