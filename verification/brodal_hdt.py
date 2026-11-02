"""Brodal et al., SODA'13, Sec. 4: construction of the hierarchical decomposition tree (HDT),
implemented literally.

  before round 1: every leaf is an L component, every internal node an I component; L -> C.
  each round:     C -> G for every downward closed C that is a child of an I;
                  an I is forking if it has >= 2 downward open children (C or I);
                  (i)   I with g >= 2 G children: floor(g/2) pairs, GG -> G  (pairing is a free choice)
                  (ii)  non-forking I with exactly one G child: IG -> C
                  (iii) maximal path of c consecutive C components: floor(c/2) consecutive pairs,
                        topmost with second, ..., CC -> C
                  all compositions of a round act on the components present at its start.
  stop when one component is left.

Every HDT node records its node set, the refined sorts it has (L, I, C1, C2, C3, G, R, in the order
in which the transformations give them) and its boundary vertices for each sort.
`choices` counts the GG steps with g >= 3 (the only freedom in the construction)."""

class Tree:
    def __init__(self, children):
        """children: dict node -> list of children; node 0 is the root"""
        self.ch = {v: list(c) for v, c in children.items()}
        self.par = {0: None}
        for v, cs in self.ch.items():
            for c in cs:
                self.par[c] = v
        self.nodes = sorted(self.par)
        for v in self.nodes:
            self.ch.setdefault(v, [])
        self.depth = {0: 0}
        stack = [0]
        while stack:
            v = stack.pop()
            for c in self.ch[v]:
                self.depth[c] = self.depth[v] + 1; stack.append(c)
        self.leaves = [v for v in self.nodes if not self.ch[v]]
        self._lca = {}

    def lca(self, a, b):
        key = (a, b) if a <= b else (b, a)
        if key in self._lca: return self._lca[key]
        x, y = a, b
        while self.depth[x] > self.depth[y]: x = self.par[x]
        while self.depth[y] > self.depth[x]: y = self.par[y]
        while x != y: x, y = self.par[x], self.par[y]
        self._lca[key] = x
        return x

    def dist(self, a, b):
        return self.depth[a] + self.depth[b] - 2 * self.depth[self.lca(a, b)]


class Comp:
    __slots__ = ('typ', 'nodes', 'parent', 'children', 'hdt')
    def __init__(self, typ, nodes, hdt):
        self.typ = typ; self.nodes = nodes; self.parent = None; self.children = []; self.hdt = hdt


class HDTNode:
    def __init__(self, nodes, kids=()):
        self.nodes = frozenset(nodes); self.kids = list(kids); self.sorts = []   # list of (sort, marks, round)


def build_hdt(T):
    rho = 0
    hdt = []
    comps = {}
    def new_hdt(nodes, kids=()):
        h = HDTNode(nodes, kids); hdt.append(h); return len(hdt) - 1
    for v in T.nodes:
        typ = 'L' if not T.ch[v] else 'I'
        comps[v] = Comp(typ, frozenset([v]), new_hdt([v]))
    for v in T.nodes:
        for c in T.ch[v]:
            comps[c].parent = comps[v]; comps[v].children.append(comps[c])
    alive = set(comps.values())

    def root_node(c):
        return min(c.nodes, key=lambda v: T.depth[v])

    def refined(c):
        if c.typ in ('L', 'I', 'G'): return c.typ
        up = rho not in c.nodes
        down = len(c.children) > 0
        return {(True, False): 'C1', (True, True): 'C2', (False, True): 'C3', (False, False): 'R'}[(up, down)]

    def marks(c):
        s = refined(c)
        if s == 'I': return {'v': next(iter(c.nodes))}
        if s == 'G': return {'T': next(iter(c.parent.nodes))}          # super root = node of the parent I
        m = {}
        if s in ('L', 'C1', 'C2'): m['T'] = T.par[root_node(c)]
        if s in ('C2', 'C3'): m['B'] = T.par[root_node(c.children[0])]
        return m

    def record(c, rnd):
        hdt[c.hdt].sorts.append((refined(c), marks(c), rnd))

    for c in alive: record(c, 0)
    for c in list(alive):                       # L -> C before round 1
        if c.typ == 'L':
            c.typ = 'C'; record(c, 0)
    rnd = 0; choices = 0
    while len(alive) > 1:
        rnd += 1
        # C -> G
        for c in list(alive):
            if c.typ == 'C' and not c.children and c.parent is not None and c.parent.typ == 'I':
                c.typ = 'G'; record(c, rnd)
        ops = []
        used = set()
        # (i) and (ii)
        for c in alive:
            if c.typ != 'I': continue
            gs = [d for d in c.children if d.typ == 'G']
            opens = [d for d in c.children if d.typ in ('C', 'I')]
            if len(gs) >= 2:
                if len(gs) >= 3: choices += 1
                for j in range(0, len(gs) - 1, 2):
                    ops.append(('GG', gs[j], gs[j + 1]))
            elif len(gs) == 1 and len(opens) <= 1:
                ops.append(('IG', c, gs[0]))
        # (iii) maximal paths of consecutive C components
        for c in alive:
            if c.typ != 'C' or (c.parent is not None and c.parent.typ == 'C'): continue
            path = [c]
            while len(path[-1].children) == 1 and path[-1].children[0].typ == 'C':
                path.append(path[-1].children[0])
            for j in range(0, len(path) - 1, 2):
                ops.append(('CC', path[j], path[j + 1]))
        for op in ops:
            for d in op[1:]:
                assert d not in used, 'overlapping compositions'
                used.add(d)
        if not ops:
            raise RuntimeError('no composition possible')
        for kind, a, b in ops:
            h = new_hdt(a.nodes | b.nodes, (a.hdt, b.hdt))
            if kind == 'GG':
                n = Comp('G', a.nodes | b.nodes, h)
                n.parent = a.parent
                a.parent.children = [d for d in a.parent.children if d is not a and d is not b] + [n]
            elif kind == 'IG':
                n = Comp('C', a.nodes | b.nodes, h)
                n.parent = a.parent
                n.children = [d for d in a.children if d is not b]
                if a.parent is not None:
                    a.parent.children = [n if d is a else d for d in a.parent.children]
            else:  # CC
                n = Comp('C', a.nodes | b.nodes, h)
                n.parent = a.parent
                n.children = list(b.children)
                if a.parent is not None:
                    a.parent.children = [n if d is a else d for d in a.parent.children]
            for d in n.children:
                d.parent = n
            alive.discard(a); alive.discard(b); alive.add(n)
            record(n, rnd)
    (last,) = alive
    return hdt, last.hdt, rnd, choices


def check_invariants(T, hdt, root):
    """sanity: node sets of the two children partition the parent's; the root is everything"""
    assert hdt[root].nodes == frozenset(T.nodes)
    for h in hdt:
        if h.kids:
            a, b = (hdt[k].nodes for k in h.kids)
            assert not (a & b) and (a | b) == h.nodes
