"""The natural counters: the counter basis of the counting scheme.

Each counter counts one kind of leaf set of a component, described by how the
set sits relative to the component's boundary (the paper's supplement,
"Counters stored by an HDT component"). On a shape with leaf set S:

  unit                       S is empty (every sort except R)
  leaf i                     one leaf of colour i
  pair i j all               C1/C3: two leaves coloured i, j
  pair i j different|same    G: lca at the super root, or below it
  pair i j F|J|L|R|D         C2: both leaves meet the spine at one vertex,
                             in different (F) or a shared (J) off-spine branch;
                             or at different vertices, i below j (R), j below
                             i (L), equal colours (D)
  triple T|B x a             a rooted triple u|vw seen from the top (T) or the
                             bottom (B); x is the colour that completes the
                             channel a with u on one side and v, w on the other
  quartet a                  four leaves resolved as channel a

Pairs are counted only for colour pairs that form a side of some channel. The
counters are checked, exactly, to span the backward space computed by the
pullback closure, and are ordered as in the paper's counter table.
"""
from __future__ import annotations

from collections import defaultdict, deque
from fractions import Fraction

import shapes as G
from shapes import LEAF, REG, U, D, Z, _has, sel_count

CHANNELS = {2: ('00|11', '01|01'), 4: ('01|23', '02|13', '03|12')}
CLASS_ORDER = {'all': 0, 'same': 1, 'different': 2,
               'J': 3, 'F': 4, 'D': 5, 'L': 6, 'R': 7}


def _leaves(frag):
    return [n for n, t in frag.kind.items() if t[0] == LEAF]


def _port(frag, kind):
    return next((n for n, t in frag.kind.items() if t[0] == kind), None)


def _color(frag, node):
    return frag.kind[node][1]


def _rooted(frag, root):
    parent, depth = {root: None}, {root: 0}
    queue = deque([root])
    while queue:
        x = queue.popleft()
        for y in sorted(frag.adj[x]):
            if y not in parent:
                parent[y], depth[y] = x, depth[x] + 1
                queue.append(y)
    return parent, depth


def _lca(parent, depth, a, b):
    while depth[a] > depth[b]:
        a = parent[a]
    while depth[b] > depth[a]:
        b = parent[b]
    while a != b:
        a, b = parent[a], parent[b]
    return a


def _channel(side_a, side_b, K):
    text = lambda side: ''.join(str(c) for c in sorted(side))
    first, second = sorted((text(side_a), text(side_b)))
    name = f'{first}|{second}'
    return name if name in CHANNELS[K] else None


def _quartet(frag, K):
    leaves = _leaves(frag)
    for x in sorted(frag.adj):
        for y in sorted(frag.adj[x]):
            if x >= y:
                continue
            side = G._side(frag, x, y)
            inside = [n for n in leaves if n in side]
            if len(inside) == 2:
                outside = [n for n in leaves if n not in side]
                channel = _channel([_color(frag, n) for n in inside],
                                   [_color(frag, n) for n in outside], K)
                return ('quartet', channel) if channel else None
    return None


def _triple(frag, root, view, K):
    parent, depth = _rooted(frag, root)
    leaves = _leaves(frag)
    for u in leaves:
        v, w = [n for n in leaves if n != u]
        vw, uv = _lca(parent, depth, v, w), _lca(parent, depth, u, v)
        if vw != uv and depth[vw] > depth[uv]:
            z, pair = _color(frag, u), sorted((_color(frag, v), _color(frag, w)))
            matches = [x for x in range(K) if _channel((x, z), pair, K)]
            if not matches:
                return None
            (x,) = matches
            return ('triple', view, x, _channel((x, z), pair, K))
    return None


def _spine_meet(frag, spine, leaf):
    """The spine vertex where the path from `leaf` meets the spine, and the path."""
    seen, queue = {leaf}, deque([(leaf, (leaf,))])
    while queue:
        x, path = queue.popleft()
        if x in spine:
            return x, path
        for y in sorted(frag.adj[x]):
            if y not in seen:
                seen.add(y)
                queue.append((y, path + (y,)))
    raise AssertionError('leaf is not connected to the spine')


def _pair(frag, sort):
    a, b = _leaves(frag)
    i, j = sorted((_color(frag, a), _color(frag, b)))
    if sort in ('C1', 'C3'):
        return ('pair', i, j, 'all')
    if sort == 'G':
        center = _port(frag, Z)
        parent, depth = _rooted(frag, center)
        return ('pair', i, j,
                'different' if _lca(parent, depth, a, b) == center else 'same')
    top, bottom = _port(frag, U), _port(frag, D)
    parent, _ = _rooted(frag, top)
    path = [bottom]
    while path[-1] != top:
        path.append(parent[path[-1]])
    position = {n: k for k, n in enumerate(reversed(path))}   # top = 0
    meet_a, path_a = _spine_meet(frag, position, a)
    meet_b, path_b = _spine_meet(frag, position, b)
    if meet_a == meet_b:
        shared = set(path_a[:-1]) & set(path_b[:-1])
        return ('pair', i, j, 'J' if shared else 'F')
    lower = a if position[meet_a] > position[meet_b] else b
    upper = b if lower == a else a
    if _color(frag, lower) == _color(frag, upper):
        return ('pair', i, i, 'D')
    return ('pair', i, j, 'R' if _color(frag, lower) == i else 'L')


def _labels(code, sort, K):
    """Natural counters that count the leaf set of shape `code` (at most two)."""
    frag = REG[code]
    size = sel_count(code)
    if sort == 'R':
        return [_quartet(frag, K)] if size == 4 else []
    if size == 0:
        return [('unit',)]
    if size == 1:
        return [('leaf', _color(frag, _leaves(frag)[0]))]
    if size == 2:
        label = _pair(frag, sort)
        side = f'{label[1]}{label[2]}'
        return [label] if any(side in c.split('|') for c in CHANNELS[K]) else []
    if size == 3:
        if sort == 'C1':
            return [_triple(frag, _port(frag, U), 'T', K)]
        if sort == 'C3':
            return [_triple(frag, _port(frag, D), 'B', K)]
        if sort == 'G':
            return [_triple(frag, _port(frag, Z), 'T', K)]
        return [_triple(frag, _port(frag, U), 'T', K),
                _triple(frag, _port(frag, D), 'B', K)]
    return [_quartet(frag, K)]


def _order(label, K):
    kind = label[0]
    if kind == 'unit':
        return (0,)
    if kind == 'leaf':
        return (1, label[1])
    if kind == 'pair':
        return (2, label[1], label[2], CLASS_ORDER[label[3]])
    if kind == 'triple':
        return (3, label[1] != 'T', label[2], CHANNELS[K].index(label[3]))
    return (4, CHANNELS[K].index(label[1]))


class CounterBasis:
    """Counters of one interface: `rows[i]` is counter `labels[i]` as a
    functional on shapes."""

    def __init__(self, labels, rows):
        self.labels = labels
        self.rows = rows

    def rank(self):
        return len(self.rows)


def _same_span(rows, closure):
    reference = G.Basis()
    for row in closure.rows:
        reference.add(dict(row))
    counters = G.Basis()
    independent = all(counters.add(dict(row)) for row in rows)
    return (independent and counters.rank() == reference.rank() and
            all(reference.in_span(row) for row in rows))


def natural_bases(K, states, closure):
    """Counter bases per interface, checked equal to the closure's spans.

    `states` maps interfaces to their shapes and `closure` to the backward
    space of each interface. The C1 interface holds sort C3 (a D port only)
    and then sort C1 (a U port); Closed is sort R.
    """
    bases = {}
    for interface in ('Closed', 'C1', 'C2', 'G'):
        table = defaultdict(dict)
        for code in sorted(states[interface], key=G.state_sort_key):
            if interface == 'C1':
                sort = 'C1' if _has(code, U) else 'C3'
            else:
                sort = {'Closed': 'R', 'C2': 'C2', 'G': 'G'}[interface]
            for label in _labels(code, sort, K):
                if label is not None:
                    key = (sort,) + label if interface == 'C1' else label
                    table[key][code] = Fraction(1)
        if interface == 'C1':
            # Production C1 vectors hold the C3 block first, then C1.
            labels = sorted(table, key=lambda l: (l[0] == 'C1', _order(l[1:], K)))
        else:
            labels = sorted(table, key=lambda l: _order(l, K))
        rows = [table[label] for label in labels]
        if not _same_span(rows, closure[interface]):
            raise AssertionError(
                f'K={K} {interface}: natural counters do not span the '
                'backward space')
        bases[interface] = CounterBasis(labels, rows)
    return bases


def coordinates(functional, basis):
    """Exact coefficients c with sum_i c_i * basis.rows[i] == functional."""
    codes = sorted({c for row in basis.rows for c in row} | set(functional),
                   key=G.state_sort_key)
    n = len(basis.rows)
    matrix = [[basis.rows[i].get(code, Fraction(0)) for i in range(n)] +
              [Fraction(functional.get(code, 0))] for code in codes]
    pivots, rank = [], 0
    for column in range(n):
        pivot = next((r for r in range(rank, len(matrix))
                      if matrix[r][column] != 0), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        lead = matrix[rank][column]
        matrix[rank] = [value / lead for value in matrix[rank]]
        for r in range(len(matrix)):
            if r != rank and matrix[r][column] != 0:
                factor = matrix[r][column]
                matrix[r] = [a - factor * b for a, b in zip(matrix[r], matrix[rank])]
        pivots.append(column)
        rank += 1
    if any(matrix[r][n] != 0 for r in range(rank, len(matrix))):
        raise AssertionError('functional is outside the counter span')
    result = [Fraction(0)] * n
    for r, column in enumerate(pivots):
        result[column] = matrix[r][n]
    return result
