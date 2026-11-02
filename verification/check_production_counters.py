"""The production counting tables are the natural scheme of natwta/.

The production generator derives its tables in the natural-counter basis
(tools/counting/natural_counters.py). This check compares them with the
natural automata in natwta/qraft{2,4}_natural_wta.json:

  1. every production operation with a counterpart in the natural signature
     (toG, sib, up_closed, up_open, close, cap, cat) has exactly the same
     update rule, term for term;
  2. the reversals map each counter to its mirror image (C1 <-> C3, top <->
     bottom views, i-below-j <-> j-below-i);
  3. the C3 compositions are the mirror images of the corresponding C1 and R
     rules.

All coefficients are 1. Run from any directory; it takes about 15 seconds.
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'tools' / 'counting'))
import shapes as G           # noqa: E402
import generate_counters as T      # noqa: E402
import natural_counters as N          # noqa: E402

COUNTERPART = {
    ('C_TO_G', ('C1',), 'G'): 'toG',
    ('IG_TF', ('G',), 'C1'): 'up_closed',
    ('IG_TT', ('G',), 'C2'): 'up_open',
    ('IG_FF', ('G',), 'Closed'): 'close',
    ('GG_TO_G', ('G', 'G'), 'G'): 'sib',
    ('CC_TO_C', ('C2', 'C1'), 'C1'): 'cap',
    ('CC_TO_C', ('C2', 'C2'), 'C2'): 'cat',
}
NATURAL_SORT = {'C1': 'C1', 'G': 'G', 'C2': 'C2', 'Closed': 'T'}


def text(label):
    return tuple(str(part) for part in label)


def mirror(label):
    """The counter a reversal maps `label` to (labels carry their sort)."""
    sort, rest = label[0], label[1:]
    sort = {'C1': 'C3', 'C3': 'C1'}.get(sort, sort)
    if rest[0] == 'triple':
        rest = (rest[0], 'B' if rest[1] == 'T' else 'T') + rest[2:]
    elif rest[0] == 'pair' and rest[3] in ('L', 'R'):
        rest = rest[:3] + ('R' if rest[3] == 'L' else 'L',)
    return (sort,) + rest


def main():
    failures = 0
    for K in (2, 4):
        artifact = T.derive_full_grammar(K, verbose=False)
        G.REG.clear()
        states, unary, gg, cc = G.enumerate_full(K)
        _, closure = G.quotient(K, ('Closed', 'C1', 'C2', 'G'), states['Closed'],
                                [t for (t,) in unary.values()], [gg, cc])
        bases = N.natural_bases(K, states, closure)

        def label(interface, index):
            raw = bases[interface].labels[index]
            if interface == 'C1':
                return raw
            return ({'Closed': 'R'}.get(interface, interface),) + raw

        natural = json.load(open(HERE / 'natwta' / f'qraft{K}_natural_wta.json'))
        natural_labels = {s: [text(l) for l in v]
                          for s, v in natural['coordinates'].items()}
        ops = [(op['name'], (op['in'],), op['out'], op['terms'])
               for op in artifact['unary']]
        ops += [(op['name'], (op['left'], op['right']), op['out'], op['terms'])
                for op in artifact['binary']]
        rules = {}
        for name, inputs, output, terms in ops:
            rule = set()
            for *indices, value in terms:
                if value != 1:
                    raise AssertionError(f'K={K} {name}: coefficient {value}')
                out, args = indices[0], indices[1:]
                rule.add((label(output, out),) +
                         tuple(label(i, a) for i, a in zip(inputs, args)))
            rules[(name, inputs, output)] = rule

        def report(what, ok):
            nonlocal failures
            failures += not ok
            print(f'K={K} {what}: {"ok" if ok else "MISMATCH"}')

        # 1. Rules shared with the natural signature.
        for key, symbol in COUNTERPART.items():
            name, inputs, output = key
            mine = {tuple(text(l[1:]) for l in term) for term in rules[key]}
            sorts = [NATURAL_SORT[s] for s in (output,) + inputs]
            theirs = set()
            for target, children, coefficient in natural['transitions'][symbol]:
                assert int(coefficient) == 1
                theirs.add((natural_labels[sorts[0]][target],) +
                           tuple(natural_labels[s][c] for s, c in zip(sorts[1:], children)))
            report(f'{name}{inputs}->{output} = natural {symbol} '
                   f'({len(mine)} terms)', mine == theirs)

        # 2. Reversals are the mirror permutations.
        for interface in ('C1', 'C2'):
            rule = rules[('REV', (interface,), interface)]
            expected = {(mirror(label(interface, i)), label(interface, i))
                        for i in range(len(bases[interface].labels))}
            report(f'REV {interface} mirrors every counter', rule == expected)

        # 3. C3 compositions mirror the C1 ones; C3.C1 -> R is self-mirrored.
        cap = rules[('CC_TO_C', ('C2', 'C1'), 'C1')]
        c3 = rules[('CC_TO_C', ('C1', 'C2'), 'C1')]
        report('C3.C2 -> C3 mirrors cap',
               c3 == {(mirror(o), mirror(b), mirror(a)) for o, a, b in cap})
        root = rules[('CC_TO_C', ('C1', 'C1'), 'Closed')]
        report('C3.C1 -> R is its own mirror',
               root == {(o, mirror(b), mirror(a)) for o, a, b in root})
    if failures:
        sys.exit(1)
    print('PRODUCTION COUNTERS ARE THE NATURAL SCHEME')


if __name__ == '__main__':
    main()
