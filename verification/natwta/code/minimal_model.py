"""Construct and evaluate the minimal automata from rational rank certificates.

The exported JSON gives every minimal counter as a rational combination of
fragment counts, all sparse unary/bilinear transitions, leaf initializers and
readouts. No counter names or tensor entries are guessed.
"""
from fractions import Fraction
import itertools
import json
from pathlib import Path
import random

import clusters as C
from certify import decode_basis
from soda_run import channels_from_certificate


class MinimalModel:
    def __init__(self, certificate):
        self.record = json.loads(Path(certificate).read_text())
        self.K = self.record['K']
        name = self.record['task']
        self.out = self.record['output_sort']
        self.channels, _ = channels_from_certificate(self.record)
        self.rows, self.pivots, self.by_code, self.counter_descriptions = {}, {}, {}, {}
        for sort, blocks in self.record['sorts'].items():
            self.rows[sort] = []
            self.pivots[sort] = []
            self.by_code[sort] = {}
            self.counter_descriptions[sort] = []
            for block in blocks:
                shapes = C.fragments(sort, block['degree'], self.K)
                for pivot, sparse in decode_basis(block):
                    i = len(self.rows[sort])
                    row = {shapes[j]: v for j, v in sparse.items()}
                    self.rows[sort].append(row)
                    self.pivots[sort].append(shapes[pivot])
                    for code, value in row.items():
                        self.by_code[sort].setdefault(code, []).append((i, value))
                    self.counter_descriptions[sort].append({
                        'degree': block['degree'], 'pivot_fragment': C.pretty(shapes[pivot]),
                        'terms': [[str(v), C.pretty(code)] for code, v in row.items()]})
        self.tensors = {}
        for op, (arity, ins, out, _) in C.OPS.items():
            entries = []
            for indices in itertools.product(*(range(len(self.pivots[s])) for s in ins)):
                codes = tuple(self.pivots[s][i] for s, i in zip(ins, indices))
                code = C.ghat(op, codes)
                for j, value in self.by_code[out].get(code, []):
                    entries.append((j, indices, value))
            self.tensors[op] = entries
        self.readouts = {}
        for name, channel in self.channels.items():
            def pattern(code):
                if self.out == 'G':
                    ((_, children),) = code
                    return (((4,), children),)
                return code
            self.readouts[name] = [channel.get(pattern(code), 0) for code in self.pivots[self.out]]

    def project(self, sort, state):
        q = [Fraction(0) for _ in self.rows[sort]]
        for code, amount in state.items():
            for j, v in self.by_code[sort].get(code, []):
                q[j] += amount * v
        return q

    def apply(self, op, vectors):
        _, ins, out, _ = C.OPS[op]
        q = [Fraction(0) for _ in self.rows[out]]
        for j, indices, value in self.tensors[op]:
            v = value
            for vector, i in zip(vectors, indices):
                v *= vector[i]
            q[j] += v
        return q

    def evaluate(self, term):
        op = term[0]
        if op == 'leaf':
            return self.project('C1', C.atom_state('leaf', term[1], term[2]))
        return self.apply(op, [self.evaluate(t) for t in term[1:]])

    def export(self):
        return {'task': self.record['task'], 'variant': self.record.get('variant'),
                'K': self.K, 'field': 'Q',
                'output_sort': self.out,
                'indexing': 'zero based; each tensor entry is [output, input_indices, coefficient]',
                'counter_basis': self.counter_descriptions,
                'transitions': {op: [[j, list(ii), str(v)] for j, ii, v in entries]
                                for op, entries in self.tensors.items()},
                'leaf_initializers': {str(c): {
                    'constant': list(map(str, self.project('C1', C.atom_state('leaf', c, 0)))),
                    'weight_coefficient': [str(b - a) for a, b in zip(
                        self.project('C1', C.atom_state('leaf', c, 0)),
                        self.project('C1', C.atom_state('leaf', c, 1)))]}
                    for c in range(self.K)},
                'readouts': self.readouts}


def main():
    audit = Path(__file__).resolve().parents[1] / 'audit'
    out = audit / 'minimal_models'
    out.mkdir(exist_ok=True)
    names = ['QRAFT2_K2', 'QRAFT4_K4', 'Ar_x3_K4', 'Br_x3_K4', 'Ar_Br_x3_K4']
    result = {}
    for name in names:
        model = MinimalModel(audit / 'certificates' / (name + '.json'))
        data = model.export()
        # Check the exported initializer, including its constant-one coordinate.
        for c in range(model.K):
            init = data['leaf_initializers'][str(c)]
            for w in (Fraction(0), Fraction(-1), Fraction(2, 3)):
                from_export = [Fraction(a) + w * Fraction(b) for a, b in
                               zip(init['constant'], init['weight_coefficient'])]
                assert from_export == model.project('C1', C.atom_state('leaf', c, w))
        rng = random.Random(127)
        for _ in range(100):
            sort = rng.choice(['T', 'C1', 'G', 'C2'])
            term, A, state = C.random_term(sort, model.K, rng,
                                          weights=lambda: rng.choice([Fraction(0), Fraction(-1), Fraction(2, 3), Fraction(2)]))
            q = model.evaluate(term)
            assert q == model.project(sort, state)
        (out / (name + '.json')).write_text(json.dumps(data, separators=(',', ':')) + '\n')
        result[name] = {'dimensions': {s: len(rows) for s, rows in model.rows.items()},
                        'random_terms': 100, 'seed': 127,
                        'tensor_nonzeros': {op: len(t) for op, t in model.tensors.items()}}
        print(name, result[name], flush=True)
    (audit / 'minimal_model_validation.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
