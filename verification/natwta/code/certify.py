"""Reproducible rational rank certificates for the completion experiment.

Only the Python standard library is needed. A certificate gives a rational
row-space upper bound and a matching finite-field lower bound. Thus each
reported finite-palette rank is proved over Q, not inferred from prime agreement.
Polynomial fits are deliberately not certified by this program.
"""
import argparse
from fractions import Fraction
from functools import reduce
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import time

import clusters as C
from soda_run import channels_for


def primitive(row):
    if not row:
        return row
    g = reduce(gcd, map(abs, row.values()))
    if row[min(row)] < 0:
        g = -g
    return {k: v // g for k, v in row.items()}


def sparse_rref(rows):
    """Integer elimination with gcd reduction, then rational back substitution."""
    echelon = {}
    seen = set()
    for dense in rows:
        key = tuple(dense)
        if key in seen:
            continue
        seen.add(key)
        row = {j: int(v) for j, v in enumerate(dense) if v}
        while row:
            p = min(row)
            if p not in echelon:
                echelon[p] = primitive(row)
                break
            b = echelon[p]
            common = gcd(abs(row[p]), abs(b[p]))
            a, c = b[p] // common, row[p] // common
            row = {j: a * row.get(j, 0) - c * b.get(j, 0)
                   for j in row.keys() | b.keys()}
            row = primitive({j: v for j, v in row.items() if v})
    reduced = {}
    for p in sorted(echelon, reverse=True):
        src = echelon[p]
        row = {j: Fraction(v, src[p]) for j, v in src.items()}
        for q in sorted(reduced):
            if q not in row:
                continue
            a = row[q]
            for j, v in reduced[q].items():
                row[j] = row.get(j, 0) - a * v
                if row[j] == 0:
                    del row[j]
        reduced[p] = row
    return [(p, reduced[p]) for p in sorted(reduced)]


def check_upper_bound(rows, basis):
    """Check each original row is exactly reconstructed using pivot columns."""
    for p, row in basis:
        assert row[p] == 1
        assert all(row.get(q, 0) == (p == q) for q, _ in basis)
    for src in rows:
        reconstructed = {}
        for p, row in basis:
            a = src[p]
            if a:
                for j, v in row.items():
                    reconstructed[j] = reconstructed.get(j, 0) + a * v
        assert all(reconstructed.get(j, 0) == a for j, a in enumerate(src))


def modular_rank_sparse(rows, p):
    echelon = {}
    for src in rows:
        row = {j: a % p for j, a in enumerate(src) if a % p}
        while row:
            q = min(row)
            if q not in echelon:
                inv = pow(row[q], -1, p)
                echelon[q] = {j: a * inv % p for j, a in row.items()}
                break
            a = row[q]
            for j, v in echelon[q].items():
                row[j] = (row.get(j, 0) - a * v) % p
                if row[j] == 0:
                    del row[j]
    return len(echelon)


def certificate(rows, ncols):
    basis = sparse_rref(rows)
    check_upper_bound(rows, basis)
    ranks = {str(p): modular_rank_sparse(rows, p) for p in (C.P1, C.P2)}
    assert max(ranks.values(), default=0) == len(basis)
    raw = json.dumps(rows, separators=(',', ':')).encode()
    return {'field': 'Q', 'rows': len(rows), 'columns': ncols,
            'rank': len(basis), 'modular_lower_bounds': ranks,
            'matrix_sha256': sha256(raw).hexdigest(),
            'basis': [{'pivot': p, 'row': [[j, str(v)] for j, v in sorted(row.items())]}
                      for p, row in basis]}


def decode_basis(cert):
    return [(b['pivot'], {j: Fraction(v) for j, v in b['row']}) for b in cert['basis']]


def evaluate_case(name, K, channels, out, dest):
    start = time.time()
    record = {'task': name, 'K': K, 'output_sort': out, 'sorts': {}}
    for sort in ('T', 'C1', 'G', 'C2'):
        blocks = []
        for m in range(C.D + 1):
            frags, rows, _ = C.completion_matrix(sort, m, K, channels, out)
            cert = certificate(rows, len(frags))
            cert['degree'] = m
            blocks.append(cert)
        record['sorts'][sort] = blocks
    (dest / f'{name}_K{K}.json').write_text(json.dumps(record, separators=(',', ':')) + '\n')
    profile = {s: [b['rank'] for b in blocks] for s, blocks in record['sorts'].items()}
    print(name, K, profile, 'seconds', round(time.time() - start, 2), flush=True)
    return {'K': K, 'output_sort': out, 'profiles': profile,
            'dimensions': {s: sum(p) for s, p in profile.items()},
            'certificate': f'certificates/{name}_K{K}.json'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--families', default='A,B,Ar,Br,Ar+Br')
    ap.add_argument('--xs', default='1,2,3,4,5,6,7')
    ap.add_argument('--skip-qraft', action='store_true')
    args = ap.parse_args()
    audit = Path(__file__).resolve().parents[1] / 'audit'
    dest = audit / 'certificates'
    dest.mkdir(parents=True, exist_ok=True)
    result = {'status': 'exact finite-palette ranks over Q', 'cases': {}}
    if not args.skip_qraft:
        for K in (2, 4):
            name = f'QRAFT{K}'
            result['cases'][name] = evaluate_case(name, K, C.qraft_channels(K), 'T', dest)
    for x in map(int, args.xs.split(',')):
        for fam in args.families.split(','):
            channels, out = channels_for(fam, x)
            name = fam.replace('+', '_') + f'_x{x}'
            result['cases'][name] = evaluate_case(name, x + 1, channels, out, dest)
            (audit / 'exact_results.json').write_text(json.dumps(result, indent=2) + '\n')
    (audit / 'exact_results.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
