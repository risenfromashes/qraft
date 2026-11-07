#!/usr/bin/env python3
"""Derive QRAFT's counter tables and write src/counting/generated/.

For K = 2 and K = 4 colours this enumerates the canonical shapes of every
component sort (shapes.py), computes the counters as an exact rational
quotient, re-expresses them in the natural counter basis
(natural_counters.py), and derives every transition (C_TO_G, IG, CLOSE,
REVERSE, GG, CC) as a sparse table with unit coefficients. Each table is
checked on every legal pair of canonical states, and the resulting counting
is checked against brute-force quartet enumeration.

For a pair of operand shapes a, b the composition satisfies
    sig_out(op(a, b)) = T_op(sig(a), sig(b)),
so with representatives whose signatures are independent,
    T_op[k, p, q] = sum_{i,j} Uinv_a[i, p] Uinv_b[j, q] sig_out(op(rep_i, rep_j))[k].

Run from anywhere:
  generate_counters.py --full                write the headers
  generate_counters.py --check-full          compare the checked-in headers
  generate_counters.py --verify-full         derive and check without writing
  generate_counters.py --check-reproducible  compare two Python hash seeds
"""
from __future__ import annotations
import argparse, hashlib, itertools, json, os, random, subprocess, sys
from collections import defaultdict
from fractions import Fraction

import shapes as G
import natural_counters
from shapes import (Basis, interface_of, sel_count, readout_rows, REG,
                             enumerate_cf, enumerate_full, quotient, reduce_and_code, Frag, U, LEAF,
                             op_c_to_g, op_gg_to_g, op_ig_to_c, _channel,
                             _four_leaf_split, brute_force, random_tree, caterpillar,
                             state_sort_key, verify_state_completeness, _has,
                             Z, I, D)

Q0 = Fraction(0)


def mat_inverse(cols):
    r = len(cols)
    A = [[cols[c][row] for c in range(r)] + [Fraction(1 if row == c else 0) for c in range(r)]
         for row in range(r)]
    for col in range(r):
        piv = next((rr for rr in range(col, r) if A[rr][col] != 0), None)
        if piv is None:
            raise ValueError("singular representative matrix")
        A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        A[col] = [x / pv for x in A[col]]
        for rr in range(r):
            if rr != col and A[rr][col] != 0:
                f = A[rr][col]
                A[rr] = [a - f * b for a, b in zip(A[rr], A[col])]
    return [[A[row][r + c] for row in range(r)] for c in range(r)]  # inverse columns


def sig_of(state, basis):
    return [row.get(state, Q0) for row in basis.rows]


# --- partial-domain tensor construction (for the full C2 grammar) ------------
# An op may be defined only on a SUBSET of an interface's states (CC_TO_C's higher
# operand needs a D-port; C_TO_G needs a U-port). Build a tensor correct on the
# op's actual domain by using domain-restricted representatives whose sig-vectors
# are independent, then extending to a full basis with standard vectors on which
# the tensor acts as 0 (those coords never occur at runtime).

def _full_basis(reps_sig, raw_cols, dim, rout):
    """reps_sig: list of dim-vectors (independent so far); raw_cols: matching
    output-sig columns (len rout). Extend with standard e_i (raw 0) to dim cols."""
    B = Basis()
    cols, raws = [], []
    for s, r in zip(reps_sig, raw_cols):
        if B.add({i: v for i, v in enumerate(s) if v != 0}):
            cols.append(s); raws.append(r)
    for i in range(dim):
        if len(cols) == dim:
            break
        if B.add({i: Fraction(1)}):
            e = [Q0] * dim; e[i] = Fraction(1)
            cols.append(e); raws.append([Q0] * rout)
    assert len(cols) == dim, f"basis extension failed ({len(cols)}/{dim})"
    return cols, raws


def build_unary_tensor_partial(sub, in_basis, out_basis):
    rin, rout = len(in_basis.rows), len(out_basis.rows)
    ordered = sorted(sub.items(), key=lambda item: state_sort_key(item[0]))
    reps_sig = [sig_of(x, in_basis) for x, _ in ordered]
    raw_cols = [sig_of(ox, out_basis) for _, ox in ordered]
    cols, raws = _full_basis(reps_sig, raw_cols, rin, rout)
    Uinv = mat_inverse(cols)               # Uinv[col][row] = INV[row][col]
    # M[k][p] = sum_i raws[i][k] * INV[i][p], INV[i][p] = Uinv[p][i]
    M = [[sum(raws[i][k] * Uinv[p][i] for i in range(rin)) for p in range(rin)]
         for k in range(rout)]
    return M


def build_binary_tensor_partial(sub, a_basis, b_basis, out_basis):
    ra, rb, rout = len(a_basis.rows), len(b_basis.rows), len(out_basis.rows)
    dom_a = sorted({a for (a, _b) in sub}, key=state_sort_key)
    dom_b = sorted({b for (_a, b) in sub}, key=state_sort_key)
    # independent reps within each operand domain, extended to full bases
    Ba = Basis(); reps_a = []
    for a in dom_a:
        if Ba.add({i: v for i, v in enumerate(sig_of(a, a_basis)) if v != 0}):
            reps_a.append(a)
    Bb = Basis(); reps_b = []
    for b in dom_b:
        if Bb.add({i: v for i, v in enumerate(sig_of(b, b_basis)) if v != 0}):
            reps_b.append(b)
    cols_a, _ = _full_basis([sig_of(a, a_basis) for a in reps_a],
                            [[Q0]] * len(reps_a), ra, 1)
    cols_b, _ = _full_basis([sig_of(b, b_basis) for b in reps_b],
                            [[Q0]] * len(reps_b), rb, 1)
    na, nb = len(reps_a), len(reps_b)
    Ua_inv = mat_inverse(cols_a)
    Ub_inv = mat_inverse(cols_b)
    # Sparse raw: only real rep pairs present in sub are nonzero.
    # raw_nz[k] = list of (i, j, value) with value = sig_out[k].
    nz = []  # (i, j, out_sig_vector)
    for i in range(na):
        for j in range(nb):
            oc = sub.get((reps_a[i], reps_b[j]))
            if oc is not None:
                nz.append((i, j, sig_of(oc, out_basis)))
    # Two-step contraction per output coord k:
    #   T[k] = INV_a^T · raw_k · INV_b,  INV_a[i][p]=Ua_inv[p][i], INV_b[j][q]=Ub_inv[q][j]
    T = [[[Q0] * rb for _ in range(ra)] for _ in range(rout)]
    for k in range(rout):
        # M1[i][q] = sum_j raw_k[i][j] * INV_b[j][q]   (sparse over nz)
        M1 = {}  # i -> [rb] row
        for (i, j, so) in nz:
            v = so[k]
            if v == 0:
                continue
            row = M1.get(i)
            if row is None:
                row = [Q0] * rb
                M1[i] = row
            ubj = Ub_inv  # Ub_inv[q][j]
            for q in range(rb):
                uq = ubj[q][j]
                if uq != 0:
                    row[q] += v * uq
        # T[k][p][q] = sum_i INV_a[i][p] * M1[i][q]
        Tk = T[k]
        for i, row in M1.items():
            uai = [Ua_inv[p][i] for p in range(ra)]  # INV_a[i][p]
            for p in range(ra):
                up = uai[p]
                if up == 0:
                    continue
                Tp = Tk[p]
                for q in range(rb):
                    rq = row[q]
                    if rq != 0:
                        Tp[q] += up * rq
    return T


def sparse_unary(M):
    return [(k, p, value) for k, row in enumerate(M)
            for p, value in enumerate(row) if value != 0]


def sparse_binary(T):
    return [(k, p, q, value) for k, plane in enumerate(T)
            for p, row in enumerate(plane)
            for q, value in enumerate(row) if value != 0]


def apply_sparse_unary(terms, rout, qin):
    result = [Q0] * rout
    for k, p, value in terms:
        if qin[p] != 0:
            result[k] += value * qin[p]
    return result


def apply_sparse_binary(terms, rout, qa, qb):
    result = [Q0] * rout
    for k, p, q, value in terms:
        if qa[p] != 0 and qb[q] != 0:
            result[k] += value * qa[p] * qb[q]
    return result


FULL_SPEC = """qraft-counters/v1
states=reduce(minimal-subtree(selected-leaves+typed-ports),suppress-nonport-degree2)
degree=selected-leaf-count<=4
interfaces=Closed,C1(U-or-D),C2(U,D),G(Z-protected)
unary=C_TO_G,IG_FF,IG_TF,IG_TT,CLOSE,REV_C1,REV_C2
binary=GG_TO_G,CC_TO_C
observables=K2:Q00_11,Q01_01;K4:Q01_23,Q02_13,Q03_12
quotient=natural-counters-spanning-exact-pullback-closure
"""


def _functional_coordinates(functional, basis):
    return natural_counters.coordinates(functional, basis)


def _coordinate_degrees(basis):
    degrees = []
    for row in basis.rows:
        row_degrees = {sel_count(code) for code, value in row.items() if value != 0}
        if len(row_degrees) != 1:
            raise AssertionError(f"quotient coordinate mixes degrees: {row_degrees}")
        degrees.append(next(iter(row_degrees)))
    return degrees


def _code_text(code):
    token = code[0]
    kind = token[0]
    if kind == LEAF:
        label = f"L{token[1]}"
    elif kind == 'EMPTY':
        label = 'E'
    else:
        label = kind
    return label + '(' + ','.join(_code_text(child) for child in code[1:]) + ')'


def _sha_lines(lines):
    payload = ''.join(line + '\n' for line in lines).encode('ascii')
    return hashlib.sha256(payload).hexdigest()


def _projection_metadata(states, basis):
    index = {code: i for i, code in enumerate(states)}
    lines = []
    terms = 0
    for row_index, row in enumerate(basis.rows):
        for code, value in sorted(row.items(), key=lambda item: state_sort_key(item[0])):
            lines.append(f"{row_index}:{index[code]}:{value.numerator}/{value.denominator}")
            terms += 1
    return _sha_lines(lines), terms


def _binary_lookup(terms):
    lookup = defaultdict(list)
    for k, p, q, value in terms:
        lookup[(p, q)].append((k, value))
    return lookup


def _apply_binary_lookup(lookup, rout, qa, qb):
    result = [Q0] * rout
    left = [(p, value) for p, value in enumerate(qa) if value != 0]
    right = [(q, value) for q, value in enumerate(qb) if value != 0]
    for p, av in left:
        for q, bv in right:
            for k, coefficient in lookup.get((p, q), ()):
                result[k] += coefficient * av * bv
    return result


def _assert_int64(values, label):
    limit = (1 << 63) - 1
    for value in values:
        if abs(value.numerator) > limit or value.denominator > limit:
            raise AssertionError(f"{label} coefficient does not fit int64: {value}")


def _assert_signed_units(values, label):
    for value in values:
        if value != 0 and (value.denominator != 1 or abs(value.numerator) != 1):
            raise AssertionError(f"{label} coefficient is not -1/0/+1: {value}")


def _jsonable(value):
    if isinstance(value, Fraction):
        return [value.numerator, value.denominator]
    if isinstance(value, dict):
        return {key: _jsonable(value[key]) for key in sorted(value)}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _artifact_digest(artifact):
    payload = json.dumps(_jsonable(artifact), sort_keys=True,
                         separators=(',', ':')).encode('ascii')
    return hashlib.sha256(payload).hexdigest()


def _full_leaf_vectors(K, basis):
    absent_frag = Frag()
    absent_frag.add(0, (U,))
    absent = G.reg(reduce_and_code(absent_frag))
    present = []
    for color in range(K):
        frag = Frag()
        frag.add(0, (U,))
        frag.add(1, (LEAF, color))
        frag.link(0, 1)
        present.append(sig_of(G.reg(reduce_and_code(frag)), basis))
    return sig_of(absent, basis), present


def _full_constant_vectors(Bq):
    center = Frag()
    center.add(0, (Z,))
    center_code = G.reg(reduce_and_code(center))
    corridor = Frag()
    corridor.add(0, (U,))
    corridor.add(1, (I,))
    corridor.add(2, (D,))
    corridor.link(0, 1)
    corridor.link(1, 2)
    corridor_code = G.reg(reduce_and_code(corridor))
    empty_code = G.reg(reduce_and_code(Frag()))
    return {
        'CENTER_EMPTY_G': sig_of(center_code, Bq['G']),
        'IDENTITY_C2': sig_of(corridor_code, Bq['C2']),
        'CLOSED_EMPTY': sig_of(empty_code, Bq['Closed']),
    }


def derive_full_grammar(K, verbose=True):
    """Derive and exhaustively verify the complete K=2/K=4 quotient algebra."""
    if verbose:
        print(f"\n=== counter tables K={K} ===")
    REG.clear()
    by, unary_tables, gg, cc = enumerate_full(K)
    ifaces = ('Closed', 'C1', 'C2', 'G')
    counts = verify_state_completeness(K, by)
    ranks, closure = quotient(K, ifaces, by['Closed'],
                              [table for (table,) in unary_tables.values()],
                              [gg, cc])
    # Coordinates are the natural counters, checked to span the closure.
    Bq = natural_counters.natural_bases(K, by, closure)
    ordered_states = {
        interface: sorted(by[interface], key=state_sort_key)
        for interface in ifaces
    }
    state_hashes = {
        interface: _sha_lines(_code_text(code) for code in ordered_states[interface])
        for interface in ifaces
    }
    projection = {
        interface: _projection_metadata(ordered_states[interface], Bq[interface])
        for interface in ifaces
    }
    degrees = {interface: _coordinate_degrees(Bq[interface]) for interface in ifaces}

    c1_orientation_rows = defaultdict(int)
    for row in Bq['C1'].rows:
        orientations = {
            U if _has(code, U) else D
            for code, value in row.items() if value != 0
        }
        if len(orientations) != 1:
            raise AssertionError(f"K={K} C1 quotient row mixes orientations")
        c1_orientation_rows[next(iter(orientations))] += 1
    if (c1_orientation_rows[U] != c1_orientation_rows[D] or
            c1_orientation_rows[U] + c1_orientation_rows[D] != ranks['C1']):
        raise AssertionError(f"K={K} C1 quotient does not split by orientation")
    # The packed kernels index C1 as the D block [0, h) then the U block.
    half = ranks['C1'] // 2
    for index, row in enumerate(Bq['C1'].rows):
        expected = D if index < half else U
        if any(_has(code, U) != (expected == U) for code in row):
            raise AssertionError(f"K={K} C1 coordinate {index} is in the wrong block")

    # Every canonical state of degree d has support only in degree-d quotient
    # coordinates. This grading is the certificate that pairs with total degree
    # greater than four are identically zero, without enumerating a million
    # deliberately omitted GxG pairs.
    for interface in ifaces:
        for code in ordered_states[interface]:
            signature = sig_of(code, Bq[interface])
            if any(value != 0 and degrees[interface][i] != sel_count(code)
                   for i, value in enumerate(signature)):
                raise AssertionError(f"signature degree mismatch {interface}")

    unary = []
    unary_checks = 0
    unary_order = ('C_TO_G', 'IG_FF', 'IG_TF', 'IG_TT', 'CLOSE', 'REV')
    for name in unary_order:
        table = unary_tables[name][0]
        buckets = defaultdict(dict)
        for source, target in sorted(table.items(), key=lambda item: state_sort_key(item[0])):
            buckets[(interface_of(source), interface_of(target))][source] = target
        for in_if, out_if in sorted(buckets):
            sub = buckets[(in_if, out_if)]
            if name in ('C_TO_G', 'CLOSE'):
                expected_domain = {
                    code for code in ordered_states['C1'] if _has(code, U)}
            elif name == 'REV':
                expected_domain = set(ordered_states[in_if])
            else:
                expected_domain = set(ordered_states['G'])
            if set(sub) != expected_domain:
                raise AssertionError(
                    f"incomplete unary domain {name} {in_if}->{out_if}: "
                    f"actual={len(sub)} expected={len(expected_domain)}")
            matrix = build_unary_tensor_partial(sub, Bq[in_if], Bq[out_if])
            terms = sparse_unary(matrix)
            _assert_int64((value for _, _, value in terms), name)
            for k, p, _ in terms:
                if degrees[out_if][k] != degrees[in_if][p]:
                    raise AssertionError(f"unary degree mismatch {name}")
            for source, target in sub.items():
                actual = apply_sparse_unary(
                    terms, ranks[out_if], sig_of(source, Bq[in_if]))
                if actual != sig_of(target, Bq[out_if]):
                    raise AssertionError(f"unary {name} {in_if}->{out_if} mismatch")
                unary_checks += 1
            unary.append({'name': name, 'in': in_if, 'out': out_if,
                          'terms': terms, 'domain_states': len(sub)})

    binary = []
    binary_checks = 0
    for name, table in (('GG_TO_G', gg), ('CC_TO_C', cc)):
        buckets = defaultdict(dict)
        for (left, right), target in sorted(
                table.items(),
                key=lambda item: (state_sort_key(item[0][0]),
                                  state_sort_key(item[0][1]))):
            key = (interface_of(left), interface_of(right), interface_of(target))
            buckets[key][(left, right)] = target
        for left_if, right_if, out_if in sorted(buckets):
            sub = buckets[(left_if, right_if, out_if)]
            if name == 'GG_TO_G':
                expected_left = ordered_states['G']
                expected_right = ordered_states['G']
            else:
                expected_left = [
                    code for code in ordered_states[left_if]
                    if left_if == 'C2' or _has(code, D)]
                expected_right = [
                    code for code in ordered_states[right_if]
                    if right_if == 'C2' or _has(code, U)]
            left_by_degree = defaultdict(list)
            right_by_degree = defaultdict(list)
            for code in expected_left:
                left_by_degree[sel_count(code)].append(code)
            for code in expected_right:
                right_by_degree[sel_count(code)].append(code)
            expected_pairs = set()
            for left_degree, left_states in left_by_degree.items():
                for right_degree in range(5 - left_degree):
                    expected_pairs.update(itertools.product(
                        left_states, right_by_degree[right_degree]))
            if set(sub) != expected_pairs:
                raise AssertionError(
                    f"incomplete binary domain {name} "
                    f"{left_if},{right_if}->{out_if}: actual={len(sub)} "
                    f"expected={len(expected_pairs)}")
            tensor = build_binary_tensor_partial(
                sub, Bq[left_if], Bq[right_if], Bq[out_if])
            terms = sparse_binary(tensor)
            _assert_int64((value for _, _, _, value in terms), name)
            for k, p, q, _ in terms:
                if degrees[out_if][k] != degrees[left_if][p] + degrees[right_if][q]:
                    raise AssertionError(f"binary degree mismatch {name}")
                if degrees[out_if][k] > 4:
                    raise AssertionError(f"untruncated tensor term {name}")
            lookup = _binary_lookup(terms)
            for (left, right), target in sub.items():
                actual = _apply_binary_lookup(
                    lookup, ranks[out_if], sig_of(left, Bq[left_if]),
                    sig_of(right, Bq[right_if]))
                if actual != sig_of(target, Bq[out_if]):
                    raise AssertionError(
                        f"binary {name} {left_if},{right_if}->{out_if} mismatch")
                binary_checks += 1
            binary.append({'name': name, 'left': left_if, 'right': right_if,
                           'out': out_if, 'terms': terms,
                           'domain_pairs': len(sub)})

    channels = (['Q00_11', 'Q01_01'] if K == 2 else
                ['Q01_23', 'Q02_13', 'Q03_12'])
    readout = [_functional_coordinates(row, Bq['Closed'])
               for row in readout_rows(K, by['Closed'])]
    leaf_absent, leaf_present = _full_leaf_vectors(K, Bq['C1'])
    constants = _full_constant_vectors(Bq)

    # The exact derivation is rational in general, but these two fixed palettes
    # happen to factor entirely over signed-unit coefficients. Keep that useful
    # runtime property as a checked theorem rather than a documented observation.
    coefficient_groups = []
    for interface in ifaces:
        coefficient_groups.extend(row.values() for row in Bq[interface].rows)
    coefficient_groups.append(
        [value for op in unary for _, _, value in op['terms']])
    coefficient_groups.append(
        [value for op in binary for _, _, _, value in op['terms']])
    coefficient_groups.extend(readout)
    coefficient_groups.append(leaf_absent)
    coefficient_groups.extend(leaf_present)
    coefficient_groups.extend(constants.values())
    _assert_signed_units(
        (value for group in coefficient_groups for value in group),
        f"K={K} full artifact")
    if any(value not in (Q0, Fraction(1))
           for vector in readout for value in vector):
        raise AssertionError(f"K={K} readout coefficient is not 0/+1")

    artifact = {
        'format': 'qraft-counters-v1',
        'K': K,
        'spec_sha256': hashlib.sha256(FULL_SPEC.encode('ascii')).hexdigest(),
        'state_counts': counts,
        'ranks': ranks,
        'coordinate_degrees': degrees,
        'state_manifest_sha256': state_hashes,
        'projection': {
            interface: {'sha256': projection[interface][0],
                        'terms': projection[interface][1]}
            for interface in ifaces
        },
        'unary': unary,
        'binary': binary,
        'channels': channels,
        'readout': readout,
        'leaf_absent': leaf_absent,
        'leaf_present': leaf_present,
        'constants': constants,
    }
    artifact['artifact_sha256'] = _artifact_digest(artifact)
    if verbose:
        print(f"  states: {counts}")
        print(f"  ranks:  {ranks}")
        print(f"  verified: {unary_checks} unary states + "
              f"{binary_checks} legal binary pairs")
        print(f"  sparse terms: unary={sum(len(op['terms']) for op in unary)} "
              f"binary={sum(len(op['terms']) for op in binary)}")
        print(f"  artifact sha256: {artifact['artifact_sha256']}")
    return artifact


def _identifier(op):
    parts = [op['name']]
    for key in ('in', 'left', 'right', 'out'):
        if key in op:
            parts.append(op[key])
    return '_'.join(parts).upper()


def _coordinate_color_masks(artifact):
    """Derive the unique required color set of every quotient coordinate."""
    K = artifact['K']
    mask_count = 1 << K
    interfaces = ('Closed', 'C1', 'C2', 'G')
    exact = [
        {interface: set() for interface in interfaces}
        for _ in range(mask_count)
    ]

    def add_vector(mask, interface, vector):
        exact[mask][interface].update(
            index for index, value in enumerate(vector) if value != 0)

    add_vector(0, 'C1', artifact['leaf_absent'])
    for color, vector in enumerate(artifact['leaf_present']):
        add_vector(1 << color, 'C1', vector)
    add_vector(0, 'G', artifact['constants']['CENTER_EMPTY_G'])
    add_vector(0, 'C2', artifact['constants']['IDENTITY_C2'])
    add_vector(0, 'Closed', artifact['constants']['CLOSED_EMPTY'])

    changed = True
    while changed:
        changed = False
        for mask in range(mask_count):
            for op in artifact['unary']:
                reached = {
                    out for out, source, value in op['terms']
                    if value != 0 and source in exact[mask][op['in']]
                }
                before = len(exact[mask][op['out']])
                exact[mask][op['out']].update(reached)
                changed |= len(exact[mask][op['out']]) != before
        for left_mask in range(mask_count):
            for right_mask in range(mask_count):
                output_mask = left_mask | right_mask
                for op in artifact['binary']:
                    reached = {
                        out for out, left, right, value in op['terms']
                        if value != 0
                        and left in exact[left_mask][op['left']]
                        and right in exact[right_mask][op['right']]
                    }
                    before = len(exact[output_mask][op['out']])
                    exact[output_mask][op['out']].update(reached)
                    changed |= len(exact[output_mask][op['out']]) != before

    coordinate_masks = {}
    for interface in interfaces:
        masks = [-1] * artifact['ranks'][interface]
        for mask in range(mask_count):
            for coordinate in exact[mask][interface]:
                if masks[coordinate] not in (-1, mask):
                    raise AssertionError(
                        f"coordinate color grading overlaps: K={K} "
                        f"{interface}[{coordinate}]")
                masks[coordinate] = mask
        if any(mask < 0 for mask in masks):
            raise AssertionError(
                f"coordinate color grading incomplete: K={K} {interface}")
        coordinate_masks[interface] = masks

    for op in artifact['unary']:
        for out, source, value in op['terms']:
            if value != 0 and (
                    coordinate_masks[op['out']][out] !=
                    coordinate_masks[op['in']][source]):
                raise AssertionError(
                    f"unary color grading mismatch: {op['name']}")
    for op in artifact['binary']:
        for out, left, right, value in op['terms']:
            if value != 0 and coordinate_masks[op['out']][out] != (
                    coordinate_masks[op['left']][left] |
                    coordinate_masks[op['right']][right]):
                raise AssertionError(
                    f"binary color grading mismatch: {op['name']}")

    active = []
    for actual_mask in range(mask_count):
        active.append({
            interface: {
                coordinate
                for coordinate, required_mask in enumerate(
                    coordinate_masks[interface])
                if required_mask & ~actual_mask == 0
            }
            for interface in interfaces
        })
    return coordinate_masks, active


def render_full_header(artifact):
    K = artifact['K']
    coordinate_masks, color_mask_support = _coordinate_color_masks(artifact)
    mask_count = 1 << K
    lines = [
        "// Generated by tools/counting/generate_counters.py --full; do not edit.",
        f"// QRAFT counter tables for K={K} colours.",
        "#pragma once",
        "#include <array>",
        "#include <cstdint>",
        "#include <string_view>",
        f"namespace qraft {{ namespace generated_counters_k{K} {{",
        "struct UnaryTerm { std::uint8_t out, in; std::int8_t sign; };",
        "struct BinaryTerm { std::uint8_t out, left, right; std::int8_t sign; };",
        "struct VectorTerm { std::uint8_t index; std::int8_t sign; };",
        f"inline constexpr int palette_K = {K};",
    ]
    for interface in ('Closed', 'C1', 'C2', 'G'):
        lines.append(f"inline constexpr int states_{interface} = "
                     f"{artifact['state_counts'][interface]};")
        lines.append(f"inline constexpr int dim_{interface} = "
                     f"{artifact['ranks'][interface]};")
        values = ','.join(str(value)
                          for value in artifact['coordinate_degrees'][interface])
        lines.append(
            f"inline constexpr std::array<std::uint8_t, dim_{interface}> "
            f"degree_{interface} = {{{values}}};")
        lines.append(f"inline constexpr std::string_view state_sha256_{interface} = "
                     f"\"{artifact['state_manifest_sha256'][interface]}\";")
        lines.append(f"inline constexpr std::string_view projection_sha256_{interface} = "
                     f"\"{artifact['projection'][interface]['sha256']}\";")
        lines.append(f"inline constexpr int projection_terms_{interface} = "
                     f"{artifact['projection'][interface]['terms']};")
        values = ','.join(str(value)
                          for value in coordinate_masks[interface])
        lines.append(
            f"inline constexpr std::array<std::uint8_t, dim_{interface}> "
            f"coordinate_color_mask_{interface} = {{{values}}};")
        active_bits = []
        for mask in range(mask_count):
            bits = sum(1 << coordinate
                       for coordinate in color_mask_support[mask][interface])
            active_bits.append(f"UINT64_C({bits})")
        lines.append(
            f"inline constexpr std::array<std::uint64_t, {mask_count}> "
            f"active_coordinate_bits_{interface} = "
            f"{{{','.join(active_bits)}}};")
        if interface == 'C1':
            half = artifact['ranks']['C1'] // 2
            low = (1 << half) - 1
            for tag, block in (('U', low << half), ('D', low)):
                packed = []
                for mask in range(mask_count):
                    bits = sum(1 << coordinate for coordinate in
                               color_mask_support[mask][interface])
                    packed.append(f"UINT64_C({bits & block})")
                lines.append(
                    f"inline constexpr std::array<std::uint64_t, {mask_count}> "
                    f"packed_coordinate_bits_C1{tag} = "
                    f"{{{','.join(packed)}}};")

    def emit_unary(name, terms, output_dim, input_interface,
                   output_interface, emit_color_mask=True, mask_variants=None):
        assert all(terms[i][0] <= terms[i + 1][0]
                   for i in range(len(terms) - 1))
        lines.append(f"inline constexpr std::array<UnaryTerm, {len(terms)}> {name} = {{{{")
        for k, p, value in terms:
            assert value.denominator == 1 and abs(value.numerator) == 1
            lines.append(f"  {{{k},{p},{value.numerator}}},")
        lines.append("}};")
        lines.append(
            f"inline void eval_{name}(const double* __restrict in, "
            "double* __restrict out) {")
        by_output = defaultdict(list)
        for term in terms:
            by_output[term[0]].append(term)
        for k in range(output_dim):
            lines.append("  {")
            lines.append("    double sum = 0.0;")
            for _, p, value in by_output[k]:
                op = "+=" if value.numerator > 0 else "-="
                lines.append(f"    sum {op} in[{p}];")
            lines.append(f"    out[{k}] = sum;")
            lines.append("  }")
        lines.append("}")
        lines.append(
            f"inline void eval_dual_{name}("
            "const double* __restrict current_in, "
            "const double* __restrict background_in, "
            "double* __restrict current_out, "
            "double* __restrict background_out) {")
        for k, group in itertools.groupby(terms, key=lambda term: term[0]):
            group = list(group)
            lines.append(f"  double current_sum_{k} = 0.0;")
            lines.append(f"  double background_sum_{k} = 0.0;")
            for _, p, value in group:
                op = "+=" if value.numerator > 0 else "-="
                lines.append(f"  current_sum_{k} {op} current_in[{p}];")
                lines.append(
                    f"  background_sum_{k} {op} background_in[{p}];")
            lines.append(f"  current_out[{k}] = current_sum_{k};")
            lines.append(f"  background_out[{k}] = background_sum_{k};")
        present_outputs = {term[0] for term in terms}
        for k in range(output_dim):
            if k not in present_outputs:
                lines.append(f"  current_out[{k}] = 0.0;")
                lines.append(f"  background_out[{k}] = 0.0;")
        lines.append("}")
        for suffix, input_allowed, output_allowed in (
                mask_variants or [('', None, None)]):
          if not emit_color_mask:
            break
          if True:
            lines.append(
                f"inline int eval_color_mask_{name}{suffix}(std::uint8_t mask, "
                "const double* in, double* out) {")
            lines.append("  (void)in;")
            lines.append("  (void)out;")
            lines.append("  switch (mask) {")
            for mask in range(mask_count):
                active_input = color_mask_support[mask][input_interface]
                active_output = color_mask_support[mask][output_interface]
                if input_allowed is not None:
                    active_input = active_input & input_allowed
                if output_allowed is not None:
                    active_output = active_output & output_allowed
                input_rank = {
                    coordinate: rank
                    for rank, coordinate in enumerate(sorted(active_input))
                }
                output_rank = {
                    coordinate: rank
                    for rank, coordinate in enumerate(sorted(active_output))
                }
                masked_terms = [term for term in terms
                                if term[1] in active_input]
                if any(term[0] not in active_output
                       for term in masked_terms):
                    raise AssertionError(
                        f"masked unary closure failed: {name} mask={mask}")
                masked_by_output = defaultdict(list)
                for term in masked_terms:
                    masked_by_output[term[0]].append(term)
                lines.append(f"    case {mask}:")
                for k in sorted(active_output):
                    lines.append("      {")
                    lines.append("        double sum = 0.0;")
                    for _, p, value in masked_by_output[k]:
                        op = "+=" if value.numerator > 0 else "-="
                        lines.append(
                            f"        sum {op} in[{input_rank[p]}];")
                    lines.append(
                        f"        out[{output_rank[k]}] = sum;")
                    lines.append("      }")
                lines.append(f"      return {len(masked_terms)};")
            lines.append("    default: return 0;")
            lines.append("  }")
            lines.append("}")

    def emit_binary(name, terms, output_dim, left_interface,
                    right_interface, output_interface, left_allowed=None,
                    right_allowed=None, output_allowed=None):
        assert all(terms[i][0] <= terms[i + 1][0]
                   for i in range(len(terms) - 1))
        lines.append(f"inline constexpr std::array<BinaryTerm, {len(terms)}> {name} = {{{{")
        for k, p, q, value in terms:
            assert value.denominator == 1 and abs(value.numerator) == 1
            lines.append(f"  {{{k},{p},{q},{value.numerator}}},")
        lines.append("}};")
        lines.append(
            f"inline void eval_{name}(const double* __restrict left, "
            "const double* __restrict right, double* __restrict out) {")
        by_output = defaultdict(list)
        for term in terms:
            by_output[term[0]].append(term)
        for k in range(output_dim):
            lines.append("  {")
            lines.append("    double sum = 0.0;")
            for _, p, q, value in by_output[k]:
                op = "+=" if value.numerator > 0 else "-="
                lines.append(
                    f"    sum {op} left[{p}] * right[{q}];")
            lines.append(f"    out[{k}] = sum;")
            lines.append("  }")
        lines.append("}")
        lines.append(
            f"inline void eval_dual_{name}("
            "const double* __restrict current_left, "
            "const double* __restrict background_left, "
            "const double* __restrict current_right, "
            "const double* __restrict background_right, "
            "double* __restrict current_out, "
            "double* __restrict background_out) {")
        for k, group in itertools.groupby(terms, key=lambda term: term[0]):
            group = list(group)
            lines.append(f"  double current_sum_{k} = 0.0;")
            lines.append(f"  double background_sum_{k} = 0.0;")
            for _, p, q, value in group:
                op = "+=" if value.numerator > 0 else "-="
                lines.append(
                    f"  current_sum_{k} {op} current_left[{p}] * "
                    f"current_right[{q}];")
                lines.append(
                    f"  background_sum_{k} {op} background_left[{p}] * "
                    f"background_right[{q}];")
            lines.append(f"  current_out[{k}] = current_sum_{k};")
            lines.append(f"  background_out[{k}] = background_sum_{k};")
        present_outputs = {term[0] for term in terms}
        for k in range(output_dim):
            if k not in present_outputs:
                lines.append(f"  current_out[{k}] = 0.0;")
                lines.append(f"  background_out[{k}] = 0.0;")
        lines.append("}")
        lines.append(
            f"inline int eval_color_mask_{name}(std::uint8_t left_mask, "
            "std::uint8_t right_mask, const double* left, "
            "const double* right, double* out) {")
        lines.append("  (void)left;")
        lines.append("  (void)right;")
        lines.append("  (void)out;")
        lines.append(
            f"  switch ((static_cast<unsigned>(left_mask) << {K}) | "
            "right_mask) {")
        for left_mask in range(mask_count):
            for right_mask in range(mask_count):
                output_mask = left_mask | right_mask
                active_left = color_mask_support[left_mask][left_interface]
                active_right = color_mask_support[right_mask][right_interface]
                active_output = color_mask_support[output_mask][output_interface]
                if left_allowed is not None:
                    active_left = active_left & left_allowed
                if right_allowed is not None:
                    active_right = active_right & right_allowed
                if output_allowed is not None:
                    active_output = active_output & output_allowed
                left_rank = {
                    coordinate: rank
                    for rank, coordinate in enumerate(sorted(active_left))
                }
                right_rank = {
                    coordinate: rank
                    for rank, coordinate in enumerate(sorted(active_right))
                }
                output_rank = {
                    coordinate: rank
                    for rank, coordinate in enumerate(sorted(active_output))
                }
                masked_terms = [
                    term for term in terms
                    if term[1] in active_left and term[2] in active_right
                ]
                if any(term[0] not in active_output
                       for term in masked_terms):
                    raise AssertionError(
                        f"masked binary closure failed: {name} "
                        f"left={left_mask} right={right_mask}")
                masked_by_output = defaultdict(list)
                for term in masked_terms:
                    masked_by_output[term[0]].append(term)
                key = (left_mask << K) | right_mask
                lines.append(f"    case {key}:")
                for k in sorted(active_output):
                    lines.append("      {")
                    lines.append("        double sum = 0.0;")
                    for _, p, q, value in masked_by_output[k]:
                        op = "+=" if value.numerator > 0 else "-="
                        lines.append(
                            f"        sum {op} left[{left_rank[p]}] * "
                            f"right[{right_rank[q]}];")
                    lines.append(
                        f"        out[{output_rank[k]}] = sum;")
                    lines.append("      }")
                lines.append(f"      return {len(masked_terms)};")
        lines.append("    default: return 0;")
        lines.append("  }")
        lines.append("}")

    def emit_vector(name, vector):
        terms = [(i, value) for i, value in enumerate(vector) if value != 0]
        lines.append(f"inline constexpr std::array<VectorTerm, {len(terms)}> {name} = {{{{")
        for i, value in terms:
            assert value.denominator == 1 and abs(value.numerator) == 1
            lines.append(f"  {{{i},{value.numerator}}},")
        lines.append("}};")

    # Production C1 is the direct sum of two orientation blocks: coordinates
    # [0, h) are the D-only sort (the draft's C3) and [h, 2h) the U-only sort
    # (the draft's C1). Every operation slot of interface C1 touches exactly
    # one block (REV maps one block onto the other), so the packed kernels
    # index each C1 slot relative to its own block; the other block is
    # identically zero there and is not stored.
    half = artifact['ranks']['C1'] // 2
    block_d = set(range(half))
    block_u = set(range(half, 2 * half))

    def c1_block(indices, what):
        indices = set(indices)
        if not indices:
            return None
        if indices <= block_u:
            return block_u
        if indices <= block_d:
            return block_d
        raise AssertionError(f"C1 slot spans both orientation blocks: {what}")

    for op in artifact['unary']:
        variants = None
        if op['in'] == 'C1' and op['out'] == 'C1':
            ud = [t for t in op['terms'] if t[1] in block_u]
            du = [t for t in op['terms'] if t[1] in block_d]
            if (not all(t[0] in block_d for t in ud) or
                    not all(t[0] in block_u for t in du) or
                    len(ud) + len(du) != len(op['terms'])):
                raise AssertionError("C1 reversal does not swap blocks")
            variants = [('_UD', block_u, block_d), ('_DU', block_d, block_u)]
        elif op['in'] == 'C1' or op['out'] == 'C1':
            in_allowed = (c1_block([t[1] for t in op['terms']], op['name'])
                          if op['in'] == 'C1' else None)
            out_allowed = (c1_block([t[0] for t in op['terms']], op['name'])
                           if op['out'] == 'C1' else None)
            variants = [('', in_allowed, out_allowed)]
        emit_unary(_identifier(op), op['terms'], artifact['ranks'][op['out']],
                   op['in'], op['out'], mask_variants=variants)
    for op in artifact['binary']:
        left_allowed = (c1_block([t[1] for t in op['terms']], op['name'])
                        if op['left'] == 'C1' else None)
        right_allowed = (c1_block([t[2] for t in op['terms']], op['name'])
                         if op['right'] == 'C1' else None)
        output_allowed = (c1_block([t[0] for t in op['terms']], op['name'])
                          if op['out'] == 'C1' else None)
        emit_binary(_identifier(op), op['terms'], artifact['ranks'][op['out']],
                    op['left'], op['right'], op['out'],
                    left_allowed, right_allowed, output_allowed)
    emit_unary('READOUT_CLOSED', sparse_unary(artifact['readout']),
               len(artifact['channels']), None, None, False)
    emit_vector('LEAF_ABSENT_C1', artifact['leaf_absent'])
    for color, vector in enumerate(artifact['leaf_present']):
        emit_vector(f'LEAF_PRESENT_{color}_C1', vector)
    for name in sorted(artifact['constants']):
        emit_vector(name, artifact['constants'][name])
    lines += [
        f"inline constexpr int n_channels = {len(artifact['channels'])};",
        f"inline constexpr std::string_view channels = \"{','.join(artifact['channels'])}\";",
        f"inline constexpr std::string_view specification_sha256 = \"{artifact['spec_sha256']}\";",
        f"inline constexpr std::string_view artifact_sha256 = \"{artifact['artifact_sha256']}\";",
        f"}} }}  // namespace generated_counters_k{K} / qraft",
        "",
    ]
    return '\n'.join(lines)


def full_header_path(repo, K):
    return os.path.join(repo, 'src', 'counting', 'generated',
                        f'counters_k{K}.hpp')


def write_or_check_full_headers(repo, artifacts, check_only):
    failed = False
    for artifact in artifacts:
        path = full_header_path(repo, artifact['K'])
        rendered = render_full_header(artifact)
        if check_only:
            try:
                with open(path, 'r', encoding='ascii') as handle:
                    existing = handle.read()
            except FileNotFoundError:
                existing = None
            if existing != rendered:
                print(f"  MISMATCH {path}")
                failed = True
            else:
                print(f"  MATCH {path}")
        else:
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, 'w', encoding='ascii', newline='\n') as handle:
                handle.write(rendered)
            print(f"  emitted {path} ({len(rendered)} bytes)")
    return not failed


def _repo_root():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--full', action='store_true', help='derive, verify, and emit full grammar')
    mode.add_argument('--verify-full', action='store_true', help='derive and verify without writing')
    mode.add_argument('--check-full', action='store_true', help='verify checked-in full headers')
    mode.add_argument('--digest-full', action='store_true', help=argparse.SUPPRESS)
    mode.add_argument('--check-reproducible', action='store_true',
                      help='compare full derivations under two Python hash seeds')
    args = parser.parse_args()

    if args.check_reproducible:
        outputs = []
        for seed in ('1', '2'):
            env = dict(os.environ)
            env['PYTHONHASHSEED'] = seed
            env['PYTHONDONTWRITEBYTECODE'] = '1'
            outputs.append(subprocess.check_output(
                [sys.executable, '-B', __file__, '--digest-full'],
                env=env, text=True))
        if outputs[0] != outputs[1]:
            print('FULL ARTIFACT REPRODUCIBILITY FAILED')
            return 1
        print(outputs[0], end='')
        print('FULL ARTIFACTS REPRODUCIBLE ACROSS PYTHONHASHSEED=1,2')
        return 0

    if args.full or args.verify_full or args.check_full or args.digest_full:
        artifacts = [derive_full_grammar(K, verbose=not args.digest_full)
                     for K in (2, 4)]
        if args.digest_full:
            for artifact in artifacts:
                print(f"K={artifact['K']} {artifact['artifact_sha256']}")
            return 0
        if args.full or args.check_full:
            if not write_or_check_full_headers(
                    _repo_root(), artifacts, check_only=args.check_full):
                return 1
        print('FULL K=2/K=4 GRAMMAR DERIVATION VERIFIED')
        return 0

    parser.print_help()
    return 2

if __name__ == '__main__':
    sys.exit(main())
