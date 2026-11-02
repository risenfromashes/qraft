#!/usr/bin/env python3
"""Run every selected small-tree case through QRAFT's C++ counting.

This driver owns only the enumeration, a strict quartet oracle and the text
protocol; the C++ runner (qraft_tests production_small_trees_protocol) runs
the HDT, Count, Contract and Extract/Contract code.

The case index is identical to check_small_trees.py: labeled binary trees,
every unique internal-edge contraction, K=2/K=4 color-composition orbit
representatives, and integer/fractional finite weights.  All shards together
partition that index exactly.
"""
from __future__ import annotations

import argparse
import itertools
import subprocess
import sys
from collections import defaultdict, deque
from fractions import Fraction

import shapes as G
import check_small_trees as S


EXPECTED_BINARY = {3: 1, 4: 3, 5: 15, 6: 105, 7: 945, 8: 10395}
EXPECTED_TOPOLOGIES = {4: 4, 5: 26, 6: 236, 7: 2752, 8: 39208}
EXPECTED_CASES = {4: 320, 5: 3544, 6: 46496, 7: 751008, 8: 14395392}


def path_nodes(adj, source, target):
    parent = {source: None}
    pending = deque([source])
    while pending and target not in parent:
        node = pending.popleft()
        for neighbor in adj[node]:
            if neighbor not in parent:
                parent[neighbor] = node
                pending.append(neighbor)
    if target not in parent:
        raise AssertionError("strict quartet oracle received a disconnected tree")
    result = set()
    node = target
    while node is not None:
        result.add(node)
        node = parent[node]
    return result


def topology_from_code(code):
    """Rebuild a deterministic adjacency from the canonical rooted code."""
    adjacency = defaultdict(set)
    position = 0
    next_internal = -1

    def parse():
        nonlocal position, next_internal
        if position >= len(code):
            raise AssertionError("truncated canonical topology code")
        if code[position] == "L":
            position += 1
            begin = position
            while position < len(code) and code[position].isdigit():
                position += 1
            if begin == position:
                raise AssertionError("canonical topology leaf has no label")
            return ("L", int(code[begin:position]))
        if code[position] != "(":
            raise AssertionError("invalid canonical topology token")
        position += 1
        node = next_internal
        next_internal -= 1
        child_count = 0
        while position < len(code) and code[position] != ")":
            child = parse()
            S.link(adjacency, node, child)
            child_count += 1
        if position >= len(code) or code[position] != ")":
            raise AssertionError("unterminated canonical topology center")
        position += 1
        if child_count == 0:
            raise AssertionError("empty canonical topology center")
        return node

    parse()
    if position != len(code):
        raise AssertionError("trailing canonical topology data")
    return adjacency


def resolved_quartet_pairings(adj, leaf_count):
    """Return the unique strict 2|2 pairing for each resolved quartet."""
    leaves = [("L", leaf) for leaf in range(leaf_count)]
    paths = {}
    for a, b in itertools.combinations(range(leaf_count), 2):
        paths[a, b] = path_nodes(adj, leaves[a], leaves[b])

    resolved = []
    pairings = ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2))
    for quartet in itertools.combinations(range(leaf_count), 4):
        for ia, ib, ic, id_ in pairings:
            a, b = sorted((quartet[ia], quartet[ib]))
            c, d = sorted((quartet[ic], quartet[id_]))
            if paths[a, b].isdisjoint(paths[c, d]):
                resolved.append((a, b, c, d))
                break
    return resolved


def weighted_pairings(resolved, weights):
    return [
        (a, b, c, d, weights[a] * weights[b] * weights[c] * weights[d])
        for a, b, c, d in resolved
    ]


def strict_quartet_counts(pairings, colors, palette):
    channels = ("Q00_11", "Q01_01") if palette == 2 else (
        "Q01_23", "Q02_13", "Q03_12")
    result = {channel: Fraction(0) for channel in channels}
    for a, b, c, d, weight in pairings:
        left = tuple(sorted((colors[a], colors[b])))
        right = tuple(sorted((colors[c], colors[d])))
        channel = G._channel(left, right, palette)
        if channel in result:
            result[channel] += weight
    values = [result[channel] for channel in channels]
    if palette == 2:
        values.append(Fraction(0))
    return values


def topology_protocol(adj, leaf_count):
    internal = sorted(node for node in adj if isinstance(node, int))
    node_id = {("L", leaf): leaf for leaf in range(leaf_count)}
    node_id.update({node: leaf_count + index
                    for index, node in enumerate(internal)})
    edges = sorted(
        (min(node_id[a], node_id[b]), max(node_id[a], node_id[b]))
        for a, neighbors in adj.items() for b in neighbors
        if node_id[a] < node_id[b])
    return len(node_id), node_id[internal[0]], edges


def background_colors(case_index, leaf_count, palette, seed):
    # Mixed immutable backgrounds are deterministic and independent of the
    # block-color representative used for the current state.
    return [int((seed + case_index + 3 * leaf + leaf * leaf) % palette)
            for leaf in range(leaf_count)]


def decimal(value):
    return format(float(value), ".17g")


def write_case(stream, *, case_index, leaf_count, palette, fractional,
               topology_code, topology, colors, weights, resolved, seed):
    all_owners = (1 << leaf_count) - 1
    selector_mask = case_index & all_owners
    nested_bits = ((case_index >> leaf_count) ^
                   (case_index * 0x9E3779B97F4A7C15)) & all_owners
    nested_mask = selector_mask & nested_bits
    if selector_mask:
        recolor_owner = (selector_mask & -selector_mask).bit_length() - 1
    else:
        recolor_owner = 0
    recolor_color = (colors[recolor_owner] + 1) % palette
    background = background_colors(case_index, leaf_count, palette, seed)

    selected_colors = [
        colors[owner] if selector_mask & (1 << owner) else background[owner]
        for owner in range(leaf_count)
    ]
    recolored_colors = list(selected_colors)
    if selector_mask:
        recolored_colors[recolor_owner] = recolor_color
    nested_colors = [
        recolored_colors[owner] if nested_mask & (1 << owner)
        else background[owner]
        for owner in range(leaf_count)
    ]

    weighted = weighted_pairings(resolved, weights)
    expected = [
        strict_quartet_counts(weighted, coloring, palette)
        for coloring in (colors, background, selected_colors,
                         recolored_colors, nested_colors)
    ]
    node_count, root, edges = topology_protocol(topology, leaf_count)
    kind = "fractional" if fractional else "integer"
    stream.write(
        f"CASE {case_index} {leaf_count} {palette} {kind} "
        f"{selector_mask} {nested_mask} {recolor_owner} {recolor_color} "
        f"{node_count} {root} {topology_code}\n")
    stream.write("COLORS " + " ".join(map(str, colors)) + "\n")
    stream.write("BACKGROUND " + " ".join(map(str, background)) + "\n")
    stream.write("WEIGHTS " + " ".join(
        f"{weight.numerator} {weight.denominator}" for weight in weights
    ) + "\n")
    stream.write("EDGES " + str(len(edges)) + " " + " ".join(
        f"{a} {b}" for a, b in edges
    ) + "\n")
    stream.write("EXPECTED " + " ".join(
        decimal(value) for counts in expected for value in counts
    ) + "\nEND\n")


def expected_selected_cases(global_cases, shard_count, shard_index):
    if shard_index >= global_cases:
        return 0
    return 1 + (global_cases - 1 - shard_index) // shard_count


def enumerate_cases(args, stream=None):
    binaries = S.binary_topologies(args.max_leaves)
    global_case = 0
    emitted = 0
    topology_counts = {}
    for leaf_count in range(4, args.max_leaves + 1):
        if len(binaries[leaf_count]) != EXPECTED_BINARY[leaf_count]:
            raise AssertionError(
                f"binary topology count n={leaf_count}: "
                f"{len(binaries[leaf_count])} != {EXPECTED_BINARY[leaf_count]}")
        topologies = S.all_contractions(binaries[leaf_count])
        topology_counts[leaf_count] = len(topologies)
        if len(topologies) != EXPECTED_TOPOLOGIES[leaf_count]:
            raise AssertionError(
                f"contracted topology count n={leaf_count}: "
                f"{len(topologies)} != {EXPECTED_TOPOLOGIES[leaf_count]}")
        for topology_code in sorted(topologies):
            # Do not retain the set-iteration-dependent representative chosen
            # by the deduplication dictionary.  The code is the canonical ABI.
            topology = topology_from_code(topology_code)
            resolved = None
            for palette in (2, 4):
                for counts in S.compositions(leaf_count, palette):
                    colors = S.block_colors(counts)
                    for fractional in (False, True):
                        case_index = global_case
                        global_case += 1
                        if case_index % args.shard_count != args.shard_index:
                            continue
                        if stream is None:
                            emitted += 1
                            continue
                        if resolved is None:
                            resolved = resolved_quartet_pairings(
                                topology, leaf_count)
                        weights = S.weights_for(
                            leaf_count, fractional, args.seed)
                        write_case(
                            stream, case_index=case_index,
                            leaf_count=leaf_count, palette=palette,
                            fractional=fractional,
                            topology_code=topology_code, topology=topology,
                            colors=colors, weights=weights,
                            resolved=resolved, seed=args.seed)
                        emitted += 1
    if global_case != EXPECTED_CASES[args.max_leaves]:
        raise AssertionError(
            f"global case count {global_case} != "
            f"{EXPECTED_CASES[args.max_leaves]}")
    wanted = expected_selected_cases(
        global_case, args.shard_count, args.shard_index)
    if emitted != wanted:
        raise AssertionError(f"emitted case count {emitted} != {wanted}")
    return global_case, emitted, topology_counts


def run(args):
    if args.max_leaves < 4 or args.max_leaves > 8:
        raise ValueError("--max-leaves must be in [4,8]")
    if args.shard_count <= 0 or not 0 <= args.shard_index < args.shard_count:
        raise ValueError("invalid deterministic shard")
    if args.count_only:
        global_cases, emitted, topology_counts = enumerate_cases(args)
        print("production-small coverage census: "
              f"seed={args.seed} shard={args.shard_index}/{args.shard_count} "
              f"max_leaves={args.max_leaves} selected={emitted}/"
              f"{global_cases} topologies={topology_counts}")
        return 0
    if not args.runner:
        raise ValueError("--runner is required unless --count-only is used")

    process = subprocess.Popen(
        [args.runner, "production_small_trees_protocol"],
        stdin=subprocess.PIPE, text=True)
    assert process.stdin is not None
    broken_pipe = False
    try:
        process.stdin.write(
            f"QRAFT_PRODUCTION_SMALL_TREES_V1 {args.seed} "
            f"{args.shard_index} {args.shard_count} {args.max_leaves}\n")
        global_cases, emitted, topology_counts = enumerate_cases(
            args, process.stdin)
        process.stdin.write(f"DONE {global_cases} {emitted}\n")
    except BrokenPipeError:
        broken_pipe = True
        global_cases = emitted = 0
        topology_counts = {}
    finally:
        try:
            process.stdin.close()
        except BrokenPipeError:
            broken_pipe = True
    return_code = process.wait()
    if return_code != 0 or broken_pipe:
        raise RuntimeError(
            f"production-small runner failed with exit code {return_code}")
    print("production-small driver completed: "
          f"seed={args.seed} shard={args.shard_index}/{args.shard_count} "
          f"max_leaves={args.max_leaves} streamed={emitted}/"
          f"{global_cases} topologies={topology_counts}")
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner")
    parser.add_argument("--max-leaves", type=int, default=6)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260713)
    parser.add_argument("--count-only", action="store_true")
    return run(parser.parse_args())


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"check_production_small_trees.py failed: {error}",
              file=sys.stderr)
        raise
