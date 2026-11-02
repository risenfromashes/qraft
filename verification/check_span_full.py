"""App. B (5): the natural counters span exactly the backward spaces of the census automaton,
computed as the closure of the gamma^a (kmw_closure_full.py), for every sort, k = 2 and 4."""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE / 'natwta'), str(HERE / 'natwta' / 'code')]
import kmw_closure_full as KF
from check_natural_closure import natural_counters, in_span
for K in (2, 4):
    sh, idx, B, dims, prof = KF.run(K, verbose=False)
    nat, ch = natural_counters(K, sh)
    for s in KF.SORTS:
        assert len(nat[s]) == dims[s], (K, s, len(nat[s]), dims[s])
        for lab, vec in nat[s]:
            assert in_span(B[s], vec), (K, s, lab)
    print(f'k={K}: natural counters span the backward spaces, dims {dims}: ok')
