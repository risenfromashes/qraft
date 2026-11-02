"""Exact minimal counter dimensions for the SODA'13 quartet families (colour 0 = outside,
colours 1..x = children of v), per sort and per degree, for a range of x; polynomial fits."""
import sys, json, time
from fractions import Fraction
import clusters as C

SORTS = ('T', 'C1', 'G', 'C2')

def dims_for(channels, K, out='T'):
    res = {}
    for s in SORTS:
        d = C.minimal_dims(s, K, channels, out=out)
        res[s] = [d[m][1] for m in range(C.D + 1)]
        res[s + '_shapes'] = [d[m][0] for m in range(C.D + 1)]
    return res

def channels_for(fam, x):
    """family spec: letters A,B,C,D,E (unrooted families, output T) joined by '+';
    or per-run rooted families Ar, Br (lower triangle, output G) and Aall, Ball (all cells, output G)."""
    chans = {}; out = None
    for t in fam.split('+'):
        this_out = 'G' if t in ('Ar', 'Br', 'Aall', 'Ball') else 'T'
        if out is not None and out != this_out:
            raise ValueError('A joint task must use one output sort; do not mix rooted and unrooted families')
        out = this_out
        if t in ('A', 'B', 'C', 'D', 'E'):
            chans[t] = C.soda_channels(x, t)
        elif t in ('Ar', 'Br'):
            chans[t] = C.brodal_run_channel(x, t[0], 'lower'); out = 'G'
        elif t in ('Aall', 'Ball'):
            chans[t] = C.brodal_run_channel(x, t[0], 'all'); out = 'G'
        else:
            raise ValueError(t)
    return {k: v for k, v in chans.items() if v}, out

def channels_from_certificate(record):
    """Resolve the actual output, including a balanced or distance variant.

    Early certificates retained an Ar_Br task tag for all three variants.
    Ignoring the variant silently selects the wrong readout and matrix rows.
    """
    name, K = record['task'], record['K']
    variant = record.get('variant')
    if variant in ('Balanced', 'Distance'):
        channels = {}
        for which in ('A', 'B'):
            channel = C.brodal_run_channel(K-1, which)
            for pattern in channel:
                form, split = C.rooted_quartet_info(pattern)
                colour_form = C.t1_type([c for side in split for c in side])
                channel[pattern] = 1 if C.ORDER[form] == C.ORDER[colour_form] else 2
            channels[which] = channel
        if variant == 'Distance':
            a, b = channels['A'], channels['B']
            channels = {'Distance': {p:2*a.get(p,0)+b.get(p,0) for p in a.keys() | b.keys()}}
        out = 'G'
    elif variant is not None:
        raise ValueError('Unknown certificate variant: ' + str(variant))
    elif name.startswith('QRAFT'):
        channels, out = C.qraft_channels(K), 'T'
    else:
        channels, out = channels_for(name.rsplit('_x',1)[0].replace('_','+'),K-1)
    assert record['output_sort'] == out
    return channels, out

def lagrange_fit(points):
    """exact polynomial through the points (x_i, y_i); returns coefficients low->high"""
    n = len(points)
    coeffs = [Fraction(0)] * n
    for i, (xi, yi) in enumerate(points):
        # basis polynomial
        num = [Fraction(1)]
        den = Fraction(1)
        for j, (xj, _) in enumerate(points):
            if j == i:
                continue
            num = polymul(num, [Fraction(-xj), Fraction(1)])
            den *= (xi - xj)
        for k in range(len(num)):
            coeffs[k] += yi * num[k] / den
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def polymul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out

def polyval(c, x):
    return sum(ci * x ** i for i, ci in enumerate(c))

def polystr(c):
    terms = []
    for i in range(len(c) - 1, -1, -1):
        if c[i] == 0:
            continue
        co = c[i]
        s = ('%s' % co) if co.denominator == 1 else ('%s/%s' % (co.numerator, co.denominator))
        if i == 0:
            terms.append(s)
        elif i == 1:
            terms.append(('' if co == 1 else s) + 'x')
        else:
            terms.append(('' if co == 1 else s) + 'x^%d' % i)
    return ' + '.join(terms) if terms else '0'

def fit_and_check(xs, ys, deg):
    """fit a polynomial of degree deg to the first deg+1 points and check the rest"""
    pts = list(zip(xs, ys))
    if len(pts) < deg + 2:
        return None, False
    c = lagrange_fit(pts[:deg + 1])
    ok = all(polyval(c, x) == y for x, y in pts[deg + 1:])
    return c, ok

if __name__ == '__main__':
    families = sys.argv[1].split(',') if len(sys.argv) > 1 else ['A', 'B', 'E', 'A+B']
    xs = [int(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1, 2, 3, 4, 5, 6, 7]
    results = {}
    for fam in families:
        results[fam] = {}
        for x in xs:
            K = x + 1
            channels, out = channels_for(fam, x)
            if not channels:
                continue
            t0 = time.time()
            r = dims_for(channels, K, out)
            results[fam][x] = r
            print('%-6s x=%d  ' % (fam, x) + '  '.join('%s=%d %s' % (s, sum(r[s]), r[s]) for s in SORTS)
                  + '   (%.1fs)' % (time.time() - t0))
            sys.stdout.flush()
        # fits
        for s in SORTS:
            xs_f = [x for x in xs if x in results[fam]]
            ys = [sum(results[fam][x][s]) for x in xs_f]
            # fit on x >= 4 (all families fully present)
            sel = [(x, y) for x, y in zip(xs_f, ys) if x >= 4]
            if len(sel) >= 4:
                c, ok = fit_and_check([p[0] for p in sel], [p[1] for p in sel], 2)
                print('   fit %-6s %-3s: %s   %s' % (fam, s, polystr(c) if c else None, 'OK' if ok else 'FAILS (not quadratic on these points)'))
                for m in range(C.D + 1):
                    ym = [results[fam][x][s][m] for x, _ in sel]
                    cm, okm = fit_and_check([p[0] for p in sel], ym, min(m, C.D - m))
                    print('        block m=%d: %s  %s   values %s' % (m, polystr(cm) if cm else None, 'OK' if okm else 'FAIL', ym))
    json.dump({fam: {str(x): v for x, v in d.items()} for fam, d in results.items()},
              open('soda_results_%s.json' % ('_'.join(families).replace('+', 'p')), 'w'), indent=1)
