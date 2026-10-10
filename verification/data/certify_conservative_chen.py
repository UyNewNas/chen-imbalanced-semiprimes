"""Exact-arithmetic Darboux certificate for a CONSERVATIVE NUMERICAL functional.

This certifies the explicit integrals below, not their number-theoretic use.
It uses integer fixed-point interval arithmetic and a proved log power-series
remainder. No floating-point computation contributes to certificate bounds.

Analytic substitutions:
* All distribution levels are the single-variable Lichtman level.
* H=h=0; S1/S2 use classical level 1/2; switching coefficient is 8.
* Upper F is exact up to s=5, and bounded by F(5) beyond 5.
* Pair lower f is replaced by 2 exp(gamma) log(s-1)/s for s>2, zero below.
* Buchstab is bounded by 4/7 in S14/S15; S11 has an explicit upper correction.
* The factors exp(gamma) cancel: the functional uses logs and rational inputs.

Run: python certify_conservative_chen.py --pair-mesh 128 --line-mesh 1000
"""
import argparse
import json
import time
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

BITS = 80
Q = 1 << BITS
LOG_TERMS = 32


def ceildiv(n, d):
    assert d > 0
    return -((-n) // d)


class I:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None, *, raw=False):
        if raw:
            self.lo, self.hi = int(lo), int(hi)
        elif hi is None:
            x = Fraction(lo)
            self.lo = (x.numerator * Q) // x.denominator
            self.hi = ceildiv(x.numerator * Q, x.denominator)
        else:
            a, b = Fraction(lo), Fraction(hi)
            self.lo = (a.numerator * Q) // a.denominator
            self.hi = ceildiv(b.numerator * Q, b.denominator)
        assert self.lo <= self.hi

    @staticmethod
    def coerce(x):
        return x if isinstance(x, I) else I(x)

    def __add__(self, other):
        other = I.coerce(other)
        return I(self.lo + other.lo, self.hi + other.hi, raw=True)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo, raw=True)

    def __sub__(self, other):
        return self + (-I.coerce(other))

    def __rsub__(self, other):
        return I.coerce(other) - self

    def __mul__(self, other):
        other = I.coerce(other)
        products = [self.lo * other.lo, self.lo * other.hi,
                    self.hi * other.lo, self.hi * other.hi]
        return I(min(products) // Q, ceildiv(max(products), Q), raw=True)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = I.coerce(other)
        assert other.lo > 0, (other.lo, other.hi)
        lows = [n * Q // d for n in (self.lo, self.hi) for d in (other.lo, other.hi)]
        highs = [ceildiv(n * Q, d) for n in (self.lo, self.hi) for d in (other.lo, other.hi)]
        return I(min(lows), max(highs), raw=True)

    def __rtruediv__(self, other):
        return I.coerce(other) / self

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        result = I(1)
        for _ in range(n):
            result = result * self
        return result

    def maximum(self, point):
        point = I.coerce(point)
        return I(max(self.lo, point.lo), max(self.hi, point.hi), raw=True)


def atanh_log_bounds(ul, uh):
    """Bounds for 2 atanh(u); 0<=u<=1/3, with exact geometric tail."""
    assert 0 <= ul <= uh < Q
    u2l = ul * ul // Q
    u2h = ceildiv(uh * uh, Q)
    pl, ph = ul, uh
    sl = sh = 0
    for j in range(LOG_TERMS):
        odd = 2 * j + 1
        sl += pl // odd
        sh += ceildiv(ph, odd)
        pl = pl * u2l // Q
        ph = ceildiv(ph * u2h, Q)
    # For j>=M, 1/(2j+1)<=1/(2M+1).
    tail = ceildiv(2 * ph * Q, (2 * LOG_TERMS + 1) * (Q - u2h))
    return 2 * sl, 2 * sh + tail


LOG2 = atanh_log_bounds(Q // 3, ceildiv(Q, 3))


@lru_cache(maxsize=250000)
def log_point(q):
    """Certified log(q / 2**BITS), q an exact positive integer."""
    assert q > 0
    if q == Q:
        return 0, 0
    k = q.bit_length() - 1 - BITS
    if k >= 0:
        rl, rh = q // (1 << k), ceildiv(q, 1 << k)
    else:
        rl = rh = q * (1 << (-k))
    ul = (rl - Q) * Q // (rh + Q)
    uh = ceildiv((rh - Q) * Q, rl + Q)
    low, high = atanh_log_bounds(ul, uh)
    if k >= 0:
        return low + k * LOG2[0], high + k * LOG2[1]
    return low + k * LOG2[1], high + k * LOG2[0]


def log(x):
    x = I.coerce(x)
    assert x.lo > 0
    return I(log_point(x.lo)[0], log_point(x.hi)[1], raw=True)


def integrate1(fn, a, b, n):
    a, b = Fraction(a), Fraction(b)
    step = (b - a) / n
    sl = sh = 0
    for j in range(n):
        box = I(a + j * step, a + (j + 1) * step)
        val = fn(box)
        sl += val.lo
        sh += val.hi
    return I(sl, sh, raw=True) * I(step)


class JTable:
    """Certified J(s)=int_2^(s-1) log(v-1)/v dv, 3<=s<=5.

    q(v)=log(v-1)/v is increasing on [2,4], because
    v/(v-1)-log(v-1) >= 4/3-log(3)>0.
    Lower/upper prefix sums and outward cell inclusion enclose each query.
    """
    def __init__(self, n):
        self.n = n
        assert log(I(3)).hi < I(Fraction(4, 3)).lo
        self.lower = [0]
        self.upper = [0]
        step = I(Fraction(2, n))
        left = I(0)
        for j in range(n):
            v = I(Fraction(2) + Fraction(2 * (j + 1), n))
            right = log(v - 1) / v
            self.lower.append(self.lower[-1] + (left * step).lo)
            self.upper.append(self.upper[-1] + (right * step).hi)
            left = right

    def __call__(self, s):
        s = I.coerce(s)
        # All calls have true s<=5 (or s=5); fixed-point outward rounding
        # can put the enclosing interval a few units above that endpoint.
        # The monotone extension J(min(s,5)) agrees on the actual domain.
        xl = max(0, min(2 * Q, s.lo - 3 * Q))
        xh = max(0, min(2 * Q, s.hi - 3 * Q))
        il = xl * self.n // (2 * Q)
        ih = ceildiv(xh * self.n, 2 * Q)
        return I(self.lower[il], self.upper[ih], raw=True)


A = Fraction(100, 1149)
B = Fraction(50, 309)
MU = Fraction(25, 128)
D = Fraction(1, 2) - 3 * A
V0 = Fraction(57, 64)


def pair_lower(a, b, triangular, higher_t, n):
    """Two-dimensional interval Darboux sum on [a,b] x [0,1]."""
    dt = (b - a) / n
    dy = Fraction(1, n)
    yboxes = [I(Fraction(j, n), Fraction(j + 1, n)) for j in range(n)]
    sl = sh = 0
    for i in range(n):
        t = I(a + i * dt, a + (i + 1) * dt)
        jac = t - A if triangular else I(B - A)
        delta0 = V0 - Fraction(5, 2) * t if higher_t else (1 - t) / 2
        for y in yboxes:
            u = A + jac * y
            delta = delta0 - u
            numerator = log((delta / I(A) - 1).maximum(1))
            val = 4 * numerator * jac / (delta * t * u)
            sl += val.lo
            sh += val.hi
    return I(sl, sh, raw=True) * I(dt * dy)


def display(x):
    """Decimals are explanatory only; exact numerator/2**BITS is decisive."""
    from decimal import Decimal, localcontext
    with localcontext() as context:
        context.prec = 35
        return {
            'lower_numerator': str(x.lo), 'upper_numerator': str(x.hi),
            'denominator': str(Q),
            'lower_decimal': str(Decimal(x.lo) / Decimal(Q)),
            'upper_decimal': str(Decimal(x.hi) / Decimal(Q)),
        }


def certify(pair_mesh, line_mesh, j_mesh):
    start = time.monotonic()
    J = JTable(j_mesh)
    s1 = Fraction(1, 2) / A
    terms = {}
    terms['S1'] = 8 * (log(I(s1 - 1)) + integrate1(lambda v: J(v) / v, 3, s1 - 1, line_mesh))
    terms['S2'] = 8 * log(I(Fraction(1, 2) / B - 1))
    print('S1/S2 completed', flush=True)

    def upper_kernel(t, delta):
        return 4 * (1 + J(delta / I(A))) / (delta * t)

    split = 1 - 10 * A
    terms['S3'] = integrate1(lambda t: 4 * (1 + J(I(5))) / (I(5 * A) * t), A, split, line_mesh)
    terms['S3'] += integrate1(lambda t: upper_kernel(t, (1 - t) / 2), split, MU, line_mesh)
    terms['S4'] = integrate1(lambda t: upper_kernel(t, V0 - Fraction(5, 2) * t), MU, Fraction(1, 4), line_mesh)
    terms['S5'] = integrate1(lambda t: upper_kernel(t, V0 - Fraction(5, 2) * t), Fraction(1, 4), Fraction(57, 224), line_mesh)
    terms['S6'] = integrate1(lambda t: 4 / (t * (Fraction(1, 2) - t)), Fraction(57, 224), Fraction(1, 3), line_mesh)
    terms['S7'] = integrate1(lambda t: upper_kernel(t, V0 - Fraction(5, 2) * t), MU, D, line_mesh)
    print('Upper prime sums completed', flush=True)
    terms['S8'] = pair_lower(A, B, True, False, pair_mesh)
    print('S8 completed', flush=True)
    terms['S9'] = pair_lower(B, MU, False, False, pair_mesh)
    print('S9 completed', flush=True)
    terms['S10'] = pair_lower(MU, D, False, True, pair_mesh)
    print('S10 completed', flush=True)

    base11 = 8 * integrate1(lambda t: log((1 - 2 * t) / t) / (t * (1 - t)), D, Fraction(1, 3), line_mesh)
    # Only an UPPER error is required. Its lower endpoint is zero.
    correction11 = I(Fraction(32, 3) * (Fraction(1, 4) - D)**2 / D**2) * log(I(1 / D - 3))
    terms['S11'] = base11 + I(0, correction11.hi, raw=True)
    terms['S12'] = 8 * integrate1(lambda t: log(2 - 3 * t) / (t * (1 - t)), A, Fraction(1, 3), line_mesh)
    terms['S13'] = 8 * integrate1(lambda t: log((1 - t - D) / I(D)) / (t * (1 - t)), B, D, line_mesh)
    terms['S14'] = Fraction(16, 7) * integrate1(
        lambda t: log(I(B) / t)**2 * log(t / I(A)) / t**2, A, B, line_mesh)
    terms['S15'] = Fraction(32, 7) * integrate1(
        lambda t: log((Fraction(1, 2) - 2 * A - t) / I(B)) / t
        * (I(1 / A) - (log(t / I(A)) + 1) / t), A, B, line_mesh)
    value = (3 * terms['S1'] + terms['S2'] + terms['S8'] + terms['S9'] + terms['S10']
             - 2 * terms['S3'] - sum(terms[f'S{i}'] for i in [4, 5, 6, 7])
             - 2 * terms['S11'] - sum(terms[f'S{i}'] for i in [12, 13, 14, 15])) / 4
    passed = 100 * value.lo >= 131 * Q
    # Independent, exact-rational logarithm upper bound for the transfer.
    # log(20/17)=2 atanh(3/37). Keep the terms u,u^3/3 and bound
    # all remaining denominators below by 5 in the positive series.
    u = Fraction(3, 37)
    log_upper = 2 * (u + u**3 / 3 + u**5 / (5 * (1 - u**2)))
    transfer_lower = Fraction(131, 100) - Fraction(40000, 4999) * log_upper
    transfer_target = Fraction(1, 125)
    transfer_passed = transfer_lower > transfer_target
    result = {
        'status': 'EXACT_INTEGER_INTERVAL_CERTIFICATE_FOR_DISPLAYED_NUMERICAL_FUNCTIONAL',
        'scope': 'Does not independently establish the number-theoretic sieve/distribution interfaces.',
        'target': '131/100', 'target_verified_by_integer_comparison': passed,
        'fixed_point_bits': BITS, 'log_series_terms': LOG_TERMS,
        'meshes': {'pairs': pair_mesh, 'lines': line_mesh, 'J': j_mesh},
        'parameters': {k: str(v) for k, v in {'a': A, 'b': B, 'mu': MU, 'd': D}.items()},
        'log2_check': display(I(*LOG2, raw=True)),
        'J5': display(J(I(5))), 'S11_correction_upper': display(correction11),
        'terms': {k: display(v) for k, v in terms.items()},
        'coefficient': display(value), 'runtime_seconds_informational': time.monotonic() - start,
        'independent_rational_transfer_certificate': {
            'scope': 'Numerical implication only; number-theoretic hypotheses are separate.',
            'chen_input_constant': '131/100', 'b': '17/20',
            'theta': '4999/10000', 'vartheta': '9999/20000',
            'atanh_u': str(u), 'log_upper_exact': str(log_upper),
            'log_upper_formula': '2*(u+u^3/3+u^5/(5*(1-u^2)))',
            'transfer_lower_exact': str(transfer_lower),
            'target': str(transfer_target),
            'strict_target_verified_by_exact_rational_comparison': transfer_passed,
        },
    }
    probe_file = Path(__file__).with_name('probe_conservative_chen.json')
    if probe_file.exists():
        probe = json.loads(probe_file.read_text(encoding='utf-8'))
        latest = probe['runs'][-1]
        result['approximate_diagnostics_not_used_for_any_certificate_bound'] = {
            'probe_file': probe_file.name, 'quadrature_order': latest['order'],
            'all_single_original_functional': latest.get('coefficient_all_single_variable_levels'),
            'all_single_and_buchstab_cap': latest.get('coefficient_all_single_and_buchstab_4_over_7'),
            'certified_elementary_functional': latest.get('coefficient_elementary_certificate_functional'),
        }
    destination = Path(__file__).with_suffix('.json')
    destination.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps({'target_verified': passed, 'coefficient': display(value), 'path': str(destination)}), flush=True)
    return passed and transfer_passed


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pair-mesh', type=int, default=128)
    parser.add_argument('--line-mesh', type=int, default=1000)
    parser.add_argument('--j-mesh', type=int, default=10000)
    args = parser.parse_args()
    if not certify(args.pair_mesh, args.line_mesh, args.j_mesh):
        raise SystemExit(2)
