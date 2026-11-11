"""First-order (in eps) plateau constant of E/Inf near P = 1/2.

For u+v = 2b+1:  r = K_b,  K_b = C(m,b) C(2m,m) / C(2m,2b)  (exact, P = 1/2 value).
For u+v = 2b:    r = 2 eps |(v-u-1)(-1)^b K_b + (-1)^v| + O(eps^2).
In the regime 2^-m << eps << 1 the second family has r -> infinity, so Psi(r) = log r^2 + log e + o(1), and
    c(m) = E/Inf - [2m - log2(1/eps) - (1/2) log2 m]
is independent of eps.  This script evaluates c(m) and the limit
    c_inf = log2(2/pi) + (3 - gamma)/2 * log2(e).
"""
import math, sys

LOG2E = 1 / math.log(2)


def lbin(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def psi_gap(x):
    """Psi(r) - log2 r^2 for r = e^x >= 1."""
    if x > 30:
        return LOG2E
    u = math.exp(-2 * x)
    return LOG2E * (math.log1p(u) / u + math.log1p(u))


def c_of_m(m):
    beta = [math.exp(lbin(m, u) - m * math.log(2)) for u in range(m + 1)]
    lK = [lbin(m, b) + lbin(2 * m, m) - lbin(2 * m, 2 * b) for b in range(m + 1)]
    tot = 0.0
    for u in range(m + 1):
        for v in range(m + 1):
            s = u + v
            wt = beta[u] * beta[v]
            if s % 2 == 1:
                x = lK[(s - 1) // 2]
                tot += wt * (2 * x * LOG2E + psi_gap(x))
            elif s >= 2:
                b = s // 2
                n = v - u - 1
                sb = 1 if b % 2 == 0 else -1
                sv = 1 if v % 2 == 0 else -1
                mag = abs(n * sb + sv * math.exp(-lK[b]))
                tot += wt * (2 * (lK[b] + math.log(mag)) * LOG2E + 2.0 + LOG2E)
            # s = 0: r_00 = O(eps) with no large factor, Psi -> 0
    return tot - (2 * m - 0.5 * math.log2(m))


c_inf = math.log2(2 / math.pi) + (3 - 0.5772156649015329) / 2 * LOG2E
print(f"c_inf = log2(2/pi) + (3-gamma)/2*log2(e) = {c_inf:.6f}")
print(f"{'m':>6} {'k':>6} {'c(m)':>10} {'(c-c_inf)*sqrt(m)':>18}")
for m in [int(a) for a in (sys.argv[1:] or "10 20 30 40 60 100 200 400".split())]:
    c = c_of_m(m)
    print(f"{m:>6} {2*m+1:>6} {c:>10.6f} {(c - c_inf) * math.sqrt(m):>18.6f}")
