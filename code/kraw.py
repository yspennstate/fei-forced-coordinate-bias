"""Closed form for E_x/Inf_x of tilde f_k (P = 1/2): sum over odd s of C(2m,s)/4^m Psi(r_s),
r_s = |[z^m](1-z)^a(1+z)^(2m-a)|, a = s-1 (even); conjectured r_s = C(a,a/2)C(2m-a,m-a/2)/C(m,a/2)."""
import sys
from math import comb, log1p, log
sys.path.insert(0, ".")
L2 = log(2)
def cm(a, m):
    return sum((-1) ** j * comb(a, j) * comb(2 * m - a, m - j) for j in range(min(a, m) + 1))
def psi_big(r):  # Psi(r) = r^2 log(1+1/r^2) + log(1+r^2), r integer >= 0
    if r == 0:
        return 0.0
    r2 = r * r
    if r2.bit_length() > 1000:
        return (1.0 / L2 + log(r2) / L2)
    return (r2 * log1p(1 / r2) + (log(r2) + log1p(1 / r2))) / L2
def closed(m):
    tot = 0.0
    for s in range(1, 2 * m + 1, 2):
        a = s - 1
        r = comb(a, a // 2) * comb(2 * m - a, m - a // 2) // comb(m, a // 2)
        tot += (comb(2 * m, s) / 4 ** m) * psi_big(r)
    return tot
if __name__ == "__main__":
    for m in range(1, 12):  # identity checks
        for a in range(0, 2 * m + 1):
            c = cm(a, m)
            if a % 2:
                assert c == 0, (m, a, c)
            else:
                b = a // 2
                num = comb(a, b) * comb(2 * m - a, m - b)
                assert num % comb(m, b) == 0 and abs(c) * comb(m, b) == num and c * (-1) ** b > 0, (m, a, c)
    print("identity ok m<=11")
    from forced_fast import ratio
    for k in (13, 41, 81):
        print(k, closed((k - 1) // 2), ratio(k, 1, 2))
    for k in (41, 81, 121, 141, 301, 1001, 4001):
        m = (k - 1) // 2
        c = closed(m)
        print(k, round(c, 6), "q=", round(c / (2 * m), 6), "offset=", round(c - m + 0.5 * log(m) / L2, 4))
