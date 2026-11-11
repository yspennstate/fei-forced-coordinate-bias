"""Same quantity as p9/code/forced_exact.py (ratio E_x/Inf_x of tilde f_k under bias P on y_1..y_m and 1-P on
y_{m+2}..y_k), with every quantity carried as an integer over a power of w (P = u/w), so no Fraction gcds.
With n = 2m-i-j and W = w^(n+1):  D*W = dn*w,  A*W = [i=j=0] W - 2 (Mn - (w-u) dn) - dn w,  dn = (-1)^j u^n,
Mn = sum_{c>=m+1} [z^c] (z-1)^(i+j) X^(m-i) Y^(m-j+1).  Weight of a coefficient pair relative to Inf_x = P^(2m):
   (u(w-u))^(i+j) C(m,i) C(m,j) w^(2m) / (w^(2(i+j)) u^(2m) W^2).  Big-int true division gives the floats.
Usage: python forced_fast.py kmin kmax t ...   (t = 0 means P = 1/2)
"""
import sys
from math import comb, log, log1p

LN2 = log(2)


def bias(t):
    if t == 0:
        return 1, 2
    return sum(comb(2 * t, j) for j in range(t, 2 * t + 1)), 4 ** t


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def l2_1p_ratio(num, den):
    """log2(1 + num/den) for nonnegative ints, den > 0."""
    r = num / den
    return log1p(r) / LN2


def ratio(k, u, w):
    m = (k - 1) // 2
    X = [w - u, u]
    Y = [u, w - u]
    Xp = [[1]]
    for _ in range(m):
        Xp.append(pmul(Xp[-1], X))
    Yp = [[1]]
    for _ in range(m + 1):
        Yp.append(pmul(Yp[-1], Y))
    Zp = [[1]]
    for _ in range(2 * m):
        Zp.append(pmul(Zp[-1], [-1, 1]))
    E = 0.0
    base_den = u ** (2 * m)
    w2m = w ** (2 * m)
    for i in range(m + 1):
        XY = None
        for j in range(m + 1):
            poly = pmul(pmul(Xp[m - i], Yp[m - j + 1]), Zp[i + j])
            Mn = sum(poly[m + 1:])
            n = 2 * m - i - j
            W = w ** (n + 1)
            dn = (-1) ** j * u ** n
            Dn = dn * w
            An = (W if i == j == 0 else 0) - 2 * (Mn - (w - u) * dn) - Dn
            a2, d2 = An * An, Dn * Dn
            wnum = (u * (w - u)) ** (i + j) * comb(m, i) * comb(m, j) * w2m
            wden = w ** (2 * (i + j)) * base_den * W * W
            if a2:
                E += (a2 * wnum / wden) * l2_1p_ratio(d2, a2)
            E += (d2 * wnum / wden) * l2_1p_ratio(a2, d2)
    return E


if __name__ == "__main__":
    kmin, kmax = int(sys.argv[1]), int(sys.argv[2])
    ts = [int(a) for a in sys.argv[3:]] or [2, 3]
    print("k    " + "".join(f"t={t:<10}" for t in ts))
    for k in range(kmin, kmax + 1, 2):
        print(f"{k:<4} " + "".join(f"{ratio(k, *bias(t)):<12.6f}" for t in ts))
        sys.stdout.flush()
