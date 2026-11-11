"""E_x/Inf_x of tilde f_k under bias P = 1/2 + eps, for eps from 0.2 down to 2^-m, at k = 2m+1."""
import math, sys
from forced_fast import ratio

c0 = 0.75*math.log2(math.e) - 0.5*math.log2(math.pi/2)
ms = [int(a) for a in sys.argv[1:]] or [10, 20, 30, 40]
print(f"{'k':>4} {'eps':>10} {'ratio':>12} {'ratio/((k-1)log2(1/P))':>24} {'ratio-(m-.5log2 m+c0)':>22}")
for m in ms:
    k = 2*m + 1
    w = 2 ** (m + 12)
    half_ref = m - 0.5*math.log2(m) + c0
    for e in [0.2, 0.1, 0.05, 0.02, 0.01, 1e-3, 1e-4, 1e-6, 2.0**(-m/2), 2.0**(-m), 2.0**(-m-6), 0.0]:
        d = int(round(e * w))
        if e > 0 and d == 0:
            continue
        u = w // 2 + d
        P = u / w
        r = ratio(k, u, w)
        print(f"{k:>4} {d/w:>10.3e} {r:>12.6f} {r/((k-1)*math.log2(1/P)):>24.6f} {r-half_ref:>22.6f}")
        sys.stdout.flush()
