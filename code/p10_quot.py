import sys
from math import log2
sys.path.insert(0, ".")
from forced_fast import ratio, bias
ks = [int(a) for a in sys.argv[1].split(",")]
ts = [int(a) for a in sys.argv[2].split(",")]
print("t  P          k    E/Inf        log2(1/Inf)  quotient")
for t in ts:
    u, w = bias(t)
    for k in ks:
        m = (k - 1) // 2
        r = ratio(k, u, w)
        L = 2 * m * (log2(w) - log2(u))
        print(f"{t}  {u}/{w:<8} {k:<4} {r:<12.6f} {L:<12.6f} {r/L:.6f}", flush=True)
