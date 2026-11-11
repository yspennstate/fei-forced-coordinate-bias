"""c = ratio - [(k-1)log2(1/P) - log2(1/eps) - (1/2)log2 m] for each row of eps_scan_out.txt (P = 1/2 + eps)."""
import math, sys
print(f"{'k':>4} {'eps':>10} {'eps*2^m':>10} {'ratio':>11} {'c':>8}")
for line in open(sys.argv[1] if len(sys.argv) > 1 else "eps_scan_out.txt", encoding="utf-8"):
    p = line.split()
    if len(p) != 5 or not p[0].isdigit():
        continue
    k, e, r = int(p[0]), float(p[1]), float(p[2])
    if e <= 0:
        continue
    m = (k - 1) // 2
    pred = (k - 1) * math.log2(1 / (0.5 + e)) - math.log2(1 / e) - 0.5 * math.log2(m)
    print(f"{k:>4} {e:>10.3e} {e * 2**m:>10.3g} {r:>11.6f} {r - pred:>8.4f}")
