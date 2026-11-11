"""Offset E/Inf - m + (1/2)log2 m of the P=1/2 single sum against c0 = (3/4)log2 e - (1/2)log2(pi/2); prints m*(offset-c0)."""
import sys
from math import log, log2, pi, e
sys.path.insert(0, ".")
from kraw import closed
c0 = 0.75 * log2(e) - 0.5 * log2(pi / 2)
print("c0 = %.10f" % c0)
print("   k      m     offset        offset-c0     m*(offset-c0)")
for k in (41, 81, 121, 301, 1001, 2001, 4001, 8001):
    m = (k - 1) // 2
    v = closed(m)
    off = v - m + 0.5 * log2(m)
    print("%5d %6d  %.10f  %.3e  %.5f" % (k, m, off, off - c0, m * (off - c0)))
