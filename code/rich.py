import math
c0 = 0.75*math.log2(math.e) - 0.5*math.log2(math.pi/2)
rows = []
for line in open("c0_check_out.txt"):
    f = line.split()
    if len(f) >= 3 and f[0].isdigit():
        k, m, off = int(f[0]), int(f[1]), float(f[2])
        rows.append((k, m, m*(off - c0)))
for (k1, m1, a1), (k2, m2, a2) in zip(rows, rows[1:]):
    lim = (m2*a2 - m1*a1)/(m2 - m1)
    print(f"k={k1}->{k2}  Richardson limit {lim:.6f}")
print(f"(5/16)log2 e = {5/16*math.log2(math.e):.6f}")
