import math
c0 = 0.75*math.log2(math.e) - 0.5*math.log2(math.pi/2)
a1 = 5/16*math.log2(math.e)
for line in open("c0_check_out.txt"):
    f = line.split()
    if len(f) >= 3 and f[0].isdigit():
        k, m, off = int(f[0]), int(f[1]), float(f[2])
        print(f"k={k:5d} m={m:5d}  m^2*(offset - c0 - a1/m) = {m*m*(off - c0 - a1/m):.4f}")
