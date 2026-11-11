"""Sum of ln n^2 over odd n, |n|<N (N even), minus half of int_{-N}^{N} ln x^2 dx; tends to 2 ln 2."""
import math
for N in (100, 1000, 10000, 100000):
    s = sum(math.log(n*n) for n in range(-N+1, N, 2))  # odd n with |n|<N, N even
    I = 4*(N*math.log(N)-N)                          # int_{-N}^{N} ln x^2 dx
    print(N, s - I/2, 2*math.log(2))
