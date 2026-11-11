# The erasure entropy of a forced coordinate under every bias

Yitzchak Shmalo

This repository contains a short paper on the Fourier entropy-influence (FEI) conjecture of Friedgut and Kalai,
together with the programs behind every numerical value in it.

Paper: [main.pdf](main.pdf) (source [main.tex](main.tex)).

## Summary

For a Boolean function f and a coordinate x, the erasure entropy E_x(f) is the conditional entropy of the event
"x is in U" given U minus x, where U is drawn from the spectral distribution. It never exceeds the binary entropy of
the influence, so E_x/Inf_x is at most log2(e/Inf_x).

F_{k,t} is the monotone function obtained from MAJ_k(y) OR (x AND [y >= y0]), with the point y0 + e_{m+1} removed
from the majority (k = 2m+1), by replacing y_1, ..., y_m with "at least t of 2t" of fresh variables and
y_{m+1}, ..., y_k with "at least t+1 of 2t". The coordinate x is the only variable lying both in a shortest minterm
and in a shortest maxterm. With P_t = Pr[Bin(2t, 1/2) >= t], the paper proves, for every fixed t >= 1,

    E_x(F_{k,t})/Inf_x(F_{k,t}) = (1 + o(1)) log2(1/Inf_x) = (1 + o(1)) (k-1) log2(1/P_t).

More generally the same holds for the biased majority under every rational bias P in (1/2, 1); the phase of the
relevant coefficients advances by arccos((2P-1)/(2P)), an irrational multiple of pi. At P = 1/2 the coefficients are
Krawtchouk values at the middle of their range, half of them vanish, and

    E_x/Inf_x = m - (1/2) log2 m + c0 + (5/16) log2(e) / m + O(1/m^2),

with c0 = (3/4) log2 e - (1/2) log2(pi/2) = 0.75627...,

from an exact single-sum formula.

## Code

All programs are Python 3 with no other dependencies.

* `code/forced_fast.py` computes E_x/Inf_x from the double sum of Lemma 5 in integer arithmetic (every quantity an
  integer over a power of the denominator of P). `python forced_fast.py kmin kmax t ...`, with t = 0 for P = 1/2.
* `code/forced_exact.py` computes the same double sum in rational arithmetic; used as a check for small k.
* `code/p10_quot.py` prints the quotient of E_x/Inf_x by log2(1/Inf_x) for given k and t. Outputs:
  `code/quot_t2.txt` to `code/quot_t4.txt` (for t = 0 see `code/kraw_out.txt`).
* `code/kraw.py` checks the Krawtchouk identity of Lemma 6 for m <= 11, compares the single sum of Theorem 3 with the
  double sum, and evaluates it up to k = 4001. Output: `code/kraw_out.txt`.
* `code/c0_check.py` compares the single sum with m - (1/2) log2 m + c0 up to k = 8001. Output: `code/c0_check_out.txt`.
* `code/rich.py` and `code/resid2.py` read that output and check the 1/m and 1/m^2 terms. Outputs:
  `code/rich_out.txt`, `code/resid2_out.txt`.

## Related papers

1. [The one-coordinate step in the Fourier entropy-influence conjecture](https://github.com/yspennstate/fourier-entropy-influence-one-coordinate)
2. [Good coordinates in the Fourier entropy-influence conjecture](https://github.com/yspennstate/fei-good-coordinates)
3. [A variance form of Han's bound on the Fourier entropy](https://github.com/yspennstate/fei-variance-han-bound)
4. [Four equivalent forms of the Fourier entropy-influence conjecture](https://github.com/yspennstate/fei-equivalent-forms)
5. [A filtration bound for the erasure entropies of Boolean functions](https://github.com/yspennstate/fei-filtration-witness)
6. [Erasure entropies of sparse linear DNFs](https://github.com/yspennstate/fei-sparse-linear-dnf)
7. [The Fourier entropy of linear monotone DNFs on few variables](https://github.com/yspennstate/fei-linear-dnf-census)
8. [Hidden coordinates and selection rules in the Fourier entropy-influence conjecture](https://github.com/yspennstate/fei-hidden-coordinates)
9. [A forced coordinate with nearly maximal erasure entropy](https://github.com/yspennstate/fei-forced-coordinate)

## License

MIT, see [LICENSE](LICENSE).
