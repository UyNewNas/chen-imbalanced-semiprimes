# Goldbach 1+1.75: imbalanced semiprime representations of even integers

Slava Chan · Frozen manuscript edition: 6 November 2026

The current edition is **Theorem (1+1.75) on the Goldbach Conjecture**, with an English manuscript of **68 pages** and a Chinese manuscript of **57 pages**.

- [English manuscript](<paper/release-20261106/Theorem (1+1.75) on the Goldbach Conjecture - EN - 20261106.pdf>)
- [Chinese manuscript](<paper/release-20261106/Theorem (1+1.75) on the Goldbach Conjecture - ZH - 20261106.pdf>)
- [Release files, source build instructions and numerical supplement](paper/release-20261106/README.md)
- [Numerical supplement on Zenodo](https://doi.org/10.5281/zenodo.23755453)
- [Lean 4 formalization project (in progress)](https://github.com/UyNewNas/goldbach-formalization)

## Main result

The manuscript proves that every sufficiently large even integer $N$ has a representation

$$
N=p+rq,\qquad p,q\text{ prime},\qquad r=1\text{ or prime},\qquad r\le q^{3/4}.
$$

For each fixed $1/2<b<1$, its quantitative bound is

$$
\liminf_{\substack{N\to\infty\\2\mid N}}\frac{D_b(N)}{H_N}
\ge c-K\log(1/b),\qquad
c=1.96176252,\quad K=6.7617368.
$$

Here $D_b(N)$ counts distinct admissible prime summands, and $H_N=C_NN/\log^2N$ uses the singular-series normalization defined in the manuscript. The threshold for $N$ is ineffective.

| Factor exponent $b$ | Eventual lower count |
| --- | --- |
| $3/4$ | $0.0165H_N$ |
| $0.74816964$ | $10^{-5}H_N$ |
| $0.74816865$ | $10^{-6}H_N$ |
| $0.9$ | $1.24889H_N$ |

## Proof and numerical materials

The proof subtracts a positive upper-sieve estimate for excluded semiprimes from a lower Chen count obtained by paired restrictions and the double sieve. The distribution argument is established for the actual upper and lower coefficient families used in the construction.

The release directory supplies the complete bilingual TeX sources, including the analytic, coefficient, sieve and numerical appendices. Appendix I defines the finite quantities and the correspondence with the supplied numerical certificates.

The numerical supplement contains rational endpoints, finite tables, source certificates, algorithm sources and numerical witnesses. Its README lists seven verification entry points, with the retained inputs distinguished for each. The release checks cover 727 manifest-bound files, 1,190,700 ordered-cover rows and 632,104 sixth-rank cube cells, as well as the finite recurrence and stated polyhedral-domain checks. The original checksum manifest binds the two PDFs, source archive and numerical archive to one release.

These finite checks have the scope stated in Appendix I. A completed formal proof of the entire manuscript and independent human peer review are not claimed. The companion Lean project has its own development status.

## AI contribution

The author directed GPT-6 Astra Ultra to investigate prime-related techniques in a collection of approximately 700 proofs and adapt promising routes to improve the earlier $(1+1.9)$ factor restriction. The contribution is declared **ai-discovered**: GPT-6 Astra Ultra carried out the mathematical discovery and substantive proof under the author's research direction. Models used subsequently assisted with manuscript polishing, editing, organization and cross-validation. The author is responsible for the mathematical content and final manuscript.

## Earlier material and review

The previous repository reading guide is preserved as [README-20261010.md](README-20261010.md). Its 82-page manuscript, proof library and older verification instructions remain available at their historical paths; they are not the current reading entry point.

Corrections and reproducibility reports are welcome. Please identify the manuscript edition and the exact page, equation, source file or verification command involved.
