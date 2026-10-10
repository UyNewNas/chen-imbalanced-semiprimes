# Imbalanced semiprime representations of even integers

**A revised modular proof with paired sieve weights and finite certificates**  
Slava Chan · 10 October 2026

[Read the manuscript](paper/SlavaChan-imbalanced-semiprimes-20261010.pdf) · [One-page reading guide](paper/SlavaChan-reading-guide-20261010.pdf) · [Proof-source index](proof-library/INDEX.html) · [Verification instructions](verification/README.txt) · [Full numerical evidence](https://github.com/UyNewNas/chen-imbalanced-semiprimes/releases/tag/v2026.10.10)

This repository publishes the revised manuscript, its source and supporting records for open mathematical review. Corrections and reproducibility reports are welcome. Please identify the version and the exact page, equation or file involved.

## Current statement

For every sufficiently large even integer $N$, the manuscript claims at least $H_N/100000$ representations

$$
N=p+rq,\qquad p,q\text{ prime},\quad r=1\text{ or prime},\quad r\le q^{0.74816964}.
$$

Here $H_N=C_NN/(\log N)^2$, where

$$
C_N=\prod_{\ell>2}\left(1-\frac{1}{(\ell-1)^2}\right)
\prod_{\substack{\ell\mid N\\\ell>2}}\frac{\ell-1}{\ell-2},
$$

and both products are over primes. Distinct prime summands $p$ are counted once. There is no extra factor of two in this normalization.

For each fixed $1/2<b<1$, the general claimed bound is

$$
\liminf_{\substack{N\to\infty\\2\mid N}}\frac{D_b(N)}{H_N}
\ge 1.96176252-6.7617368\log(1/b).
$$

| Choice of $b$ | Eventual lower count |
| --- | --- |
| **0.74816964** | $H_N/100000$ |
| 0.74816865 | $H_N/1000000$ |
| 0.9 | $1.24889H_N$ |

No effective threshold for $N$ is supplied. These are the manuscript's mathematical claims for review. Finite numerical checks do not establish the correctness of the full analytic proof. Complete machine formalization and independent human peer review are not claimed.

## What this revision changes

The editorial update adds Yuta Suzuki's dated sieve notes as reference [10], gives a citation key for the appendix references, and links the repository and release from the reproduction section. The manuscript remains 82 pages; its mathematical statements and frozen proof modules are unchanged. A [clean LaTeX upload package](paper/SlavaChan-arxiv-source-20261010.zip) and [suggested arXiv Comments text](paper/arxiv-comments.txt) are supplied from the same source.

The selected constants are exactly $K=8452171/1250000=6.7617368$ and $c=49044063/25000000=1.96176252$. This snapshot uses the adopted same-mask v2 refinement, frozen in `corpus-renewal-improvement-assembly-20261010.json`. Its threshold is $7/3$; it replaces the earlier $12/5$ candidate and the gains from those alternatives are not added.

The main construction now includes a restricted moving-residue distribution argument for the finite coefficient families actually used. It retains the oscillatory phase and accounts for its loss. The lower count uses the full eleven-term weight, its pointwise proof, and paired deletion of the relevant positive and negative domains. The upper estimate uses the same positive majorant for its density and progression estimate.

The newly adopted refinement keeps the same 17,572 paired cells, the same $K$, and the same finite tables. It proves $\omega(s)\ge 11/20$ for $s\ge7/3$, adds a certified $0.000558583$ of gross mass, and replaces the earlier paid saved-estimator loss by $0.002334885$. After two additional reserves, the admitted lower coefficient rises to $1.96176252$. Appendices K–L give this refinement; the older constants in the preserved baseline modules are superseded by it.

The old branch that used the external coefficient **1.9728** is no longer an input or a current corollary. This revision does not restore an unrestricted distribution theorem or claim that all results of any cited paper are invalid. Historical supporting records are retained only with explicit supersession notes.

Five source passages print `17772` for the pair-mask union; the correct sum is `6091 + 11481 = 17572`. [ERRATA.txt](proof-library/ERRATA.txt) identifies each passage. The selected cell lists and numerical endpoints are unchanged.

## Read the proof

The new **82-page manuscript** contains a typeset main argument, references, and twelve complete consolidation and refinement modules in their original ASCII mathematical notation. It replaces the former 308-page assembled draft as the current reading entry point.

1. Start with the one-page guide and Section 1 for the result and normalization.
2. Sections 2–5 give the restricted distribution argument, upper sieve and lower-weight construction. The analytic details are in Appendices A–H.
3. Section 6 and Appendix I specify the selected finite inequalities and reproduction scope.
4. Section 7 and Appendix J give the ordered bad-semiprime mass and exact final transfer. Section 8 describes comparison limits and substantive AI assistance.

The [proof-source library](proof-library/README.txt) contains **151 texts** with an index and manifest. Its index distinguishes the current analytic framework, the adopted v2 refinement, and superseded numerical candidates. It also preserves referenced historical derivations, including an old manuscript whose constants are superseded. This library is a reading dependency collection, not a proof-assistant dependency certificate.

## Reproduce the finite checks

Use Python **3.11 or later**. The two release runners use only the standard library. From the repository root:

```bash
python verification/run_verification.py --work-dir replay-quick
```

The output directory must be new. Omitting `--work-dir` uses and removes a fresh temporary directory.

The compact run validates its manifest, reconstructs the baseline upper/lower ledgers from 49 numerical inputs, checks the reserves, and replays 106 actions for 430 simultaneous rounds from the seed. It also regenerates the adopted refinement on all 17,572 paired cells with two interval implementations, compares the saved-floor estimates, and checks the new budgets and logarithm margins. All 400 final entries are compared with the frozen result. A successful run reports `PASS_COMPACT_VERIFICATION`.

For the full evidence audit, download and extract [the numerical evidence archive](https://github.com/UyNewNas/chen-imbalanced-semiprimes/releases/download/v2026.10.10/chen-numerical-evidence-20261010.zip). From its extracted directory:

```bash
python run_full_verification.py --work-dir full-replay
```

The full runner checks the public manifest, audits the inherited file bindings, and repeats the selected baseline and refinement calculations. A successful run reports `PASS_PUBLIC_FULL_BINDINGS_LEDGER_AND_430_REPLAY`.

Both documented runs were repeated successfully for this release; see [verification-results.json](verification/verification-results.json). The full archive contains the 1,129 baseline numerical binding files, the entry records, and the additional adopted refinement evidence. **Neither command regenerates every historical quadrature or proves the analytic arguments.** Saved interval endpoints remain inherited inputs. The archive includes their generators and witnesses for further review; some historical generators require additional packages and should be run in a fresh copy.

The final primary count has a strict rational surplus of more than $1.1497\times10^{-8}$ after subtracting $10^{-5}$; the secondary surplus exceeds $6.4162\times10^{-8}$. These margins are checked with rational arithmetic and an explicit logarithm-series remainder.

## Files and provenance

| Location | Contents |
| --- | --- |
| `paper/` | Current manuscript and one-page reading guide. |
| `paper/source/` | Complete LaTeX sources and the twelve included proof modules. |
| `proof-library/` | 151 supporting proof texts, index, bibliography map, manifest and errata. |
| `verification/` | Compact runnable package, full-archive runner, results and manifests. |
| Release assets | Full selected numerical evidence archive and SHA-256 checksums. |

Public copies omit private correspondence and third-party paper files. Eleven source-library texts have minimal edits to remove private-communication attribution or local paths. A small set of files in the full numerical archive has documented nonmathematical normalizations, including local paths; its manifest lists the exact affected files. Original and public hashes are recorded; mathematical coefficients, formulas and certificate arrays are unchanged. The full replay accepts only the listed original-hash/public-hash substitutions and reports them separately. The unchanged historical runner is supplied as a reference and is not expected to accept the normalized files without the adapter.

The repository is maintained as a single public snapshot commit under **Slava Chan**. The release tag and checksum manifests identify this edition.

## Rebuild the PDFs

Use a pdfLaTeX installation with the packages declared in the two source files:

```bash
cd paper/source
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error reading-guide.tex
```

The outputs are `main.pdf` and `reading-guide.pdf`. The guide's page references correspond to the supplied edition. Some `.tex` files in `proof-library/source/` are historical fragments rather than standalone documents; they are not required to compile `paper/source/main.tex`.

## Review and authorship

Issues and proposed corrections are welcome. Explain the specific step at issue and, for computational reports, provide the command, version and output. Reports of a checked step should state the scope and method of the check.

OpenAI ChatGPT/Codex provided substantive assistance with mathematical exploration, candidate proofs, parameter searches, checking, program development and writing. Separate implementations and AI-assisted internal reviews are not independent human peer review. The human author remains responsible for the work. Section 8 gives the full disclosure.

An explicit license for redistribution and adaptation has not yet been selected. This release preserves the supplied materials and their existing notices.

## 中文说明

本版正式采用 $b=0.74816964$，对应每个充分大偶数至少 $H_N/100000$ 个表示；另列 $b=0.74816865$、计数 $H_N/1000000$ 的次要推论。公开署名统一为 **Slava Chan**。

新版包含 82 页英文手稿、1 页阅读指南、151 份证明文本、可直接运行的小验证包，以及 Release 中的完整数值证据包。旧版依赖外部 $1.9728$ 的分支不再作为当前结论。原始文本中的掩码数量笔误另列勘误。

本次已重跑基线精确账本、430 轮递推、新增同掩码改进的两套区间计算及文件绑定核验。没有重新生成全部历史积分，也没有完成形式化验证或独立的人类同行审查。欢迎针对具体证明步骤、解析条件和数值证书提出问题。
