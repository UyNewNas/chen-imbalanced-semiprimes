SELECTED NUMERICAL VERIFICATION -- RENEWAL V2, 10 OCTOBER 2026

Frozen current selection
------------------------
K = 6.7617368 = 8452171/1250000
c = 1.96176252 = 49044063/25000000
b = 0.74816964 = 18704241/25000000, count coefficient 1/100000.
Secondary b = 0.74816865 = 14963373/20000000, coefficient 1/1000000.
Selected record: corpus-renewal-improvement-assembly-20261010.json.
The earlier b=0.74821909 result remains the explicitly identified baseline.
The selected refinement is threshold 7/3 v2; the earlier 12/5 version is replaced,
not added. No new mask cells or external distribution theorem are adopted.

Quick run
---------
Python 3.11 or newer; standard library only. Tested with Python 3.12.14 on Windows.

    python -B run_verification.py

The default uses fresh temporary directories and removes them afterward. To retain
all detailed results in a NEW directory:

    python -B run_verification.py --work verification-results

The command first validates the compact manifest, then replays the baseline's
22 upper saving terms, lower ledger, reserves and 106-action recurrence from the
original seed for 430 rounds. Historical arrays at rounds 248/312 are comparison
objects, never recurrence starting conditions. The final array SHA-256 is
4243db3238c28e8ac7f22d848db98e7872ad3c134fc4d697058583248f46f52d.

It then ACTUALLY recalculates the adopted v2 refinement on all 17,572 cells:
6,091 old E0 cells and 11,481 E1 cells, including 227 zero-crossing E1 cells.
The author engine replays the literal Q80 saved floors and directed Q160 gross
mass intervals. An independent decimal-integer/Mercator engine recomputes the
same mask and loss upper proxies. A third program checks every row, geometry,
mask, saved-floor domination, both reserves, the feedback slope and final scalar
inequalities. Their mathematical fields must match the frozen v2 certificates.
A separate 96-term rational atanh check also verifies both selected b values and
failure at the preceding eight-decimal value. No floating-point decision is made.

The v2 identity is additional gross 0.000558583 plus recovered saved-floor loss
0.000037209 = raw increment 0.000595792. Each of two reserves receives gain/8;
the improvement after both reserves is 0.000446844. Only additional gross/K is
subtracted from the preceding feedback slope. Old gains and v1 gains are not added.
The first staged full delta replay took 53.84 seconds (author 5.80 s, independent
engine 46.84 s). Hardware-dependent whole quick/full timings are in the results.

Scope and adapters
------------------
The baseline higher-dimensional quadratures and primitive enclosures remain
inherited inputs. The 430 recurrence and selected scalar ledgers are freshly run;
the NEW v2 local cell/floor calculations are freshly run as described above.
This is not a fresh reproof of retained analytic distribution or compiler modules,
and is not a proof-assistant certificate. Matching a hash or stored PASS label is
not itself a mathematical proof of an inherited endpoint.

replay_selected_ledger.py preserves the original baseline formulae but stops before
its full historical audit. replay_refinement_author.py preserves the complete v2
numerical body while omitting mutable research-state and historical-source guards.
replay_refinement_review.py uses the frozen baseline assembly in place of mutable
research-state and leaves historical hash auditing to the full wrapper. The original
independent-v2 program runs unchanged. All three mathematical outputs are compared
with their frozen records; metadata differences are explicitly excluded. The compact
manifest verifies supplied bytes. Original programs are retained for review.

Full public evidence archive
----------------------------
Extract chen-numerical-evidence-20261010.zip and run its
run_full_verification.py. It verifies every public file, the 1,689 baseline binding
comparisons and 143 selected refinement/source comparisons, then actually performs
the same baseline ledger/430 replay AND new v2 author/independent/row replay.
The 143 comparisons include the original author's 53 source-dependency comparisons.
It does not recursively treat literature-scan metadata as mathematical premises.
No third-party papers, research-state dump, private correspondence or attachment
transcripts are bundled. Older v1 records are provenance only, never another gain.

public_normalizations.json documents every nonmathematical privacy normalization,
with original/public hashes and exact JSON fields or text spans. Mathematical data,
formulae and executable calculation logic are unchanged. The public full adapter
accepts only manifested file/original-hash/public-hash triples and reports actual
exceptions separately. Original frozen programs cannot be described as having run
unchanged on privacy-normalized source bytes. Wider analytic proof modules are
provided separately in the publication's proof library.

The historical 17,772-cell typo is corrected in the adopted v2 record: the union
is 6,091 + 11,481 = 17,572. Historical immutable source records retain their spelling;
see the publication ERRATA. This does not alter masks or numerical endpoints.

中文说明
--------
本版本采用7/3门槛的v2改进。快速命令会真正重算17572格的新增质量和保存下界，
再用独立区间算法和逐行检查交叉核验；此前的大型积分仍为继承输入。原430轮
递推和基线账本继续单独复现。主参数b=0.74816964，次参数b=0.74816865；v1与v2
收益不会重复累计。所有输出在新目录生成，不修改冻结输入。
