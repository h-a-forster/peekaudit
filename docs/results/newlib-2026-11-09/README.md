# Three more engines: tea-tasting, expectation, gsdesign (2026-11-09)

A study outside the full study's grid. peekaudit's own simulators (`simulate.generate`,
same DGPs, seeds and arm counts as the study), `stats.clopper_pearson`, the one-sided exact
binomial test and Holm were called against three more libraries, without mirrors: every
number is a real library call. REPORT.md, "Three more engines", has the findings; this
folder holds the scripts and the results they were drawn from.

| Library | Version | Claim | Result |
|---|---|---|---|
| tea-tasting | 2.0.0 | Fixed-horizon tests, no sequential claim | The default `Proportion` confidence interval misses its level at any sample size when the rates differ: 0.89 (absolute) and 0.78 (relative) at a nominal 0.95 |
| expectation | 0.6.1 | Anytime-valid sequential tests and confidence sequences | 9 of 51 cells violate after Holm, all already reported in the maintainer's #154 (default mean test 27.3% [23.8, 30.9] FPR at a nominal 5%). New: the proportion test raises while all observations are 0 or all 1, and the README's first import fails |
| gsdesign | 0.1.5 | Group sequential densities (numerical core) | Clean: per-look crossing probabilities within 8e-7 of `peekaudit.gst` for K up to 50 |

## tea-tasting 2.0.0

`Mean` (Welch t, z, Student t) and `Proportion` (`method="norm"`, pooled and unpooled)
through `analyze_aggregates`: 5 metrics x N 1,000 and 20,000 x allocation 50/50 and 20/80,
A/A and A/B, one look at N_max, R = 4,000; 14 daily looks with R = 1,000. 184 cells, 0
errors, 0 NaN (`tea-tasting/report.txt`).

- No sequential claim, so peeking inflation is informational: 20.6-22.9% on normal data with
  14 looks, as the exact 22.0% predicts.
- A/A, one look: 4 of 76 cells violate after Holm, all `Mean` Welch and z with 20% to B at
  N = 1,000 (1% conversions: 14.2% [13.2, 15.3]). The same small-arm behaviour peekaudit
  found in spotify-confidence and its own reference z-test: a property of the method.
- A/B coverage: 23 of 152 cells below 95% after Holm. `student_t` with unequal arms is what
  `equal_var=True` asks for; Welch and z on small skewed arms fade with N.
- **Defect: the default `Proportion` interval.** `equal_var=True` (the default) pools the
  variance for the confidence interval as well as for the p-value, so the miss does not
  shrink with n (`repro_tea.py`, 4,000 replicates, seed 123): 0.10 -> 0.30 at 8,000/2,000,
  absolute interval covers 0.888 (large-sample limit 0.885), 0.891 at 80,000/20,000; 0.30
  -> 0.10, relative interval 0.782. An independent 18-cell run (`pk_tea_prop_ci.py`,
  `prop_ci.out`) gives the same, while every `equal_var=False` control lies in
  0.948-0.957. Filed as tea-tasting#327.

## expectation 0.6.1

`SequentialTesting` mean (default and known variance) and proportion, `ParallelSequentialTest`,
`KSampleSequentialTest` A/A and both confidence sequences: 51 cells, Holm over 51
(`expectation/report.txt`). Each `SequentialTesting` cell had a 45 s budget, so R varies by
cell; the per-cell R is in `cells.jsonl`.

- 9 of 51 cells violate, all duplicates of the maintainer's open #154: the default mean
  test with empirical variance, one observation at a time, has 27.3% FPR at a nominal 5% on
  N(0, 1) (known variance: 2.1%); `ParallelSequentialTest` with empirical variance, SD 3:
  28.6%; `EmpiricalBernsteinConfidenceSequence` with a range other than width 1 misses in
  300 of 300 runs.
- **New:** `update()` raises `ValueError` while all observations so far are 0 or all are 1
  (every one-observation first update failed), and the README's
  `from expectation.seqtest import SequentialTesting` raises `ImportError`
  (`repro_exp.py`, `repro_exp.out`). Filed together as expectation#155.

The licence is GPL-3.0 plus `LicenseRef-AI-Training-Prohibited`. Running it for an audit is
not training, and no expectation code is copied here.

## gsdesign 0.1.5

The whole public API (`gridpts`, `h1`, `hupdate`) against `peekaudit.gst`, exact normal
values and one Monte Carlo (`pk_gs.py`, deterministic, `gsdesign/out.txt`). It computes
densities and makes no peeking claim, so the FPR battery does not apply. Per-look crossing
probabilities are within 8e-7 of `peekaudit.gst` for every K up to 50 with constant,
O'Brien-Fleming and Pocock bounds; K = 2, 5, 10 at z = 1.96 give 0.0831, 0.1417, 0.1933
(published 0.083, 0.142, 0.193). Nothing to report.

## Files

| Path | Contents |
|---|---|
| `scripts/pk_tea.py`, `pk_exp.py`, `pk_exp2.py` | Study drivers; one JSON line per finished cell, skipped on rerun |
| `scripts/tea_report.py`, `exp_report.py` | Holm and the tables in `tea-tasting/report.txt` and `expectation/report.txt` |
| `scripts/repro_tea.py`, `pk_tea_prop_ci.py`, `repro_exp.py` | Standalone checks of the two defects; outputs beside the results |
| `scripts/pk_gs.py` | gsdesign against `peekaudit.gst` |
| `scripts/mkvenv.sh` | A fresh wheels-only venv per library under the repository's `tmp/newlib/` |

## Reproducing

A fixed-size subset of these checks (tea-tasting A/A FPR and the `Proportion` interval;
expectation's mean and proportion tests and the README import) runs weekly against each new
release in the living audit (`peekaudit.living`, `scripts/living.py`, AGENTS.md). The scripts
below are kept as they were run.

```sh
sh docs/results/newlib-2026-11-09/scripts/mkvenv.sh tea-tasting 3.12 tea-tasting==2.0.0
tmp/newlib/tea-tasting/venv/Scripts/python -I docs/results/newlib-2026-11-09/scripts/pk_tea.py
uv run python -I docs/results/newlib-2026-11-09/scripts/tea_report.py docs/results/newlib-2026-11-09/tea-tasting/cells.jsonl
uv run python -I docs/results/newlib-2026-11-09/scripts/exp_report.py
```

The drivers write into this folder and skip cells already recorded; move the committed
`cells*.jsonl` aside to rerun from scratch. The report scripts need only this repository and
reproduce the committed `report.txt` files from the committed cells. The study process died
three times under shared-machine load (exit 127, exit 139, and a numpy `TypeError` that did
not recur with the same seed) and resumed from its checkpoint each time; the causes were not
established.
