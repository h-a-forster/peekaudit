# peekaudit report

Do the sequential tests in two open-source A/B testing engines keep their false-positive
guarantees under repeated looks? This report summarises the full-size study. The generated
reports hold every cell:

- full study (headline): [docs/results/full/REPORT.md](docs/results/full/REPORT.md),
  per-cell counts in [docs/results/full/summary.json](docs/results/full/summary.json);
- quick study (cross-check): [docs/results/quick/REPORT.md](docs/results/quick/REPORT.md);
- quick vs full comparison:
  [docs/results/full/quick_vs_full.md](docs/results/full/quick_vs_full.md);
- mirror vs engine on randomised inputs: [docs/results/verify.json](docs/results/verify.json).

All rates are percentages with exact 95% Clopper-Pearson intervals. "Violates" means the
one-sided exact binomial test of FPR <= 5% (or coverage >= 95%) rejects after Holm
adjustment at family-wise level 0.05 within the method's family (DESIGN 6.2). Calibration
findings describe statistical behaviour on the stated inputs; they are not claims of
software defects. Engine behaviours that may be worth raising with the maintainers are
listed under [Engine behaviours](#engine-behaviours), each with a standalone reproduction.

## What was run

| | full (headline) | quick (cross-check) |
|---|---|---|
| replicates per A/A unit / A/B unit | 20,000 / 10,000 | 4,000 / 2,000 |
| units x methods x schedules | 52 x 11 x 4 (2,224 cells) | same |
| master seed, chunk size | 20261023, 250 | same |
| config hash | `dde79da82e86a243` | `82e303ecf6547d58` |
| code | git `e4e8dcc` (all 52 units) | git `fa12471` (all 52 units) |
| wall time / CPU time | 53.2 min / 211.8 CPU-min, 4 workers | 50.3 min / 200.6 CPU-min, 4 workers |
| failed / unverified units | 0 / 0 | 0 / 0 |
| engine checks in the loop | first 20 replicates per unit and schedule (spotify GST on 140 looks: 4) | 10 (2) |
| sanity gates (DESIGN 5.3) | 40 of 40 pass | 40 of 40 pass |
| sensitivity (both arms >= 150) | 28 A/A units, 5.7 min, 3 workers | 28 units, 1.6 min |

Versions (both studies): Python 3.11.9, gbstats 0.8.0, spotify-confidence 4.1.0, numpy
1.26.4, pandas 1.5.3, scipy 1.17.1, statsmodels 0.15.0. The code changes between `fa12471`
and `e4e8dcc` touch reporting, provenance and validation only; the quick-vs-full check
below confirms that the simulated outcomes agree.

Metrics: normal, Bernoulli 10% and 1%, lognormal revenue (sigma 1.5), Pareto revenue with
tail index 2.5 and 1.5 (the latter has infinite variance and forms a separate stress
family), and a ratio metric. Planned sizes N_max = 1,000 and 20,000; 50% or 20% of traffic
to B. Schedules: one look at N_max (`fixed1`), 14 daily looks, 140 looks (`stream140`), and
21 looks running to 1.5 x N_max (`overrun21`). Full definitions: [docs/DESIGN.md](docs/DESIGN.md).

**Detection limit.** With 20,000 A/A replicates the smallest FPR a cell can flag at the
strictest Holm step is 5.43-5.52% for the sequential methods (quick: about 6.2%).
"Consistent" therefore means "no excess above about 5.5% detected", not "calibrated".

## Headline results (full study)

The most important pattern: many of the largest excesses arise at the first looks of small
experiments, where one arm holds a handful of users. With N_max = 1,000 and 20% of traffic
to B, the first of 140 looks has 6 users in A and 1 in B. A zero-variance arm gives a
standard error of 0, or a z-statistic rests on a few observations, and the engines report
significance on 1 to 30 users in B. The engine libraries apply no minimum sample size
(GrowthBook's product layer hides results below 150 units; that layer is not audited). The
"guarded" rates below count only rejections at looks where both arms hold at least 150
units (same data and analysis; a sensitivity analysis, not a model of any product). Rare
conversions are a separate problem: with 1% conversions the arms are large but the events
are few, and the guard changes nothing.

### gbstats 0.8.0, sequential test (asymptotic confidence sequence)

The vendor documentation says results can be checked "as many times as you like"; it states
no minimum sample size or tail condition.

- **Relative effect, default tuning N* = 5,000:** violates in 14 of 96 primary cells and 6 of
  16 stress cells. It holds on normal data, 10% conversions and the ratio metric. The
  violations are lognormal revenue with 20% to B (4 cells), 1% conversions (9 cells) and one
  Pareto tail-2.5 cell (5.6% [5.3, 5.9], N_max = 20,000, 20% to B, 140 looks).
  - Lognormal, N_max = 1,000, 20% to B, 140 looks: 16.7% [16.2, 17.2]; 93% of the false
    positives first occur with fewer than 150 users in an arm; guarded 3.7% [3.5, 4.0].
  - Lognormal, N_max = 20,000, 20% to B, 140 looks: 11.9% [11.5, 12.4]; guarded still
    9.0% [8.6, 9.4]. Here the arms are large and the excess comes from skew.
  - 1% conversions: from 6.3% [6.0, 6.7] (N_max = 20,000, 50/50, daily looks) to 14.1%
    [13.7, 14.6] (same, 140 looks). The guard changes none of them.
  - Guarded, 10 of 96 primary cells still violate (the nine 1% conversion cells and
    lognormal at N_max = 20,000); stress family guarded: 3 of 16.
  - Stress (Pareto tail 1.5): up to 14.3% [13.8, 14.8] (N_max = 20,000, 20% to B, 140 looks).
- **Absolute effect:** no violation in the primary family (0 of 96). One stress cell
  violates: 5.7% [5.3, 6.0] (Pareto tail 1.5, N_max = 20,000, 20% to B, 140 looks; guarded
  5.5% [5.2, 5.9]). The relative-scale excesses on 1% conversions therefore come with the
  relative (delta-method) interval, not with the confidence-sequence boundary as such.
- **N* tuned to the planned size:** worse on small experiments: 22 of 96 primary cells, up
  to 37.2% [36.6, 37.9] (lognormal, N_max = 1,000, 20% to B, 140 looks), where 98% of false
  positives start with fewer than 150 users in an arm; guarded 5.7% [5.3, 6.0], still a
  violation. Guarded, 9 of 96 primary cells violate (7 with 1% conversions, 2 lognormal at
  N_max = 1,000 with 20% to B). Stress family: 8 of 16 (guarded 5). One coverage violation:
  94.0% [93.5, 94.5] (lognormal, N_max = 1,000, 20% to B, 140 looks).
- Coverage at stopping is consistent with 95% in all 96 claimed A/B cells for the default
  relative and the absolute variants.

### spotify-confidence 4.1.0, group sequential test (Lan-DeMets power family, rho = 2)

The claim covers schedules within the planned sample size (`fixed1`, `daily14`, `stream140`).

- Violates in 21 of 60 primary cells and 6 of 12 stress cells, all with 20% of traffic to B
  except one: normal data, N_max = 1,000, 50/50, 140 looks, 5.9% [5.5, 6.2] (guarded 4.9%
  [4.6, 5.2], not flagged).
- Three primary violations are single-look failures, not sequential ones. At one look the
  GST is a single z-test at N_max with 200 users in B, which over-rejects on 1% conversions
  (13.3% [12.9, 13.8]), lognormal revenue (7.1% [6.7, 7.4]) and Pareto tail-2.5 revenue
  (6.1% [5.8, 6.4]), exactly like the fixed-horizon z-test. For 1% conversions, B has zero
  conversions with probability 0.99^200 = 13.4%, which accounts for essentially the whole
  rate ([scripts/zero_count_example.py](scripts/zero_count_example.py)).
- The largest sequential excess is 20.5% [19.9, 21.1] (1% conversions, N_max = 1,000, 20% to
  B, daily looks; guarded 16.8% [16.3, 17.3]). On 10% conversions at N_max = 1,000 with 140
  looks it is 20.4% [19.8, 20.9]: 44% of the false positives occur at the very first look (6
  users in A, 1 in B) and 90% at looks with fewer than 150 in an arm; guarded 5.1% [4.8, 5.4],
  not flagged.
- Guarded, 14 of 60 primary cells still violate (3 of them single-look); stress 6 of 12.
- Coverage at stopping violates in 8 of 60 cells, worst 81.7% [80.9, 82.4] (10% conversions,
  N_max = 1,000, 20% to B, 140 looks).
- The independently implemented reference GST (Lan-DeMets O'Brien-Fleming) shows the same
  pattern: 17 of 72 primary cells (17 guarded), 6 of 12 stress. This is how a z-based group
  sequential test behaves on small, skewed or rare-event arms, not an implementation defect.
- Past the planned sample size (informational; the library makes no claim there) it reaches
  7.3-8.2% on normal data, for example 8.2% [7.8, 8.5] at N_max = 1,000, 20% to B.

### Controls and contrasts

- **Peeking at fixed-horizon tests** (informational): the naive z-test and spotify's z-test
  reach 21.4-23.7% with 14 daily looks on normal data and 39.4-40.4% with 140 looks at
  N_max = 20,000, matching the exact repeated-significance values 22.0% and 39.9%. At
  N_max = 1,000, 140 looks give 46%, because the early looks have very few users. The gbstats
  t-test behaves the same (21.2-21.7%, 37.3-39.3%).
- **gbstats Bayesian "chance to win >= 95%"** is exactly a two-sided level-0.10 z-test on the
  relative scale (single-look FPR 9.8%, 99.9% interval [9.1, 10.5]) and reaches 37-39% with
  daily looks and 61-62% with 140 looks on normal data.
- **Reference anytime-valid methods:** mSPRT violates in 11 of 96 primary cells (guarded 4,
  all 1% conversions); the asymptotic CS with N* = N_max in 3 of 96 (guarded 1: 8.5%
  [8.1, 8.9], 1% conversions, N_max = 1,000, 20% to B, overrun).
- All 40 sanity gates pass: positive controls against exact crossing probabilities,
  negative controls, fixed-horizon and power calibration, in-loop engine agreement and two
  mirror self-consistency checks (code paths that share arithmetic, so not independent
  evidence of engine correctness).

### Engine mirror verification

Engine rates come from vectorised mirrors of the engines' arithmetic. In the full run the
real engines were called on 1,126,400 replicate-looks inside the loop (first 20 replicates
of every unit and schedule): 0 decision, tie, no-result or interval mismatches. The
randomised grid (`peekaudit verify`, including empty and single-unit arms and extreme
variances) also has 0 mismatches.

## Cross-check: quick vs full

The quick study's replicates are the first 4,000 (A/A) and 2,000 (A/B) replicates of the
full study: both use the same master seed and chunk size, and the per-replicate caches are
identical over the shared replicates (checked on `normal-n1000-b50-aa`, all 220 arrays). The
two studies are therefore nested, not independent. `peekaudit compare` compares each quick
rate with the full study's remaining 16,000 (8,000) replicates, which are independent of it,
using a two-sided Fisher exact test per cell and Holm across cells:

- **FPR / power:** 2,224 cells compared, 0 differ after Holm (131 raw p < 0.05, in line with
  chance; smallest Holm-adjusted p = 1).
- **Coverage at stopping:** 2,224 cells, 0 differ after Holm (104 raw p < 0.05).
- **Verdicts:** all 132 quick violations are violations in the full study. 43 cells that were
  consistent in quick are violations in full (20 FPR, 23 coverage). Each full estimate lies
  within the quick 95% interval, the largest FPR among them is 6.7% and the lowest coverage
  92.7%: these are excesses between the two detection limits (about 6.2% and 5.5%), found
  because the full study has five times the replicates. No verdict moved the other way.

The full study agrees with the quick study within Monte Carlo error and sharpens it. For
example, gbstats sequential, lognormal, N_max = 1,000, 20% to B, 140 looks: quick 15.8%
[14.7, 17.0], full 16.7% [16.2, 17.2].

## Engine behaviours

Each item has a standalone reproduction that imports only the engine (plus numpy, scipy and
pandas), so it cannot come from peekaudit's adapters. Draft issue texts are in
[docs/upstream/](docs/upstream/README.md); none has been filed.

1. **spotify-confidence 4.1.0, exception:** `ZTest.difference` raises `IndexError` when the
   control group has a single unit with value 0
   ([scripts/spotify_indexerror_repro.py](scripts/spotify_indexerror_repro.py)). Degenerate
   input, outside the study grid.
2. **spotify-confidence 4.1.0, bounds after an overrun:** once the cumulative sample in a
   sequential call exceeds the planned size, every earlier bound is recomputed as if the
   plan had been the largest sample in the call. With a plan of 800 and looks every 100 up
   to 1,200, look 8 has adjusted alpha 0.0311 when analysed at the time and 0.0128 in a
   re-analysis with all 12 looks; look 10: 0.0297 vs 0.0198
   ([scripts/spotify_overrun_repro.py](scripts/spotify_overrun_repro.py), deterministic; the
   re-analysis values equal an independent solver's bounds at fractions k/12). Past results
   are not reproducible after an overrun.
3. **gbstats 0.8.0, sequential FPR with few events:** relative effect, 1% conversions,
   20,000 users 50/50, 140 looks: 15.2% [13.0, 17.6] false positives in a direct 1,000-
   replicate engine simulation (study: 14.1% [13.7, 14.6]); the absolute effect gives 1.5%
   [0.8, 2.5] (study: 2.1% [1.9, 2.3])
   ([scripts/gbstats_rare_conversion_repro.py](scripts/gbstats_rare_conversion_repro.py)).
   The draft also covers the tiny- and skewed-arm results above and asks for a documented
   minimum sample size.
4. **gbstats packaging:** gbstats 0.8.0 on PyPI (2024-06-20) is not the code on GrowthBook's
   main branch, which also reports version 0.8.0. Based on a source comparison; to be
   re-checked before filing.

The spotify GST excesses on small, skewed or rare-event arms are not drafted: the reference
GST shows the same pattern, so they are not engine defects.

## Limitations

- Results apply to gbstats 0.8.0 and spotify-confidence 4.1.0 as published on PyPI and to
  the simulated conditions. GrowthBook's hosted code is newer than the PyPI release.
- Product-layer behaviour (minimum sample sizes, data-quality checks) is not audited.
- 140 looks understate truly continuous monitoring; more looks can only add rejections.
- Excesses below about 5.5% FPR (or coverage above about 94.5%) can go undetected.
- Engine rates come from mirrors verified on a subsample of replicates and on randomised
  grids; the two identity gates are self-consistency checks only.
- The engine fixtures from the gbstats sdist and the spotify repository are not reproduced
  (only the scaffold's gbstats fixture is); see DESIGN section 13.

## Reproduce

```sh
uv sync --locked
uv run python scripts/run_study.py --config full --out runs/full --workers 3  # ~55-70 min
uv run peekaudit sensitivity runs/full --workers 3                             # ~6 min
uv run peekaudit report runs/full --publish docs/results/full
uv run peekaudit compare docs/results/quick docs/results/full --out docs/results/full/quick_vs_full
```

The run checkpoints one file per unit and resumes where it stopped; the wrapper retries up
to three times. On a loaded machine the full run takes longer than 45 minutes; it can be
split into shorter sessions with `--units` (for example `--units 'normal-*'`), all writing
to the same output directory.
