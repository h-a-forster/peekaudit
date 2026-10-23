# peekaudit design

Status: v0.1 design, written before implementation. Implemented so far: `peekaudit.stats`,
`peekaudit.gst` and the CLI `versions` command (see "Milestones").

## 1. Question and claims under test

Experimentation platforms say you can check sequential results as often as you like. This
project tests that claim against the shipped open-source engine code, in conditions vendor
studies leave out: skewed and heavy-tailed revenue, ratio metrics, unequal allocation, small
samples, very frequent looks, and running past the planned sample size.

| Engine (PyPI, licence) | Method audited | What the vendor documents |
|---|---|---|
| `gbstats` 0.8.0 (MIT), GrowthBook's stats engine | `SequentialTwoSidedTTest`: asymptotic confidence sequence (Waudby-Smith et al.), tuning N* = 5000 | You can look at results "as many times as you like" while keeping false positives below the expected rate ([docs](https://docs.growthbook.io/statistics/sequential)). The docs give no minimum sample size, burn-in or tail condition. |
| `gbstats` 0.8.0 | `TwoSidedTTest` peeked (sequential testing is off by default in GrowthBook) | Fixed-horizon test. No sequential claim. Control. |
| `gbstats` 0.8.0 | `EffectBayesianABTest`, default improper prior, "chance to win" >= 95% or <= 5% | The docs warn that stopping early can still inflate false positives. No sequential claim. Contrast case. |
| `spotify-confidence` 4.1.0 (Apache-2.0) | `ZTest` group sequential test (Lan-DeMets power-family spending, rho = 2, planned sample size) | Spotify's comparison: GST holds 5% when the planned sample size is right and inflates when you oversample ([blog](https://engineering.atspotify.com/2023/03/choosing-sequential-testing-framework-comparisons-and-discussions)). It used normal data with known variance, 500 per group, and no skewed data. |
| `spotify-confidence` 4.1.0 | `ZTest` fixed-horizon, peeked | No sequential claim. Control. |

The GrowthBook product defaults come from `packages/shared/src/constants.ts` at commit
`de93c5d6e1`: Bayesian engine, confidence level 0.95, sequential testing disabled, tuning
parameter 5000, improper prior, minimum sample size 150. The minimum sample size is applied by
the product layer, not by gbstats; see section 4.3.

**PyPI version and hosted version.** `gbstats` 0.8.0 was released to PyPI on 2024-06-20. The
GrowthBook main branch has changed since then (for example, post-stratification is on by
default) but still reports version 0.8.0. This audit covers the PyPI release. It does not
cover the code GrowthBook runs today. The report must say so.

## 2. Architecture

```
scenarios.py    metric generators, true moments, effect calibration, look schedules, grids
simulate.py     chunked, seeded generation of cumulative sufficient statistics (LookStats)
methods/
  base.py       Method protocol, MethodOutput, capabilities and claims, ScalarMethod helper
  reference.py  naive z, mSPRT, asymptotic CS, Lan-DeMets OBF GST (uses gst.py)
  gbstats_adapter.py   engine-in-the-loop calls + vectorised mirror, lazy import
  spotify_adapter.py   engine-in-the-loop calls + vectorised mirror, lazy import
  registry.py   name -> factory ("gbstats.seq.rel", "ref.acs", ...)
gst.py          independent Lan-DeMets solver (recursive integration)        [implemented]
stats.py        Clopper-Pearson, exact binomial test, Holm, mean CI           [implemented]
verify.py       engine vs mirror agreement on a replicate subsample
runner.py       work units, process pool (<= 4 workers), checkpoints, progress log
conformance.py  batteries, check_conformance(), assert_conforms()
report.py       tables (Markdown/JSON), Holm families, figures (matplotlib, optional)
cli.py          peekaudit run | status | report | conform | verify | versions
```

Data flow: `scenarios` defines a unit, which is (metric, N_max, allocation, hypothesis).
`simulate` produces cumulative per-arm sufficient statistics for R replicates at every
fine-grid look. Each method returns, per replicate and look, a reject flag (signed) and an
interval, or "no result". The `runner` derives each schedule from the fine grid, applies
the stopping rule, aggregates and checkpoints. `report` applies the exact inference and
Holm, and writes REPORT.md and the figures.

### 2.1 Core types

```python
@dataclass(frozen=True)
class LookStats:  # arrays of shape (R, K) per arm, cumulative at each look
    n_a: np.ndarray
    sum_a: np.ndarray
    sumsq_a: np.ndarray
    n_b: np.ndarray
    sum_b: np.ndarray
    sumsq_b: np.ndarray
    # ratio metrics only (numerator m, denominator d, per unit):
    dsum_a: np.ndarray | None
    dsumsq_a: ...
    mdsum_a: ...  # and the same for b
    binary: bool  # Bernoulli outcomes (engines switch variance formulas on this)


@dataclass(frozen=True)
class MethodOutput:  # arrays of shape (R, K)
    reject: np.ndarray  # int8: +1 reject with B > A, -1 reject with B < A, 0 no reject
    lower: np.ndarray  # interval at that look (NaN when the method gives none)
    upper: np.ndarray
    no_result: np.ndarray  # bool: the engine returned its "cannot analyse" default output
    scale: Literal["absolute", "relative"]
```

`Method` protocol: `name`, `scale`, `supports(metric) -> bool`, a
`claims(schedule, n_max) -> bool` predicate (section 6.3), `evaluate(stats, schedule_ctx) ->
MethodOutput` (vectorised), and, for engine adapters only, `evaluate_engine(stats, rows,
schedule_ctx)`, which calls the real engine for each listed replicate. `ScalarMethod` wraps
a per-look callback `decide(n_a, sum_a, sumsq_a, n_b, sum_b, sumsq_b) -> (reject, lo, hi)`,
so third parties can add an engine without vectorising it.

## 3. Simulation design

### 3.1 Metrics (generators and true moments)

| id | per-unit outcome | mean | variance | notes |
|---|---|---|---|---|
| `normal` | N(1, 1) | 1 | 1 | Spotify's vendor setting; oracle cell |
| `bern_p10` | Bernoulli(0.10) | p | p(1-p) | |
| `bern_p01` | Bernoulli(0.01) | p | p(1-p) | rare conversions; many zero-count early looks |
| `lognorm` | LogNormal(mu=0, sigma=1.5) | e^{mu+sigma^2/2} = 3.080 | (e^{sigma^2}-1) e^{2mu+sigma^2} = 80.6 | skewness 33.5 |
| `pareto_a25` | classical Pareto, x_m = 1, tail a = 2.5 | a x_m/(a-1) = 1.667 | a x_m^2/((a-1)^2(a-2)) = 2.222 | finite variance, infinite third moment |
| `pareto_a15` | classical Pareto, x_m = 1, a = 1.5 | 3 | infinite | stress: outside every method's assumptions; A/A only |
| `ratio` | sessions S = 1 + Poisson(1.5), clicks C given S ~ Binomial(S, 0.2); metric sum C / sum S | q = 0.2 | linearised q(1-q)/E[S] = 0.064 | delta-method ratio; gbstats and references only |

`numpy.random.Generator.pareto(a)` draws Lomax (Pareto II), so classical Pareto =
`x_m * (1 + rng.pareto(a))`. Generators must be validated against the closed-form moments
(tests with 10^6 draws, tolerance from the CLT; for `pareto_a25`, check the tail index with
a Hill estimator, because the sample variance converges slowly).

Exact shortcuts that keep the per-interval sufficient statistics exactly distributed:
Bernoulli uses sum ~ Binomial(m, p) per interval. Normal uses sum ~ N(m mu, m sigma^2) and
sumsq = sum^2/m + sigma^2 chi^2_{m-1}, which are independent. Everything else is drawn per
unit.

### 3.2 Effects (A/B units)

Treatment is a multiplicative lift r on the control distribution (B = (1 + r) A for
continuous metrics; p_B = p(1 + r) for Bernoulli; q_B = q(1 + r) for the ratio clicks).
This makes the relative effect exactly r and the absolute effect r * mean. For each
(metric, N_max, allocation), r solves |Delta(r)| = (z_{0.975} + z_{0.80}) *
sqrt(var_A / n_A + var_B(r) / n_B) at N_max with the true variances (brentq). In words, a
fixed-horizon two-sided z-test at N_max has 80% power. `pareto_a15` has no A/B unit.
Coverage targets are the exact true effects: r for relative engines, r * mean for absolute
ones.

### 3.3 Allocation and looks

- Total sample sizes N_max in {1,000, 20,000}. Allocation to B pi in {0.5, 0.2}.
- Arm counts are deterministic, as in blocked randomisation and Spotify's simulation:
  n_B,k = round(pi N_k), n_A,k = N_k - n_B,k. Random sample ratio noise (SRM) is out of
  scope. This also makes GST information fractions identical across replicates.
- Fine grid: 210 looks, N_k = round(k N_max / 140) for k = 1..210, so data runs to
  1.5 N_max. Schedules are views of the fine grid:
  - `fixed1`: one look at k = 140 (N_max). Fixed-horizon baseline.
  - `daily14`: k = 10, 20, ..., 140. "Daily for two weeks".
  - `stream140`: k = 1..140. The first look has N_max/140 units (7 at N_max = 1,000). This
    is a discrete stand-in for continuous monitoring, so it under-estimates the
    continuous-time error rate of AVI methods. The report states this.
  - `overrun21`: k = 10, 20, ..., 210. Daily looks that continue to 1.5 N_max while GSTs
    still plan for N_max. Spotify says oversampling breaks GST guarantees; this measures by
    how much in the shipped code.
- Stopping rule: stop at the first look that rejects. Otherwise the run ends at the
  schedule's last look. A/A: false positive = any rejection. A/B: power = rejection with
  the correct sign; type S rate = rejection with the wrong sign; stopping N = N at the stop.

### 3.4 Units, replicates, seeds

- Unit = (metric, N_max, pi, hypothesis): 7 x 2 x 2 x 2 - 4 (no `pareto_a15` A/B) = 52
  units. Within a unit, every method sees the same data (common random numbers). Methods
  are compared on paired outcomes.
- Default study ("quick", this machine): R_AA = 4,000, R_AB = 2,000 replicates, engine
  verification on the first k = 10 replicates per unit (spotify `stream140`: k = 2; see
  section 5.3).
- Full study (documented, not run by default): R_AA = 20,000, R_AB = 10,000, k = 100
  (`peekaudit run --config full`).
- Seeds: master seed 20261023. Unit seed = `SeedSequence(master, spawn_key=(h(unit_id),
  chunk))` where h is the first 8 bytes of SHA-256 of the unit id. Adding or reordering
  units does not change any other unit's stream. Replicates are generated in fixed chunks
  of 250 (memory: 250 x 30,000 float64 = 60 MB per array at N_max = 20,000), so results
  do not depend on worker count.

### 3.5 Grid size and compute budget

- Methods: 11 (section 4). Evaluations per unit: 11 methods x 4 schedules, minus methods
  that do not support the metric.
- Measured costs (this machine, Python 3.11):
  - gbstats `SequentialTwoSidedTTest` ~0.18 ms per call, `TwoSidedTTest` ~0.8 ms,
    `EffectBayesianABTest` ~5.7 ms (its risk computation uses `truncnorm`).
  - spotify-confidence `ZTest.difference` with sequential bounds: ~0.4 s for 14 looks and
    ~8.5 s for 140 looks per replicate path.
  - Vectorised mirrors: microseconds per replicate-look.
- Estimate for the quick study with 4 workers: data generation ~3 CPU-min, vectorised
  evaluation ~5 CPU-min, engine verification ~35-45 CPU-min, about 15-20 min wall in
  total, within the 45-minute target. M5 measures a one-unit pilot first and adjusts k
  (not R) if the projection exceeds 45 minutes. Full study: estimate 3-5 h wall. Confirm
  with the pilot projection before launching.

## 4. Methods: exact specifications

Notation per look: n, S = sum, Q = sum of squares per arm; xbar = S/n;
s^2 = (Q - S^2/n)/(n - 1); Bernoulli MLE variance v = xbar (1 - xbar).
Absolute effect D = xbar_B - xbar_A with V = s_A^2/n_A + s_B^2/n_B.
Relative effect (gbstats) D_rel = (xbar_B - xbar_A)/xbar_A with
V_rel = v_B/(xbar_A^2 n_B) + v_A xbar_B^2/(xbar_A^4 n_A) (delta method;
`frequentist_variance`). All tests are two-sided with alpha = 0.05.

### 4.1 gbstats 0.8.0 (engine code paths: `gbstats/frequentist/tests.py`, `gbstats/bayesian/tests.py`, `gbstats/models/statistics.py`)

Statistics: Bernoulli uses `ProportionStatistic(n, sum)` (variance p(1-p), MLE). Continuous
metrics use `SampleMeanStatistic(n, sum, sum_squares)` (n - 1 variance; 0 when n <= 1).
Ratio uses `RatioStatistic(n, m_statistic, d_statistic, m_d_sum_of_products)` with
delta-method variance m_var/d_mean^2 - 2 cov m_mean/d_mean^3 + m_mean^2 d_var/d_mean^4.

Default output ("no result": p = 1, CI [0, 0], chance to win 0.5): the control mean is 0, or
either variance is <= 0. Bayes also returns no result when either n = 0. Mirrors must
reproduce these branches exactly, and they get dedicated branch tests.

1. `gbstats.seq.rel` (default), `gbstats.seq.rel.tuned` (N* = N_max), `gbstats.seq.abs`
   (absolute, N* = 5000). With N = n_A + n_B:
   rho = sqrt((-2 ln a + ln(1 - 2 ln a)) / N*) (arXiv:2103.06476v7, eq. 161),
   half-width = sqrt(V N) * sqrt(2 (N rho^2 + 1) ln(sqrt(N rho^2 + 1)/a) / (N rho)^2)
   (eq. 9), p = min(1, sqrt(N rho^2 + 1) exp(-rho^2 N D^2 / (2 V (N rho^2 + 1)))) (eq. 155).
   Reject iff p < alpha. This is exactly "the CS excludes 0" (tested:
   `tests/test_engines_env.py::test_gbstats_sequential_ci_and_pvalue_are_dual`).
2. `gbstats.ttest.rel` (fixed horizon, peeked): Welch df
   nu = (v_A/n_A + v_B/n_B)^2 / ((v_A/n_A)^2/(n_A-1) + (v_B/n_B)^2/(n_B-1)) computed from
   raw variances, also for relative effects as in the engine. p = 2(1 - F_t(|D_rel|/sqrt(V_rel); nu)).
   CI D_rel +/- t_{0.975,nu} sqrt(V_rel).
3. `gbstats.bayes.ctw95` (improper prior, relative): posterior N(D_rel, V_rel), chance to
   win = Phi(D_rel / sqrt(V_rel)), reject iff chance to win >= 0.95 or <= 0.05, interval =
   95% credible interval. **Identity:** this rule is exactly a two-sided z-test at level
   0.10 on the relative scale. The mirror is tested for this identity, and the report says
   so.

### 4.2 spotify-confidence 4.1.0 (`analysis/frequentist/confidence_computers/z_test_computer.py`, `sequential_bound_solver.py`)

Input: one DataFrame row per (look, group) with cumulative `sum`, `sum_squares` and `n`.
`ZTest(..., categorical_group_columns="group", ordinal_group_column="look")`, then
`difference("A", "B", groupby="look", final_expected_sample_size_column=..., verbose=True)`.
Variance: if `sum_squares == sum` on every row (binary), v = Q/n - xbar^2. Otherwise
(Q - S^2/n)/(n - 1). SE = sqrt(v_A/n_A + v_B/n_B). Absolute difference.

4. `spotify.gst`: information fraction t_k = N_k / max(N_max, max_{j<=k} N_j) with
   N_k = n_A,k + n_B,k. Bounds come from `sequential_bounds(t, alpha=0.05, sides=2)` =
   `bounds(t, alpha, rho=2, ztrun=8, sides, max_nints=1000)`, a port of R ldbounds with
   power-family spending 0.025 t^2 per side. adjusted alpha_k = 2(1 - Phi(zb_k)),
   adjusted CI D +/- zb_k SE, decision = engine `is_significant`.
   - Look-by-look semantics: at look k the analyst passes looks 1..k. Verified: one call
     on the full path gives identical rows when no total exceeds the planned size. Past the
     planned size the engine rescales every earlier fraction by the current total, so a
     full-path call is NOT equivalent (measured on a 12-look path planned for 8: look-10
     adjusted alpha 0.0198 full-path vs 0.0297 look-by-look; M4 turns this into a
     regression test). In-the-loop verification of `overrun21` therefore calls
     look-by-look, and the mirror recomputes bounds per overrun look with
     t_j = N_j / N_k, j <= k.
   - Mirror: reject iff |D|/SE > zb_k, with zb from the engine's own `sequential_bounds`,
     cached per schedule (the bounds do not depend on the data, given deterministic N_k).
     Independent cross-check: `peekaudit.gst.lan_demets_bounds(t, spend_power(0.025, 2))`
     agrees with the engine to 2e-3 (tested for K = 5, 10).
5. `spotify.z` (fixed horizon, peeked): same input without `final_expected_sample_size`;
   p = 2(1 - Phi(|D|/SE)), reject iff p < 0.05, CI D +/- z_{0.975} SE.
   Edge cases seen: arms with n = 1 give SE NaN or 0 and `is_significant` False. A zero
   denominator raises `ValueError`; the adapter records that as no result.

### 4.3 Product-layer guard (optional variant, M3)

GrowthBook hides results until the baseline reaches `DEFAULT_MIN_SAMPLE_SIZE = 150`. The
exact quantity it is compared against (units, conversions or metric total) must be confirmed
from the GrowthBook source before it is implemented. If confirmed, add
`gbstats.seq.rel.guard150`: no decision before the guard passes. If it cannot be confirmed
from source, drop the variant and say so in the report.

### 4.4 Reference implementations (absolute scale, n - 1 variances)

6. `ref.naive_z`: reject iff |D|/sqrt(V) > z_{0.975}; CI D +/- z_{0.975} sqrt(V).
7. `ref.msprt` (Johari, Koomen, Pekelis and Walsh 2022, normal mixture, plug-in variance):
   Lambda_k = sqrt(V/(V + tau^2)) exp(tau^2 D^2 / (2 V (V + tau^2))). Reject iff
   Lambda_k >= 1/alpha. Always-valid CI D +/- sqrt(V (V + tau^2)/tau^2 * (2 ln(1/alpha) +
   ln((V + tau^2)/V))). tau = the unit's design effect r * mean (section 3.2), the same for
   the A/A and A/B units of a cell.
8. `ref.acs` (Waudby-Smith et al. 2024): formula 1 on the absolute scale with N* = N_max.
   **Identity:** on non-binary metrics this equals `gbstats.seq.abs` run with N* = N_max
   (the same formula and variance). For Bernoulli the variances differ (MLE vs n - 1).
   Tested as an exact mirror identity.
9. `ref.gst_obf`: Lan-DeMets O'Brien-Fleming-type spending (gsDesign `sfLDOF`
   convention, 0.025 per side), t_k = N_k/N_max. Spending is capped at t = 1, so there are
   no rejections after the planned maximum, which is the correct protocol when a study
   overruns. Reject iff |D|/sqrt(V) >= b_k; repeated CI D +/- b_k sqrt(V).

Known mathematical relations, used as oracles:
- `ref.acs` is `ref.msprt` with tau^2 = rho^2 * N V (data-dependent tau).
- For Gaussian data with known variance, the mixture martingale behind msprt and ACS is
  continuous in continuous time, so P(ever reject) = alpha exactly. Discrete looks give
  less than alpha, approaching alpha as looks get denser.

## 5. Validation strategy (before any audit number is trusted)

### 5.1 Oracles and reference values (pinned)

| Check | Reference value | Where |
|---|---|---|
| Naive repeated z-tests, K equal looks, nominal 5% | K = 2: 0.083, 3: 0.107, 4: 0.126, 5: 0.142, 10: 0.193, 20: 0.248 (Armitage et al. 1969) | `tests/test_gst.py` (done) |
| Pocock constant, two-sided 0.05 | K = 2..5: 2.178, 2.289, 2.361, 2.413 (Jennison and Turnbull 2000, Table 2.1) | done |
| O'Brien-Fleming constant | K = 2..5: 1.977, 2.004, 2.024, 2.040 (Table 2.3) | done |
| LD OBF-type bounds, K = 5, 0.025 per side | 4.8769, 3.3570, 2.6803, 2.2898, 2.0310 (ldbounds/gsDesign) | done |
| Power-family rho = 2, K = 5 | 3.0902, 2.7141, 2.4728, 2.2799, 2.1140 (spotify-confidence 4.1.0; independently reproduced) | done |
| Power-family rho = 2, K = 10 | 3.4808, 3.1523, 2.9508, 2.7935, 2.6611, 2.5450, 2.4406, 2.3451, 2.2566, 2.1737 | engine cross-check (done) |
| gbstats fixtures | `test_sequential_test_runs_prop`: expected 0.50386, CI [0.40098, 0.60675], sd 0.03386, plus the other cases in gbstats' own `tests/` (sdist) | adapter tests (one done) |
| gbstats CS / p-value duality | reject iff CS excludes 0 | done |
| Clopper-Pearson, Holm | closed forms; statsmodels `proportion_confint(method="beta")`, `multipletests(method="holm")` | `tests/test_stats.py` (done) |

### 5.2 Adapter validation (gate before audit, M3/M4)

1. Each adapter reproduces the engine's own fixtures: all gbstats 0.8.0 frequentist,
   sequential and Bayesian fixtures from its sdist `tests/`, and for spotify-confidence a
   small subset of `tests/frequentist/test_ztest.py` sequential cases at commit `0f12dd39af`.
   Numerical test vectors are copied with attribution; no engine code is vendored.
2. Mirror equals engine: on a randomised grid of 2,000 sufficient-statistic inputs per
   method, including degenerate inputs (n in {0, 1, 2}, zero sums, all-ones Bernoulli,
   zero control mean, huge variances, unequal n), decisions must agree exactly. Intervals
   must agree to rtol 1e-9 / atol 1e-12, and no-result flags must agree. A tie is
   |statistic - threshold| < 1e-9 relative; ties are counted, not failed.
3. In-the-loop agreement in every study unit: the first k replicates are re-run through
   the real engine look by look (spotify: full-path calls where equivalent, look-by-look
   past the planned size). Any decision mismatch fails the unit, and the unit is not
   reported. The counts are reported per unit.
4. Bug policy: an apparent engine defect is reported only with (a) a minimal standalone
   reproduction that calls the engine directly, without peekaudit, (b) the mirror and the
   paper formula both disagreeing with the engine, and (c) a second independent check that
   it is not an adapter mistake (wrong statistic type, scale, alpha, sides or input
   semantics). Calibration failures under conditions outside a method's assumptions are
   findings about the method, not bugs.

### 5.3 Monte Carlo sanity gates (must pass before REPORT.md is written)

- Positive control: `ref.naive_z` on `normal`, N_max = 20,000, pi = 0.5. FPR for `daily14`
  and `stream140` must contain the exact value `gst.crossing_probability(t, z_{0.975})` with
  the realised fractions t inside its 99.9% CP interval.
- Negative controls: on the same cells, `ref.gst_obf` FPR has 0.05 inside its 99.9% CP
  interval on `fixed1`/`daily14`/`stream140`, and `ref.acs`/`ref.msprt` have CP lower bound
  <= 0.05.
- Fixed-horizon calibration: every method's `fixed1` A/A FPR is consistent with its
  fixed-look level (0.05; `gbstats.bayes.ctw95`: 0.10 by the identity) on `normal`,
  N_max = 20,000.
- Power calibration: `ref.naive_z` `fixed1` A/B power has 0.80 inside its 99.9% CP
  interval on `normal` and `bern_p10` at N_max = 20,000 (validates section 3.2).
- Identities: `gbstats.bayes.ctw95` equals the level-0.10 relative z-test per replicate.
  `ref.acs` equals `gbstats.seq.abs` (N* = N_max) per replicate on non-binary metrics.
- Determinism: rerunning one unit with a different worker count or chunk order gives
  byte-identical unit JSON.

### 5.4 Monte Carlo size and power

FPR violation test: one-sided exact binomial test of H0 FPR <= 0.05, Holm across the
m <= 84 primary A/A cells of a method. Critical rate and power (computed exactly):

| R | Holm level | critical FPR | power at 6.0% | 6.5% | 7.0% | 8.0% |
|---|---|---|---|---|---|---|
| 4,000 (quick) | 0.05/84 | 6.17% | 0.33 | 0.81 | 0.98 | 1.00 |
| 20,000 (full) | 0.05/84 | 5.51% | 0.998 | 1.00 | 1.00 | 1.00 |

95% CP half-widths: about +/-0.7 pp at 5% (R = 4,000), +/-0.3 pp (R = 20,000), and
+/-1.8 pp at 80% power (R_AB = 2,000). The quick study can reliably flag inflation to 7%
or more. Smaller excesses need the full study. The report must state this detection limit.
`scripts/reference_values.py` regenerates this table and the reference values in 5.1.

## 6. Inference on the simulation output and reporting

### 6.1 Per cell (method x unit x schedule)

- FPR (A/A) with a 95% CP interval, plus a one-sided exact p-value for FPR > alpha.
- Coverage at stopping (A/B; also A/A for the true effect 0): the share of replicates whose
  interval at the stopping look contains the true effect, with a CP interval. One-sided
  exact p-value for misses > alpha. Reported over replicates with an interval; the
  no-result share is reported next to it.
- Power, type S rate (CP intervals), mean stopping N / N_max (t-interval), and the
  stopping-look histogram.
- Engine verification counts (k, mismatches, ties).

### 6.2 Families and multiple testing

Holm (strong FWER 0.05) per method and per question:
- FPR, primary family: cells within the method's claim (6.3) on finite-variance metrics.
- FPR, stress family: `pareto_a15` cells, reported separately.
- Coverage family: A/B cells within claim.
Classification: **violation** (Holm-adjusted p < 0.05), **consistent** (not rejected;
report the CP upper bound, e.g. "FPR <= 6.1%"), **conservative** (CP upper < 0.025,
informational). Cells outside a method's claim (peeked fixed-horizon tests, Bayes, GST past
N_max) get **informational** rates with intervals and no pass/fail.

### 6.3 Claims per method

| Method | Claimed regime |
|---|---|
| gbstats.seq.* , ref.acs, ref.msprt | every schedule, every finite-variance metric, any n (the vendor states no minimum) |
| spotify.gst, ref.gst_obf | schedules that stay within N_max (`fixed1`, `daily14`, `stream140`); `overrun21` is informational |
| gbstats.ttest.rel, spotify.z, ref.naive_z | `fixed1` only |
| gbstats.bayes.ctw95 | none (informational everywhere) |

### 6.4 Headline deliverable

`docs/results/<study>/`: `summary.json` (every cell, with counts, so all numbers can be
recomputed), and `REPORT.md` with:
1. A cross-engine table of realised FPR with 95% CP intervals per method x metric x N_max x
   allocation x schedule, with Holm classification.
2. Coverage at stopping, power, type S and mean stopping N for A/B cells.
3. Figures: FPR by schedule per method (dots with CP intervals, faceted by metric, nominal
   line), power vs mean stopping N, and FPR vs N_max for the asymptotic methods.
4. "What was run": R, seeds, versions (`peekaudit versions`), git SHA, wall time, the MC
   detection limit, and every deviation from this design.
5. A findings list in the form "engine X: FPR Y% [CP interval] on metric Z at N_max = n,
   schedule s (Holm-adjusted p = ...)". No finding without its interval and conditions.

## 7. Runner, checkpointing, logging

- Output root `runs/<study>/` (gitignored): `manifest.json` (config, config hash, versions,
  git SHA, start time), `units/<unit_id>.json` (aggregates; written to a temp file, then
  `os.replace`), `cache/<unit_id>.npz` (per-replicate stop look, sign and interval at stop,
  for debugging and re-analysis), and `progress.log`. The log has one line per finished
  unit: done/total, elapsed, rate, ETA.
- Resume: a unit is skipped if its JSON exists with a matching config hash and peekaudit
  version. `--force` recomputes. Units are independent, so a crash loses at most one unit
  per worker.
- Pool: `concurrent.futures.ProcessPoolExecutor(max_workers=min(4, os.cpu_count()))`, with
  the cap enforced in code. BLAS threads are pinned to 1 per worker (`OMP_NUM_THREADS`,
  `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` set before numpy import in workers).
- Resume-on-failure: `scripts/run_study.py` wraps `peekaudit run` in up to 3 attempts. The
  run is idempotent because of the checkpoints.
- `peekaudit status runs/<study>` prints done/total, failed units and the ETA from the log.
- Published results: `peekaudit report` copies `summary.json`, `REPORT.md` and the figures
  (PNG, small) into `docs/results/<study>/`. Raw `runs/` and `cache/` are never committed.

## 8. Conformance suite (reusable deliverable)

- `peekaudit.conformance.check_conformance(method, battery="quick") -> ConformanceReport`
  and `assert_conforms(method, battery="quick")`, which raises `AssertionError` with the
  failing cells. Any project can call it from its own pytest suite. A dedicated pytest
  plugin is not in v0.1; the function covers the same use.
- CLI: `peekaudit conform --method <registry name | module:factory> --battery
  smoke|quick|standard [--json out.json]`. Exit code 0 = no violation, 1 = violation,
  2 = error.
- Batteries (A/A only, the method's claimed schedules):
  - `smoke`: `normal`, `bern_p10`; N_max = 2,000; `daily14`; R = 1,000; under 30 s. Used in
    CI.
  - `quick`: adds `lognorm`, `pareto_a25`, N_max = 1,000 and 20,000, `stream140`;
    R = 2,000; a few minutes.
  - `standard`: the A/A half of the study grid at R = 4,000.
- The verdict rules are those of section 6.2. The report lists the battery's detection
  limit (5.4), so "pass" is never read as "proven calibrated".

## 9. CLI

```
peekaudit run [--config quick|full|PATH] [--out runs/NAME] [--workers N<=4] [--units GLOB] [--force]
peekaudit status runs/NAME
peekaudit report runs/NAME [--publish docs/results/NAME]
peekaudit verify [--method NAME] [--n 2000]          # mirror-vs-engine randomised grid
peekaudit conform --method NAME --battery smoke|quick|standard [--json PATH]
peekaudit versions
```

## 10. Scope decisions (v0.1)

Kept: the headline cross-engine table, all metrics in the original scope (normal added as
the oracle cell), unequal allocation, small n, daily vs stream looks, coverage at stopping,
power and expected N, exact CIs with Holm, checkpointing, and the conformance function and
CLI.

Cut or deferred, with reasons:
- Random (Bernoulli) assignment and SRM noise: deterministic counts make GST fractions
  exact and keep the engines out of empty-arm failure modes that are artefacts of the
  simulation. Deferred.
- A dedicated pytest plugin: `assert_conforms` covers it. Deferred.
- gbstats from the GrowthBook main branch (what the product runs today): the version
  string is unchanged and there is no PyPI artefact. Deferred; the report states the
  limitation.
- confseq and other e-value libraries: out of scope per the original plan. confseq 0.0.11
  (2023) has a C++ extension and no recent wheels.
- CUPED/regression adjustment, multi-arm, novelty drift, closed engines: out of scope.
- One-sided tests: the engines' defaults are two-sided. Deferred.
- Zero-inflated revenue (conversion x lognormal): a stretch goal if the budget allows.
- No upstream issues and no contact with maintainers until owner review.

## 11. Dependencies and environment

- Core: numpy, scipy. Engines (optional extra `engines`, pinned): `gbstats==0.8.0` (MIT),
  `spotify-confidence==4.1.0` (Apache-2.0), plus `bokeh<3.10` because chartify 5.0.1, a
  spotify-confidence dependency, fails to import with bokeh 3.10. Engines are imported,
  never vendored.
- gbstats 0.8.0 pins numpy<2 and pandas<2. These have no wheels for Python 3.12+, so
  `[tool.uv] override-dependencies` lifts them on 3.12+ only. Python 3.11 keeps the declared
  pins and is the canonical study environment (`.python-version`). gbstats' own 43 tests pass
  on Python 3.13 with numpy 2.x and pandas 2.x. Installing the `engines` extra with pip on
  3.12+ fails without the same overrides; the README documents this.
- CI runs ubuntu 3.11 (as shipped), ubuntu 3.13 and windows 3.13 (overridden), plus lint.
  Triggers: `pull_request` and `workflow_dispatch` only.

## 12. Milestones (each ends with `uv run pytest`, `uv run ruff check .` and `uv run ruff format --check .` green, a CHANGELOG entry and an AGENTS.md status update)

| # | Deliverable | Done-check |
|---|---|---|
| M0 | Scaffold, `stats.py`, `gst.py`, engine environment tests, this design | done (this commit) |
| M1 | `scenarios.py` + `simulate.py`: metrics, true moments, effect calibration, schedules, seeded chunked generation | moment tests (10^6 draws, CLT tolerances; Hill tail index for Pareto); exact-shortcut generators match per-unit generators in distribution (KS on 10^4 sums); calibrated r gives analytic power 0.80; same seed gives byte-identical arrays for any chunk order |
| M2 | `methods/base.py`, `methods/reference.py`, registry | closed-form unit tests per formula; msprt/ACS identity; GST uses `gst.py`; slow MC: naive-z FPR matches `crossing_probability`; ACS/msprt <= alpha on Gaussian known-variance data |
| M3 | `gbstats_adapter.py` (seq rel/rel.tuned/abs, ttest, bayes; guard variant if confirmed) | all gbstats sdist fixtures reproduced via the adapter; mirror equals engine on the 2,000-input randomised grid incl. degenerate branches; Bayes identity |
| M4 | `spotify_adapter.py` (gst, z) | engine fixtures; full-path equals look-by-look within N_max; overrun look-by-look semantics mirrored; bounds cache equals engine; randomised grid agreement |
| M5 | `runner.py`, `verify.py`, CLI `run`/`status`/`verify`; one-unit pilot | kill-and-resume test (unit JSON unchanged, completed units skipped); worker cap enforced; pilot timing projected for quick and full and recorded in AGENTS.md |
| M6 | Quick study run, `report.py`, figures | section 5.3 gates all pass; `docs/results/quick/` with summary.json, REPORT.md and figures; wall time < 45 min |
| M7 | `conformance.py` + `peekaudit conform`, smoke battery in CI | `ref.gst_obf` and `ref.acs` pass the smoke battery; `ref.naive_z` claimed for all schedules fails it (positive control); runs in < 30 s |
| M8 | README findings, CHANGELOG 0.1.0, handoff | every number in README traceable to summary.json; findings phrased per section 6.4; limitations listed |

## 13. References

- Armitage, P., McPherson, C. K. and Rowe, B. C. (1969). Repeated significance tests on
  accumulating data. JRSS A 132(2), 235-244.
- Howard, S. R., Ramdas, A., McAuliffe, J. and Sekhon, J. (2021). Time-uniform, nonparametric,
  nonasymptotic confidence sequences. Annals of Statistics 49(2), 1055-1080.
- Jennison, C. and Turnbull, B. W. (2000). Group Sequential Methods with Applications to
  Clinical Trials. Chapman and Hall/CRC.
- Johari, R., Koomen, P., Pekelis, L. and Walsh, D. (2022). Always valid inference:
  continuous monitoring of A/B tests. Operations Research 70(3), 1806-1821.
  arXiv:1512.04922.
- Kim, K. and DeMets, D. L. (1987). Design and analysis of group sequential tests based on
  the type I error spending rate function. Biometrika 74(1), 149-154.
- Lan, K. K. G. and DeMets, D. L. (1983). Discrete sequential boundaries for clinical
  trials. Biometrika 70(3), 659-663.
- O'Brien, P. C. and Fleming, T. R. (1979). A multiple testing procedure for clinical
  trials. Biometrics 35(3), 549-556.
- Pocock, S. J. (1977). Group sequential methods in the design and analysis of clinical
  trials. Biometrika 64(2), 191-199.
- Waudby-Smith, I., Arbour, D., Sinha, R., Kennedy, E. H. and Ramdas, A. (2024). Time-uniform
  central limit theory and asymptotic confidence sequences. Annals of Statistics 52(6).
  arXiv:2103.06476 (v7 equation numbers as cited in gbstats).
- Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian
  Journal of Statistics 6(2), 65-70.
- Clopper, C. J. and Pearson, E. S. (1934). The use of confidence or fiducial limits
  illustrated in the case of the binomial. Biometrika 26(4), 404-413.
- GrowthBook sequential testing docs: https://docs.growthbook.io/statistics/sequential
- Spotify, "Choosing a sequential testing framework" (2023):
  https://engineering.atspotify.com/2023/03/choosing-sequential-testing-framework-comparisons-and-discussions
