# Upstream issue drafts

Draft issue texts for engine behaviours found while building peekaudit. **None of them has
been filed.** Each draft has a standalone reproduction in `scripts/` that imports only the
engine (plus numpy, scipy and pandas), so the behaviour cannot come from peekaudit's
adapters.

| Draft | Engine | Kind | Reproduction |
|---|---|---|---|
| [spotify-confidence-indexerror-single-unit-control.md](spotify-confidence-indexerror-single-unit-control.md) | spotify-confidence 4.1.0 | exception on a degenerate input | `scripts/spotify_indexerror_repro.py` |
| [spotify-confidence-sequential-bounds-after-overrun.md](spotify-confidence-sequential-bounds-after-overrun.md) | spotify-confidence 4.1.0 | earlier looks' bounds change after an overrun | `scripts/spotify_overrun_repro.py` |
| [spotify-confidence-significance-single-unit-binary-arm.md](spotify-confidence-significance-single-unit-binary-arm.md) | spotify-confidence 4.1.0 | significant result from a one-unit binary group (SE 0 or one-sided) | `scripts/spotify_single_unit_repro.py` |
| [gbstats-sequential-relative-rare-events.md](gbstats-sequential-relative-rare-events.md) | gbstats 0.8.0 | sequential FPR above alpha (relative effect, few events; tiny or skewed arms) | `scripts/gbstats_rare_conversion_repro.py` |
| [gbstats-pypi-release-behind-main.md](gbstats-pypi-release-behind-main.md) | gbstats 0.8.0 | packaging: PyPI release differs from main at the same version | none (source comparison; re-check before filing) |

Not drafted, on purpose:

- spotify-confidence group sequential test exceeding 5% on small, skewed or rare-event arms,
  apart from the one-unit binary arm above. A reference GST with the same spending function
  (power family, rho = 2) and sample variances violates in exactly the same cells (21 of 60
  primary, 6 of 12 stress; `docs/results/full/gst_reference.md`), so the rest is a property
  of z-based group sequential tests on such data, not an engine defect.
- Over-rejection of fixed-horizon z- and t-tests with 1% conversions and 200 users in an arm
  (`scripts/zero_count_example.py`): a known property of the Wald test.

Before filing any draft:

1. Re-run the reproduction against the engine's latest release and main branch; update the
   version lines or drop the draft if the behaviour has changed.
2. Search the tracker for existing issues and link them instead of duplicating.
3. Keep the tone factual; the study numbers are in `docs/results/full/`.
