# Upstream reports

Issues filed from peekaudit's findings, with their state read from GitHub on 2026-11-13.
The issue texts are in [upstream/](upstream/README.md); REPORT.md, "Engine behaviours" and
"Three more engines", has the evidence.

| Project | Issue | Finding | Filed | State on 2026-11-13 | Outcome |
|---|---|---|---|---|---|
| spotify/confidence | [#118](https://github.com/spotify/confidence/issues/118) | `ZTest.difference` raises `IndexError` on a one-unit control group (behaviour 1) | 2026-11-07 | open | No reply |
| growthbook | [#7642](https://github.com/growthbook/growthbook/issues/7642) | gbstats 0.8.0 on PyPI is not the code on `main`, which reports the same version (behaviour 4) | 2026-11-07 | open | No reply; PyPI still has 0.8.0 from 2024-06-20 |
| growthbook | [#7649](https://github.com/growthbook/growthbook/issues/7649) | Sequential A/A false-positive rate above alpha for relative effects with few events or small arms (behaviour 3) | 2026-11-07 | open | No reply |
| spotify/confidence | [#119](https://github.com/spotify/confidence/issues/119) | `import spotify_confidence` fails on a fresh Python 3.12+ install: chartify 5.0.1 breaks on bokeh 3.10 (why peekaudit pins `bokeh<3.10`) | 2026-11-09 | open | No reply |
| tea-tasting | [#327](https://github.com/e10v/tea-tasting/issues/327) | Default `Proportion` confidence interval pools the variance and misses its level at any n | 2026-11-09 | open | No reply |
| expectation | [#155](https://github.com/jakorostami/expectation/issues/155) | Proportion test raises while all observations are 0 or all 1; README import fails | 2026-11-09 | open | No reply |
| spotify/confidence | [#110](https://github.com/spotify/confidence/issues/110#issuecomment-6833819535) (comment on an existing issue) | Bounds of earlier looks change after an overrun (behaviour 2), answering the question asked in #110 | 2026-11-12 | open | No reply |

Not filed: behaviour 5 (spotify-confidence significance from a one-unit binary group) was
dropped on 2026-11-12 as a degenerate input close to the open #118. spotify-confidence has
had no release or commit since 4.1.0 (2026-02-26), and its maintainers have not replied on
issues since 2025-04.
