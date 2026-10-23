# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). This
project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- Project scaffold: uv + hatchling, `src/` layout, CI (pull requests and manual dispatch),
  repository hygiene files.
- `docs/DESIGN.md`: architecture, exact method specifications, validation strategy, study
  design and milestones.
- `peekaudit.stats`: Clopper-Pearson intervals, one-sided exact binomial test, Holm
  adjustment, t-interval for means.
- `peekaudit.gst`: independent Lan-DeMets group sequential solver (recursive numerical
  integration) with O'Brien-Fleming-type, power-family and Pocock-type spending, validated
  against Armitage et al. (1969), Jennison and Turnbull (2000) and ldbounds/gsDesign values.
- Engine environment checks: pinned gbstats 0.8.0 and spotify-confidence 4.1.0, a gbstats
  fixture, the gbstats confidence-sequence/p-value duality, and spotify-confidence bounds
  against the reference solver.
- CLI: `peekaudit --version`, `peekaudit versions`.
