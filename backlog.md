# Backlog

## In progress

## Waiting on me

- Pull request from `polish-2026-10-30` opened, CI green on all four jobs (CI has never run on this repository), then merged to `main`.
- `v0.1.1` tagged on the merged commit, with `date-released` in `CITATION.cff` matching the tag date (`v0.1.0` was never tagged; optional tag at `6790d89`).
- Private vulnerability reporting enabled in the repository settings, because `SECURITY.md` sends reporters to it and it is off.
- Decision on which of the five `docs/upstream/` drafts to file, after the re-check steps in `docs/upstream/README.md`.

## Decided

## Proposed

- Sensitivity column for GrowthBook's actual product rule (larger arm's metric total >= 150, both arms non-zero), so the headline rare-conversion and revenue rates can be read against what the product shows.
- gbstats sdist and spotify-confidence repository fixtures reproduced through the adapters (DESIGN 5.2 item 1; needs a download).
- The gbstats code on GrowthBook main audited as a second engine, since the audited PyPI release is older than what the product runs.
- CUPED study (`runs/cuped`, `docs/results/cuped`); plan in `AGENTS.md`.
