# Security policy

Security fixes are made for the latest release.

## Reporting

Report vulnerabilities privately through the repository's **Security** tab ("Report a
vulnerability"). Do not open a public issue. Include the peekaudit version, the operating
system, the impact, and steps that reproduce it.

## Scope

peekaudit runs simulations locally and imports the engines you install. It makes no network
calls at run time.

In scope:

- peekaudit reading, writing or deleting files outside the output directory you chose (for
  example through unit ids, study names or `--units` patterns);
- `peekaudit conform --method module:factory` importing code other than the module named;
- reports that execute script from data in a run directory.

Out of scope:

- code in a module you pass to `--method` (it runs with your permissions by design);
- vulnerabilities in the audited engines or their dependencies (report those to their
  maintainers);
- statistical findings about an engine (see CONTRIBUTING.md for how those are handled).
