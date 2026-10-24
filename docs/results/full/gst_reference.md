# spotify.gst against like-for-like reference GSTs (A/A units)

Study `full`, config hash `dde79da82e86a243`, 20,000 replicates per unit, master seed 20261023. Same data for all three methods. `ref.gst_pow2` uses spotify's spending function (power family, rho = 2) with reference variances (n - 1; no result for an arm with n <= 1 or zero variance). FPR in % with 95% Clopper-Pearson intervals. McNemar: exact test of spotify.gst vs ref.gst_pow2 on the same replicates, Holm-adjusted over all 72 cells. 'ref n/a' = spotify rejections at looks where ref.gst_pow2 gives no result; 'SE 0' = spotify rejections with a zero standard error. Rows: any method >= 6.0%. Run counts reproduced: yes, every cell.

Violations (FPR > 5%, one-sided binomial test, Holm per method and family; primary / stress cells): spotify.gst 21 of 60 / 6 of 12; ref.gst_pow2 21 of 60 / 6 of 12; ref.gst_obf 17 of 60 / 6 of 12.

| unit | schedule | look-1 arms (A, B) | spotify.gst | ref.gst_pow2 | ref.gst_obf | spotify at look 1 | spotify only / pow2 only | McNemar Holm p | ref n/a | SE 0 |
|---|---|---|---|---|---|---|---|---|---|---|
| bern_p01-n1000-b20-aa | daily14 | 57, 14 | 20.5 [19.9, 21.1] | 20.5 [19.9, 21.1] | 19.7 [19.1, 20.2] | 0% | 0 / 0 | 1 | 0 | 0 |
| bern_p01-n1000-b20-aa | fixed1 | 800, 200 | 13.3 [12.9, 13.8] | 13.3 [12.9, 13.8] | 13.3 [12.9, 13.8] | 100% | 1 / 0 | 1 | 0 | 0 |
| bern_p01-n1000-b20-aa | stream140 | 6, 1 | 19.7 [19.2, 20.3] | 18.8 [18.2, 19.3] | 19.6 [19.0, 20.1] | 5% | 192 / 0 | 2.3e-56 | 192 | 180 |
| bern_p01-n20000-b20-aa | daily14 | 1143, 286 | 7.9 [7.5, 8.2] | 7.9 [7.5, 8.2] | 5.7 [5.4, 6.1] | 19% | 0 / 0 | 1 | 0 | 0 |
| bern_p01-n20000-b20-aa | stream140 | 114, 29 | 9.3 [8.9, 9.7] | 9.2 [8.8, 9.7] | 5.9 [5.6, 6.2] | 0% | 2 / 0 | 1 | 0 | 0 |
| bern_p10-n1000-b20-aa | daily14 | 57, 14 | 9.9 [9.5, 10.4] | 9.6 [9.2, 10.0] | 6.7 [6.3, 7.0] | 5% | 69 / 0 | 2.4e-19 | 0 | 0 |
| bern_p10-n1000-b20-aa | stream140 | 6, 1 | 20.4 [19.8, 20.9] | 11.5 [11.0, 11.9] | 7.0 [6.7, 7.4] | 44% | 1780 / 0 | < 1e-300 | 1812 | 1079 |
| lognorm-n1000-b20-aa | daily14 | 57, 14 | 10.0 [9.6, 10.4] | 10.0 [9.6, 10.4] | 8.7 [8.3, 9.1] | 2% | 0 / 0 | 1 | 0 | 0 |
| lognorm-n1000-b20-aa | fixed1 | 800, 200 | 7.1 [6.7, 7.4] | 7.1 [6.7, 7.4] | 7.1 [6.7, 7.4] | 100% | 0 / 0 | 1 | 0 | 0 |
| lognorm-n1000-b20-aa | stream140 | 6, 1 | 10.4 [9.9, 10.8] | 10.4 [9.9, 10.8] | 8.9 [8.5, 9.3] | 0% | 2 / 2 | 1 | 0 | 0 |
| lognorm-n20000-b20-aa | daily14 | 1143, 286 | 6.7 [6.3, 7.0] | 6.7 [6.3, 7.0] | 6.0 [5.7, 6.4] | 3% | 1 / 0 | 1 | 0 | 0 |
| lognorm-n20000-b20-aa | stream140 | 114, 29 | 6.9 [6.5, 7.2] | 6.9 [6.5, 7.2] | 6.1 [5.8, 6.4] | 0% | 1 / 0 | 1 | 0 | 0 |
| normal-n1000-b20-aa | stream140 | 6, 1 | 6.5 [6.2, 6.9] | 6.5 [6.2, 6.9] | 5.2 [4.9, 5.5] | 0% | 0 / 0 | 1 | 0 | 0 |
| pareto_a15-n1000-b20-aa | daily14 | 57, 14 | 10.0 [9.6, 10.4] | 10.0 [9.6, 10.4] | 9.1 [8.7, 9.5] | 2% | 0 / 0 | 1 | 0 | 0 |
| pareto_a15-n1000-b20-aa | fixed1 | 800, 200 | 7.6 [7.2, 7.9] | 7.6 [7.2, 7.9] | 7.6 [7.2, 7.9] | 100% | 0 / 0 | 1 | 0 | 0 |
| pareto_a15-n1000-b20-aa | stream140 | 6, 1 | 10.3 [9.9, 10.7] | 10.3 [9.9, 10.7] | 8.9 [8.5, 9.3] | 0% | 0 / 1 | 1 | 0 | 0 |
| pareto_a15-n20000-b20-aa | daily14 | 1143, 286 | 8.1 [7.7, 8.5] | 8.1 [7.7, 8.5] | 7.4 [7.1, 7.8] | 1% | 0 / 0 | 1 | 0 | 0 |
| pareto_a15-n20000-b20-aa | fixed1 | 16000, 4000 | 6.5 [6.2, 6.9] | 6.5 [6.2, 6.9] | 6.5 [6.2, 6.9] | 100% | 0 / 0 | 1 | 0 | 0 |
| pareto_a15-n20000-b20-aa | stream140 | 114, 29 | 8.4 [8.0, 8.8] | 8.4 [8.0, 8.8] | 7.4 [7.1, 7.8] | 0% | 1 / 0 | 1 | 0 | 0 |
| pareto_a25-n1000-b20-aa | daily14 | 57, 14 | 8.5 [8.1, 8.9] | 8.5 [8.1, 8.9] | 7.4 [7.0, 7.8] | 3% | 0 / 0 | 1 | 0 | 0 |
| pareto_a25-n1000-b20-aa | fixed1 | 800, 200 | 6.1 [5.8, 6.4] | 6.1 [5.8, 6.4] | 6.1 [5.8, 6.4] | 100% | 0 / 0 | 1 | 0 | 0 |
| pareto_a25-n1000-b20-aa | stream140 | 6, 1 | 9.3 [8.9, 9.7] | 9.3 [8.9, 9.7] | 7.6 [7.2, 8.0] | 0% | 1 / 0 | 1 | 0 | 0 |
| pareto_a25-n20000-b20-aa | daily14 | 1143, 286 | 6.3 [5.9, 6.6] | 6.3 [5.9, 6.6] | 5.8 [5.5, 6.1] | 2% | 0 / 0 | 1 | 0 | 0 |
| pareto_a25-n20000-b20-aa | stream140 | 114, 29 | 6.5 [6.1, 6.8] | 6.5 [6.1, 6.8] | 5.8 [5.5, 6.2] | 0% | 0 / 0 | 1 | 0 | 0 |
