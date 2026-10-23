# Cross-check: quick vs full

Config hashes `82e303ecf6547d58` (quick) and `dde79da82e86a243` (full); 2224 cells in both.

The two studies share the master seed and chunk size, so the quick replicates are the first replicates of full (nested, not independent). Each quick rate is therefore compared with the replicates of full not in quick (counts differenced), which is independent of it.

Two-sided Fisher exact test per cell, Holm-adjusted within each rate kind. Under unchanged, deterministic code every null is true, so rejections would point to a code or seeding difference.

| rate | cells compared | raw p < 0.05 | Holm p < 0.05 | smallest Holm p |
|---|---|---|---|---|
| FPR / power | 2,224 | 131 | 0 | 1 |
| coverage at stopping | 2,224 | 104 | 0 | 1 |

95% CP intervals of the two studies overlap in 2,224 of 2,224 cells for FPR / power and in 2,222 for coverage (the intervals are not independent when nested; shown for orientation).

Verdicts on claimed cells (DESIGN 6.2):

| quick | full | cells |
|---|---|---|
| conservative | conservative | 382 |
| conservative | consistent | 1 |
| consistent | conservative | 11 |
| consistent | consistent | 907 |
| consistent | violation | 43 |
| violation | violation | 132 |

Cells whose violation status differs (43):

| method | cell | family | quick | full |
|---|---|---|---|---|
| gbstats.seq.rel | pareto_a25-n20000-b20-aa stream140 | fpr_primary | 5.2 [4.6, 6.0] consistent | 5.6 [5.3, 5.9] violation |
| gbstats.seq.rel.tuned | bern_p01-n20000-b20-aa daily14 | fpr_primary | 5.2 [4.5, 5.9] consistent | 5.8 [5.5, 6.1] violation |
| gbstats.seq.rel.tuned | bern_p01-n20000-b20-aa overrun21 | fpr_primary | 5.8 [5.1, 6.5] consistent | 6.3 [6.0, 6.7] violation |
| gbstats.seq.rel.tuned | lognorm-n1000-b20-ab stream140 | coverage | 93.7 [92.5, 94.7] consistent | 94.0 [93.5, 94.5] violation |
| gbstats.ttest.rel | bern_p01-n20000-b20-aa fixed1 | fpr_primary | 5.6 [4.9, 6.4] consistent | 5.6 [5.3, 6.0] violation |
| gbstats.ttest.rel | bern_p10-n1000-b20-aa fixed1 | fpr_primary | 5.5 [4.8, 6.3] consistent | 5.8 [5.5, 6.1] violation |
| gbstats.ttest.rel | bern_p10-n1000-b50-aa fixed1 | fpr_primary | 5.4 [4.7, 6.1] consistent | 5.5 [5.2, 5.8] violation |
| gbstats.ttest.rel | bern_p10-n20000-b20-ab fixed1 | coverage | 93.7 [92.5, 94.7] consistent | 94.2 [93.7, 94.6] violation |
| gbstats.ttest.rel | lognorm-n20000-b20-ab fixed1 | coverage | 94.2 [93.1, 95.2] consistent | 94.3 [93.8, 94.7] violation |
| gbstats.ttest.rel | pareto_a25-n1000-b20-ab fixed1 | coverage | 93.9 [92.8, 94.9] consistent | 93.0 [92.5, 93.5] violation |
| gbstats.ttest.rel | pareto_a25-n20000-b20-ab fixed1 | coverage | 94.7 [93.6, 95.6] consistent | 94.3 [93.9, 94.8] violation |
| ref.acs | lognorm-n1000-b20-aa stream140 | fpr_primary | 5.3 [4.6, 6.1] consistent | 5.5 [5.2, 5.9] violation |
| ref.gst_obf | bern_p01-n1000-b20-ab fixed1 | coverage | 93.5 [92.3, 94.5] consistent | 92.7 [92.1, 93.2] violation |
| ref.gst_obf | bern_p01-n20000-b20-aa daily14 | fpr_primary | 5.7 [5.0, 6.4] consistent | 5.7 [5.4, 6.1] violation |
| ref.gst_obf | bern_p01-n20000-b20-aa stream140 | fpr_primary | 5.7 [5.0, 6.5] consistent | 5.9 [5.6, 6.2] violation |
| ref.gst_obf | bern_p10-n1000-b20-aa daily14 | fpr_primary | 6.0 [5.3, 6.8] consistent | 6.7 [6.3, 7.0] violation |
| ref.gst_obf | lognorm-n1000-b50-ab fixed1 | coverage | 93.7 [92.5, 94.7] consistent | 94.3 [93.8, 94.7] violation |
| ref.gst_obf | pareto_a25-n1000-b20-aa fixed1 | fpr_primary | 6.1 [5.4, 6.9] consistent | 6.1 [5.8, 6.4] violation |
| ref.gst_obf | pareto_a25-n1000-b20-ab fixed1 | coverage | 93.8 [92.6, 94.8] consistent | 92.9 [92.4, 93.4] violation |
| ref.gst_obf | pareto_a25-n20000-b20-aa daily14 | fpr_primary | 5.6 [4.9, 6.4] consistent | 5.8 [5.5, 6.1] violation |
| ref.gst_obf | pareto_a25-n20000-b20-aa stream140 | fpr_primary | 5.6 [4.9, 6.4] consistent | 5.8 [5.5, 6.2] violation |
| ref.msprt | lognorm-n20000-b20-aa stream140 | fpr_primary | 5.5 [4.8, 6.3] consistent | 5.6 [5.3, 6.0] violation |
| ref.naive_z | bern_p10-n1000-b20-aa fixed1 | fpr_primary | 5.5 [4.8, 6.3] consistent | 5.5 [5.1, 5.8] violation |
| ref.naive_z | bern_p10-n20000-b20-ab fixed1 | coverage | 93.9 [92.8, 94.9] consistent | 94.3 [93.8, 94.7] violation |
| ref.naive_z | lognorm-n1000-b50-ab fixed1 | coverage | 93.7 [92.5, 94.7] consistent | 94.3 [93.8, 94.7] violation |
| ref.naive_z | lognorm-n20000-b20-ab fixed1 | coverage | 94.5 [93.4, 95.4] consistent | 94.3 [93.9, 94.8] violation |
| ref.naive_z | pareto_a25-n1000-b20-ab fixed1 | coverage | 93.8 [92.6, 94.8] consistent | 92.9 [92.4, 93.4] violation |
| ref.naive_z | pareto_a25-n20000-b20-ab fixed1 | coverage | 94.8 [93.8, 95.8] consistent | 94.4 [93.9, 94.8] violation |
| spotify.gst | bern_p01-n1000-b20-ab fixed1 | coverage | 93.5 [92.3, 94.5] consistent | 92.7 [92.1, 93.2] violation |
| spotify.gst | bern_p10-n20000-b20-aa stream140 | fpr_primary | 6.0 [5.3, 6.8] consistent | 5.7 [5.3, 6.0] violation |
| spotify.gst | bern_p10-n20000-b20-ab fixed1 | coverage | 93.9 [92.8, 94.9] consistent | 94.3 [93.8, 94.7] violation |
| spotify.gst | lognorm-n1000-b20-ab daily14 | coverage | 93.5 [92.3, 94.5] consistent | 94.3 [93.8, 94.7] violation |
| spotify.gst | lognorm-n1000-b50-ab fixed1 | coverage | 93.7 [92.5, 94.7] consistent | 94.3 [93.8, 94.7] violation |
| spotify.gst | normal-n1000-b20-aa daily14 | fpr_primary | 5.1 [4.5, 5.9] consistent | 5.5 [5.2, 5.8] violation |
| spotify.gst | normal-n1000-b50-aa stream140 | fpr_primary | 5.8 [5.1, 6.5] consistent | 5.9 [5.5, 6.2] violation |
| spotify.gst | pareto_a25-n1000-b20-ab fixed1 | coverage | 93.8 [92.6, 94.8] consistent | 92.9 [92.4, 93.4] violation |
| spotify.gst | pareto_a25-n20000-b20-aa daily14 | fpr_primary | 5.9 [5.2, 6.7] consistent | 6.3 [5.9, 6.6] violation |
| spotify.z | bern_p10-n1000-b20-aa fixed1 | fpr_primary | 5.5 [4.8, 6.3] consistent | 5.5 [5.2, 5.8] violation |
| spotify.z | bern_p10-n20000-b20-ab fixed1 | coverage | 93.9 [92.8, 94.9] consistent | 94.3 [93.8, 94.7] violation |
| spotify.z | lognorm-n1000-b50-ab fixed1 | coverage | 93.7 [92.5, 94.7] consistent | 94.3 [93.8, 94.7] violation |
| spotify.z | lognorm-n20000-b20-ab fixed1 | coverage | 94.5 [93.4, 95.4] consistent | 94.3 [93.9, 94.8] violation |
| spotify.z | pareto_a25-n1000-b20-ab fixed1 | coverage | 93.8 [92.6, 94.8] consistent | 92.9 [92.4, 93.4] violation |
| spotify.z | pareto_a25-n20000-b20-ab fixed1 | coverage | 94.8 [93.8, 95.8] consistent | 94.4 [93.9, 94.8] violation |
