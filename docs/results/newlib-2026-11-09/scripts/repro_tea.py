import numpy as np
from scipy.stats import norm
from tea_tasting.aggr import Aggregates
from tea_tasting.metrics.proportion import Proportion

m = Proportion("x")  # defaults: method="auto" (norm at n >= 1000), equal_var=True
rng = np.random.default_rng(123)
R = 4000
for pc, pt, nc, nt in [(0.10, 0.30, 8000, 2000), (0.10, 0.30, 80000, 20000), (0.30, 0.10, 8000, 2000), (0.05, 0.10, 18000, 2000)]:
    sc, st = rng.binomial(nc, pc, R), rng.binomial(nt, pt, R)
    abs_hit = rel_hit = 0
    for a, b in zip(sc, st):
        r = m.analyze_aggregates(Aggregates(count_=nc, mean_={"x": a / nc}), Aggregates(count_=nt, mean_={"x": b / nt}))
        abs_hit += r.effect_size_ci_lower <= pt - pc <= r.effect_size_ci_upper
        rel_hit += r.rel_effect_size_ci_lower <= pt / pc - 1 <= r.rel_effect_size_ci_upper
    pp = (pc * nc + pt * nt) / (nc + nt)
    ratio = np.sqrt(pp * (1 - pp) * (1 / nc + 1 / nt) / (pc * (1 - pc) / nc + pt * (1 - pt) / nt))
    limit = 2 * norm.cdf(norm.ppf(0.975) * ratio) - 1
    print(f"{pc} -> {pt}, n {nc}/{nt}: absolute CI covers {abs_hit / R:.3f} (large-sample limit {limit:.3f}), relative CI covers {rel_hit / R:.3f}, nominal 0.950")
