"""tea-tasting Proportion CI with default equal_var=True: coverage of the difference when the arms really differ.
Independent seed; real engine calls. Run: <venv python> -I pk_tea_prop_ci.py"""
import sys, warnings, json
warnings.filterwarnings("ignore")
import numpy as np
from scipy.stats import beta
from tea_tasting.aggr import Aggregates
from tea_tasting.metrics.proportion import Proportion
def cp(k, n): return (beta.ppf(.025, k, n-k+1) if k else 0.0, beta.ppf(.975, k+1, n-k) if k < n else 1.0)
ms = {"pooled(default)": Proportion("x", method="norm"), "unpooled": Proportion("x", method="norm", equal_var=False)}
def ag(n, s):
    m = s / n; return Aggregates(count_=n, mean_={"x": m}, var_={"x": m*(1-m)*n/(n-1)})
R = 4000; rng = np.random.default_rng(20261109)
grid = [(0.10,0.20,8000,2000),(0.10,0.30,8000,2000),(0.05,0.10,18000,2000),(0.30,0.60,18000,2000),(0.10,0.30,80000,20000),(0.05,0.06,8000,2000),(0.10,0.12,8000,2000),(0.10,0.30,5000,5000),(0.30,0.10,8000,2000)]
print("pA pB nA nB method  cov_abs [95% CP] cov_rel  | exact_se_ratio pooled/true")
for pa, pb, na, nb in grid:
    sa = rng.binomial(na, pa, R); sb = rng.binomial(nb, pb, R)
    pp = (pa*na+pb*nb)/(na+nb); ratio = np.sqrt(pp*(1-pp)*(1/na+1/nb) / (pa*(1-pa)/na + pb*(1-pb)/nb))
    for nm, m in ms.items():
        ca = cr = 0
        for a, b in zip(sa, sb):
            r = m.analyze_aggregates(ag(na, float(a)), ag(nb, float(b)))
            ca += r.effect_size_ci_lower <= pb-pa <= r.effect_size_ci_upper
            cr += r.rel_effect_size_ci_lower <= pb/pa-1 <= r.rel_effect_size_ci_upper
        l, u = cp(ca, R)
        print(f"{pa} {pb} {na} {nb} {nm:16s} {ca/R:.4f} [{l:.4f},{u:.4f}] rel {cr/R:.4f} | {ratio:.3f}", flush=True)
