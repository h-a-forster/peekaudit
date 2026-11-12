"""expectation 0.6.1: parallel engine, k-sample Bernoulli, confidence sequences under the null / true mean.
Run: <venv python> -I pk_exp2.py. Appends to expectation/cells2.jsonl, skip on rerun."""
import sys, os, json, time, warnings
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "..", "..", "..", "src"))  # this repository's peekaudit
warnings.filterwarnings("ignore")
import numpy as np
from peekaudit.stats import clopper_pearson, binom_sf_test
from expectation.par_seqtest import ParallelSequentialTest, ParallelTestConfig
from expectation.ksample.ksample_test import KSampleSequentialTest
from expectation.ksample.config import KSampleConfig
from expectation.confseq.confidencesequence import ConfidenceSequence, EmpiricalBernsteinConfidenceSequence
from expectation.confseq.confidenceconfig import ConfidenceSequenceConfig, EmpiricalBernsteinConfig
D = os.path.join(BASE, "expectation"); OUT = os.path.join(D, "cells2.jsonl")
done = {json.loads(l)["id"] for l in open(OUT)} if os.path.exists(OUT) else set()
def emit(id_, row, k, n, kind, **kw):
    rec = dict(id=id_, row=row, k=int(k), n=int(n), kind=kind, rate=k / n, ci=clopper_pearson(int(k), int(n)), p=binom_sf_test(int(k), int(n), 0.05), **kw)
    with open(OUT, "a") as f: f.write(json.dumps(rec) + "\n")
    print(id_, f"{k}/{n}={k/n:.4f}", [round(x, 4) for x in rec["ci"]], kw, flush=True)
# 1 parallel engine, documented config, N(0,1) with variance=1 (model correct), T=500
HOR = 500
for alt, extra, sd, tag in (("greater", {}, 1.0, "documented"), ("two_sided", {}, 1.0, "two_sided_known"), ("greater", dict(variance_mode="empirical"), 1.0, "empirical_var"), ("greater", dict(variance_mode="empirical"), 3.0, "empirical_var_sd3")):
    id_ = f"parallel/{tag}"
    if id_ in done: continue
    NT = 100_000; rng = np.random.default_rng(1)
    kw = dict(null_values=0.0, variance=1.0) if "variance_mode" not in extra else dict(null_values=0.0)
    e = ParallelSequentialTest(config=ParallelTestConfig(n_tests=NT, alpha=0.05, alternative=alt, combiner="empirically_adaptive", **extra), **kw)
    for t in range(HOR): e.step(rng.standard_normal(NT) * sd)
    emit(id_, "parallel", e.rejected().sum(), NT, "fpr", horizon=HOR)
# 2 k-sample Bernoulli, A/A, default config, peek each block, 10 per arm per block
for p in (0.5, 0.1, 0.01):
    for bs in (1, 10):
        id_ = f"ksample/default/p{p}/block{bs}"
        if id_ in done: continue
        R = 200; T = 200 if bs == 10 else 500; rej = 0; t1 = time.time()
        for rep in range(R):
            rng = np.random.default_rng([2, int(p * 100), bs, rep]); test = KSampleSequentialTest(KSampleConfig(k=2)); r_any = False
            for t in range(T):
                res = test.update({0: (rng.random(bs) < p).astype(int), 1: (rng.random(bs) < p).astype(int)})
                if res.reject_null: r_any = True; break
            rej += r_any
        emit(id_, "ksample", rej, R, "fpr", secs=round(time.time() - t1))
# 3 confidence sequences: simultaneous coverage to T
R = 300; T = 500
for name, lb, ub, gen, mu in (("eb[0,1]/beta(7,3)", 0.0, 1.0, lambda g, n: g.beta(7, 3, n), 0.7), ("eb[0,1]/uniform", 0.0, 1.0, lambda g, n: g.uniform(0, 1, n), 0.5),
                              ("eb[0,2]/uniform", 0.0, 2.0, lambda g, n: g.uniform(0, 2, n), 1.0), ("eb[0,10]/uniform", 0.0, 10.0, lambda g, n: g.uniform(0, 10, n), 5.0), ("eb[-1,1]/uniform", -1.0, 1.0, lambda g, n: g.uniform(-1, 1, n), 0.0)):
    id_ = f"cs/{name}"
    if id_ in done: continue
    miss = 0; first = []; t1 = time.time(); errs = 0; ran = 0
    for rep in range(R):
        g = np.random.default_rng([3, rep, int(ub * 10)]); x = gen(g, T)
        cs = EmpiricalBernsteinConfidenceSequence(config=EmpiricalBernsteinConfig(alpha=0.05, lower_bound=lb, upper_bound=ub))
        bad = None
        try:
            for i in range(0, T, 5):
                res = cs.update(x[i:i + 5])
                if not (res.lower <= mu <= res.upper): bad = i + 5; break
        except RuntimeError as ex:
            errs += 1; continue
        ran += 1; miss += bad is not None
        if bad: first.append(bad)
    emit(id_, "cs_empirical_bernstein", miss, max(ran, 1), "miscover", errs=errs, first_median=float(np.median(first)) if first else None, last=[res.lower, res.upper], secs=round(time.time() - t1))
for name, gen, mu in (("base/normal", lambda g, n: g.standard_normal(n), 0.0), ("base/normal_sd3", lambda g, n: 3 * g.standard_normal(n), 0.0), ("base/bern50", lambda g, n: (g.random(n) < .5).astype(float), 0.5)):
    id_ = f"cs/{name}"
    if id_ in done: continue
    miss = 0; t1 = time.time(); errs = 0; ran = 0
    for rep in range(R):
        g = np.random.default_rng([4, rep]); x = gen(g, T); cs = ConfidenceSequence(config=ConfidenceSequenceConfig(alpha=0.05)); bad = False
        try:
            for i in range(0, T, 5):
                res = cs.update(x[i:i + 5])
                if not (res.lower <= mu <= res.upper): bad = True; break
        except RuntimeError:
            errs += 1; continue
        ran += 1; miss += bad
    emit(id_, "cs_base", miss, max(ran, 1), "miscover", errs=errs, secs=round(time.time() - t1))
