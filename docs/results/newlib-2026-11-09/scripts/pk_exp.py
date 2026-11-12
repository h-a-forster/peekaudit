"""peekaudit's standard (FPR <= alpha with exact CP and Holm) on `expectation` public API, streaming one observation at a time.
A/A = data at the null value. Run: <venv python> -I pk_exp.py [reps]. One JSON line per finished cell, skipped on rerun."""
import sys, os, json, time, warnings
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "..", "..", "..", "src"))  # this repository's peekaudit
warnings.filterwarnings("ignore")
import numpy as np
from peekaudit.stats import clopper_pearson, binom_sf_test
from expectation.seqtest.sequential_e_testing import SequentialTesting
REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 45.0  # seconds per cell; rep seeds are fixed, so a rerun with --reps=R reproduces the first R
D = os.path.join(BASE, "expectation"); OUT = os.path.join(D, "cells.jsonl"); LOG = os.path.join(D, "progress.log")
done = {tuple(json.loads(l)["id"]) for l in open(OUT)} if os.path.exists(OUT) else set()
def log(s):
    with open(LOG, "a") as f: f.write(time.strftime("%H:%M:%S ") + s + "\n")
HOR = 500
def dgp(name, rng, n):
    return {"normal": lambda: rng.standard_normal(n), "normal_sd3": lambda: 3 * rng.standard_normal(n), "normal_sd01": lambda: 0.1 * rng.standard_normal(n),
            "uniform": lambda: rng.uniform(-1, 1, n), "t3": lambda: rng.standard_t(3, n),
            "bern50": lambda: (rng.random(n) < .5).astype(float), "bern10": lambda: (rng.random(n) < .1).astype(float), "bern01": lambda: (rng.random(n) < .01).astype(float),
            }[name]()
CELLS = []
for d in ("normal", "normal_sd3", "normal_sd01", "uniform", "t3"):
    for alt in ("two_sided", "greater"):
        CELLS.append(("mean/default", d, alt, dict(test_type="mean", null_value=0.0, alternative=alt), 0))
CELLS += [("mean/known_var", "normal", "two_sided", dict(test_type="mean", null_value=0.0, alternative="two_sided", known_variance=1.0), 0),
          ("mean/known_var", "normal", "greater", dict(test_type="mean", null_value=0.0, alternative="greater", known_variance=1.0), 0)]
for d, p in (("bern50", .5), ("bern10", .1), ("bern01", .01)):
    for alt in ("two_sided", "greater", "less"):
        CELLS.append(("proportion/default", d, alt, dict(test_type="proportion", null_value=p, alternative=alt), 0))
t0 = time.time()
for ci, (row, d, alt, kw, _) in enumerate(CELLS):
    for batch in (1, 10):
        cid = (row, d, alt, batch)
        if cid in done: continue
        t1 = time.time(); rej = 0; first = []; err = {}; cov = 0; ncov = 0
        ran = 0
        for rep in range(REPS):
            if time.time() - t1 > BUDGET: break
            ran += 1
            rng = np.random.default_rng([20261109, ci, batch, rep])
            x = dgp(d, rng, HOR); T = SequentialTesting(**kw); r_any = False; cov_ok = True; tt = None
            try:
                for i in range(0, HOR, batch):
                    res = T.update(x[i:i + batch])
                    if res.reject_null and not r_any: r_any = True; tt = i + batch
                    cb = res.confidence_bounds
                    if cb is not None and kw["test_type"] == "mean" and not (cb[0] <= 0.0 <= cb[1]): cov_ok = False
                    if cb is not None and kw["test_type"] == "proportion" and not (cb[0] <= kw["null_value"] <= cb[1]): cov_ok = False
            except Exception as e:
                k = type(e).__name__ + ":" + str(e)[:80]; err[k] = err.get(k, 0) + 1; continue
            rej += r_any; cov += cov_ok; ncov += 1
            if tt: first.append(tt)
        n = ran - sum(err.values())
        rec = dict(id=cid, row=row, dgp=d, alt=alt, batch=batch, reps=n, horizon=HOR, k=rej, rate=rej / max(n, 1), ci=clopper_pearson(rej, n) if n else None, p=binom_sf_test(rej, n, .05) if n else None,
                   cs_covered=cov / max(ncov, 1), median_first=float(np.median(first)) if first else None, err=err, secs=round(time.time() - t1, 1))
        with open(OUT, "a") as f: f.write(json.dumps(rec) + "\n")
        log(f"{row} {d} {alt} batch={batch} fpr={rec['rate']:.3f} cs_cov={rec['cs_covered']:.3f} {rec['secs']}s err={err}")
log(f"done {time.time()-t0:.0f}s")
