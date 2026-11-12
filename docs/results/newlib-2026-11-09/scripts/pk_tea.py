"""peekaudit on tea-tasting: fixed-horizon FPR and coverage, peeked FPR (informational). Real engine calls, no mirror.
Run: <venv python> -I pk_tea.py [reps]. One JSON line per finished cell in tea-tasting/cells.jsonl, skipped on rerun."""
import sys, os, json, time, warnings
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "..", "..", "..", "src"))  # this repository's peekaudit
warnings.filterwarnings("ignore")
import numpy as np
from peekaudit.scenarios import Unit, SCHEDULES, METRICS
from peekaudit.simulate import generate
from peekaudit.stats import clopper_pearson, binom_sf_test
from tea_tasting.aggr import Aggregates
from tea_tasting.metrics.mean import Mean
from tea_tasting.metrics.proportion import Proportion
REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
D = os.path.join(BASE, "tea-tasting"); OUT = os.path.join(D, "cells.jsonl"); LOG = os.path.join(D, "progress.log")
done = {tuple(json.loads(l)["id"]) for l in open(OUT)} if os.path.exists(OUT) else set()
def log(s):
    with open(LOG, "a") as f: f.write(time.strftime("%H:%M:%S ") + s + "\n")
def aggs(n, s, q):
    mean = s / n; var = (q - n * mean * mean) / (n - 1) if n > 1 else 0.0
    return Aggregates(count_=n, mean_={"x": mean}, var_={"x": max(var, 0.0)})
def bin_aggs(n, s):
    mean = s / n
    return Aggregates(count_=n, mean_={"x": mean}, var_={"x": mean * (1 - mean) * n / (n - 1)})
M = {"welch_t": Mean("x"), "z": Mean("x", use_t=False), "student_t": Mean("x", equal_var=True)}
P = {"prop_norm_pooled": Proportion("x", method="norm"), "prop_norm_unpooled": Proportion("x", method="norm", equal_var=False)}
def one(name, nA, sA, qA, nB, sB, qB):
    if name in M: r = M[name].analyze_aggregates(aggs(nA, sA, qA), aggs(nB, sB, qB))
    else: r = P[name].analyze_aggregates(bin_aggs(nA, sA), bin_aggs(nB, sB))
    return r
metrics = ("normal", "bern_p10", "bern_p01", "lognorm", "pareto_a25")
t0 = time.time()
for hyp in ("aa", "ab"):
  for mid in metrics:
    for nmax in (1000, 20000):
      for pi in (0.5, 0.2):
        unit = Unit(mid, nmax, pi, hyp)
        names = list(M) + (list(P) if METRICS[mid].binary else [])
        scheds = ["fixed1"] + (["daily14"] if hyp == "aa" and mid in ("normal", "bern_p10") else [])
        todo = [(n, s) for n in names for s in scheds if (unit.id, n, s) not in done]
        if not todo: continue
        data = generate(unit, REPS)
        for name, sched in todo:
            cols = SCHEDULES[sched]; R = REPS if sched == "fixed1" else min(REPS, 1000)
            rej = np.zeros(R, bool); cov_abs = np.zeros(R, float); cov_rel = np.zeros(R, float); nan = 0; err = 0; t1 = time.time(); errs = {}
            for i in range(R):
                for j in cols:
                    try:
                        r = one(name, int(data.n_a[j]), float(data.sum_a[i, j]), float(data.sumsq_a[i, j]), int(data.n_b[j]), float(data.sum_b[i, j]), float(data.sumsq_b[i, j]))
                    except Exception as e:
                        err += 1; errs[type(e).__name__ + ":" + str(e)[:60]] = errs.get(type(e).__name__ + ":" + str(e)[:60], 0) + 1; continue
                    if not np.isfinite(r.pvalue): nan += 1; continue
                    if r.pvalue < 0.05: rej[i] = True
                if sched == "fixed1" and hyp == "ab":
                    try:
                        cov_abs[i] = r.effect_size_ci_lower <= unit.true_effect("absolute") <= r.effect_size_ci_upper
                        cov_rel[i] = r.rel_effect_size_ci_lower <= unit.true_effect("relative") <= r.rel_effect_size_ci_upper
                    except Exception: pass
            k = int(rej.sum())
            rec = dict(id=(unit.id, name, sched), unit=unit.id, method=name, schedule=sched, reps=R, k=k, rate=k / R, ci=clopper_pearson(k, R), p=binom_sf_test(k, R, 0.05),
                       nan=nan, err=err, errs=errs, cov_abs=float(cov_abs.mean()) if hyp == "ab" else None, cov_rel=float(cov_rel.mean()) if hyp == "ab" else None, secs=round(time.time() - t1, 1))
            with open(OUT, "a") as f: f.write(json.dumps(rec) + "\n")
            log(f"{unit.id} {name} {sched} rate={k/R:.4f} cov={rec['cov_abs']} rel={rec['cov_rel']} err={err} nan={nan} {rec['secs']}s")
log(f"done {time.time()-t0:.0f}s")
