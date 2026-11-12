"""expectation 0.6.1 tiny standalone checks. Run: <venv python> -I repro_exp.py"""
import warnings, traceback; warnings.filterwarnings("ignore")
import numpy as np
from expectation.seqtest.sequential_e_testing import SequentialTesting
print("1 README import:")
try:
    from expectation.seqtest import SequentialTesting as S; print("  ok")
except Exception as e: print(" ", type(e).__name__, e)
print("2 proportion one observation:")
for p, alt in ((0.5, "two_sided"), (0.1, "greater")):
    T = SequentialTesting(test_type="proportion", null_value=p, alternative=alt)
    try:
        r = T.update([1.0]); print("  ok", r.e_value)
    except Exception as e:
        print(" ", p, alt, type(e).__name__, str(e)[:100]); tb = traceback.extract_tb(e.__traceback__); print("   at", [f"{f.name}:{f.lineno}" for f in tb][-4:])
T = SequentialTesting(test_type="proportion", null_value=0.5)
for i, b in enumerate([[1,0,1,1,0,0,1,0,1,1], [0,0,1,1,0,1,1,0,0,1]]):
    try: print("  batch ok", T.update(np.array(b, float)).e_value)
    except Exception as e: print("  batch ERR", i, str(e)[:80])
print("3 mean, default, one at a time, N(0,1), 500 obs, 400 reps; known_variance=1 as control")
for kw, name in ((dict(), "default"), (dict(known_variance=1.0), "known_var")):
    rej = 0; first = []
    for rep in range(400):
        rng = np.random.default_rng([99, rep]); T = SequentialTesting(test_type="mean", null_value=0.0, **kw); hit = None
        for i, x in enumerate(rng.standard_normal(500)):
            if T.update(float(x)).reject_null: hit = i + 1; break
        rej += hit is not None
        if hit: first.append(hit)
    from scipy.stats import beta
    print(f"  {name}: FPR {rej/400:.3f} [{beta.ppf(.025,rej,401-rej) if rej else 0:.3f},{beta.ppf(.975,rej+1,400-rej):.3f}] median first rejection n={np.median(first) if first else None}")
print("4 first rejection anatomy: x=[0.1, 0.11] (sd small)")
T = SequentialTesting(test_type="mean", null_value=0.0)
for x in (0.1, 0.11, 0.09):
    r = T.update(x); print("  ", x, "e-process", round(r.e_process.cumulative_value, 3), "reject", r.reject_null)
