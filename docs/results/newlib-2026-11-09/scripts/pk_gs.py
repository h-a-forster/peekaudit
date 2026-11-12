"""gsdesign (grid + density recursion) against peekaudit.gst (validated on published tables) and exact values.
Deterministic. Run: <venv python> -I pk_gs.py"""
import os, sys, warnings, math
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE, "..", "..", "..", "src"))  # this repository's peekaudit
import numpy as np
from scipy.stats import norm
from peekaudit.gst import crossing_probability, lan_demets_bounds, spend_obf, spend_pocock
from gsdesign import gridpts, h1, hupdate
def cross(bounds, infos, theta=0.0, r=18):
    """P(first |Z_k|>=b_k) using gsdesign's recursion: the upper crossing via Simpson on [b, inf) is the documented use."""
    out = []; prev = None
    for k, (b, I) in enumerate(zip(bounds, infos)):
        args = dict(r=r, theta=theta, info=I)
        if prev is None:
            up = h1(a=b, **args); lo = h1(b=-b, **args); cont = h1(a=-b, b=b, **args)
        else:
            up = hupdate(r, theta, I, b, np.inf, theta, infos[k - 1], prev); lo = hupdate(r, theta, I, -np.inf, -b, theta, infos[k - 1], prev); cont = hupdate(r, theta, I, -b, b, theta, infos[k - 1], prev)
        out.append(up[2].sum() + lo[2].sum()); prev = cont
    return np.array(out)
def rel(a, b): return abs(a - b) / max(abs(b), 1e-300)
print("gridpts: weights integrate the density"); 
for r in (2, 5, 18, 50):
    for mu, a, b in ((0, -np.inf, np.inf), (1.3, -np.inf, np.inf), (0, 1.96, np.inf), (0, -1.96, 1.96), (2, -1, 0.5), (0, 8, np.inf)):
        z, w = gridpts(r, mu, a, b); val = (w * norm.pdf(z - mu)).sum(); ex = norm.cdf(b - mu) - norm.cdf(a - mu)
        print(f" r={r} mu={mu} [{a},{b}] sum={val:.12g} exact={ex:.12g} relerr={rel(val, ex):.2e} sorted={bool(np.all(np.diff(z) > 0))} npts={z.size}")
print("h1 + hupdate vs peekaudit.gst crossing_probability (ngrid 2001)")
worst = 0
for K in (2, 3, 5, 10, 20, 50):
    t = np.arange(1, K + 1) / K
    for label, b in (("z1.96", np.full(K, 1.96)), ("obf", lan_demets_bounds(t, spend_obf(0.025))), ("pocock", lan_demets_bounds(t, spend_pocock(0.025)))):
        ref = crossing_probability(t, b, ngrid=2001); got = cross(list(b), list(t * 100.0))
        e = float(np.max(np.abs(got - ref))); worst = max(worst, e)
        print(f" K={K} {label}: total gsdesign {got.sum():.8f} ref {ref.sum():.8f} max|diff| {e:.2e}")
print("worst", worst)
print("exact checks: peeked z 1.96 K=2 exact 0.083 (Armitage), K=5 0.142, K=10 0.193")
for K, ex in ((2, .083), (5, .142), (10, .193)):
    t = np.arange(1, K + 1) / K; print(K, cross([1.96] * K, list(t * 100.0)).sum(), ex)
print("power check: K=1 theta: P(Z>=1.96) with drift matches normal")
for th in (0.1, 0.3):
    z = h1(r=18, theta=th, info=100, a=1.96); print(th, z[2].sum(), norm.sf(1.96 - th * 10))
print("two looks with drift vs Monte Carlo seed 1, 4e6")
rng = np.random.default_rng(1); I = np.array([50., 100.]); th = .2; n = 4_000_000
s1 = rng.normal(th * I[0], np.sqrt(I[0]), n); s2 = s1 + rng.normal(th * (I[1] - I[0]), np.sqrt(I[1] - I[0]), n)
z1, z2 = s1 / np.sqrt(I[0]), s2 / np.sqrt(I[1]); hit = (z1 >= 2.5) | ((np.abs(z1) < 2.5) & (z2 >= 2.0)); pmc = hit.mean(); se = math.sqrt(pmc * (1 - pmc) / n)
a = h1(18, th, I[0], a=2.5)[2].sum(); c = h1(18, th, I[0], a=-2.5, b=2.5); b = hupdate(18, th, I[1], 2.0, np.inf, th, I[0], c)[2].sum()
print(f" recursion {a + b:.6f} mc {pmc:.6f} +- {1.96 * se:.6f}")
print("edge input: r<2, a>=b, info<=0, info order, gm1 bad, a=b tiny, b below grid")
for f, args in ((gridpts, dict(r=1)), (gridpts, dict(r=18, a=1, b=1)), (h1, dict(info=0)), (hupdate, dict(r=18, theta=0, info=1, a=-1, b=1, theta_prev=0, info_prev=1, gm1=h1()))):
    try: f(**args); print(" no error", f.__name__, args if f is gridpts else "")
    except Exception as e: print(" ", type(e).__name__, str(e)[:60])
for a, b in ((20, np.inf), (-np.inf, -20), (0, 1e-9), (5.0, 5.0000001)):
    try:
        z, w = gridpts(18, 0, a, b); val = (w * norm.pdf(z)).sum(); ex = norm.cdf(b) - norm.cdf(a); print(f" [{a},{b}] npts={z.size} sum={val:.6g} exact={ex:.6g}")
    except Exception as e: print(" ERR", a, b, type(e).__name__, str(e)[:80])
