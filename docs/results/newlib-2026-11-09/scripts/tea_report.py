import os, sys,json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "..", "..", "..", "src"))  # this repository's peekaudit
from peekaudit.stats import holm
import numpy as np
from scipy.stats import binom
rs=[json.loads(l) for l in open(sys.argv[1])]
aa=[r for r in rs if r['unit'].endswith('-aa') and r['schedule']=='fixed1']
adj=holm([r['p'] for r in aa])
print('fixed1 A/A cells',len(aa),'violate after Holm',sum(a<.05 for a in adj))
for r,a in zip(aa,adj):
    if a<.05: print(' ',r['unit'],r['method'],f"{100*r['rate']:.1f}% [{100*r['ci'][0]:.1f},{100*r['ci'][1]:.1f}] holm {a:.2g} nan={r['nan']}")
print('peeked daily14 (informational)')
for r in rs:
    if r['schedule']=='daily14': print(' ',r['unit'],r['method'],f"{100*r['rate']:.1f}% [{100*r['ci'][0]:.1f},{100*r['ci'][1]:.1f}]",r['reps'])
ab=[r for r in rs if r['unit'].endswith('-ab')]
from scipy.stats import binomtest
ps=[];L=[]
for r in ab:
    for k in ('cov_abs','cov_rel'):
        n=r['reps'];c=round(r[k]*n);p=binom.cdf(c,n,.95);ps.append(p);L.append((r['unit'],r['method'],k,r[k]))
adj=holm(ps);print('A/B coverage cells',len(L),'below 95% after Holm',sum(a<.05 for a in adj))
for l,a in zip(L,adj):
    if a<.05: print(' ',l[0],l[1],l[2],f"{l[3]:.4f} holm {a:.2g}")
print('errs',{r['unit']+r['method']:r['errs'] for r in rs if r['err']}, 'nan cells',[(r['unit'],r['method'],r['schedule'],r['nan']) for r in rs if r['nan']][:20])
