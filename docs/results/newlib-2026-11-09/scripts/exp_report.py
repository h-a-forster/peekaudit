import os, sys,json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "..", "..", "..", "src"))  # this repository's peekaudit
from peekaudit.stats import holm
E = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'expectation')
a=[json.loads(l) for l in open(os.path.join(E, 'cells.jsonl'))]; b=[json.loads(l) for l in open(os.path.join(E, 'cells2.jsonl'))]
rows=[]
for r in a:
    if r['reps']>0: rows.append(('%s %s %s b%d'%(r['row'],r['dgp'],r['alt'],r['batch']),r['k'],r['reps'],r['p'],r['ci']))
for r in b: rows.append((r['id'],r['k'],r['n'],r['p'],r['ci']))
adj=holm([x[3] for x in rows]); print(len(rows),'cells; violating after Holm',sum(x<.05 for x in adj))
for x,h in zip(rows,adj):
    if h<.05: print(' ',x[0],f"{x[1]}/{x[2]}={x[1]/x[2]:.3f} [{x[4][0]:.3f},{x[4][1]:.3f}] holm {h:.2g}")
print('errors by row:')
for r in a:
    if sum(r['err'].values()): print(' ',r['row'],r['dgp'],r['alt'],r['batch'],sum(r['err'].values()),'/',r['reps']+sum(r['err'].values()),list(r['err'])[0][:50])
