"""Numerical check of Proposition 9: global search (differential evolution) of the ratio
wid[z]_SH / wid[z]_E over operands with n symbols. Run from the repository root."""
import sys, math; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.dirname(__import__('os').path.abspath(__file__))))
import numpy as np
from scipy.optimize import differential_evolution
from zcaa.zcaa import _uv_range, _sorted_generators, _uv_range_sorted
def ratio(p,n):
    x0,y0,xr,yr=p[0],p[1],abs(p[2]),abs(p[3]); a=p[4:4+n]; b=p[4+n:4+2*n]
    rlo,rhi=_uv_range(list(zip(a,b))+[(xr,0.),(0.,yr)]); A=abs(y0)*xr+abs(x0)*yr
    sh=rhi-rlo+2*A
    sg=_sorted_generators(zip(a,b)); lo,hi=math.inf,-math.inf
    for s in (-1,1):
        for t in (-1,1):
            l,h=_uv_range_sorted(sg,xr*s,yr*t); k=x0*yr*t+y0*xr*s; lo=min(lo,l+k); hi=max(hi,h+k)
    ex=hi-lo
    return -(sh/ex) if ex>1e-9 else 0
for n in (1,2,3,5):
    best=0
    for seed in range(6):
        r=differential_evolution(ratio,[(-3,3)]*(4+2*n),args=(n,),seed=seed,maxiter=600,popsize=30,tol=1e-12,polish=True)
        best=max(best,-r.fun)
        if -r.fun>1.5+1e-9: print("  >1.5 found:",-r.fun,np.round(r.x,4))
    print(f"n={n}: max ratio {best:.6f}")
