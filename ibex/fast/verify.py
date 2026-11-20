from fractions import Fraction as F
from functools import cmp_to_key
import sys
def rng(gens, c0, c1):
    g=[(-u,-v) if (u<0 or (u==0 and v<0)) else (u,v) for u,v in gens if u or v]
    if not g: return c0*c1, c0*c1
    g.sort(key=cmp_to_key(lambda a,b: -1 if a[0]*b[1]-a[1]*b[0]>0 else (1 if a[0]*b[1]-a[1]*b[0]<0 else 0)))
    pu=c0-sum(t[0] for t in g); pv=c1-sum(t[1] for t in g); lo=hi=pu*pv
    for s in (2,-2):
        for gu,gv in g:
            du,dv=s*gu,s*gv; qu,qv=pu+du,pv+dv; lo=min(lo,qu*qv); hi=max(hi,qu*qv)
            A=du*dv; B=pu*dv+pv*du
            if A!=0 and 0< -B/(2*A) <1:
                val=pu*pv-B*B/(4*A); lo=min(lo,val); hi=max(hi,val)
            pu,pv=qu,qv
    return lo,hi
bad=0; badI=0; slack=[]
for line in open(sys.argv[1]):
    L,Rr=line.split("|"); t=L.split(); e,n=int(t[0]),int(t[1])
    x0,ex,y0,ey=[F(float.fromhex(z)) for z in t[2:6]]
    v=[F(float.fromhex(z)) for z in t[6:6+2*n]]; gens=list(zip(v[0::2],v[1::2]))
    fl,fh,il,ih=[F(float.fromhex(z)) for z in Rr.split()]
    if e==1 and min(abs(x0)*ey,abs(y0)*ex)>0:
        lo,hi=F(10**9),F(-10**9)
        for s in (-1,1):
            for tt in (-1,1):
                a,b=rng(gens,ex*s,ey*tt); k=x0*ey*tt+y0*ex*s; lo=min(lo,a+k); hi=max(hi,b+k)
    else:
        a,b=rng(gens+[(ex,F(0)),(F(0),ey)],F(0),F(0)); l=abs(y0)*ex+abs(x0)*ey; lo,hi=a-l,b+l
    if fl>lo or fh<hi: bad+=1
    if il>lo or ih<hi: badI+=1
    w=hi-lo
    if w>0: slack.append(float(((lo-fl)+(fh-hi))/w))
print("casos",len(slack),"violaciones fast",bad,"violaciones interval",badI,"holgura relativa max",max(slack),"mediana",sorted(slack)[len(slack)//2])
