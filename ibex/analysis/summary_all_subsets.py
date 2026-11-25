import re,glob,os
exec(open(__import__("os").path.join(__import__("os").path.dirname(__file__),"summary_all_by_set.py")).read().split("groups=")[0])
st={}
for f in glob.glob("/home/user/raa-affine-arithmetic/ibex/results/sdiag/*.txt"):
    t=open(f).read(); n=os.path.basename(f)[:-4]
    m=dict(re.findall(r"(\w+)=(\d+)",t))
    st[n]={k:int(v) for k,v in m.items()} if m else None
none=[n for n in st if st[n] is None]; print("sin stats:",len(none),none[:10])
def frac(n,k): s=st.get(n); return s[k]/s["productos"] if s and s["productos"] else 0
for thr in [0, 0.01, 0.1]:
  for key in ["s_mas_ancho","s_mas_ancho_10"]:
    sel=[i for i in D if all(v in D[i] for v in V) and frac(i,key)>thr]
    A=[i for i in sel if all(ok(D[i][v]) for v in V)]
    solved=[sum(ok(D[i][v]) for i in sel) for v in V]
    rc=lambda v: gm([max(D[i][v]["c"],1)/max(D[i]["splus"]["c"],1) for i in A])
    rt=lambda v: gm([max(D[i][v]["t"],1e-2)/max(D[i]["splus"]["t"],1e-2) for i in A])
    es=[max(D[i]["exact"]["c"],1)/max(D[i]["sh"]["c"],1) for i in A]
    print(f"{key} frac>{thr}: n={len(sel)} solved={solved} common={len(A)} nodes S={rc('af2'):.3f} SH={rc('sh'):.3f} E={rc('exact'):.3f} cpu SH={rt('sh'):.3f} E={rt('exact'):.3f} E/SH={gm(es):.3f} {sum(x<.99 for x in es)}/{sum(x>1.01 for x in es)}")
sel=[i for i in D if frac(i,"s_mas_ancho")>0]
import collections; print(collections.Counter(man[i] for i in sel))
cp=[i for i in D if st.get(i) and st[i]["acoplados"]>0]; print("con productos acoplados:",len(cp))
for thr in [0,0.01,0.1]:
    sel=[i for i in D if frac(i,"acoplados")>thr]
    A=[i for i in sel if all(ok(D[i][v]) for v in V)]
    es=[max(D[i]["exact"]["c"],1)/max(D[i]["sh"]["c"],1) for i in A]
    et=[max(D[i]["exact"]["t"],1e-2)/max(D[i]["sh"]["t"],1e-2) for i in A]
    print(f"acoplados>{thr}: n={len(sel)} common={len(A)} E/SH nodes={gm(es):.3f} cpu={gm(et):.3f} {sum(x<.99 for x in es)}/{sum(x>1.01 for x in es)} solved SH/E={sum(ok(D[i]['sh']) for i in sel)}/{sum(ok(D[i]['exact']) for i in sel)}")
