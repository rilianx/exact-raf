import csv, glob, math, collections
RES="/home/user/ibex-local/results/interleaved_all"
man={r["instance"][:-4]:r["set"] for r in csv.DictReader(open("/home/user/x2taylor/benchs/optim/all/manifest.csv"))}
D=collections.defaultdict(dict)
for f in glob.glob(RES+"/rows/*.csv"):
    for r in csv.reader(open(f)):
        st,inst=r[0].split("/")
        if st!="all": continue
        D[inst][r[1].split("_")[1]]=dict(st=r[2],c=int(r[3] or 0),t=float(r[6] or "nan"))
V=["af2","splus","sh","exact"]
gm=lambda x: math.exp(sum(map(math.log,x))/len(x)) if x else float("nan")
ok=lambda r:r["st"]=="SUCCESS"
groups=collections.defaultdict(list)
for i in D:
    if all(v in D[i] for v in V): groups[man.get(i,"?")].append(i); groups["ALL"].append(i)
for g,inst in sorted(groups.items()):
    solved=[sum(ok(D[i][v]) for i in inst) for v in V]
    A=[i for i in inst if all(ok(D[i][v]) for v in V)]
    rc=lambda v: gm([max(D[i][v]["c"],1)/max(D[i]["splus"]["c"],1) for i in A])
    rt=lambda v: gm([max(D[i][v]["t"],1e-2)/max(D[i]["splus"]["t"],1e-2) for i in A])
    es=[max(D[i]["exact"]["c"],1)/max(D[i]["sh"]["c"],1) for i in A]
    print(f"{g:9s} n={len(inst):3d} solved S/S+/SH/E={solved} common={len(A)} nodes S={rc('af2'):.3f} SH={rc('sh'):.3f} E={rc('exact'):.3f} cpu SH={rt('sh'):.3f} E={rt('exact'):.3f} E/SH={gm(es):.3f} {sum(x<.99 for x in es)}/{sum(x>1.01 for x in es)}")
print("instancias con 4 variantes:",len(groups["ALL"]),"de",len(D))
for i in groups["ALL"]:
    s=[ok(D[i][v]) for v in V]
    if len(set(s))>1: print("diff", i, man.get(i), [(v,D[i][v]["st"],D[i][v]["c"],round(D[i][v]["t"],1)) for v in V])
print("---- |E/SH| > 1.1 o < 0.9")
for i in groups["ALL"]:
    if ok(D[i]["sh"]) and ok(D[i]["exact"]):
        r=max(D[i]["exact"]["c"],1)/max(D[i]["sh"]["c"],1)
        if r<0.9 or r>1.1: print(f"{i:22s} {man.get(i):10s} S+={D[i]['splus']['c']} SH={D[i]['sh']['c']} E={D[i]['exact']['c']} r={r:.3f}")
