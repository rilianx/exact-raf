"""Task files for run_interleaved.sh (Table 4 of the paper).

    python ibex/make_tasks.py COCONUT_DIR      # writes ibex/tasks_synthetic.txt and ibex/tasks_coconut.txt

COCONUT_DIR holds the 298 instances of the optimization benchmark distributed with Ibex
(ibex-lib/benchs/optim: easy, medium, hard, others, unsolved, blowup and coconutbenchmark-library2,
listed with their original set in ibex/results/coconut_manifest.csv), copied into one folder.
Variants: af2 (S), splus (S+), sh (SH), exact (E); synq and synt contain no squares, so that S+ = S
is not run on them. Time limits: 300 s for synq, 60 s for the others."""
import csv, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
V4, V3 = ["af2", "splus", "sh", "exact"], ["af2", "sh", "exact"]
lines = []
for s, tl, vs in (("synq", 300, V3), ("synt", 60, V3), ("synp", 60, V4)):
    d = os.path.join(ROOT, "benchmarks", "ibex", s)
    for f in sorted(x for x in os.listdir(d) if x.endswith(".bch")):
        lines += [f"{s} {os.path.join(d, f)} {tl} {v}" for v in vs]
open(os.path.join(HERE, "tasks_synthetic.txt"), "w").write("\n".join(lines) + "\n")
if len(sys.argv) > 1:
    cdir = os.path.abspath(sys.argv[1])
    names = [r["instance"] for r in csv.DictReader(open(os.path.join(HERE, "results", "coconut_manifest.csv")))]
    missing = [n for n in names if not os.path.exists(os.path.join(cdir, n))]
    if missing: sys.exit(f"{len(missing)} instances missing in {cdir}, e.g. {missing[:3]}")
    lines = [f"all {os.path.join(cdir, n)} 60 {v}" for n in names for v in V4]
    open(os.path.join(HERE, "tasks_coconut.txt"), "w").write("\n".join(lines) + "\n")
print("written", "ibex/tasks_synthetic.txt" + (" and ibex/tasks_coconut.txt" if len(sys.argv) > 1 else ""))
