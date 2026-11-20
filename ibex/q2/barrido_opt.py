"""Barrido de ibexopt-raak con la relajacion afin (--lr=art) y varias aritmeticas, en paralelo.

Uso: python3 barrido_opt.py <ibexopt-raak> <salida.csv> <tiempo_limite_s> <trabajos> <variantes,separadas> <bch>...
Escribe una fila por (instancia, variante): estado (SUCCESS si se alcanza la precision, TIME si se agota el tiempo),
celdas, f_lb, f_ub, cpu_s, y las estadisticas de QUAD2 si las hay. Reanuda: omite pares ya presentes en el CSV.
"""
import csv, os, re, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

exe, out_csv, tlim, jobs, variants, files = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), sys.argv[5].split(","), sys.argv[6:]
FIELDS = ["instancia", "conjunto", "variante", "estado", "celdas", "f_lb", "f_ub", "cpu_s", "quad2"]
done = set()
if os.path.exists(out_csv):
    done = {(r["instancia"], r["variante"]) for r in csv.DictReader(open(out_csv))}
else:
    with open(out_csv, "w", newline="") as fh: csv.DictWriter(fh, fieldnames=FIELDS).writeheader()

def run(task):
    f, v = task
    cmd = [exe, f, "--lr=art", f"--arith={v}", "-t", str(tlim)]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=tlim * 3 + 60)
        out = p.stdout + p.stderr
    except subprocess.TimeoutExpired:
        out = "TIMEOUT-EXTERNO"
    r = {"instancia": Path(f).stem, "conjunto": Path(f).parent.name, "variante": v}
    r["estado"] = "SUCCESS" if "[passed]" in out and "time limit" not in out else ("TIME" if "time limit" in out else "ERROR")
    m = re.search(r"f\* in\s*\[([^,]+),([^\]]+)\]", out); r["f_lb"], r["f_ub"] = (m.group(1), m.group(2)) if m else ("", "")
    m = re.search(r"number of cells:\s*(\d+)", out); r["celdas"] = m.group(1) if m else ""
    m = re.search(r"cpu time used:\s*([0-9.eE+-]+)s", out); r["cpu_s"] = m.group(1) if m else ""
    m = re.search(r"QUAD2 stats: (.*)", out); r["quad2"] = m.group(1) if m else ""
    if r["estado"] == "ERROR" and "infeasible" in out.lower(): r["estado"] = "INFEASIBLE"
    return r

tasks = [(f, v) for f in files for v in variants if (Path(f).stem, v) not in done]
print(len(tasks), "corridas", flush=True)
with ThreadPoolExecutor(max_workers=jobs) as ex:
    for r in ex.map(run, tasks):
        with open(out_csv, "a", newline="") as fh: csv.DictWriter(fh, fieldnames=FIELDS).writerow(r)
        print(r["instancia"], r["variante"], r["estado"], r["celdas"], r["cpu_s"], flush=True)
