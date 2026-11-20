"""Barrido de ibexsolve-raak (--ctc=default-xn+affine: contractor por omision con la relajacion afin en lugar de
X-Newton) con varias aritmeticas, en paralelo.
Uso: python3 barrido_solve.py <ibexsolve-raak> <salida.csv> <tiempo_s> <trabajos> <variantes> <bch>...
Fila por (instancia, variante) con la linea RAAK-SOLVE (estado, celdas, tiempo, soluciones, ...). Reanuda."""
import csv, os, re, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
exe, out_csv, tlim, jobs, variants, files = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), sys.argv[5].split(","), sys.argv[6:]
FIELDS = ["instancia", "variante", "estado", "celdas", "tiempo", "solution", "boundary", "unknown", "pending", "q2_evals", "q2_failed", "q2_unsupported"]
done = {(r["instancia"], r["variante"]) for r in csv.DictReader(open(out_csv))} if os.path.exists(out_csv) else set()
if not done:
    with open(out_csv, "w", newline="") as fh: csv.DictWriter(fh, fieldnames=FIELDS).writeheader()
def run(task):
    f, v = task
    try:
        out = subprocess.run([exe, f, "--ctc=default-xn+affine", f"--arith={v}", "-t", str(tlim), "-q"],
                             capture_output=True, text=True, timeout=tlim * 3 + 60).stdout
    except subprocess.TimeoutExpired:
        out = ""
    m = re.search(r"RAAK-SOLVE (.*)", out)
    kv = dict(re.findall(r"(\w+)=(\S+)", m.group(1))) if m else {"status": "ERROR"}
    return {"instancia": Path(f).stem, "variante": v, "estado": kv.get("status"), "celdas": kv.get("cells", ""),
            "tiempo": kv.get("time", ""), "solution": kv.get("solution", ""), "boundary": kv.get("boundary", ""),
            "unknown": kv.get("unknown", ""), "pending": kv.get("pending", ""), "q2_evals": kv.get("q2_evals", ""),
            "q2_failed": kv.get("q2_failed", ""), "q2_unsupported": kv.get("q2_unsupported", "")}
tasks = [(f, v) for f in files for v in variants if (Path(f).stem, v) not in done]
print(len(tasks), "corridas", flush=True)
with ThreadPoolExecutor(max_workers=jobs) as ex:
    for r in ex.map(run, tasks):
        with open(out_csv, "a", newline="") as fh: csv.DictWriter(fh, fieldnames=FIELDS).writerow(r)
        print(r["instancia"], r["variante"], r["estado"], r["celdas"], r["tiempo"], flush=True)
