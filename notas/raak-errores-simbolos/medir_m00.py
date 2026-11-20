"""M00 = fAF2 de Ibex por omision (formula (5) de Skalna-Hladik), en las mismas cajas que medir_m3.py, con el
tiempo de evaluacion. Uso: python3 medir_m00.py <salida.csv> <bch>..."""
import csv, sys, time, zlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import aa_sym as A, medir_cajas as MC  # noqa: E402
rows = []
for f in sys.argv[2:]:
    box0, consts, trees = A.read_bch(f)
    for w, d, box in MC.walks(box0, seed=zlib.crc32(Path(f).stem.encode())):
        for name, tree in trees:
            try:
                t = time.perf_counter(); lo, hi = A.enclose(tree, box, consts, "M00"); t = time.perf_counter() - t
            except (ValueError, OverflowError, ZeroDivisionError):
                continue
            rows.append({"instancia": Path(f).stem, "expr": name, "caminata": w, "prof": d, "M00": hi - lo, "t_M00": t})
with open(sys.argv[1], "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
print(len(rows))
