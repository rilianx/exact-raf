"""Genera las formas homogeneizadas H de M3 centrado (mismas cajas que medir_m3c.py), corre prodk en C++
(prodk_driver) y en Python (prodk_fast) y compara anchos y tiempos. Python solo con m <= 40 (como en las
mediciones); C++ con todas.
Uso: python3 comparar.py <prodk_driver> <salida.csv> <bch>..."""
import csv, subprocess, sys, time, zlib
from pathlib import Path
import numpy as np
H_ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(H_ / "raak-errores-simbolos")); sys.path.insert(0, str(H_ / "raak-qplib"))
import aa_sym as A, aa_quad as AQ, medir_cajas as MC, prodk_fast as PF  # noqa: E402

drv, out_csv, files = sys.argv[1], sys.argv[2], sys.argv[3:]
forms = []
for f in files:
    box0, consts, trees = A.read_bch(f)
    for w, d, box in MC.walks(box0, seed=zlib.crc32(Path(f).stem.encode())):
        for name, tree in trees:
            try: F = AQ.EvalQ(box, consts, centered=True).evq(tree)
            except (ValueError, OverflowError, ZeroDivisionError): continue
            if F.Q: forms.append(({"instancia": Path(f).stem, "expr": name, "caminata": w, "prof": d}, AQ.homog(F)))
print(len(forms), "formas", flush=True)
inp = "\n".join(f"{H.shape[0]} " + " ".join(f"{x:.17g}" for x in H.ravel()) for _, H in forms) + "\n"
out = subprocess.run([drv, "5"], input=inp, capture_output=True, text=True, check=True).stdout.strip().splitlines()
rows = []
for (meta, H), line in zip(forms, out):
    L, U, nb, t = line.split()
    r = {**meta, "m": H.shape[0], "cpp_ancho": float(U) - float(L), "cpp_bloques": int(nb), "cpp_us": float(t)}
    if H.shape[0] <= 40:
        t0 = time.perf_counter(); l2, u2, nb2 = PF.prodk_interval(H); r["py_us"] = 1e6 * (time.perf_counter() - t0)
        r["py_ancho"] = u2 - l2; r["py_bloques"] = nb2
    rows.append(r)
keys = list(dict.fromkeys(k for r in rows for k in r))
with open(out_csv, "w", newline="") as fh:
    wr = csv.DictWriter(fh, fieldnames=keys); wr.writeheader(); wr.writerows(rows)
