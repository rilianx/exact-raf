"""Mide el ancho del encierro de cada objetivo/restriccion con M0, M1 y M2 (aa_sym.py) en la caja raiz y en cajas
internas, y lo compara con el rango interior (muestreo + busqueda local).

Cajas internas: caminatas aleatorias de biseccion desde la raiz (se biseca la variable de mayor ancho, como el
bisector por omision de Ibex con precision relativa, y se elige un hijo al azar); se guardan las cajas a las
profundidades DEPTHS. Varias caminatas por instancia.

Uso: python3 medir_cajas.py <salida.csv> <bch>...
"""
import csv, math, sys, time, zlib
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import aa_sym as A  # noqa: E402

DEPTHS = (0, 5, 10, 20)
WALKS = 4
MODES = ("M0", "M1", "M2")


def walks(box, seed):
    rng = np.random.default_rng(seed)
    out = [(0, 0, dict(box))]
    for w in range(WALKS):
        b = dict(box)
        for d in range(1, max(DEPTHS) + 1):
            v = max(b, key=lambda k: (b[k][1] - b[k][0]) / max(1.0, abs(b[k][0]) + abs(b[k][1])))
            lo, hi = b[v]; mid = (lo + hi) / 2
            b[v] = (lo, mid) if rng.random() < 0.5 else (mid, hi)
            if d in DEPTHS: out.append((w + 1, d, dict(b)))
    return out


def main():
    out_csv, files = sys.argv[1], sys.argv[2:]
    rows = []
    for f in files:
        t0 = time.time()
        box0, consts, trees = A.read_bch(f)
        for w, d, box in walks(box0, seed=zlib.crc32(Path(f).stem.encode())):
            for name, tree in trees:
                r = {"instancia": Path(f).stem, "n": len(box0), "expr": name, "caminata": w, "prof": d}
                try:
                    for m in MODES:
                        lo, hi = A.enclose(tree, box, consts, m)
                        r[m] = hi - lo
                    lo, hi = A.inner_range(tree, box, consts)
                    r["INNER"] = hi - lo
                except (ValueError, OverflowError, ZeroDivisionError) as exc:
                    continue
                if not all(math.isfinite(r[k]) for k in MODES + ("INNER",)): continue
                for m in MODES:  # validez: ningun encierro por debajo del rango interior
                    assert r[m] >= r["INNER"] * (1 - 1e-6) - 1e-9 * max(1, abs(r["INNER"])), (f, name, d, m, r[m], r["INNER"])
                rows.append(r)
        print(f"{Path(f).stem:20s} {time.time() - t0:6.1f}s", flush=True)
        with open(out_csv, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)


if __name__ == "__main__":
    main()
