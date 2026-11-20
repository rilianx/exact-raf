"""Presupuesto de error en M1: el radio del encierro final = sum |coef| de las variables + sum |coef| de cada
simbolo de error; los simbolos llevan la etiqueta de la operacion que los creo. Se reporta la fraccion del radio
que aporta cada tipo de operacion, en las mismas cajas que medir_cajas.py.

Uso: python3 atribuir.py <salida.csv> <bch>...
"""
import csv, math, sys, time
from collections import defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import aa_sym as A  # noqa: E402
import medir_cajas as MC  # noqa: E402

TAGS = ["mul_lin", "mul_comp", "sqr_lin", "sqr_comp", "recip", "pow", "exp", "log", "sqrt", "sincos", "abs", "otro"]


def main():
    out_csv, files = sys.argv[1], sys.argv[2:]
    rows = []
    for f in files:
        t0 = time.time()
        box0, consts, trees = A.read_bch(f)
        for w, d, box in MC.walks(box0, seed=hash(Path(f).stem) % 2**32):
            for name, tree in trees:
                try:
                    F = A.Evaluator(box, consts, "M1").ev(tree)
                except (ValueError, OverflowError, ZeroDivisionError):
                    continue
                rad = F.rad()
                if not (math.isfinite(rad) and rad > 0): continue
                parts = defaultdict(float)
                for k, v in F.lin.items():
                    parts[k.split(":")[1] if k.startswith("#") else "lineal"] += abs(v)
                r = {"instancia": Path(f).stem, "expr": name, "caminata": w, "prof": d, "radio": rad}
                for t in ["lineal"] + TAGS: r[t] = parts.get(t, 0.0) / rad
                rows.append(r)
        print(f"{Path(f).stem:20s} {time.time() - t0:6.1f}s", flush=True)
        with open(out_csv, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)


if __name__ == "__main__":
    main()
