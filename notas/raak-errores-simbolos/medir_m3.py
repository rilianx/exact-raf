"""Una sola pasada, en las mismas cajas (semilla determinista): M0, M1, M2 (aa_sym), M3a y M3b (aa_quad) y el
rango interior; ademas la atribucion del radio de M1 y de M3 por tipo de operacion.

Uso: python3 medir_m3.py <salida.csv> <bch>...
"""
import csv, math, sys, time, zlib
from collections import defaultdict
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import aa_sym as A    # noqa: E402
import aa_quad as AQ  # noqa: E402
import medir_cajas as MC  # noqa: E402

TAGS = ["lineal", "mul_lin", "mul_comp", "sqr_lin", "sqr_comp", "recip", "pow", "exp", "log", "sqrt", "sincos", "abs", "otro"]


def attrib(lin, rad, pref):
    parts = defaultdict(float)
    for k, v in lin.items(): parts[k.split(":")[1] if k.startswith("#") else "lineal"] += abs(v)
    return {f"{pref}_{t}": parts.get(t, 0.0) / rad for t in TAGS}


def main():
    out_csv, files = sys.argv[1], sys.argv[2:]
    rows = []
    for f in files:
        t0 = time.time()
        box0, consts, trees = A.read_bch(f)
        for w, d, box in MC.walks(box0, seed=zlib.crc32(Path(f).stem.encode())):
            for name, tree in trees:
                r = {"instancia": Path(f).stem, "n": len(box0), "expr": name, "caminata": w, "prof": d}
                try:
                    for m in ("M0", "M1", "M2"):
                        lo, hi = A.enclose(tree, box, consts, m); r[m] = hi - lo
                    F1 = A.Evaluator(box, consts, "M1").ev(tree)
                    E3 = AQ.EvalQ(box, consts); F3 = E3.evq(tree)
                    t1 = time.time(); r["M3a"], r["M3b"] = AQ.final_bounds(F3); r["t_M3b"] = time.time() - t1
                    r["M3_terminos_Q"] = len(F3.Q)
                    lo, hi = A.inner_range(tree, box, consts); r["INNER"] = hi - lo
                except (ValueError, OverflowError, ZeroDivisionError):
                    continue
                if not all(math.isfinite(r[k]) for k in ("M0", "M1", "M2", "M3a", "M3b", "INNER")): continue
                for m in ("M0", "M1", "M2", "M3a", "M3b"):
                    assert r[m] >= r["INNER"] * (1 - 1e-6) - 1e-9 * max(1, abs(r["INNER"])), (f, name, d, m, r[m], r["INNER"])
                if F1.rad() > 0: r.update(attrib(F1.lin, F1.rad(), "a1"))
                rows.append(r)
        print(f"{Path(f).stem:20s} {time.time() - t0:6.1f}s", flush=True)
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(out_csv, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=keys); wr.writeheader(); wr.writerows(rows)


if __name__ == "__main__":
    main()
