"""M3 original y M3 centrado (aa_quad.EvalQ centered=False/True), con los encierros finales por separado
(tt, prodk, kn6, sdp), en las mismas cajas que medir_m3.py. Se une despues con cajas_m3.csv (M0, M1, M2, INNER).
Uso: python3 medir_m3c.py <driver knomial> <salida.csv> <bch>..."""
import csv, math, sys, time, zlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import aa_sym as A, aa_quad as AQ, medir_cajas as MC  # noqa: E402


def main():
    driver, out_csv, files = sys.argv[1], sys.argv[2], sys.argv[3:]
    rows = []
    for f in files:
        t0 = time.time()
        box0, consts, trees = A.read_bch(f)
        for w, d, box in MC.walks(box0, seed=zlib.crc32(Path(f).stem.encode())):
            for name, tree in trees:
                r = {"instancia": Path(f).stem, "expr": name, "caminata": w, "prof": d}
                try:
                    for tag, cen in (("o", False), ("c", True)):
                        F = AQ.EvalQ(box, consts, centered=cen).evq(tree)
                        t1 = time.time()
                        for k, v in AQ.final_all(F, driver).items(): r[f"M3{tag}_{k}"] = v
                        r[f"t_M3{tag}"] = time.time() - t1
                except (ValueError, OverflowError, ZeroDivisionError):
                    continue
                rows.append(r)
        print(f"{Path(f).stem:20s} {time.time() - t0:6.1f}s", flush=True)
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(out_csv, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=keys); wr.writeheader(); wr.writerows(rows)


if __name__ == "__main__":
    main()
