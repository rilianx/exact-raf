"""Vuelca las formas homogeneizadas H de M3 centrado (mismas cajas que las mediciones) a un archivo de texto
(una por linea, formato de prodk_driver) y su metadato a un CSV.
Uso: python3 volcar_formas.py <formas.txt> <meta.csv> <bch>..."""
import csv, sys, zlib
from pathlib import Path
H_ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(H_ / "raak-errores-simbolos")); sys.path.insert(0, str(H_ / "raak-qplib"))
import aa_sym as A, aa_quad as AQ, medir_cajas as MC  # noqa: E402
out, meta, files = sys.argv[1], sys.argv[2], sys.argv[3:]
with open(out, "w") as fo, open(meta, "w", newline="") as fm:
    w = csv.writer(fm); w.writerow(["instancia", "expr", "caminata", "prof", "m"])
    for f in files:
        box0, consts, trees = A.read_bch(f)
        for wk, d, box in MC.walks(box0, seed=zlib.crc32(Path(f).stem.encode())):
            for name, tree in trees:
                try: F = AQ.EvalQ(box, consts, centered=True).evq(tree)
                except (ValueError, OverflowError, ZeroDivisionError): continue
                if not F.Q: continue
                H = AQ.homog(F)
                fo.write(f"{H.shape[0]} " + " ".join(f"{x:.17g}" for x in H.ravel()) + "\n")
                w.writerow([Path(f).stem, name, wk, d, H.shape[0]])
