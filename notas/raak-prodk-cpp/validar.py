"""Valida que el ancho de prodk (C++) nunca quede por debajo del rango exacto de z^T H z (enumeracion de caras KKT,
m <= 12). Regenera las formas como comparar.py.
Uso: python3 validar.py <prodk_driver> <bch>..."""
import subprocess, sys, zlib
from pathlib import Path
H_ = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(H_ / "raak-errores-simbolos")); sys.path.insert(0, str(H_ / "raak-qplib"))
import aa_sym as A, aa_quad as AQ, medir_cajas as MC, medir_q as S  # noqa: E402
drv, files = sys.argv[1], sys.argv[2:]
Hs = []
for f in files:
    box0, consts, trees = A.read_bch(f)
    for w, d, box in MC.walks(box0, seed=zlib.crc32(Path(f).stem.encode())):
        for name, tree in trees:
            try: F = AQ.EvalQ(box, consts, centered=True).evq(tree)
            except (ValueError, OverflowError, ZeroDivisionError): continue
            if F.Q:
                H = AQ.homog(F)
                if H.shape[0] <= 12: Hs.append((Path(f).stem, name, d, H))
inp = "\n".join(f"{H.shape[0]} " + " ".join(f"{x:.17g}" for x in H.ravel()) for *_, H in Hs) + "\n"
out = subprocess.run([drv, "1"], input=inp, capture_output=True, text=True, check=True).stdout.strip().splitlines()
bad, ratio = 0, []
for (i, e, d, H), line in zip(Hs, out):
    L, U = map(float, line.split()[:2])
    lo, hi = S.quad_range_exact(H)
    sc = max(1.0, abs(lo), abs(hi))
    if L > lo + 1e-9 * sc or U < hi - 1e-9 * sc:
        bad += 1
        if bad <= 5: print("INVALIDO", i, e, d, H.shape[0], (L, U), (lo, hi))
    if hi - lo > 0: ratio.append((U - L) / (hi - lo))
import statistics as st
print(f"formas validadas (m <= 12): {len(Hs)}; invalidas: {bad}; ancho C++ / exacto: mediana {st.median(ratio):.3f}, p90 {sorted(ratio)[int(.9*len(ratio))]:.3f}")
