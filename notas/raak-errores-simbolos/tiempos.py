"""Tiempos por expresion y caja (Python, salvo el k-nomial en C++): evaluacion de M0, M1, M2 y M3 centrado, y los
encierros finales de M3 por separado: tt (barato), prodk (bloques producto, solo la reduccion), kn6 (k-nomial,
C++, medido en lote para no contar el arranque del proceso) y sdp (Shor con SCS). Mismas cajas que las mediciones
(caminata 1, profundidades 0 y 10).
Uso: python3 tiempos.py <driver> <salida.csv> <bch>..."""
import csv, subprocess, sys, time, zlib
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import aa_sym as A, aa_quad as AQ, medir_cajas as MC  # noqa: E402
import prodk_fast as PF  # noqa: E402
import medir_q as S      # noqa: E402

def tick(f, *a):
    t = time.perf_counter(); r = f(*a); return r, time.perf_counter() - t

def main():
    driver, out_csv, files = sys.argv[1], sys.argv[2], sys.argv[3:]
    rows, Hs = [], []
    for f in files:
        box0, consts, trees = A.read_bch(f)
        for w, d, box in MC.walks(box0, seed=zlib.crc32(Path(f).stem.encode())):
            if w > 1 or d not in (0, 10): continue
            for name, tree in trees:
                r = {"instancia": Path(f).stem, "expr": name, "prof": d}
                try:
                    for m in ("M0", "M1", "M2"):
                        _, r[f"t_{m}"] = tick(A.enclose, tree, box, consts, m)
                    F, r["t_M3c_eval"] = tick(lambda: AQ.EvalQ(box, consts, centered=True).evq(tree))
                except (ValueError, OverflowError, ZeroDivisionError):
                    continue
                _, r["t_tt"] = tick(F.rng)
                r["terminos"] = len(F.lin) + len(F.Q)
                if F.Q:
                    H, r["t_homog"] = tick(AQ.homog, F)
                    r["m"] = H.shape[0]
                    if H.shape[0] <= 40: _, r["t_prodk"] = tick(PF.prodk_interval, H)
                    if H.shape[0] <= 60:
                        try: _, r["t_sdp"] = tick(S.sdp_width, H)
                        except Exception: r["t_sdp"] = float("nan")   # SCS fallo (en la medicion se usa el respaldo)
                    Hs.append(H)
                rows.append(r)
        print(Path(f).stem, flush=True)
    # k-nomial en lote: tiempo total del driver / numero de formas
    lines = [f"Q {H.shape[0]} 6 " + " ".join(f"{x:.17g}" for x in H.ravel()) for H in Hs]
    t = time.perf_counter(); subprocess.run([driver], input="\n".join(lines) + "\n", capture_output=True, text=True, check=True)
    t_kn = (time.perf_counter() - t) / max(1, len(Hs))
    for r in rows:
        if "m" in r: r["t_kn6_medio_lote"] = t_kn
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(out_csv, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys); wr.writeheader(); wr.writerows(rows)

if __name__ == "__main__":
    main()
