"""Resumen de synq_raiz.csv: medias geometricas de anchos relativos a SH (suma por producto) y al exacto."""
import csv, sys
from math import exp, log
rows = list(csv.DictReader(open(sys.argv[1] if len(sys.argv) > 1 else "synq_raiz.csv")))
M = ["S", "SH", "RAAp0", "RAAp6", "EIG", "Qk0", "Qk2", "Qk4", "Qk6", "Qk8", "EXACT"]
gm = lambda xs: exp(sum(map(log, xs)) / len(xs))
def tabla(grupos, ref):
    print(f"\n| grupo | inst | " + " | ".join(M) + " |")
    print("|---|---|" + "---|" * len(M))
    for nom, rs in grupos:
        print(f"| {nom} | {len(rs)} | " + " | ".join(f"{gm([float(r[m]) / float(r[ref]) for r in rs]):.3f}" for m in M) + " |")
gr = [("todas", rows)] + [(f"n = {n}", [r for r in rows if r["n"] == n]) for n in sorted({r["n"] for r in rows}, key=int)] \
   + [(f"d = {d}", [r for r in rows if r["densidad"] == d]) for d in sorted({r["densidad"] for r in rows})] \
   + [(f"convexa = {c}", [r for r in rows if r["convexa"] == c]) for c in ("True", "False")]
print("Ancho / ancho de SH (suma por producto); media geometrica"); tabla(gr, "SH")
print("\nAncho / ancho exacto; media geometrica"); tabla(gr, "EXACT")
# ganancia marginal del k-nomial sobre la Q completa termino a termino
for k in (2, 4, 6, 8):
    xs = [float(r[f"Qk{k}"]) / float(r["Qk0"]) for r in rows]
    print(f"Qk{k}/Qk0: media geom {gm(xs):.3f}, min {min(xs):.3f}, max {max(xs):.3f}; "
          f"brecha cerrada (Qk0-Qk{k})/(Qk0-EXACT) = {gm([max(1e-9,(float(r['Qk0'])-float(r[f'Qk{k}']))/(float(r['Qk0'])-float(r['EXACT']))) for r in rows if float(r['Qk0'])-float(r['EXACT'])>1e-6]):.3f} (media geom)")
print("instancias con Qk6 >= 0,9 SH:", sum(float(r["Qk6"]) >= 0.9 * float(r["SH"]) for r in rows), "de", len(rows))
print("instancias con Qk0 > SH:", [r["instancia"] for r in rows if float(r["Qk0"]) > float(r["SH"])])
