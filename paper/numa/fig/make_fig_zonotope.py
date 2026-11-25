"""Figura del Lema 2, del Algoritmo 1 y del Ejemplo 1: (a) zonotopo Z de los generadores ordenados
por angulo; (b) curvas de nivel de (u+delta1)(v+delta2) con delta = (1, 1) y sus extremos en el borde.
    python make_fig_zonotope.py          # zonotope.pdf (dos paneles lado a lado)
    python make_fig_zonotope.py kindle   # zonotope_kindle.pdf (paneles uno sobre otro)"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

KINDLE = len(sys.argv) > 1 and sys.argv[1] == "kindle"
g = np.array([[2.0, -1.0], [1.0, 0.0], [2.0, 1.0]])   # Ejemplo 1, ya en u > 0 y ordenados por angulo
p0 = -g.sum(axis=0)
V = [p0]
for k in range(3): V.append(V[-1] + 2 * g[k])
for k in range(3): V.append(V[-1] - 2 * g[k])
V = np.array(V)                                         # p0..p6 (p6 = p0)
d = np.array([1.0, 1.0])                                # delta para (s, t) = (1, 1)

plt.rcParams.update({"font.size": 12 if KINDLE else 9, "font.family": "serif", "mathtext.fontset": "stix"})
fig, ax = plt.subplots(2, 1, figsize=(3.4, 5.2)) if KINDLE else plt.subplots(1, 2, figsize=(6.8, 2.6))
XL, YL = (-6.0, 6.0), (-2.9, 3.2)

# (a) borde de Z
a = ax[0]
a.axhline(0, color="0.75", lw=0.6, zorder=0); a.axvline(0, color="0.75", lw=0.6, zorder=0)
a.fill(V[:, 0], V[:, 1], color="0.88", zorder=1)
a.plot(V[:, 0], V[:, 1], "k-", lw=1.0, zorder=2)
cols = ["0.0", "0.38", "0.6"]          # blanco y negro: los generadores se distinguen por el tono de gris
styles = ["-", "-", "-"]
for k in range(3):
    a.plot([V[k][0], V[k + 1][0]], [V[k][1], V[k + 1][1]], color=cols[k], ls=styles[k], lw=4.0, solid_capstyle="butt", zorder=3)
    a.annotate("", xy=g[k], xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=cols[k], lw=1.8, mutation_scale=10), zorder=4)
a.plot(V[:-1, 0], V[:-1, 1], "ko", ms=3, zorder=5)
off = {0: (-0.95, -0.15), 1: (-0.55, -0.75), 2: (0.05, -0.75), 3: (0.25, -0.15), 4: (0.1, 0.25), 5: (-0.6, 0.25)}
for k in range(6):
    a.annotate(f"$p_{k}$", V[k] + np.array(off[k]))
for k, o in enumerate([(0.15, -0.75), (0.25, -0.15), (-1.3, 0.3)]):
    a.annotate(f"$\\tilde g_{{({k+1})}}$", g[k] + np.array(o), color=cols[k])
a.set_title("(a) boundary of $Z$", fontsize=9 if not KINDLE else 11)
a.set_xlabel("$u=e^Tx$"); a.set_ylabel("$v=e^Ty$")
a.set_aspect("equal"); a.set_xlim(*XL); a.set_ylim(YL[0], 3.0 if KINDLE else 5.2)

# (b) curvas de nivel de (u + d1)(v + d2) sobre Z
b = ax[1]
uu, vv = np.meshgrid(np.linspace(*XL, 500), np.linspace(YL[0], 5.2, 400))
PH = (uu + d[0]) * (vv + d[1])
b.contour(uu, vv, PH, levels=[-12, -10, -8, -6, -4, -2, 2, 4, 6, 8, 10, 12], colors="0.6", linewidths=0.6)
b.fill(V[:, 0], V[:, 1], color="0.88", zorder=0)
b.plot(V[:, 0], V[:, 1], "k-", lw=1.0)
pts = np.vstack([V[k] + t * (V[k + 1] - V[k]) for k in range(6) for t in np.linspace(0, 1, 4001)])
val = (pts[:, 0] + d[0]) * (pts[:, 1] + d[1])
imax, imin = val.argmax(), val.argmin()
b.plot(*pts[imax], "^", color="black", ms=7, label=r"max $8$ (stationary point)")
b.plot(*pts[imin], "v", mfc="white", mec="black", mew=1.4, ms=7, label=r"min $-9/2$ (stationary point)")
b.plot(-d[0], -d[1], "kx", ms=7, label=r"saddle $(-\delta_1,-\delta_2)$")
if KINDLE:
    b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), fontsize=9, frameon=False)
else:
    b.legend(loc="upper left", fontsize=6.5, frameon=True, framealpha=0.95)
b.set_title(r"(b) level curves of $\varphi_\delta$, $\delta=(1,1)$" if KINDLE else r"(b) level curves of $(u+\delta_1)(v+\delta_2)$, $\delta=(1,1)$", fontsize=9 if not KINDLE else 11)
b.set_xlabel("$u$"); b.set_aspect("equal"); b.set_xlim(*XL); b.set_ylim(YL[0], 3.0 if KINDLE else 5.2)
fig.tight_layout()
name = "zonotope_kindle" if KINDLE else "zonotope"
fig.savefig(name + ".pdf"); fig.savefig(name + ".png", dpi=150)
print("max", pts[imax], val.max(), "min", pts[imin], val.min())
