"""Figura del Lema 2 y del Algoritmo 1: (a) zonotopo plano y su borde a partir de los generadores
ordenados por angulo; (b) curvas de nivel de (u+c1)(v+c2) y sus extremos en el borde."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

g = np.array([[0.6, -0.5], [1.0, 0.3], [0.4, 0.8]])          # ya en u>0 y ordenados por angulo
p0 = -g.sum(axis=0)
V = [p0]
for k in range(3): V.append(V[-1] + 2 * g[k])
for k in range(3): V.append(V[-1] - 2 * g[k])
V = np.array(V)                                                # p0..p6 (p6 = p0)

plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "stix"})
import sys
KINDLE = len(sys.argv) > 1 and sys.argv[1] == "kindle"
plt.rcParams.update({"font.size": 12} if KINDLE else {})
fig, ax = plt.subplots(2, 1, figsize=(3.2, 5.6)) if KINDLE else plt.subplots(1, 2, figsize=(6.6, 3.0))

# (a)
a = ax[0]
a.fill(V[:, 0], V[:, 1], color="0.88", zorder=0)
a.plot(V[:, 0], V[:, 1], "k-", lw=1.2)
a.plot(V[:-1, 0], V[:-1, 1], "ko", ms=3)
for k in range(6):
    off = {0: (-0.32, -0.12), 1: (0.05, -0.22), 2: (0.08, -0.05), 3: (0.05, 0.08), 4: (-0.25, 0.1), 5: (-0.32, 0.0)}[k]
    a.annotate(f"$p_{k}$", V[k] + off)
cols = ["C0", "C1", "C2"]
for k in range(3):
    a.annotate("", xy=g[k], xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=cols[k], lw=1.4))
    a.annotate(f"$\\tilde g_{{({k+1})}}$", g[k] * 1.08 + np.array([0.04, -0.06 if k == 0 else 0.02]), color=cols[k])
    a.annotate("", xy=V[k + 1], xytext=V[k], arrowprops=dict(arrowstyle="-", color=cols[k], lw=2.4))
a.plot(0, 0, "k+", ms=6)
a.set_title("(a) boundary of $Z$", fontsize=9)
a.set_xlabel("$u=e^Tx$"); a.set_ylabel("$v=e^Ty$")
a.set_aspect("equal"); a.set_xlim(-2.4, 2.4); a.set_ylim(-2.0, 2.9)

# (b) shifted product over Z: phi_c(u,v) = (u+c1)(v+c2), c = (0.5, -0.4)
c = np.array([0.5, -0.4])
b = ax[1]
uu, vv = np.meshgrid(np.linspace(-2.4, 2.4, 400), np.linspace(-2.0, 3.4, 400))
PH = (uu + c[0]) * (vv + c[1])
b.contour(uu, vv, PH, levels=np.linspace(-4, 4, 17), colors="0.6", linewidths=0.6)
b.fill(V[:, 0], V[:, 1], color="0.88", zorder=0)
b.plot(V[:, 0], V[:, 1], "k-", lw=1.2)
# extremos sobre el borde (muestreo fino de las aristas)
pts = np.vstack([V[k] + t * (V[k + 1] - V[k]) for k in range(6) for t in np.linspace(0, 1, 2001)])
val = (pts[:, 0] + c[0]) * (pts[:, 1] + c[1])
imax, imin = val.argmax(), val.argmin()
b.plot(*pts[imax], "^", color="C3", ms=7, label="max (stationary point of an edge)")
b.plot(*pts[imin], "v", color="C0", ms=7, label="min (vertex)")
b.plot(-c[0], -c[1], "kx", ms=7, label="saddle $(-c_1,-c_2)$")
b.legend(loc="upper left", fontsize=7.5 if KINDLE else 6.5, frameon=True, framealpha=0.95)
b.set_title("(b) level curves of $(u+c_1)(v+c_2)$", fontsize=9)
b.set_xlabel("$u$"); b.set_aspect("equal"); b.set_xlim(-2.4, 2.4); b.set_ylim(-2.0, 3.4 if KINDLE else 2.9)
fig.tight_layout()
name = "zonotope_kindle" if KINDLE else "zonotope"
fig.savefig(name + ".pdf"); fig.savefig(name + ".png", dpi=150)
print("max", pts[imax], val.max(), "min", pts[imin], val.min())
