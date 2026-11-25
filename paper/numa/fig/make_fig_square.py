"""Figura: por que el cuadrado de Chebyshev no es exacto en una forma reducida cuando el error domina.
Ejemplo x = 1 + 2 eps_x (n = 0): t = x en [-1, 3].
(a) banda de Chebyshev de t^2 con pendiente 2 x0 = 2: optima entre las bandas de esa pendiente.
(b) la forma reducida no puede conservar la pendiente (t depende solo del simbolo de error), asi que la banda se
    reemplaza por su proyeccion [-3, 9]; la banda horizontal exacta es [0, 9].
    python make_fig_square.py [kindle]"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
KINDLE = len(sys.argv) > 1 and sys.argv[1] == "kindle"
plt.rcParams.update({"font.size": 12 if KINDLE else 9, "font.family": "serif", "mathtext.fontset": "stix"})
fig, ax = plt.subplots(2, 1, figsize=(3.4, 5.4)) if KINDLE else plt.subplots(1, 2, figsize=(6.8, 2.7))
t = np.linspace(-1, 3, 200)
for a in ax:
    a.set_xlim(-1.4, 3.4); a.set_ylim(-4.5, 10.5); a.axhline(0, color="0.8", lw=0.6); a.axvline(0, color="0.8", lw=0.6)
    a.set_xlabel(r"$t=\hat x=1+2\varepsilon_x$")
# (a) Chebyshev
a = ax[0]
a.fill_between(t, 2 * t - 1, 2 * t + 3, color="0.85", lw=0)
a.plot(t, 2 * t - 1, "k--", lw=1); a.plot(t, 2 * t + 3, "k--", lw=1)
a.plot(t, t ** 2, "k-", lw=1.8)
a.set_title(r"(a) $t^2\in 2t-1+[0,4]$" if KINDLE else r"(a) Chebyshev: $t^2\in 2t-1+[0,4]$", fontsize=9 if not KINDLE else 11)
a.set_ylabel(r"$t^2$")
# (b) bandas horizontales: la forma reducida solo conserva el centro (n = 0)
b = ax[1]
b.axhspan(0, 9, color="0.75", lw=0)
b.axhspan(-3, 0, facecolor="white", edgecolor="0.35", hatch="////", lw=0)
b.axhline(9, color="k", ls="--", lw=1); b.axhline(-3, color="k", ls="--", lw=1)
b.plot(t, t ** 2, "k-", lw=1.8)
fs = 8 if not KINDLE else 9.5
b.text(-1.3, 6.2, "exact: $[0,9]$", fontsize=fs, bbox=dict(facecolor="white", edgecolor="none", pad=1))
b.text(-1.3, -2.2, "lost by Chebyshev", fontsize=fs, bbox=dict(facecolor="white", edgecolor="none", pad=1))
b.set_title(r"(b) without the slope: $[-3,9]$" if KINDLE else r"(b) the slope cannot be kept: $\hat x^2\in[-3,9]$", fontsize=9 if not KINDLE else 11)
fig.tight_layout()
name = "square_cheb_kindle" if KINDLE else "square_cheb"
fig.savefig(name + ".pdf"); fig.savefig(name + ".png", dpi=150)
