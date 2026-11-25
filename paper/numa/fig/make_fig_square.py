"""Figura 2: el cuadrado de x = 1 + eps_1 + 2 eps_x (x0 = 1, ||x||_1 = 1, xr = 2).
(a) imagen clasica: en el plano (t, t^2), t = x in [-2, 4], la banda de Chebyshev de pendiente 2 es la mas angosta.
(b) la forma reducida expresa x^2 en funcion de u = e^T x solamente (el simbolo de error no es una coordenada del
    resultado): en el plano (u, x^2) los valores posibles forman una region (para cada u, t recorre un intervalo de
    ancho 2 xr). La banda de Chebyshev, 1 + 2u + [-4, 13], no es la mas angosta de pendiente 2 que contiene la
    region; la exacta es 1 + 2u + [-3, 13].
    python make_fig_square.py [kindle]"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
KINDLE = len(sys.argv) > 1 and sys.argv[1] == "kindle"
plt.rcParams.update({"font.size": 12 if KINDLE else 9, "font.family": "serif", "mathtext.fontset": "stix"})
fig, ax = plt.subplots(2, 1, figsize=(3.4, 6.0)) if KINDLE else plt.subplots(1, 2, figsize=(6.8, 2.9), sharey=True)
x0, S, xr = 1.0, 1.0, 2.0
fs = 9.5 if KINDLE else 8
ts = 11 if KINDLE else 9

# (a) plano (t, t^2)
a = ax[0]
t = np.linspace(-2, 4, 300)
a.fill_between(t, 2 * t - 1, 2 * t + 8, color="0.9", lw=0)
a.plot(t, 2 * t - 1, "k--", lw=1); a.plot(t, 2 * t + 8, "k--", lw=1)
a.plot(t, t ** 2, "k-", lw=1.8)
a.set_xlabel(r"$t=\hat x=1+\varepsilon_1+2\varepsilon_x\in[-2,4]$")
a.set_ylabel(r"$\hat x^2$")
a.set_title(r"(a) in $(t,t^2)$" if KINDLE else r"(a) in $(t,t^2)$: Chebyshev band $2t-1+[0,9]$", fontsize=ts)
a.set_ylim(-6, 18)

# (b) plano (u, x^2): region de valores posibles
b = ax[1]
u = np.linspace(-S, S, 300)
lo = np.array([min(((x0 + uu + xr * e) ** 2) for e in np.linspace(-1, 1, 801)) for uu in u])
hi = (x0 + u + xr) ** 2
b.fill_between(u, lo, hi, color="0.75", lw=0)
b.plot(u, x0 ** 2 + 2 * x0 * u - 4, "k--", lw=1)                    # Chebyshev, inferior
b.plot(u, x0 ** 2 + 2 * x0 * u + 13, "k-", lw=1)                    # superior (igual en ambos)
b.plot(u, x0 ** 2 + 2 * x0 * u - 3, "k-", lw=1.4)                   # exacta, inferior
b.plot(S, 0, "o", color="k", ms=4); b.plot(S, 16, "o", color="k", ms=4)
b.text(-0.6, 6, "possible values\nof $\\hat x^2$", fontsize=fs)
b.text(-0.95, 14.4, r"$1+2u+13$ (both)", fontsize=fs)
b.text(0.05, -4.9, r"Chebyshev: $1+2u-4$", fontsize=fs)
b.text(-0.95, -1.2, r"exact: $1+2u-3$", fontsize=fs)
b.set_xlabel(r"$u=e^Tx=\varepsilon_1\in[-1,1]$")
b.set_title(r"(b) in $(u,\hat x^2)$" if KINDLE else r"(b) in $(u,\hat x^2)$: a region, not a curve", fontsize=ts)
b.set_xlim(-1.05, 1.05); b.set_ylim(-6, 18)
if KINDLE: b.set_ylabel(r"$\hat x^2$")
fig.tight_layout()
name = "square_cheb_kindle" if KINDLE else "square_cheb"
fig.savefig(name + ".pdf"); fig.savefig(name + ".png", dpi=150)
