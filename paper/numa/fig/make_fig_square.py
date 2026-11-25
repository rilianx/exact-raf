"""Figura 2: N_2 = 2 x0 xr xi + (u + xr xi)^2 en funcion de u = e^T x in [-S, S], una curva por valor del
simbolo de error xi. La forma reducida debe encerrar todas las curvas con un solo intervalo (su rango);
la cota inferior de Chebyshev, -2|x0| xr, solo se alcanza si la curva de xi = -sign(x0) tiene su minimo
dentro de [-S, S], es decir, si xr <= S. Ejemplo x = 1 + eps_1 + xr eps_x (x0 = 1, S = 1).
    python make_fig_square.py [kindle]"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
KINDLE = len(sys.argv) > 1 and sys.argv[1] == "kindle"
plt.rcParams.update({"font.size": 12 if KINDLE else 9, "font.family": "serif", "mathtext.fontset": "stix"})
fig, ax = plt.subplots(2, 1, figsize=(3.4, 5.6)) if KINDLE else plt.subplots(1, 2, figsize=(6.8, 2.9), sharey=True)
x0, S = 1.0, 1.0
u = np.linspace(-S, S, 300)
for a, xr, tag in [(ax[0], 0.5, "a"), (ax[1], 2.0, "b")]:
    xis = np.linspace(-1, 1, 9)
    allv = np.array([2 * x0 * xr * xi + (u + xr * xi) ** 2 for xi in xis])
    a.fill_between(u, allv.min(0), allv.max(0), color="0.88", lw=0)
    for xi, v in zip(xis, allv):
        a.plot(u, v, color="0.6", lw=0.7)
    a.plot(u, allv[0], "k-", lw=1.8)                                    # xi = -1
    cheb_lo = -2 * abs(x0) * xr
    exact_lo = allv.min()
    a.axhline(cheb_lo, color="k", ls="--", lw=1)
    a.axhline(exact_lo, color="k", ls=":", lw=1.2) if exact_lo > cheb_lo + 1e-9 else None
    a.set_xlim(-1.15, 2.45)
    a.set_xlabel(r"$u=e^Tx$")
    fs = 8 if not KINDLE else 9.5
    if tag == "a":
        a.set_title(r"(a) $x_r=0.5\leq\|x\|_1=1$: exact", fontsize=9 if not KINDLE else 11)
        a.text(1.08, cheb_lo - 3.0, "Chebyshev\n$=\\min N_2=-1$", fontsize=fs)
        a.set_ylabel(r"$N_2$")
    else:
        a.set_title(r"(b) $x_r=2>\|x\|_1=1$: not exact", fontsize=9 if not KINDLE else 11)
        a.text(-1.1, cheb_lo - 1.3, r"Chebyshev $=-4$", fontsize=fs)
        a.text(-1.1, exact_lo + 0.45, r"$\min N_2=-3$", fontsize=fs)
    ue = np.linspace(S, xr + 0.35, 100)                                  # prolongacion de la curva xi = -1 fuera de [-S, S]
    if xr > S:
        a.plot(ue, -2 * x0 * xr + (ue - xr) ** 2, "k--", lw=1.2)
        a.plot(xr, cheb_lo, "o", mfc="white", mec="k", ms=6)
        a.text(1.08, cheb_lo + 2.4, "vertex outside\nthe range of $u$", fontsize=fs)
    else:
        a.plot(xr, cheb_lo, "o", color="k", ms=5)
    a.axvspan(S, 3, color="0.97", lw=0, zorder=0)
    a.axvline(S, color="0.5", lw=0.6); a.axvline(-S, color="0.5", lw=0.6)
ax[0].set_ylim(-6.2, 13.5); ax[1].set_ylim(-6.2, 13.5)
fig.tight_layout()
name = "square_cheb_kindle" if KINDLE else "square_cheb"
fig.savefig(name + ".pdf"); fig.savefig(name + ".png", dpi=150)
