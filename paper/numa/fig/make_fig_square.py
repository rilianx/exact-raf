"""Figura 2: el cuadrado de x = 1 + eps_1 + 2 eps_x (x0 = 1, ||x||_1 = 1, xr = 2), t = x in [-2, 4].
(a) imagen clasica de Chebyshev de t^2 en una variable: banda de pendiente 2 x0 = 2, t^2 in 2t - 1 + [0, 9].
(b) con el error hay dos variables, u = e^T x y eps_x: N_2(u, eps_x) = 4 eps_x + (u + 2 eps_x)^2 sobre el
    rectangulo [-1,1]^2; la cota de Chebyshev (-4) es su minimo en (u, eps_x) = (2, -1), fuera del rectangulo;
    el minimo en el rectangulo es -3, en (1, -1).
(c) cortes: una parabola en u por cada valor de eps_x; la forma reducida las encierra con un solo intervalo.
    python make_fig_square.py [kindle]"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
KINDLE = len(sys.argv) > 1 and sys.argv[1] == "kindle"
plt.rcParams.update({"font.size": 12 if KINDLE else 8.5, "font.family": "serif", "mathtext.fontset": "stix"})
fig, ax = plt.subplots(3, 1, figsize=(3.4, 8.4)) if KINDLE else plt.subplots(1, 3, figsize=(7.0, 2.6))
x0, S, xr = 1.0, 1.0, 2.0
fs = 9.5 if KINDLE else 7.5
ts = 11 if KINDLE else 8.5
N2 = lambda u, e: 2 * x0 * xr * e + (u + xr * e) ** 2

# (a) Chebyshev en una variable
a = ax[0]
t = np.linspace(-2, 4, 300)
a.fill_between(t, 2 * t - 1, 2 * t + 8, color="0.88", lw=0)
a.plot(t, 2 * t - 1, "k--", lw=1); a.plot(t, 2 * t + 8, "k--", lw=1)
a.plot(t, t ** 2, "k-", lw=1.8)
a.set_xlabel(r"$t=\hat x\in[-2,4]$"); a.set_ylabel(r"$t^2$")
a.set_title(r"(a) one variable: $t^2\in 2t-1+[0,9]$", fontsize=ts)

# (b) dos variables: curvas de nivel de N_2 sobre (u, eps_x)
b = ax[1]
uu, ee = np.meshgrid(np.linspace(-1.3, 2.6, 400), np.linspace(-1.35, 1.35, 300))
cs = b.contour(uu, ee, N2(uu, ee), levels=[-3.5, -3, -1, 2, 6], colors="0.45", linewidths=0.7)
b.clabel(cs, fmt="%g", fontsize=fs - 1, manual=[(1.6, 0.3), (-0.3, -0.75), (0.1, 0.55), (1.55, -0.55), (1.95, -0.8)])
b.add_patch(plt.Rectangle((-S, -1), 2 * S, 2, fill=True, facecolor="0.9", edgecolor="k", lw=1.2, zorder=0))
b.plot(xr, -1, "o", mfc="white", mec="k", ms=6, zorder=5)
b.plot(S, -1, "o", color="k", ms=5, zorder=5)
b.annotate(r"$-4$ (Chebyshev)", xy=(xr, -1), xytext=(1.25, -0.25), fontsize=fs, arrowprops=dict(arrowstyle="->", lw=0.7))
b.annotate(r"$\min N_2=-3$", xy=(S, -1), xytext=(-1.2, -0.62), fontsize=fs, arrowprops=dict(arrowstyle="->", lw=0.7))
b.plot(S, 1, "^", color="k", ms=6, zorder=5)
b.annotate(r"$\max N_2=13$", xy=(S, 1), xytext=(1.2, 1.05), fontsize=fs, arrowprops=dict(arrowstyle="->", lw=0.7))
b.set_xlim(-1.3, 2.6); b.set_ylim(-1.35, 1.35)
b.set_xlabel(r"$u=e^Tx$"); b.set_ylabel(r"$\varepsilon_x$")
b.set_title(r"(b) two variables: $N_2(u,\varepsilon_x)$", fontsize=ts)

# (c) cortes: parabolas en u para cada eps_x
c = ax[2]
u = np.linspace(-S, S, 300)
xis = np.linspace(-1, 1, 9)
allv = np.array([N2(u, xi) for xi in xis])
c.fill_between(u, allv.min(0), allv.max(0), color="0.88", lw=0)
for v in allv: c.plot(u, v, color="0.6", lw=0.7)
c.plot(u, allv[0], "k-", lw=1.8)
ue = np.linspace(S, xr + 0.3, 100)
c.plot(ue, N2(ue, -1), "k--", lw=1.1)
c.plot(xr, -2 * x0 * xr, "o", mfc="white", mec="k", ms=6)
c.plot(S, N2(S, -1), "o", color="k", ms=5)
c.axhline(-4, color="k", ls="--", lw=0.8); c.axhline(-3, color="k", ls=":", lw=1.1)
c.axvline(S, color="0.5", lw=0.6); c.axvline(-S, color="0.5", lw=0.6)
c.text(-1.0, -2.6, r"$\min N_2=-3$", fontsize=fs); c.text(-1.0, -5.6, r"Chebyshev $=-4$", fontsize=fs)
hi = N2(S, 1)
c.plot(S, hi, "^", color="k", ms=6)
c.axhline(hi, color="k", ls="-", lw=0.6)
c.text(-1.0, hi + 0.6, r"$\max N_2=13$", fontsize=fs)
# rangos: exacto [-3, 13] y Chebyshev [-4, 13]
c.plot([1.55, 1.55], [-3, hi], color="k", lw=2.2); c.plot([1.5, 1.6], [-3, -3], "k-", lw=1.2); c.plot([1.5, 1.6], [hi, hi], "k-", lw=1.2)
c.plot([2.05, 2.05], [-4, hi], color="k", lw=1.2, ls="--"); c.plot([2.0, 2.1], [-4, -4], "k-", lw=1.2); c.plot([2.0, 2.1], [hi, hi], "k-", lw=1.2)
c.text(1.3, 5.5, r"$\mathrm{rng}\,N_2$" + "\n" + r"$=[-3,13]$", fontsize=fs - 0.5, ha="right")
c.text(2.12, 5.5, "Cheb.\n$[-4,13]$", fontsize=fs - 0.5)
c.set_xlim(-1.15, 2.75); c.set_ylim(-6.2, 16.5)
c.set_xlabel(r"$u=e^Tx$"); c.set_ylabel(r"$N_2$")
c.set_title(r"(c) slices at fixed $\varepsilon_x$", fontsize=ts)
fig.tight_layout()
name = "square_cheb_kindle" if KINDLE else "square_cheb"
fig.savefig(name + ".pdf"); fig.savefig(name + ".png", dpi=150)
