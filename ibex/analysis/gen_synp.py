"""Familia synp: f(x) = sum_k q_k(x) r_k(x) + lineal, con q_k, r_k cuadraticas con constante.
Ambos factores tienen termino de error (por los cuadrados y productos) y centro distinto de 0,
que es el caso en que la multiplicacion exacta puede ganar a Skalna-Hladik (Corolario 11)."""
import random, sys, os
out = sys.argv[1]
def lin(rng, n, k):
    idx = rng.sample(range(1, n+1), k)
    return " + ".join(f"({rng.gauss(0,1):.3f})*x{i}" for i in idx)
def quad(rng, n):
    i, j = rng.sample(range(1, n+1), 2)
    c = rng.choice([-1, 1]) * rng.uniform(0.5, 2.0)
    return f"({c:.3f} + {lin(rng,n,2)} + ({rng.gauss(0,1):.3f})*x{i}^2 + ({rng.gauss(0,1):.3f})*x{i}*x{j})"
for n in (3, 5, 8):
    for s in range(1, 7):
        rng = random.Random(1000*n + s)
        K = n
        obj = " + ".join(f"{quad(rng,n)}*{quad(rng,n)}" for _ in range(K)) + " + " + lin(rng, n, n)
        cons = f"{lin(rng,n,n)} <= {rng.uniform(0.2,0.8):.3f};"
        var = "\n".join(f"x{i} in [-1, 1];" for i in range(1, n+1))
        open(os.path.join(out, f"synp_n{n:02d}_s{s}.bch"), "w").write(
            f"variables\n{var}\n\nminimize\n{obj};\n\nconstraints\n{cons}\n\nend\n")
