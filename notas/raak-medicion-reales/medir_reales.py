"""Medicion en la raiz de instancias reales: sumas de productos dentro de cada expresion.

Para cada objetivo y restriccion se recorre el arbol y se buscan los nodos 'suma' (combinaciones lineales maximales:
+, -, producto o division por constante) con al menos 2 terminos cuadraticos. Un termino es cuadratico si es un
producto de dos subexpresiones no constantes, A*B, un cuadrado, A^2, o una division A/B (= A * (1/B) en fAF2);
A y B pueden ser compuestas (por ejemplo sin(x)*y, x*ln(x) o x/(x+y)).

En la caja raiz cada factor se evalua en aritmetica afin: A = a0 + a^T eps + ea xi, con eps los simbolos de las
variables y xi un simbolo de error propio de la subexpresion A (el mismo si A aparece varias veces). Con
zeta = (eps, xi's), el residuo cuadratico de la suma es
    N = sum_k c_k (u_k^T zeta)(v_k^T zeta),
y lo demas (centro y terminos lineales, incluidos a0*eb*eta) es igual para todos los metodos. Se mide el ancho
del encierro de N con:

  S      suma por producto de la formula de fAF2
  SH     suma por producto del rango exacto de (u^T zeta)(v^T zeta)  (Skalna-Hladik; linea base)
  Qk0..  Q agregada de la suma: termino a termino (k=0) y k-nomial SDSOS (k=2,4,6) de ibex_KNomial.h
  EIG    cota espectral de la Q agregada
  EXACT  rango exacto de zeta^T Q zeta en [-1,1]^m si m <= 12; si no, una cota inferior (busqueda local)

Uso: python3 medir_reales.py <driver> <salida_nodos.csv> <bch>...
"""
import ast, csv, itertools, math, re, subprocess, sys
from pathlib import Path
import numpy as np

KS = (0, 2, 4, 6)
M_EXACT = 12


# ---------------------------------------------------------------- lectura de Minibex
def read_bch(path):
    txt = re.sub(r"//[^\n]*", "", Path(path).read_text())
    sec = re.split(r"^\s*(constants|variables|minimize|constraints|end)\b", txt, flags=re.M)
    parts = dict(zip(sec[1::2], sec[2::2]))
    consts = {}
    for m in re.finditer(r"(\w+)\s*=\s*([^;]+);", parts.get("constants", "")):
        consts[m.group(1)] = float(eval(m.group(2).replace("^", "**")))
    box = {}
    for m in re.finditer(r"(\w+)\s+in\s*\[([^,\]]+),([^\]]+)\]", parts["variables"]):
        box[m.group(1)] = (float(eval(m.group(2), {}, consts)), float(eval(m.group(3), {}, consts)))
    exprs = []
    obj = parts.get("minimize", "").strip().rstrip(";").strip()
    if obj:
        exprs.append(("obj", obj))
    for i, c in enumerate(s for s in parts.get("constraints", "").split(";") if s.strip()):
        lhs, rhs = re.split(r"<=|>=|=", c, maxsplit=1)
        exprs.append((f"c{i + 1}", f"({lhs}) - ({rhs})"))
    trees = []
    for name, e in exprs:
        e = " ".join(e.split()).replace("^", "**")
        e = re.sub(r"\bln\s*\(", "log(", e)
        trees.append((name, ast.parse(e, mode="eval").body))
    return box, consts, trees


# ---------------------------------------------------------------- aritmetica afin de los factores
class AF:
    """a0 + sum lin[v]*eps_v + err*xi (eps_v en [-1,1], ya escalado por el radio)."""
    def __init__(self, c, lin=None, err=0.0):
        self.c, self.lin, self.err = float(c), dict(lin or {}), float(err)
    def rad(self): return sum(map(abs, self.lin.values())) + self.err
    def rng(self): r = self.rad(); return self.c - r, self.c + r
    def add(self, o, s=1.0):
        lin = dict(self.lin)
        for k, v in o.lin.items(): lin[k] = lin.get(k, 0.0) + s * v
        return AF(self.c + s * o.c, lin, self.err + o.err)
    def scale(self, a): return AF(a * self.c, {k: a * v for k, v in self.lin.items()}, abs(a) * self.err)


def zono_vertices(u, v):
    g = [(a, b) if (b > 0 or (b == 0 and a > 0)) else (-a, -b) for a, b in zip(u, v) if a or b]
    if not g:
        return [(0.0, 0.0)]
    g.sort(key=lambda t: math.atan2(t[1], t[0]))
    p = np.array([-sum(t[0] for t in g), -sum(t[1] for t in g)])
    V = [p.copy()]
    for t in g: p = p + 2 * np.array(t); V.append(p.copy())
    for t in g: p = p - 2 * np.array(t); V.append(p.copy())
    return V


def product_range(u, v):
    """Rango exacto de (u^T z)(v^T z) en [-1,1]^m: x*y sobre el zonotopo plano (u^T z, v^T z)."""
    V = zono_vertices(u, v)
    lo, hi = math.inf, -math.inf
    for p, q in zip(V, V[1:] + V[:1]):
        d = q - p
        ts = [0.0, 1.0]
        den = 2 * d[0] * d[1]
        if den != 0:
            t = -(p[0] * d[1] + p[1] * d[0]) / den
            if 0 < t < 1: ts.append(t)
        for t in ts:
            x, y = p + t * d
            lo, hi = min(lo, x * y), max(hi, x * y)
    return lo, hi


def af_mul(x, y):
    keys = sorted(set(x.lin) | set(y.lin))
    u = [x.lin.get(k, 0.0) for k in keys] + [x.err, 0.0]
    v = [y.lin.get(k, 0.0) for k in keys] + [0.0, y.err]
    lo, hi = product_range(u, v)
    lin = {k: y.c * x.lin.get(k, 0.0) + x.c * y.lin.get(k, 0.0) for k in keys}
    return AF(x.c * y.c + (lo + hi) / 2, lin, abs(y.c) * x.err + abs(x.c) * y.err + (hi - lo) / 2)


def af_unary_minrange(x, f, df_end, lo, hi):
    """Min-range con pendiente alpha = df_end; f - alpha*t monotona en [lo, hi]."""
    a = df_end
    d1, d2 = f(lo) - a * lo, f(hi) - a * hi
    d1, d2 = min(d1, d2), max(d1, d2)
    r = x.scale(a)
    return AF(r.c + (d1 + d2) / 2, r.lin, r.err + (d2 - d1) / 2)


def af_secant(x, f, npts=4001):
    """Linealizacion por secante para funciones no monotonas o sin forma cerrada (sin, cos, x^k, abs):
    pendiente de la secante y residuo f - alpha t acotado por muestreo denso (medicion, no rigurosa).
    Conserva la parte lineal en eps, que es la que genera los terminos bilineales en un producto."""
    lo, hi = x.rng()
    if hi - lo <= 1e-12: return AF(f(x.c))
    a = (f(hi) - f(lo)) / (hi - lo)
    ts = np.linspace(lo, hi, npts)
    d = np.array([f(t) - a * t for t in ts])
    pad = 0.01 * (d.max() - d.min())  # margen por el muestreo
    d1, d2 = d.min() - pad, d.max() + pad
    r = x.scale(a)
    return AF(r.c + (d1 + d2) / 2, r.lin, r.err + (d2 - d1) / 2)


def af_interval(lo, hi):
    return AF((lo + hi) / 2, {}, (hi - lo) / 2)


def af_unary(name, x):
    lo, hi = x.rng()
    if name == "exp":
        return af_unary_minrange(x, math.exp, math.exp(lo), lo, hi) if hi - lo > 1e-12 else AF(math.exp(x.c))
    if name in ("log", "sqrt", "recip"):
        if lo <= 0 and name != "recip":
            lo = max(lo, 1e-300)
        f = {"log": math.log, "sqrt": math.sqrt, "recip": lambda t: 1 / t}[name]
        if hi - lo <= 1e-12: return AF(f(x.c))
        if name == "recip":
            if lo <= 0 <= hi: raise ValueError("division por intervalo que contiene 0")
            a = -1 / (lo * hi)  # pendiente de la secante: f - a t es monotona? usar intervalo seguro
            ts = [lo, hi, -math.sqrt(-1 / a) if lo < 0 else math.sqrt(-1 / a)]
            ds = [f(t) - a * t for t in ts if lo <= t <= hi]
            r = x.scale(a)
            return AF(r.c + (min(ds) + max(ds)) / 2, r.lin, r.err + (max(ds) - min(ds)) / 2)
        df = {"log": 1 / hi, "sqrt": 0.5 / math.sqrt(hi)}[name]
        return af_unary_minrange(x, f, df, lo, hi)
    if name in ("sin", "cos"):
        return af_secant(x, getattr(math, name))
    if name == "abs":
        if lo >= 0: return x
        if hi <= 0: return x.scale(-1)
        return af_secant(x, abs)
    raise ValueError(name)


def af_pow(x, k):
    if k == 2: return af_mul(x, x)
    return af_secant(x, lambda t: t ** k)


# ---------------------------------------------------------------- arbol
class Ev:
    def __init__(self, box, consts):
        self.box, self.consts, self.memo = box, consts, {}
        self.mid = {v: (a + b) / 2 for v, (a, b) in box.items()}
        self.r = {v: (b - a) / 2 for v, (a, b) in box.items()}

    def is_const(self, n):
        return not any(isinstance(s, ast.Name) and s.id in self.box for s in ast.walk(n))

    def const(self, n):
        return float(eval(compile(ast.Expression(n), "<c>", "eval"), {"exp": math.exp, "log": math.log, "sqrt": math.sqrt, "abs": abs, "sin": math.sin, "cos": math.cos}, self.consts))

    def af(self, n):
        key = ast.dump(n)
        if key in self.memo: return self.memo[key]
        if self.is_const(n): r = AF(self.const(n))
        elif isinstance(n, ast.Name): r = AF(self.mid[n.id], {n.id: self.r[n.id]})
        elif isinstance(n, ast.UnaryOp): r = self.af(n.operand).scale(-1 if isinstance(n.op, ast.USub) else 1)
        elif isinstance(n, ast.BinOp):
            if isinstance(n.op, (ast.Add, ast.Sub)):
                r = self.af(n.left).add(self.af(n.right), 1 if isinstance(n.op, ast.Add) else -1)
            elif isinstance(n.op, ast.Mult): r = af_mul(self.af(n.left), self.af(n.right))
            elif isinstance(n.op, ast.Div):
                d = self.af(n.right)
                r = self.af(n.left).scale(1 / d.c) if not d.lin and d.err == 0 else af_mul(self.af(n.left), af_unary("recip", d))
            elif isinstance(n.op, ast.Pow):
                k = self.const(n.right)
                if k == int(k) and k > 0: r = af_pow(self.af(n.left), int(k))
                else: r = af_unary("exp", af_mul(af_unary("log", self.af(n.left)), AF(k)))
            else: raise ValueError(ast.dump(n))
        elif isinstance(n, ast.Call): r = af_unary(n.func.id, self.af(n.args[0]))
        else: raise ValueError(ast.dump(n))
        self.memo[key] = r
        return r

    # combinacion lineal maximal: lista de (coef, termino)
    def lincomb(self, n, c=1.0, out=None):
        out = [] if out is None else out
        if self.is_const(n): out.append((c * self.const(n), None)); return out
        if isinstance(n, ast.BinOp) and isinstance(n.op, (ast.Add, ast.Sub)):
            self.lincomb(n.left, c, out); self.lincomb(n.right, c if isinstance(n.op, ast.Add) else -c, out)
        elif isinstance(n, ast.UnaryOp):
            self.lincomb(n.operand, -c if isinstance(n.op, ast.USub) else c, out)
        elif isinstance(n, ast.BinOp) and isinstance(n.op, ast.Mult) and self.is_const(n.left):
            self.lincomb(n.right, c * self.const(n.left), out)
        elif isinstance(n, ast.BinOp) and isinstance(n.op, ast.Mult) and self.is_const(n.right):
            self.lincomb(n.left, c * self.const(n.right), out)
        elif isinstance(n, ast.BinOp) and isinstance(n.op, ast.Div) and self.is_const(n.right):
            self.lincomb(n.left, c / self.const(n.right), out)
        else: out.append((c, n))
        return out

    def quad_factors(self, t):
        if isinstance(t, ast.BinOp) and isinstance(t.op, ast.Mult): return t.left, t.right
        if isinstance(t, ast.BinOp) and isinstance(t.op, ast.Div):  # fAF2: A/B = A * (1/B)
            return t.left, ast.Call(func=ast.Name("recip", ast.Load()), args=[t.right], keywords=[])
        if isinstance(t, ast.BinOp) and isinstance(t.op, ast.Pow) and self.is_const(t.right) and self.const(t.right) == 2:
            return t.left, t.left
        return None

    def sum_nodes(self, n, out):
        """Todos los nodos suma maximales del arbol, con sus terminos."""
        terms = [(c, t) for c, t in self.lincomb(n) if t is not None]
        if len(terms) > 1 or (terms and terms[0][1] is not n):
            out.append(terms)
            children = [t for _, t in terms]
        else:
            children = list(ast.iter_child_nodes(n))
        for t in children:
            for ch in ast.iter_child_nodes(t):
                if isinstance(ch, ast.expr) and not self.is_const(ch): self.sum_nodes(ch, out)
        return out


# ---------------------------------------------------------------- cotas
def quad_range_exact(Q):
    n = Q.shape[0]
    lo, hi = math.inf, -math.inf
    for mask in range(1 << n):
        F = [i for i in range(n) if mask >> i & 1]
        B = [i for i in range(n) if not mask >> i & 1]
        XB = np.array(list(itertools.product((-1.0, 1.0), repeat=len(B)))).T if B else np.zeros((0, 1))
        X = np.zeros((n, XB.shape[1])); X[B] = XB
        if F:
            QFF, QFB = Q[np.ix_(F, F)], Q[np.ix_(F, B)]
            rhs = -QFB @ XB if B else np.zeros((len(F), 1))
            XF = np.linalg.lstsq(QFF, rhs, rcond=None)[0]
            sc = max(1.0, np.abs(Q).max())
            ok = (np.abs(QFF @ XF - rhs).max(axis=0) <= 1e-9 * sc) & (np.abs(XF).max(axis=0) <= 1 + 1e-12)
            X[F] = XF; X = X[:, ok]
            if X.shape[1] == 0: continue
        vals = np.einsum("ij,ik,kj->j", X, Q, X)
        lo, hi = min(lo, vals.min()), max(hi, vals.max())
    return lo, hi


def quad_range_inner(Q, starts=300, seed=0):
    """Cota interior del rango (busqueda local por coordenadas desde vertices aleatorios)."""
    rng = np.random.default_rng(seed)
    n = Q.shape[0]
    best = [math.inf, -math.inf]
    for sgn in (1, -1):
        for _ in range(starts):
            x = rng.choice([-1.0, 1.0], size=n)
            for _ in range(50):
                old = x.copy()
                for i in range(n):
                    g = Q[i] @ x - Q[i, i] * x[i]  # f = Qii xi^2 + 2 g xi + ...
                    cand = [-1.0, 1.0]
                    if sgn * Q[i, i] < 0: cand.append(float(np.clip(-g / Q[i, i], -1, 1)))
                    x[i] = max(cand, key=lambda t: sgn * (Q[i, i] * t * t + 2 * g * t))
                if np.array_equal(old, x): break
            val = x @ Q @ x
            best[0], best[1] = min(best[0], val), max(best[1], val)
    return best[0], best[1]


def run_driver(driver, lines):
    out = subprocess.run([driver], input="\n".join(lines) + "\n", capture_output=True, text=True, check=True)
    return [tuple(map(float, l.split())) for l in out.stdout.strip().splitlines()]


def measure_node(ev, terms, driver):
    quads = []
    for c, t in terms:
        f = ev.quad_factors(t)
        if f: quads.append((c, f[0], f[1]))
    if len(quads) < 2: return None
    syms = {}
    vecs = []
    for c, A, B in quads:
        uv = []
        for X in (A, B):
            fx = ev.af(X)
            d = dict(fx.lin)
            if fx.err > 0: d["E:" + ast.unparse(X)] = fx.err
            uv.append(d)
            for k in d: syms.setdefault(k, len(syms))
        vecs.append((c, uv[0], uv[1]))
    m = len(syms)
    U = np.zeros((len(vecs), m)); V = np.zeros((len(vecs), m))
    for i, (c, du, dv) in enumerate(vecs):
        for k, val in du.items(): U[i, syms[k]] = c * val
        for k, val in dv.items(): V[i, syms[k]] = val
    Q = sum((np.outer(u, v) + np.outer(v, u)) / 2 for u, v in zip(U, V))
    # se comparten simbolos si algun simbolo aparece en dos terminos distintos
    sup = [set(np.nonzero(u)[0]) | set(np.nonzero(v)[0]) for u, v in zip(U, V)]
    shared = any(sup[i] & sup[j] for i in range(len(sup)) for j in range(i + 1, len(sup)))
    r = {"K": len(vecs), "m": m, "comparte": shared, "rango_Q": int(np.linalg.matrix_rank(Q, tol=1e-9 * max(1, np.abs(Q).max())))}
    r["S"] = sum(2 * (np.abs(u).sum() * np.abs(v).sum() - 0.5 * np.abs(u * v).sum()) for u, v in zip(U, V))
    r["SH"] = sum(np.subtract(*product_range(u, v)[::-1]) for u, v in zip(U, V))
    res = run_driver(driver, [f"Q {m} {k} " + " ".join(f"{x:.17g}" for x in Q.ravel()) for k in KS])
    for k, (L, Uu, nb, kmax) in zip(KS, res):
        r[f"Qk{k}"] = Uu - L
    ev_ = np.linalg.eigvalsh(Q)
    r["EIG"] = m * (max(ev_[-1], 0) - min(ev_[0], 0))
    if m <= M_EXACT:
        lo, hi = quad_range_exact(Q); r["exacto"] = True
    else:
        lo, hi = quad_range_inner(Q); r["exacto"] = False
    r["EXACT"] = hi - lo
    for key in ("S", "SH", "EIG") + tuple(f"Qk{k}" for k in KS):
        assert r[key] >= r["EXACT"] * (1 - 1e-9) - 1e-9, (key, r[key], r["EXACT"])
    return r


def main():
    driver, out_csv, files = sys.argv[1], sys.argv[2], sys.argv[3:]
    rows = []
    for f in files:
        box, consts, trees = read_bch(f)
        ev = Ev(box, consts)
        n_nodes = 0
        for name, tree in trees:
            for terms in ev.sum_nodes(tree, []):
                try:
                    r = measure_node(ev, terms, driver)
                except (ValueError, OverflowError, ZeroDivisionError) as exc:
                    print(f"  {Path(f).stem} {name}: omitido ({exc})", flush=True); continue
                if r is None: continue
                n_nodes += 1
                r = {"instancia": Path(f).stem, "n": len(box), "expr": name, **r}
                rows.append(r)
        print(f"{Path(f).stem:20s} n={len(box):3d} expresiones={len(trees):3d} sumas con >=2 productos={n_nodes}", flush=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
