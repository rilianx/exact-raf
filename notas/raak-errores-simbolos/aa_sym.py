"""Evaluador afin de expresiones Minibex con tres modos, para medir el efecto de guardar los errores como simbolos.

  M0  como fAF2 hoy: forma = c + sum a_i eps_i + e*xi con e escalar (los errores de distintas subexpresiones son
      independientes y se acumulan en un solo termino); productos con SH (rango exacto de la parte cuadratica,
      con los errores de los factores como generadores propios) y los cruces centro*error sumados en valor absoluto.
  M1  errores como simbolos: cada operacion no lineal (producto, funcion, potencia) crea un simbolo nuevo para su
      error de aproximacion y los simbolos se propagan linealmente (como la aritmetica afin estandar AF1/Affine3).
      Productos con SH sobre TODOS los simbolos; los cruces centro*error quedan lineales y exactos.
  M2  M1 y ademas cada suma de >= 2 productos se acota en conjunto: la parte cuadratica agregada
      sum_k c_k (u_k^T z)(v_k^T z) se encierra con la menor de {suma de SH, bloques producto, termino a termino}.

Las subexpresiones repetidas se evaluan una vez por expresion (como el DAG de Ibex), asi que en M1/M2 comparten
sus simbolos. Funciones: exp, log/ln, sqrt, 1/x (min-range), sin, cos, abs, x^k (secante con residuo muestreado;
medicion, no rigurosa). Las tres variantes usan las mismas linealizaciones unarias.
"""
import ast, itertools, math, re
from pathlib import Path
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "raak-qplib"))
import prodk_fast as PF  # noqa: E402


# ------------------------------------------------------------------ lectura (Minibex con vectores y constantes)
def read_bch(path):
    txt = re.sub(r"//[^\n]*", "", Path(path).read_text())
    sec = re.split(r"^\s*(constants|variables|minimize|constraints|end)\b", txt, flags=re.M | re.I)
    parts = {}
    for k, v in zip(sec[1::2], sec[2::2]): parts[k.lower()] = v
    consts = {}
    for m in re.finditer(r"(\w+)\s*(?:=|in)\s*([^;\[]+);", parts.get("constants", "")):
        consts[m.group(1)] = float(eval(m.group(2).replace("^", "**"), {}, consts))
    box = {}
    for m in re.finditer(r"(\w+)(?:\[(\d+)\])?\s+in\s*\[([^,\]]+),([^\]]+)\]", parts["variables"]):
        lo, hi = float(eval(m.group(3), {}, consts)), float(eval(m.group(4), {}, consts))
        if m.group(2):
            for i in range(1, int(m.group(2)) + 1): box[f"{m.group(1)}_{i}"] = (lo, hi)
        else:
            box[m.group(1)] = (lo, hi)
    vecs = {re.sub(r"_\d+$", "", v) for v in box if re.search(r"_\d+$", v)}
    exprs = []
    obj = parts.get("minimize", "").strip().rstrip(";").strip()
    if obj: exprs.append(("obj", obj))
    for i, c in enumerate(s for s in parts.get("constraints", "").split(";") if s.strip()):
        lhs, rhs = re.split(r"<=|>=|=", c, maxsplit=1)
        exprs.append((f"c{i + 1}", f"({lhs}) - ({rhs})"))
    trees = []
    for name, e in exprs:
        e = " ".join(e.split()).replace("^", "**")
        for v in vecs: e = re.sub(rf"\b{v}\s*[\(\[]\s*(\d+)\s*[\)\]]", rf"{v}_\1", e)
        e = re.sub(r"\bln\s*\(", "log(", e)
        trees.append((name, ast.parse(e, mode="eval").body))
    return box, consts, trees


# ------------------------------------------------------------------ formas afines
class AF:
    __slots__ = ("c", "lin", "err")
    def __init__(self, c, lin=None, err=0.0):
        self.c, self.lin, self.err = float(c), dict(lin or {}), float(err)
    def rad(self): return sum(map(abs, self.lin.values())) + self.err
    def rng(self): r = self.rad(); return self.c - r, self.c + r
    def add(self, o, s=1.0):
        lin = dict(self.lin)
        for k, v in o.lin.items(): lin[k] = lin.get(k, 0.0) + s * v
        return AF(self.c + s * o.c, lin, self.err + o.err)
    def scale(self, a): return AF(a * self.c, {k: a * v for k, v in self.lin.items()}, abs(a) * self.err)


def product_range(u, v):
    w = PF.sh_width(u, v)  # ancho; tambien hace falta el centro
    V = PF.zono_vertices(u, v)
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


class Evaluator:
    def __init__(self, box, consts, mode):
        self.box, self.consts, self.mode = box, consts, mode
        self.mid = {v: (a + b) / 2 for v, (a, b) in box.items()}
        self.r = {v: (b - a) / 2 for v, (a, b) in box.items()}
        self.memo, self.nsym = {}, 0

    # -- simbolos de error
    def new_err(self, form, e, tag="otro"):
        """Agrega un error de magnitud e >= 0 a la forma: escalar (M0) o simbolo nuevo (M1, M2). En M1/M2 el
        simbolo lleva la etiqueta de la operacion que lo creo ('#e<k>:<tag>'), para el presupuesto de error."""
        if e <= 0: return form
        if self.mode == "M0": return AF(form.c, form.lin, form.err + e)
        self.nsym += 1
        lin = dict(form.lin); lin[f"#e{self.nsym}:{tag}"] = e
        return AF(form.c, lin, form.err)

    def is_const(self, n):
        return not any(isinstance(s, ast.Name) and s.id in self.box for s in ast.walk(n))

    def const(self, n):
        env = {"exp": math.exp, "log": math.log, "sqrt": math.sqrt, "abs": abs, "sin": math.sin, "cos": math.cos}
        return float(eval(compile(ast.Expression(n), "<c>", "eval"), env, self.consts))

    # -- operaciones
    def mul(self, x, y):
        keys = sorted(set(x.lin) | set(y.lin))
        u = [x.lin.get(k, 0.0) for k in keys] + [x.err, 0.0]
        v = [y.lin.get(k, 0.0) for k in keys] + [0.0, y.err]
        lo, hi = product_range(u, v)
        lin = {k: y.c * x.lin.get(k, 0.0) + x.c * y.lin.get(k, 0.0) for k in keys}
        f = AF(x.c * y.c + (lo + hi) / 2, lin, abs(y.c) * x.err + abs(x.c) * y.err)
        haserr = lambda z: z.err > 0 or any(k.startswith("#") for k in z.lin)
        tag = ("sqr" if x is y else "mul") + ("_comp" if haserr(x) or haserr(y) else "_lin")
        return self.new_err(f, (hi - lo) / 2, tag)

    def lin_approx(self, x, a, d1, d2, tag="otro"):
        r = x.scale(a)
        return self.new_err(AF(r.c + (d1 + d2) / 2, r.lin, r.err), (d2 - d1) / 2, tag)

    def unary(self, name, x):
        lo, hi = x.rng()
        if hi - lo <= 1e-15 * max(1, abs(lo)):
            f = {"exp": math.exp, "log": math.log, "sqrt": math.sqrt, "recip": lambda t: 1 / t,
                 "sin": math.sin, "cos": math.cos, "abs": abs}[name]
            return AF(f(x.c))
        if name == "exp":
            a = math.exp(lo); return self.lin_approx(x, a, *sorted((math.exp(lo) - a * lo, math.exp(hi) - a * hi)), tag="exp")
        if name in ("log", "sqrt"):
            if hi <= 0: raise ValueError(f"{name} de intervalo no positivo")
            lo = max(lo, 1e-300)
            f, a = (math.log, 1 / hi) if name == "log" else (math.sqrt, 0.5 / math.sqrt(hi))
            return self.lin_approx(x, a, *sorted((f(lo) - a * lo, f(hi) - a * hi)), tag=name)
        if name == "recip":
            if lo <= 0 <= hi: raise ValueError("division por intervalo que contiene 0")
            a = -1 / (lo * hi)
            t0 = math.copysign(math.sqrt(lo * hi), lo)
            ds = [1 / t - a * t for t in (lo, hi, t0) if lo <= t <= hi]
            return self.lin_approx(x, a, min(ds), max(ds), tag="recip")
        if name == "abs":
            if lo >= 0: return x
            if hi <= 0: return x.scale(-1)
        f = {"sin": math.sin, "cos": math.cos, "abs": abs}.get(name)
        if f is None: raise ValueError(name)
        return self.secant(x, f, lo, hi, tag="sincos" if name in ("sin", "cos") else name)

    def secant(self, x, f, lo, hi, npts=2001, tag="otro"):
        a = (f(hi) - f(lo)) / (hi - lo)
        ts = np.linspace(lo, hi, npts)
        d = np.array([f(t) for t in ts]) - a * ts
        pad = 0.01 * (d.max() - d.min()) + 1e-15
        return self.lin_approx(x, a, d.min() - pad, d.max() + pad, tag)

    def power(self, x, k):
        if k == 2: return self.mul(x, x)
        lo, hi = x.rng()
        if hi - lo <= 0: return AF(x.c ** k)
        return self.secant(x, lambda t: t ** k, lo, hi, tag="pow")

    # -- combinaciones lineales y sumas de productos
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
        if isinstance(t, ast.BinOp) and isinstance(t.op, ast.Mult): return t.left, t.right, False
        if isinstance(t, ast.BinOp) and isinstance(t.op, ast.Div): return t.left, t.right, True
        if isinstance(t, ast.BinOp) and isinstance(t.op, ast.Pow) and self.is_const(t.right) and self.const(t.right) == 2:
            return t.left, t.left, False
        return None

    def sum_of_products(self, terms):
        """M2: parte lineal exacta de cada producto y parte cuadratica agregada con una sola cota."""
        const, acc, quads = 0.0, AF(0.0), []
        for c, t in terms:
            if t is None: const += c; continue
            qf = self.quad_factors(t) if t is not None else None
            if qf is None: acc = acc.add(self.ev(t).scale(c)); continue
            A, B, div = qf
            x = self.ev(A); y = self.ev(B)
            if div: y = self.unary("recip", y)
            quads.append((c, x, y))
        if len(quads) < 2:
            for c, x, y in quads: acc = acc.add(self.mul(x, y).scale(c))
            return acc.add(AF(const))
        keys = sorted(set().union(*[set(x.lin) | set(y.lin) for _, x, y in quads]))
        m = len(keys) + 2 * len(quads)                  # + un generador por error escalar (solo M0 los tiene)
        U, V, sh_sum = [], [], 0.0
        for k, (c, x, y) in enumerate(quads):
            u = np.zeros(m); v = np.zeros(m)
            for j, s in enumerate(keys):
                u[j] = c * x.lin.get(s, 0.0); v[j] = y.lin.get(s, 0.0)
            u[len(keys) + 2 * k] = c * x.err; v[len(keys) + 2 * k + 1] = y.err
            U.append(u); V.append(v)
            lin = {s: c * (y.c * x.lin.get(s, 0.0) + x.c * y.lin.get(s, 0.0)) for s in keys}
            acc = acc.add(AF(c * x.c * y.c, lin, abs(c) * (abs(y.c) * x.err + abs(x.c) * y.err)))
        ranges = [product_range(u, v) for u, v in zip(U, V)]
        L, H = sum(r[0] for r in ranges), sum(r[1] for r in ranges)          # suma de SH
        nz = [j for j in range(m) if any(u[j] or v[j] for u, v in zip(U, V))]
        Q = sum((np.outer(u, v) + np.outer(v, u)) / 2 for u, v in zip(U, V))[np.ix_(nz, nz)]
        l2, h2 = PF.tt_interval(Q)                                           # termino a termino
        L, H = max(L, l2), min(H, h2)
        if len(nz) <= 40:                                                    # bloques producto
            l3, h3, _ = PF.prodk_interval(Q)
            L, H = max(L, l3), min(H, h3)
        best_c, best_r = (L + H) / 2, max(0.0, (H - L) / 2)
        return self.new_err(acc.add(AF(const + best_c)), best_r, "sumprod")

    def ev(self, n):
        key = ast.dump(n)
        if key in self.memo: return self.memo[key]
        if self.is_const(n): r = AF(self.const(n))
        elif isinstance(n, ast.Name): r = AF(self.mid[n.id], {n.id: self.r[n.id]})
        elif self.mode == "M2" and isinstance(n, (ast.BinOp, ast.UnaryOp)) and \
                (isinstance(n, ast.UnaryOp) or isinstance(n.op, (ast.Add, ast.Sub))):
            r = self.sum_of_products(self.lincomb(n))
        elif isinstance(n, ast.UnaryOp): r = self.ev(n.operand).scale(-1 if isinstance(n.op, ast.USub) else 1)
        elif isinstance(n, ast.BinOp):
            if isinstance(n.op, (ast.Add, ast.Sub)):
                r = self.ev(n.left).add(self.ev(n.right), 1 if isinstance(n.op, ast.Add) else -1)
            elif isinstance(n.op, ast.Mult):
                if self.is_const(n.left): r = self.ev(n.right).scale(self.const(n.left))
                elif self.is_const(n.right): r = self.ev(n.left).scale(self.const(n.right))
                else: r = self.mul(self.ev(n.left), self.ev(n.right))
            elif isinstance(n.op, ast.Div):
                if self.is_const(n.right): r = self.ev(n.left).scale(1 / self.const(n.right))
                else: r = self.mul(self.ev(n.left), self.unary("recip", self.ev(n.right)))
            elif isinstance(n.op, ast.Pow):
                k = self.const(n.right)
                if k == int(k) and k > 0: r = self.power(self.ev(n.left), int(k))
                else: r = self.unary("exp", self.ev(n.left).scale(0) if False else self.mul(self.unary("log", self.ev(n.left)), AF(k)))
            else: raise ValueError(ast.dump(n))
        elif isinstance(n, ast.Call): r = self.unary(n.func.id, self.ev(n.args[0]))
        else: raise ValueError(ast.dump(n))
        self.memo[key] = r
        return r


def enclose(tree, box, consts, mode):
    return Evaluator(box, consts, mode).ev(tree).rng()


# ------------------------------------------------------------------ rango interior (referencia)
def inner_range(tree, box, consts, samples=400, seed=0):
    names = list(box)
    code = compile(ast.Expression(tree), "<f>", "eval")
    env = {"exp": np.exp, "log": np.log, "sqrt": np.sqrt, "abs": np.abs, "sin": np.sin, "cos": np.cos}
    lo = np.array([box[v][0] for v in names]); hi = np.array([box[v][1] for v in names])
    rng = np.random.default_rng(seed)
    def f(X):
        with np.errstate(all="ignore"):
            val = eval(code, env, {**consts, **{v: X[:, i] for i, v in enumerate(names)}})
        return np.broadcast_to(np.asarray(val, dtype=float), (X.shape[0],))
    X = lo + (hi - lo) * rng.random((samples, len(names)))
    corners = lo + (hi - lo) * rng.integers(0, 2, size=(samples, len(names)))
    X = np.vstack([X, corners, (lo + hi)[None, :] / 2])
    out = []
    for sgn in (1.0, -1.0):
        y = sgn * f(X)
        y[~np.isfinite(y)] = -np.inf
        idx = np.argsort(-y)[:5]
        best = y[idx[0]]
        for i in idx:                           # busqueda por patron desde los 5 mejores
            x = X[i].copy(); fx = y[i]; step = (hi - lo) / 4
            for _ in range(60):
                cand = np.clip(x + step * rng.choice([-1, 0, 1], size=(2 * len(names) + 4, len(names))), lo, hi)
                fc = sgn * f(cand); fc[~np.isfinite(fc)] = -np.inf
                j = int(np.argmax(fc))
                if fc[j] > fx: x, fx = cand[j], fc[j]
                else: step = step / 2
            best = max(best, fx)
        out.append(sgn * best)
    return out[1], out[0]
