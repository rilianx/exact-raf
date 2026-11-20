"""M3: aritmetica de formas CUADRATICAS  F = c + a^T z + z^T Q z  (+ simbolos de error nuevos, como M1).

z = (eps de las variables, simbolos de error). Respecto de M1 (aa_sym.py), la diferencia es que el producto
conserva todo el grado 2:
    X*Y = cx cy + (cx b + cy a)^T z + z^T (cx Qy + cy Qx + sym(a b^T)) z  +  R3,
    R3  = (a^T z)(z^T Qy z) + (b^T z)(z^T Qx z) + (z^T Qx z)(z^T Qy z)   (grado 3 y 4)
y solo R3 se acota (con intervalos) y pasa a un simbolo de error nuevo. Las funciones usan la misma linealizacion
que M1 (min-range o secante) sobre el rango de la forma completa: f(F) ~ alpha F + beta +- delta, asi que la
parte cuadratica se conserva escalada.

Rango de una forma (para las funciones y para el encierro final 'M3a'): para cada simbolo i, a_i z_i + Q_ii z_i^2
tiene rango exacto en [-1,1]; las entradas fuera de la diagonal, termino a termino. Es O(terminos).
Encierro final 'M3b': la menor de M3a, bloques producto (PRODK) y SDP de Shor sobre la forma homogeneizada
[[Q, a/2], [a^T/2, 0]] en [-1,1]^(m+1), que tiene el mismo rango que la forma (es lineal en la variable extra t,
y t = -1 equivale a cambiar z por -z).
"""
import ast, math
from pathlib import Path
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "raak-qplib"))
import aa_sym as AS  # noqa: E402
import prodk_fast as PF  # noqa: E402


class QF:
    __slots__ = ("c", "lin", "Q")
    def __init__(self, c, lin=None, Q=None):
        self.c, self.lin, self.Q = float(c), dict(lin or {}), dict(Q or {})   # Q: {(i, j): v} con i <= j

    def add(self, o, s=1.0):
        lin = dict(self.lin)
        for k, v in o.lin.items(): lin[k] = lin.get(k, 0.0) + s * v
        Q = dict(self.Q)
        for k, v in o.Q.items(): Q[k] = Q.get(k, 0.0) + s * v
        return QF(self.c + s * o.c, lin, Q)

    def scale(self, a):
        return QF(a * self.c, {k: a * v for k, v in self.lin.items()}, {k: a * v for k, v in self.Q.items()})

    def quad_interval(self):
        """Rango (cota) de z^T Q z solo: diagonal [min(0,q), max(0,q)], fuera de la diagonal +-|q| (Q ya simetrica
        guardada como i<j con el coeficiente total del monomio z_i z_j)."""
        lo = hi = 0.0
        for (i, j), v in self.Q.items():
            if i == j: lo += min(0.0, v); hi += max(0.0, v)
            else: lo -= abs(v); hi += abs(v)
        return lo, hi

    def rng(self):
        """c + sum_i rango(a_i z + Q_ii z^2) + sum_{i<j} +-|Q_ij|."""
        lo = hi = self.c
        diag = {i: v for (i, j), v in self.Q.items() if i == j}
        for k in set(self.lin) | set(diag):
            a, q = self.lin.get(k, 0.0), diag.get(k, 0.0)
            vals = [q + a, q - a]
            if q != 0 and abs(a / (2 * q)) < 1: vals.append(-a * a / (4 * q))
            lo += min(vals); hi += max(vals)
        for (i, j), v in self.Q.items():
            if i != j: lo -= abs(v); hi += abs(v)
        return lo, hi

    def lin_abs(self): return sum(map(abs, self.lin.values()))


def key(i, j): return (i, j) if i <= j else (j, i)


class EvalQ(AS.Evaluator):
    """Reutiliza la estructura de aa_sym.Evaluator (memo, lincomb, constantes) con formas cuadraticas.
    centered=True: en el producto, cada parte cuadratica q = z^T Q z se escribe q = m + qhat con m el centro de su
    rango [lo, hi] y qhat en [-r, r]; m pasa a la constante (exacta) y el resto de grado 3 y 4 se acota con qhat:
        X*Y = Cx Cy - Cx my - Cy mx + (Cx b + Cy a)^T z + z^T (Cx Qy + Cy Qx + sym(a b^T)) z
              + [ (a^T z) qhat_y + qhat_x (b^T z) + qhat_x qhat_y ],     Cx = cx + mx, Cy = cy + my,
    y para X*X el resto es 2 (a^T z) qhat_x + qhat_x^2, con qhat_x^2 en [0, r_x^2]."""
    def __init__(self, box, consts, centered=False):
        super().__init__(box, consts, "M1")
        self.centered = centered

    def new_errQ(self, F, e, tag):
        if e <= 0: return F
        self.nsym += 1
        lin = dict(F.lin); lin[f"#e{self.nsym}:{tag}"] = e
        return QF(F.c, lin, F.Q)

    def mulQ(self, X, Y):
        same = X is Y
        qx, qy = X.quad_interval(), Y.quad_interval()
        if self.centered:
            mx, my = (qx[0] + qx[1]) / 2, (qy[0] + qy[1]) / 2
            rx, ry = (qx[1] - qx[0]) / 2, (qy[1] - qy[0]) / 2
        else:
            mx = my = 0.0
        Cx, Cy = X.c + mx, Y.c + my
        lin = {k: Cy * X.lin.get(k, 0.0) + Cx * Y.lin.get(k, 0.0) for k in set(X.lin) | set(Y.lin)}
        Q = {}
        for k, v in X.Q.items(): Q[k] = Q.get(k, 0.0) + Cy * v
        for k, v in Y.Q.items(): Q[k] = Q.get(k, 0.0) + Cx * v
        for i, a in X.lin.items():                      # (a^T z)(b^T z) completo, en grado 2
            for j, b in Y.lin.items():
                kk = key(i, j); Q[kk] = Q.get(kk, 0.0) + a * b
        sa, sb = X.lin_abs(), Y.lin_abs()
        iv = lambda a, b: (min(a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1]), max(a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1]))
        if self.centered:
            if same:
                r = [(-2 * sa * rx, 2 * sa * rx), (0.0, rx * rx)]
            else:
                r = [(-sa * ry, sa * ry), (-sb * rx, sb * rx), (-rx * ry, rx * ry)]
            const = Cx * Cy - Cx * my - Cy * mx
        else:                                           # version original: rangos completos de q
            r = [iv((-sa, sa), qy), iv((-sb, sb), qx), iv(qx, qy)]
            const = X.c * Y.c
        lo, hi = sum(t[0] for t in r), sum(t[1] for t in r)
        F = QF(const + (lo + hi) / 2, lin, Q)
        haserr = lambda Z: bool(Z.Q) or any(k.startswith("#") for k in Z.lin)
        tag = ("sqr" if same else "mul") + ("_comp" if haserr(X) or haserr(Y) else "_lin")
        return self.new_errQ(F, (hi - lo) / 2, tag)

    def unaryQ(self, name, X):
        lo, hi = X.rng()
        if hi - lo <= 1e-15 * max(1, abs(lo)):
            f = {"exp": math.exp, "log": math.log, "sqrt": math.sqrt, "recip": lambda t: 1 / t,
                 "sin": math.sin, "cos": math.cos, "abs": abs}[name]
            return QF(f(X.c))
        # misma linealizacion que M1: se obtiene (alpha, d1, d2) aplicando aa_sym a una forma afin con ese rango
        probe = AS.AF((lo + hi) / 2, {"_t": (hi - lo) / 2})
        if name == "abs" and lo >= 0: return X
        if name == "abs" and hi <= 0: return X.scale(-1)
        F1 = AS.Evaluator.unary(self, name, probe) if name != "pow" else None
        alpha = F1.lin.get("_t", 0.0) / ((hi - lo) / 2)
        delta = sum(abs(v) for k, v in F1.lin.items() if k != "_t") + F1.err
        beta = F1.c - alpha * probe.c
        tagk = [k for k in F1.lin if k.startswith("#")]
        tag = tagk[0].split(":")[1] if tagk else name
        R = X.scale(alpha); R = QF(R.c + beta, R.lin, R.Q)
        return self.new_errQ(R, delta, tag)

    def powQ(self, X, k):
        if k == 2: return self.mulQ(X, X)
        lo, hi = X.rng()
        if hi - lo <= 0: return QF(X.c ** k)
        probe = AS.AF((lo + hi) / 2, {"_t": (hi - lo) / 2})
        F1 = AS.Evaluator.secant(self, probe, lambda t: t ** k, lo, hi, tag="pow")
        alpha = F1.lin.get("_t", 0.0) / ((hi - lo) / 2)
        delta = sum(abs(v) for kk, v in F1.lin.items() if kk != "_t") + F1.err
        R = X.scale(alpha); R = QF(R.c + F1.c - alpha * probe.c, R.lin, R.Q)
        return self.new_errQ(R, delta, "pow")

    def evq(self, n):
        k = ast.dump(n)
        if k in self.memo: return self.memo[k]
        if self.is_const(n): r = QF(self.const(n))
        elif isinstance(n, ast.Name): r = QF(self.mid[n.id], {n.id: self.r[n.id]})
        elif isinstance(n, ast.UnaryOp): r = self.evq(n.operand).scale(-1 if isinstance(n.op, ast.USub) else 1)
        elif isinstance(n, ast.BinOp):
            if isinstance(n.op, (ast.Add, ast.Sub)):
                r = self.evq(n.left).add(self.evq(n.right), 1 if isinstance(n.op, ast.Add) else -1)
            elif isinstance(n.op, ast.Mult):
                if self.is_const(n.left): r = self.evq(n.right).scale(self.const(n.left))
                elif self.is_const(n.right): r = self.evq(n.left).scale(self.const(n.right))
                else:
                    a = self.evq(n.left); b = a if ast.dump(n.left) == ast.dump(n.right) else self.evq(n.right)
                    r = self.mulQ(a, b)
            elif isinstance(n.op, ast.Div):
                if self.is_const(n.right): r = self.evq(n.left).scale(1 / self.const(n.right))
                else: r = self.mulQ(self.evq(n.left), self.unaryQ("recip", self.evq(n.right)))
            elif isinstance(n.op, ast.Pow):
                kk = self.const(n.right)
                if kk == int(kk) and kk > 0: r = self.powQ(self.evq(n.left), int(kk))
                else: r = self.unaryQ("exp", self.mulQ(self.unaryQ("log", self.evq(n.left)), QF(kk)))
            else: raise ValueError(ast.dump(n))
        elif isinstance(n, ast.Call): r = self.unaryQ(n.func.id, self.evq(n.args[0]))
        else: raise ValueError(ast.dump(n))
        self.memo[k] = r
        return r


def homog(F):
    syms = sorted(set(F.lin) | {i for k in F.Q for i in k})
    idx = {s: t for t, s in enumerate(syms)}
    m = len(syms)
    H = np.zeros((m + 1, m + 1))
    for (i, j), v in F.Q.items():
        a, b = idx[i], idx[j]
        if a == b: H[a, a] += v
        else: H[a, b] += v / 2; H[b, a] += v / 2
    for s, v in F.lin.items():
        H[idx[s], m] += v / 2; H[m, idx[s]] += v / 2
    return H


def final_all(F, driver=None, prodk_max=40, sdp_max=60):
    """Anchos de los encierros finales de F: 'tt' (barato), 'prodk' (min con bloques producto), 'kn6' (min con el
    k-nomial k=6 de ibex_KNomial.h, via el driver en C++), 'sdp' (min con Shor). NaN si no se calculo."""
    lo, hi = F.rng(); wa = hi - lo
    out = {"tt": wa, "prodk": wa, "kn6": wa, "sdp": wa}
    if not F.Q: return out
    H = homog(F); m1 = H.shape[0]
    if m1 <= prodk_max:
        L, U, _ = PF.prodk_interval(H); out["prodk"] = min(wa, U - L)
    if driver is not None:
        import subprocess
        line = f"Q {m1} 6 " + " ".join(f"{x:.17g}" for x in H.ravel())
        res = subprocess.run([driver], input=line + "\n", capture_output=True, text=True, check=True).stdout.split()
        out["kn6"] = min(wa, float(res[1]) - float(res[0]))
    if m1 <= sdp_max:
        try:
            import medir_q as S
            out["sdp"] = min(wa, S.sdp_width(H))
        except Exception:
            pass
    return out


def final_bounds(F, strong=True, sdp_max=60):
    """(M3a, M3b): anchos del encierro final de la forma cuadratica F."""
    lo, hi = F.rng()
    wa = hi - lo
    if not strong or not F.Q: return wa, wa
    syms = sorted(set(F.lin) | {i for k in F.Q for i in k})
    idx = {s: t for t, s in enumerate(syms)}
    m = len(syms)
    H = np.zeros((m + 1, m + 1))
    for (i, j), v in F.Q.items():
        a, b = idx[i], idx[j]
        if a == b: H[a, a] += v
        else: H[a, b] += v / 2; H[b, a] += v / 2
    for s, v in F.lin.items():
        H[idx[s], m] += v / 2; H[m, idx[s]] += v / 2
    best = wa
    if m + 1 <= 40:
        L, U, _ = PF.prodk_interval(H); best = min(best, U - L)
    if m + 1 <= sdp_max:
        try:
            import medir_q as S          # ../raak-qplib: Shor con SCS y reparacion (cota valida)
            best = min(best, S.sdp_width(H))
        except Exception:
            pass
    return wa, best
