"""Reduccion avida por bloques producto (PRODK), version incremental para n grande.

Misma idea y mismos candidatos que ../raak-pares/medir_prodk.py, pero:
  - la ganancia de un candidato P sobre el soporte S es local: TT(R_S) - TT(R_S - P_S) - SH(P), porque P solo
    toca las entradas de S y TT es separable por entradas;
  - cada soporte guarda su mejor candidato y, tras extraer un bloque en S, solo se recalculan los soportes que
    intersectan S;
  - soportes: pares con entrada no nula (o ambas diagonales no nulas) y, si m <= max_triples_m, triples con al
    menos 2 aristas; el soporte completo si m <= full_max.

ancho(PRODK) = sum_p SH(u_p, v_p) + TT(R_final), valido para cualquier eleccion de bloques.
"""
import heapq, itertools, math
import numpy as np

TS = (1.0, 0.75, 0.5, 0.25)


def tt(R):
    d = np.abs(np.diag(R)).sum()
    return d + 2 * (np.abs(R).sum() - d)


def zono_vertices(u, v):
    g = [(a, b) if (b > 0 or (b == 0 and a > 0)) else (-a, -b) for a, b in zip(u, v) if a or b]
    if not g: return [np.zeros(2)]
    g.sort(key=lambda t: math.atan2(t[1], t[0]))
    p = -np.array([sum(t[0] for t in g), sum(t[1] for t in g)])
    V = [p.copy()]
    for t in g: p = p + 2 * np.array(t); V.append(p.copy())
    for t in g: p = p - 2 * np.array(t); V.append(p.copy())
    return V


def sh_width(u, v):
    lo, hi = sh_range(u, v)
    return hi - lo


def sh_range(u, v):
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


def local_candidates(RS):
    """Candidatos sobre la submatriz RS (coordenadas locales): listas de pares (u, v)."""
    out = []
    lam, V = np.linalg.eigh(RS)
    sc = max(1e-300, np.abs(lam).max()) * 1e-12
    ip, ineg = int(np.argmax(lam)), int(np.argmin(lam))
    if lam[ip] > sc and lam[ineg] < -sc:
        a, b = np.sqrt(lam[ip]) * V[:, ip], np.sqrt(-lam[ineg]) * V[:, ineg]
        out.append([(a + b, a - b)])
    k = int(np.argmax(np.abs(lam)))
    if abs(lam[k]) > sc:
        a = np.sqrt(abs(lam[k])) * V[:, k]
        out.append([(a, np.sign(lam[k]) * a)])
    s = RS.shape[0]
    for t in range(s):
        v = 2 * RS[t].copy(); v[t] = RS[t, t]
        if np.abs(v).sum() > 0:
            u = np.zeros(s); u[t] = 1.0
            out.append([(u, v)])
    return out


def sym(pairs, s):
    P = np.zeros((s, s))
    for u, v in pairs: P += (np.outer(u, v) + np.outer(v, u)) / 2
    return P


def best_on(R, S):
    """Mejor candidato en el soporte S: (ganancia, pares en coordenadas locales)."""
    RS = R[np.ix_(S, S)]
    base = tt(RS)
    if base == 0: return 0.0, None
    best, gain = None, 0.0
    for pairs in local_candidates(RS):
        for t in TS:
            ps = [(t * u, v) for u, v in pairs]
            g = base - tt(RS - sym(ps, len(S))) - sum(sh_width(u, v) for u, v in ps)
            if g > gain: gain, best = g, ps
    return gain, best


def supports(Q, max_triples_m=40, full_max=12):
    m = Q.shape[0]
    nz = np.abs(Q) > 0
    S = [(i, j) for i, j in itertools.combinations(range(m), 2) if nz[i, j] or (nz[i, i] and nz[j, j])]
    if m <= max_triples_m:
        S += [(i, j, k) for i, j, k in itertools.combinations(range(m), 3) if nz[i, j] + nz[i, k] + nz[j, k] >= 2]
    if m <= full_max: S.append(tuple(range(m)))
    return S


def tt_interval(R):
    d = np.diag(R)
    off = np.abs(R).sum() - np.abs(d).sum()
    return np.minimum(d, 0).sum() - off, np.maximum(d, 0).sum() + off


def prodk(Q, max_iter=100000, **kw):
    """Devuelve (ancho, numero de bloques). prodk_interval devuelve ademas el intervalo [L, U]."""
    L, U, nb = prodk_interval(Q, max_iter, **kw)
    return U - L, nb


def prodk_interval(Q, max_iter=100000, **kw):
    m = Q.shape[0]
    R = Q.copy()
    sup = supports(Q, **kw)
    by_var = [[] for _ in range(m)]
    for s, S in enumerate(sup):
        for i in S: by_var[i].append(s)
    ver = [0] * len(sup)
    heap = []
    cache = {}
    def push(s):
        g, ps = best_on(R, list(sup[s]))
        cache[s] = ps
        if ps is not None: heapq.heappush(heap, (-g, s, ver[s]))
    for s in range(len(sup)): push(s)
    total_sh, nblocks, Lsh, Ush = 0.0, 0, 0.0, 0.0
    scale = max(1.0, tt(Q)) * 1e-12
    for _ in range(max_iter):
        while heap and heap[0][2] != ver[heap[0][1]]: heapq.heappop(heap)
        if not heap or -heap[0][0] <= scale: break
        _, s, _ = heapq.heappop(heap)
        S = list(sup[s]); ps = cache[s]
        R[np.ix_(S, S)] -= sym(ps, len(S))
        for u, v in ps:
            lo_, hi_ = sh_range(u, v); Lsh += lo_; Ush += hi_
        nblocks += len(ps)
        touched = set(t for i in S for t in by_var[i])
        for t in touched:
            ver[t] += 1; push(t)
    Lr, Ur = tt_interval(R)
    return Lsh + Lr, Ush + Ur, nblocks
