"""SVG de los diagramas de la presentacion (hexagono Z, Z + delta y curvas de nivel de uv)."""
import numpy as np
BLUE, ORANGE, GREEN, INK = "#2F6EA5", "#A0561A", "#3C8D5A", "#1D2B3A"
G_SORTED = [(2, -1), (1, 0), (2, 1)]
COLS = [BLUE, ORANGE, GREEN]

def verts(delta=(0, 0)):
    g = np.array(G_SORTED, float); V = [np.array(delta, float) - g.sum(0)]
    for k in range(3): V.append(V[-1] + 2 * g[k])
    for k in range(3): V.append(V[-1] - 2 * g[k])
    return V

class Map:
    def __init__(self, X0, Y0, s=64, W=800, H=420): self.X0, self.Y0, self.s, self.W, self.H = X0, Y0, s, W, H
    def __call__(self, u, v): return (self.X0 + self.s * u, self.Y0 - self.s * v)

def axes(M):
    x0, y0 = M(0, 0)
    return (f'<line x1="0" y1="{y0:.1f}" x2="{M.W}" y2="{y0:.1f}" stroke="#B9B3A6" stroke-width="2"/>'
            f'<line x1="{x0:.1f}" y1="0" x2="{x0:.1f}" y2="{M.H}" stroke="#B9B3A6" stroke-width="2"/>')

def poly(M, V, fill="#E2DCCF", stroke=INK, sw=3):
    pts = " ".join(f"{M(*p)[0]:.1f},{M(*p)[1]:.1f}" for p in V[:-1])
    return f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>'

def seg(M, a, b, color, sw=8, dash=None):
    (x1, y1), (x2, y2) = M(*a), M(*b)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" stroke-width="{sw}" stroke-linecap="round"{d}/>'

def dot(M, p, r=8, color=INK, ring=False):
    x, y = M(*p); extra = ' stroke="#F7F5F0" stroke-width="4"' if ring else ""
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{color}"{extra}/>'

def levels(M, cs=(-8, -6, -4, -2, 2, 4, 6, 8)):
    """Curvas uv = c (coordenadas del dibujo), recortadas al rectangulo."""
    out = []
    umin, umax = (0 - M.X0) / M.s, (M.W - M.X0) / M.s
    vmin, vmax = (M.Y0 - M.H) / M.s, M.Y0 / M.s
    for c in cs:
        for sgn in (1, -1):
            us = np.linspace(0.05, max(abs(umin), abs(umax)), 600) * sgn
            pts = [(u, c / u) for u in us if umin <= u <= umax and vmin - 1 <= c / u <= vmax + 1]
            if len(pts) < 2: continue
            d = " ".join(("M" if i == 0 else "L") + f"{M(*p)[0]:.1f} {M(*p)[1]:.1f}" for i, p in enumerate(pts))
            dash = ' stroke-dasharray="8 6"' if c < 0 else ""
            out.append(f'<path d="{d}" fill="none" stroke="#9AA3AD" stroke-width="2"{dash} clip-path="url(#clip)"/>')
    return "".join(out)

def svg(M, body, label):
    return (f'<svg aria-label="{label}" width="{M.W}" height="{M.H}" viewBox="0 0 {M.W} {M.H}">'
            f'<defs><clipPath id="clip"><rect x="0" y="0" width="{M.W}" height="{M.H}"/></clipPath>'
            f'<marker id="m0" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto"><path d="M0 0 L10 5 L0 10 z" fill="{BLUE}"/></marker>'
            f'<marker id="m1" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto"><path d="M0 0 L10 5 L0 10 z" fill="{ORANGE}"/></marker>'
            f'<marker id="m2" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto"><path d="M0 0 L10 5 L0 10 z" fill="{GREEN}"/></marker>'
            f'</defs>{body}</svg>')

def arrow(M, g, k):
    (x1, y1), (x2, y2) = M(0, 0), M(*g)
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{COLS[k]}" stroke-width="5" marker-end="url(#m{k})"/>'

def lab(M, u, v, text, cls="", dx=0, dy=0, w=160):
    x, y = M(u, v)
    c = f' class="lab {cls}"' if cls else ' class="lab"'
    return f'<p{c} style="left:{x + dx:.0f}px; top:{y + dy:.0f}px; width:{w}px;">{text}</p>'
