import os
from svg import *
OUT = "../deck/project/slides"
def write(name, html): open(os.path.join(OUT, name + ".html"), "w").write(html)
DIAG = '<div class="diagram" style="width:800px; height:420px;">'

# ---------- idea-borde: Z con los generadores originales ----------
M = Map(400, 210)
V = verts()
body = axes(M) + poly(M, V) + "".join(arrow(M, g, k) for k, g in enumerate([(2, -1), (1, 0), (2, 1)]))
labs = (lab(M, 2, -1, "g₃ = (2, −1)", "b", 8, 0, 200) + lab(M, 1, 0, "g₁ = (1, 0)", "o", -10, -40, 200) +
        lab(M, 2, 1, "g₂ = (2, 1)", "g", 8, -34, 200) + lab(M, 6, 0, "u", "", -30, 8, 40) + lab(M, 0, 3, "v", "", 10, 0, 40))
write("idea-borde", f'''<section id="idea-borde">
<p class="eyebrow">Idea 2 · Lema 2</p>
<h2>(u, v) recorre un zonotopo y uv tiene sus extremos en el borde</h2>
<div class="cols" style="align-items:center;">
<div class="col">
<p>El par (u, v) = (eᵀx, eᵀy) recorre Z = {{ Σ εᵢ gᵢ }}, un zonotopo plano con generadores gᵢ = (xᵢ, yᵢ).</p>
<p>El producto (u + δ₁)(v + δ₂) tiene un solo punto crítico, una silla: no hay máximo ni mínimo en el interior de Z.</p>
<p class="lead">Basta recorrer el borde de Z.</p>
<p class="note">En el ejemplo hay tres generadores, así que Z es un hexágono.</p>
</div>
{DIAG}
{svg(M, body, "Zonotopo Z del ejemplo: hexágono con vértices (−5,0), (−1,−2), (1,−2), (5,0), (1,2) y (−1,2), y sus tres generadores dibujados desde el origen")}
{labs}
</div>
</div>
<aside>La matriz hessiana de uv es [[0,1],[1,0]], indefinida: el único punto crítico es una silla. Por eso el máximo y el mínimo sobre un conjunto compacto como Z están en su borde. Con tres generadores, el borde tiene seis aristas, dos por generador.</aside>
</section>
''')

# ---------- paso-orden: borde con vertices ----------
body = axes(M) + poly(M, V, fill="#F7F5F0", sw=0)
for k in range(3):
    body += seg(M, V[k], V[k + 1], COLS[k]) + seg(M, V[3 + k], V[4 + k], COLS[k], dash="16 12")
body += "".join(dot(M, p) for p in V[:-1])
offs = [(-48, -14), (-30, 12), (6, 12), (14, -14), (6, -40), (-40, -40)]
labs = "".join(lab(M, *V[k], f"p{'₀₁₂₃₄₅'[k]}", "", *offs[k], 60) for k in range(6))
labs += lab(M, -3, -1, "+2g̃₍₁₎", "b", -120, 0, 140) + lab(M, 0, -2, "+2g̃₍₂₎", "o", -50, 16, 140) + lab(M, 3, -1, "+2g̃₍₃₎", "g", 20, 0, 140)
labs += lab(M, 3, 1, "−2g̃₍₁₎", "b", 20, -30, 140) + lab(M, 0, 2, "−2g̃₍₂₎", "o", -50, -46, 140) + lab(M, -3, 1, "−2g̃₍₃₎", "g", -130, -30, 140)
write("paso-orden", f'''<section id="paso-orden" class="alt">
<p class="eyebrow">Líneas 1–2 del algoritmo</p>
<h2>Ordenar los generadores por ángulo y armar el borde</h2>
<div class="cols" style="align-items:center;">
<div class="col">
<p class="step">Línea 1</p>
<p>Los tres ya están en el semiplano u &gt; 0. Ordenados por ángulo: g̃₍₁₎ = (2, −1), g̃₍₂₎ = (1, 0), g̃₍₃₎ = (2, 1).</p>
<p class="step">Línea 2</p>
<p>σ = (5, 0). El borde parte en p₀ = −σ, suma 2g̃₍₁₎, 2g̃₍₂₎, 2g̃₍₃₎ y después los resta en el mismo orden:</p>
<p class="formula" style="font-size:28px;">(−5, 0) → (−1, −2) → (1, −2) → (5, 0) → (1, 2) → (−1, 2) → (−5, 0)</p>
</div>
{DIAG}
{svg(M, body, "Borde de Z: seis aristas; las tres primeras suman dos veces cada generador ordenado y las tres últimas los restan")}
{labs}
</div>
</div>
<aside>El orden por ángulo es lo que garantiza que el camino recorra el borde y no se meta al interior. Esto se hace una sola vez: los cuatro recorridos usan el mismo polígono, trasladado por δ.</aside>
</section>
''')

# ---------- paso-esquina: Z + delta con curvas de nivel ----------
M2 = Map(336, 274)
V2 = verts((1, 1))
body = axes(M2) + levels(M2) + poly(M2, V2, fill="#E2DCCF") + "".join(dot(M2, p) for p in V2[:-1])
uvs = [float(p[0] * p[1]) + 0.0 for p in V2[:-1]]
uvs = [0.0 if abs(x) < 1e-12 else x for x in uvs]
offs = [(-70, -10), (-40, 14), (8, 14), (14, -14), (6, -42), (-70, -42)]
labs = "".join(lab(M2, *V2[k], f"{uvs[k] + 0.0:g}".replace("-", "−"), "", *offs[k], 70) for k in range(6))
labs += lab(M2, 7, 0, "u", "", -30, 6, 40) + lab(M2, 0, 4, "v", "", 10, 0, 40)
write("paso-esquina", f'''<section id="paso-esquina">
<p class="eyebrow">Líneas 3–6 · una esquina</p>
<h2>Esquina (s, t) = (1, 1): trasladar el borde y mirar uv</h2>
<div class="cols" style="align-items:center;">
<div class="col">
<p class="formula" style="font-size:34px;">δ = (xᵣs, yᵣt) = (1, 1)<br>κ = x₀yᵣt + y₀xᵣs = 2<br>p = δ − σ = (−4, 1)</p>
<p>El hexágono Z + δ sobre las curvas de nivel de uv (hipérbolas; punteadas, las negativas). El número en cada vértice es uv.</p>
<p class="note">En los vértices, uv va de −4 a 6. ¿Hay valores mayores o menores dentro de una arista?</p>
</div>
{DIAG}
{svg(M2, body, "Hexágono Z trasladado por (1,1), con vértices (−4,1), (0,−1), (2,−1), (6,1), (2,3) y (0,3), sobre las curvas de nivel de uv")}
{labs}
</div>
</div>
<aside>Las curvas de nivel de uv son hipérbolas con centro en el origen, que es el punto silla. Una arista puede cortar varias de ellas: si es tangente a una dentro del tramo, ahí hay un extremo de uv sobre la arista que no está en sus vértices.</aside>
</section>
''')

# ---------- paso-arista: arista 4 (max) y arista 6 (min) ----------
body = axes(M2) + levels(M2) + poly(M2, V2, fill="#F7F5F0")
body += seg(M2, V2[3], V2[4], "#C0392B", sw=9) + seg(M2, V2[5], V2[6], BLUE, sw=9)
body += dot(M2, V2[3]) + dot(M2, V2[4]) + dot(M2, (4, 2), r=13, color="#C0392B", ring=True) + dot(M2, (-3, 1.5), r=13, color=BLUE, ring=True)
labs = (lab(M2, 6, 1, "p", "", 12, -10, 40) + lab(M2, 2, 3, "q", "", -10, -44, 40) +
        lab(M2, 4, 2, "máx: uv = 8", "o", 16, -46, 220) + lab(M2, -3, 1.5, "mín: uv = −9/2", "b", -20, 22, 240) +
        lab(M2, 7, 0, "u", "", -30, 6, 40) + lab(M2, 0, 4, "v", "", 10, 0, 40))
write("paso-arista", f'''<section id="paso-arista" class="alt">
<p class="eyebrow">Líneas 7–14 · una arista</p>
<h2>En cada arista, uv es una parábola en τ</h2>
<div class="cols" style="align-items:center;">
<div class="col">
<p>Cuarta arista: p = (6, 1), d = −2g̃₍₁₎ = (−4, 2), q = (2, 3). Sus puntos son p + τd, τ ∈ [0, 1]:</p>
<p class="formula" style="font-size:32px;">uv = P + Bτ + Aτ² = 6 + 8τ − 8τ²<br><span style="color:#5A6672;">P = pᵤpᵥ = 6,&nbsp; B = pᵤdᵥ + pᵥdᵤ = 8,&nbsp; A = dᵤdᵥ = −8</span></p>
<p>τ* = −B/(2A) = ½ ∈ (0, 1): en (4, 2), uv = P − B²/(4A) = 8, aunque en ambos extremos vale 6.</p>
<p class="note">La sexta arista da el mínimo, −9/2, con τ* = ¾ en (−3, 3/2).</p>
</div>
{DIAG}
{svg(M2, body, "Hexágono trasladado con la cuarta arista resaltada y su punto estacionario (4,2), donde uv vale 8, y la sexta arista con su punto estacionario (−3,1.5), donde uv vale −9/2")}
{labs}
</div>
</div>
<aside>La arista es tangente a la hipérbola uv = 8 en (4, 2): ahí está el máximo, en el medio de la arista. Si solo se miraran los vértices, se tendría máximo 6 y mínimo −4, y el intervalo sería incorrecto.</aside>
</section>
''')
print("ok")
