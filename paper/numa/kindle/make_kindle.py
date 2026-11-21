"""Version para Kindle de exact_raf.tex: pagina de 100 x 135 mm, tablas fijas ([H]) junto al texto y
ajustadas al ancho; la tabla del solver se parte en dos bloques. Uso: python make_kindle.py && pdflatex ..."""
import os, re
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "exact_raf.tex")
RES = os.path.join(HERE, "..", "..", "..", "results")

def fit(tex):
    tex = re.sub(r"\\begin\{tabular\}", r"\\begin{adjustbox}{max width=\\linewidth}\\begin{tabular}", tex)
    tex = tex.replace(r"\end{tabular}", r"\end{tabular}\end{adjustbox}")
    return tex.replace(r"\begin{table}[t]", r"\begin{table}[H]")

s = open(SRC).read()
s = s.replace(r"\documentclass[11pt]{article}", r"\documentclass[10pt]{article}")
s = s.replace(r"\usepackage[a4paper,margin=27mm]{geometry}",
              r"\usepackage[paperwidth=100mm,paperheight=135mm,left=4mm,right=4mm,top=5mm,bottom=7mm,footskip=4mm]{geometry}"
              "\n\\usepackage{float}\n\\usepackage{adjustbox}\n\\sloppy\n\\setlength{\\emergencystretch}{3em}")
s = s.replace(r"\begin{algorithm}[t]", r"\begin{algorithm}[H]")
s = fit(s)
s = re.sub(r"\\input\{\.\./\.\./results/([a-z_]+)\.tex\}", lambda m: r"\input{" + m.group(1) + "_k.tex}", s)

# ecuaciones largas: partidas en varias lineas (solo en la version Kindle)
BREAKS = [
 (r"R(e,\ex,\ey) = (e^Tx + x_r\ex)(e^Ty + y_r\ey), \qquad Q(e) = (e^Tx)(e^Ty) = \sum_i x_iy_i\eps_i^2 + \sum_{i<j}(x_iy_j+x_jy_i)\eps_i\eps_j .",
  r"\begin{gathered}R(e,\ex,\ey) = (e^Tx + x_r\ex)(e^Ty + y_r\ey),\\ Q(e) = (e^Tx)(e^Ty) = \sum_i x_iy_i\eps_i^2 + \sum_{i<j}(x_iy_j+x_jy_i)\eps_i\eps_j .\end{gathered}"),
 (r"\hat x = \ell_x(v)+\rho_x,\qquad \ell_x(v)=x_0+\sum_i\frac{x_i}{w_i}(v_i-m_i),\qquad \rho_x\in[-x_r,x_r],",
  r"\begin{gathered}\hat x = \ell_x(v)+\rho_x,\qquad \rho_x\in[-x_r,x_r],\\ \ell_x(v)=x_0+\sum_i\frac{x_i}{w_i}(v_i-m_i),\end{gathered}"),
 (r"N=\ell_x(m)\,\rho_y+\ell_y(m)\,\rho_x+\big(\ell_x(v)-\ell_x(m)+\rho_x\big)\big(\ell_y(v)-\ell_y(m)+\rho_y\big)",
  r"\begin{aligned}N={}&\ell_x(m)\,\rho_y+\ell_y(m)\,\rho_x\\&+\big(\ell_x(v)-\ell_x(m)+\rho_x\big)\big(\ell_y(v)-\ell_y(m)+\rho_y\big)\end{aligned}"),
 (r"N(e,s,t) = x_0y_rt+y_0x_rs + (u + x_rs)(v + y_rt),\qquad (u,v)=(e^Tx,\,e^Ty)\in Z,",
  r"\begin{gathered}N(e,s,t) = x_0y_rt+y_0x_rs + (u + x_rs)(v + y_rt),\\ (u,v)=(e^Tx,\,e^Ty)\in Z,\end{gathered}"),
 (r"p_0=-\sum_k\tilde g_{(k)},\quad p_k = p_{k-1}+2\tilde g_{(k)}\ (k=1,\dots,m),\quad p_{m+k} = p_{m+k-1}-2\tilde g_{(k)}\ (k=1,\dots,m).",
  r"\begin{gathered}p_0=-\sum_k\tilde g_{(k)},\\ p_k = p_{k-1}+2\tilde g_{(k)},\quad p_{m+k} = p_{m+k-1}-2\tilde g_{(k)}\quad (k=1,\dots,m).\end{gathered}"),
 (r"(\|x\|_1+x_r)(\|y\|_1+y_r) - \tfrac12|x|^T|y| = \tfrac12|x|^T|y|+\sum_{i\ne j}|x_iy_j| + x_r\|y\|_1+y_r\|x\|_1+x_ry_r .",
  r"\begin{aligned}&(\|x\|_1+x_r)(\|y\|_1+y_r) - \tfrac12|x|^T|y|\\ &\quad= \tfrac12|x|^T|y|+\sum_{i\ne j}|x_iy_j| + x_r\|y\|_1+y_r\|x\|_1+x_ry_r .\end{aligned}"),
 (r"\max N\ \ge\ R(p)+|L(p)|,\qquad \min N\ \le\ R(p)-|L(p)|\qquad\text{for every } p\in\U^{n+2}.",
  r"\begin{gathered}\max N\ \ge\ R(p)+|L(p)|,\qquad \min N\ \le\ R(p)-|L(p)|\\ \text{for every } p\in\U^{n+2}.\end{gathered}"),
 (r"x_0y_0 + e^T(y_0x+x_0y) + k_1\eps_\pm + k_2\eps_+ + k_3\eps_-,\quad k_2=\sum_{x_iy_i>0}x_iy_i,\ k_3=\sum_{x_iy_i<0}|x_iy_i|,",
  r"\begin{gathered}x_0y_0 + e^T(y_0x+x_0y) + k_1\eps_\pm + k_2\eps_+ + k_3\eps_-,\\ k_2=\sum_{x_iy_i>0}x_iy_i,\quad k_3=\sum_{x_iy_i<0}|x_iy_i|,\end{gathered}"),
]
for a, b in BREAKS:
    assert a in s, a[:40]
    s = s.replace(a, b)
s = s.replace(r"\bibliography{bibliography}", r"\bibliography{../bibliography}")
open(os.path.join(HERE, "exact_raf_kindle.tex"), "w").write(s)

for name in ["table_rules_random", "table_rules_expressions", "table_cost"]:
    open(os.path.join(HERE, name + "_k.tex"), "w").write(fit(open(os.path.join(RES, name + ".tex")).read()))

# tabla del solver en dos bloques
rows = []
for line in open(os.path.join(RES, "table_solver.tex")):
    if line.count("&") == 12:
        c = [x.strip() for x in line.replace(r"\\", "").split("&")]
        if c[0] != "Set": rows.append(c)
t1 = [r"\begin{tabular}{lrrrrr}", r"\toprule", r" & & & \multicolumn{3}{c}{Solved} \\", r"\cmidrule(lr){4-6}",
      r"Set & Inst. & $t_{\max}$ & S & SH & E \\", r"\midrule"] + \
     [" & ".join(c[0:6]) + r" \\" for c in rows] + [r"\bottomrule", r"\end{tabular}"]
t2 = [r"\begin{tabular}{lrrrrrrr}", r"\toprule",
      r" & \multicolumn{3}{c}{Nodes / S} & \multicolumn{2}{c}{CPU / S} & \multicolumn{2}{c}{E vs SH} \\",
      r"\cmidrule(lr){2-4}\cmidrule(lr){5-6}\cmidrule(lr){7-8}",
      r"Set & SH & E & SH$^\dagger$ & SH & E & nodes & $<$/$>$ \\", r"\midrule"] + \
     [" & ".join([c[0]] + c[6:13]) + r" \\" for c in rows] + [r"\bottomrule", r"\end{tabular}"]
open(os.path.join(HERE, "table_solver_k.tex"), "w").write(fit("\n".join(t1)) + "\n\n\\medskip\n" + fit("\n".join(t2)) + "\n")
print("ok")
