"""Diff para Kindle (compilar con lualatex): texto nuevo con fondo gris (visible en blanco y negro),
eliminado en rojo tachado. Usa lua-ul, que resalta texto con formulas y citas y corta lineas."""
import sys
s = open(sys.argv[1]).read()
s = s.replace(r"\providecommand{\DIFadd}[1]{{\protect\color{blue} \sf #1}}", r"\providecommand{\DIFadd}[1]{\highLight[hlgray]{#1}}")
s = s.replace(r"\providecommand{\DIFdel}[1]{{\protect\color{red} \scriptsize #1}}", r"\providecommand{\DIFdel}[1]{{\protect\color{red}\strikeThrough{#1}}}")
s = s.replace(r"\begin{document}", "\\usepackage{luacolor}\n\\usepackage{lua-ul}\n\\definecolor{hlgray}{gray}{0.78}\n\\begin{document}", 1)
s = s.replace(r"\usepackage{mathptmx}", "\\usepackage[T1]{fontenc}\n\\usepackage{mathptmx}", 1)
open(sys.argv[2], "w").write(s)
