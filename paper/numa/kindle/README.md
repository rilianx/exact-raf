# Versión Kindle

    python make_kindle.py v4 AAAA-MM-DD           # exact_raf_kindle.tex (100 x 135 mm, tablas fijas junto al texto)
    pdflatex exact_raf_kindle && bibtex exact_raf_kindle && pdflatex exact_raf_kindle && pdflatex exact_raf_kindle
    latexdiff --type=CFONT --math-markup=coarse kindle_vN.tex exact_raf_kindle.tex > diff_raw.tex
    python mkdiff.py diff_raw.tex exact_raf_kindle_vM_cambios.tex    # nuevo con fondo gris, borrado en rojo tachado
    lualatex ... (tres pasadas y bibtex)

`kindle_vN.tex` guarda la fuente de la última versión enviada, para hacer el diff solo contra ella.
