"""Builds submission/ (flat folder required by the journal: no subfolders, figures Fig1/Fig2, tables inlined)
and submission.zip.   python make_submission.py   (then compile in submission/ with pdflatex+bibtex)"""
import re, os, shutil, zipfile
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
os.makedirs("submission", exist_ok=True)
s = open("exact_raf.tex").read()
s = s.replace("{fig/zonotope.pdf}", "{Fig1.pdf}").replace("{fig/square_cheb.pdf}", "{Fig2.pdf}")
s = re.sub(r"\\input\{(\.\./\.\./results/[^}]+)\}",
           lambda m: "% " + os.path.basename(m.group(1)) + "\n" + open(os.path.normpath(m.group(1))).read().rstrip("\n"), s)
assert "\\input" not in s and "fig/" not in s
open("submission/exact_raf.tex", "w").write(s)
shutil.copy("fig/zonotope.pdf", "submission/Fig1.pdf"); shutil.copy("fig/square_cheb.pdf", "submission/Fig2.pdf")
shutil.copy("bibliography.bib", "submission/bibliography.bib")
