# Anonymous LaTeX Source Package

This directory contains the anonymous manuscript source for CJL/RCL review.

Compile with:

```sh
pdflatex main.tex
biber main
pdflatex main.tex
pdflatex main.tex
```

The included `main.bbl` is provided so the file can also be checked with `pdflatex main.tex` alone. The bibliography file `refs.bib` contains only the entries cited by the anonymous manuscript. The identified title page/accompanying document is separate: `../CJL_TITLE_PAGE.docx`.
