HOPPET documentation
=====================

The latest official documentation for HOPPET can be found 
at https://hoppet-code.github.io/hoppet/, which links to
the manual and Python interface documentation.

This directory contains the source for the documentation,
which can be built locally if desired. There are two parts:

- [manual/](manual/): The HOPPET manual in LaTeX format, including the v2 update
  note. Run `latexmk -pdf HOPPET-doc.tex` in that directory to generate
  the PDF from the source. You will need a good LaTeX installation (e.g.
  texlive) for this to work.

- [python-auto/](python-auto/): skeleton for automatic generation of documentation for
  the Python interface. Run `pip install -r requirements.txt; make html`
  in that directory to generate the HTML documentation from the source.
  You will need a working Python environment, with pip able to install
  the required packages. 
