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

- [doxygen-fortran/](doxygen-fortran/): Doxygen configuration file for generating
  documentation from the Fortran source code. From the top-level directory,
  run `doxygen docs/doxygen-fortran/Doxyfile` to generate the documentation, which 
  will be placed in `docs/doxygen-fortran-output/html/`. Note that
  much of the documentation will remain incomplete until hoppet has replaced
  "!!" comments with proper "!>" Doxygen comments.