---
project: HOPPET
author: Frederic Dreyer, Alexander Karlberg, Paolon Nason, Juan Rojo, Gavin Salam and Giulia Zanderighi
preprocess: true
src_dir: ../../src
docmark: <
predocmark: !
exclude_dir: ../../src/inc
---

HOPPET is a code for evolution and manipulation of parton distribution
functions. This documentation is for the Fortran interface.

Key types:
  * [[grid_def]]
  * a single flavour is held in a double-precision pointer array,
    typically allocated with [[AllogGridQuant]]
  * [[grid_conv]]: 
