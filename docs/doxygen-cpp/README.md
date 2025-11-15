# Documentation for the HOPPET C++ interface

## Introduction
Hoppet has two interfaces: 

* the streamlined interface: a set of C-style functions that wrap the
  basic functionality of hoppet (evolution, tabulation, access) of
  parton distribution functions. All names are came-case, such as
  hoppetStart(...)

* the general or object-oriented interface: a set of C++ classes that
  encapsulate the main objects of hoppet (grids, pdfs, splitting
  functions, etc). These allow finer-grained manipulation of PDFs,
  splitting functions, etc. All objects are in the `hoppet` namespace.


## General (object-oriented) interface

The C++ interface provides a series of objects that wrap the underlying
hoppet's underlying Fortran types. The main objects are listed below, 
grouped by category.

For each class `hoppet::classname` (which owns underlying data), there is
a corresponding C++ wrapper class `hoppet::classname_view`, which is a
lightweight non-owning view of the underlying data. The view classes
are useful to avoid unnecessary copies when returning objects. 

### Grid representation

| class | brief description |
|-------|-------------------|
| hoppet::grid_def  | An object that defines a grid in \f$y=\ln(1/x)\f$ for pdfs     |

The hoppet::grid_def object is used to define the grid on which
parton distribution functions (pdfs) and splitting functions are
represented. Many objects will take "views" of the grid_def object
and the grid_def object must therefore be kept alive as long as
those objects are in use.

### Parton distribution function (pdf) objects

| class | brief description |
|-------|-------------------|
| hoppet::pdf_flav  | a single flavour of parton distribution function (pdf) defined on a grid in \f$y\f$ (alias for hoppet::grid_quant) |
| hoppet::pdf       | An object that holds all flavours of parton distribution functions (pdfs) defined on a grid  (alias for hoppet::grid_quant_2d)|
| hoppet::pdf_table | An object that holds a tabulation of parton distribution functions (pdfs) in \f$y\f$ and \f$Q\f$ |
| hoppet::pdf_table_array | An array of hoppet::pdf_table objects |

### Splitting function and related objects

| class | brief description |
|-------|-------------------|
| hoppet::split_fn  | a grid representation of a splitting function object (alias for hoppet::grid_conv)|
| hoppet::split_mat | a (sparse) matrix of QCD splitting functions  |
| hoppet::mass_threshold_mat | a (sparse) matrix of mass threshold terms |
| hoppet::dglap_holder | an object that holds all the splitting functions and mass threshold terms needed for QCD DGLAP evolution |

### Other classes

| class | brief description |
|-------|-------------------|
| hoppet::running_coupling | an object that represents the running coupling \f$\alpha_s(Q)\f$ |

## Object-oriented access to streamlined-interface objects

See the hoppet::sl namespace for C++ classes that wrap the streamlined interface objects.