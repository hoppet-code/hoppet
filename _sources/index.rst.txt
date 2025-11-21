.. HOPPET documentation master file, created by
   sphinx-quickstart on Mon Nov 10 10:28:55 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HOPPET documentation
====================

This is the documentation for HOPPET's Python interface.

HOPPET is a Fortran package for carrying out DGLAP evolution and other
common manipulations of parton distribution functions (PDFs). It also
includes C++ and Python interfaces to many of the key features.

Within HOPPET, PDFs are represented on a grid in x-space so as to
avoid limitations on the functional form of input
distributions. Good speed and accuracy are obtained through the
representation of splitting functions in terms of their convolution
with a set of piecewise polynomial basis functions, and Runge-Kutta
techniques are used for the evolution in Q.

Key features include:

* unpolarized evolution in the MSbar scheme to N3LO including heavy-quark thresholds
* QED evolution
* massless structure function evaluation up to N3LO
* polarized evolution up to NLO
* simple access (in Fortran 90) to the objects representing splitting
  function and PDFs, making it possible for a user to extend the
  facilities already provided.


For more complete explanations of how HOPPET works, as well as its
underlying Fortran interface, see the full manual in :doc:`PDF format <_generated/github_links>`
for details.


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   api
   _generated/hoppet_constants
   _generated/github_links
   examples.md
