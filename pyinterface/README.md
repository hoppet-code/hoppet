# HOPPET: A Higher Order Perturbative Parton Evolution Toolkit

HOPPET is a package for carrying out DGLAP evolution and other common
manipulations of parton distribution functions (PDFs) in particle
physics. The core package is written in Fortran. This Python interface
provides access to a key subset of the functionality, the part known 
in the documentation as the "streamlined" interface.

The full underlying code and documentation can be found at
https://github.com/hoppet-code/hoppet.

## Example usage

```Python
import hoppet as hp
import numpy as np

def main():
    # choose the underlying spacing in y=log(1/x): 
    # smaller spacing gives higher accuracy but slower evolution
    dy = 0.1  
    # choose NNLO evolution
    nloop = 3 
    # Start hoppet
    hp.Start(dy, nloop)
    
    asQ0 = 0.35
    Q0 = np.sqrt(2.0)
    # Do the evolution starting from a standard benchmark initial condition 
    hp.Evolve(asQ0, Q0, nloop, 1.0, hp.BenchmarkPDFunpol, Q0)

    # Evaluate the PDFs at some x values and print them
    xvals = [1e-5,1e-4,1e-3,1e-2,0.1,0.3,0.5,0.7,0.9]
    Q = 100.0

    print('')
    print('           Evaluating PDFs at Q =',Q, ' GeV')
    print('    x      u-ubar      d-dbar    2(ubr+dbr)    c+cbar       gluon')
    for ix in range(9):
        pdf_array = hp.Eval(xvals[ix], Q)
        print('{:7.1E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E}'.format(
            xvals[ix],
            pdf_array[6 + 2] - pdf_array[6 - 2], 
            pdf_array[6 + 1] - pdf_array[6 - 1], 
            2 * (pdf_array[6 - 1] + pdf_array[6 - 2]),
            pdf_array[6 - 4] + pdf_array[6 + 4],
            pdf_array[6 + 0]
        ))

    hp.DeleteAll()
```   
For more examples take a look at
[example_py](https://github.com/hoppet-code/hoppet/tree/master/example_py).
The above example is essentially identical to
[tabulation_example.py](https://github.com/hoppet-code/hoppet/tree/master/example_py/tabulation_example.py)
and prints the output of a typical benchmark PDF.


## Citation policy

If you use this program in a scientific publication we ask that you cite

G.P. Salam, J. Rojo, 'A Higher Order Perturbative Parton Evolution Toolkit (HOPPET)', 
Comput. Phys. Commun. 180 (2009) 120-156, [arXiv:0804.3755](https://arxiv.org/abs/0804.3755)

and                                                       

A. Karlberg, P. Nason, G.P. Salam, G. Zanderighi & F. Dreyer [arXiv:2509.XXXXX](https://arxiv.org/abs/2509.XXXXX). 

## Version numbering

A version number `X.Y.Z.a` means that the interface uses hoppet version
`X.Y.Z`. The `.a` part of the version number gets incremented when the
Python pip package and/or wheels are updated, but the underlying hoppet
code stays the same. 