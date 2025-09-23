#!/usr/bin/env python3

import hoppet as hp
import numpy as np

# This is the PDF subroutine. Notice that the signature is different
# from the one in the C++ code to allow for the python interface, that
# handles pointers in a more complicated way. In the actual example we
# use the internal BenchmarkPDFunpol which is identical

#def hera_lhc(x, Q):
#    N_g = 1.7
#    N_ls = 0.387975
#    N_uv = 5.107200
#    N_dv = 3.064320
#    N_db = N_ls / 2
#
#    uv = N_uv * pow(x, 0.8) * pow((1 - x), 3)
#    dv = N_dv * pow(x, 0.8) * pow((1 - x), 4)
#    dbar = N_db * pow(x, -0.1) * pow(1 - x, 6)
#    ubar = dbar * (1 - x)
#
#    pdf = np.zeros(13) # Initialise the array with zeros
#
#    pdf[ 0+6] = N_g * pow(x, -0.1) * pow(1 - x, 5)
#    pdf[-3+6] = 0.2 * (dbar + ubar)
#    pdf[ 3+6] = pdf[-3 + 6]
#    pdf[ 2+6] = uv + ubar
#    pdf[-2+6] = ubar
#    pdf[ 1+6] = dv + dbar
#    pdf[-1+6] = dbar
#    
## Overkill
#    pdf[ 4+6] = 0.0
#    pdf[ 5+6] = 0.0
#    pdf[ 6+6] = 0.0
#    pdf[-4+6] = 0.0
#    pdf[-5+6] = 0.0
#    pdf[-6+6] = 0.0
#
#    return pdf

def main():
    dy = 0.1    
    nloop = 3
    # Start hoppet
    hp.Start(dy, nloop)
    
    asQ0 = 0.35
    Q0 = np.sqrt(2.0)
    # Do the evolution. 
    hp.Evolve(asQ0, Q0, nloop, 1.0, hp.BenchmarkPDFunpol, Q0)

    # Evaluate the PDFs at some x values and print them
    xvals = [1e-5,1e-4,1e-3,1e-2,0.1,0.3,0.5,0.7,0.9]
    Q = 100.0

    print("")
    print("           Evaluating PDFs at Q =",Q, " GeV")
    print("    x      u-ubar      d-dbar    2(ubr+dbr)    c+cbar       gluon")
    for ix in range(9):
        pdf_array = hp.Eval(xvals[ix], Q)
        print("{:7.1E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E}".format(
            xvals[ix],
            pdf_array[6 + 2] - pdf_array[6 - 2], 
            pdf_array[6 + 1] - pdf_array[6 - 1], 
            2 * (pdf_array[6 - 1] + pdf_array[6 - 2]),
            pdf_array[6 - 4] + pdf_array[6 + 4],
            pdf_array[6 + 0]
        ))

    #hp.WriteLHAPDFGrid("test_python",0)

    hp.DeleteAll()

if __name__ == "__main__":
    main()
