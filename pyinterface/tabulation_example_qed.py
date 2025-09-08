#!/usr/bin/env python3

import hoppet as hp
import numpy as np

# This is the PDF subroutine. Notice that the signature is different
# from the one in the C++ code to allow for the python interface, that
# handles pointers in a more complicated way
def hera_lhc(x, Q):
    N_g = 1.7
    N_ls = 0.387975
    N_uv = 5.107200
    N_dv = 3.064320
    N_db = N_ls / 2

    gluon = N_g * pow(x, -0.1) * pow(1 - x, 5)
    uv = N_uv * pow(x, 0.8) * pow((1 - x), 3)
    dv = N_dv * pow(x, 0.8) * pow((1 - x), 4)
    dbar = N_db * pow(x, -0.1) * pow(1 - x, 6)
    ubar = dbar * (1 - x)

    pdf = np.zeros(18) # Initialise the array with zeros

    pdf[ 0+6] = 0.99 * gluon
    pdf[-3+6] = 0.2 * (dbar + ubar)
    pdf[ 3+6] = pdf[-3 + 6]
    pdf[ 2+6] = uv + 0.9 * ubar
    pdf[-2+6] = 0.9 * ubar
    pdf[ 1+6] = dv + 0.9 * dbar
    pdf[-1+6] = 0.9 * dbar

    # qed part
    pdf[ 8+6] = gluon / 100.0 + ubar / 10.0  # photon 
    pdf[ 9+6] = dbar / 10.0  # electron
    pdf[ 10+6] = ubar / 10.0 # muon
    pdf[ 11+6] = dbar / 10.0 # tau
    
# Overkill
    pdf[ 4+6] = 0.0
    pdf[ 5+6] = 0.0
    pdf[ 6+6] = 0.0
    pdf[-4+6] = 0.0
    pdf[-5+6] = 0.0
    pdf[-6+6] = 0.0

    return pdf

def main():
    dy = 0.1    
    nloop = 3
    
    # set QED
    use_qcd_qed  = True
    use_Plq_nnlo = False
    hp.SetQED(True, use_qcd_qed, use_Plq_nnlo)
    
    # Start hoppet
    hp.StartExtended(12.0, dy, 1.0, 28000.0, dy/4.0, nloop, -6, hp.factscheme_MSbar)

    # Set heavy flavour scheme
    mc = 1.414213563
    mb = 4.5
    mt = 175.0
    hp.SetVFN(mc, mb, mt)

    
    asQ0 = 0.35
    Q0 = np.sqrt(2.0)
    # Do the evolution. 
    hp.Evolve(asQ0, Q0, nloop, 1.0, hera_lhc, Q0)

    # Evaluate the PDFs at some x values and print them
    xvals = [1e-5,1e-4,1e-3,1e-2,0.1,0.3,0.5,0.7,0.9]
    Q = 100.0
    
    print("")
    print("           Evaluating PDFs at Q =",Q, " GeV")
    print("    x      u-ubar      d-dbar    2(ubr+dbr)    c+cbar       gluon      photon     e+ + e-    mu+ + mu-   tau+ + tau-")
    for ix in range(9):
        pdf_array = hp.Eval(xvals[ix], Q)
        print("{:7.1E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E}".format(
            xvals[ix],
            pdf_array[6 + 2] - pdf_array[6 - 2], 
            pdf_array[6 + 1] - pdf_array[6 - 1], 
            2 * (pdf_array[6 - 1] + pdf_array[6 - 2]),
            pdf_array[6 - 4] + pdf_array[6 + 4],
            pdf_array[6 + 0], pdf_array[6 + 8], pdf_array[6 + 9], pdf_array[6 + 10], pdf_array[6 + 11]
        ))

    #hp.WriteLHAPDFGrid("test_python",0)

    hp.DeleteAll()

if __name__ == "__main__":
    main()
