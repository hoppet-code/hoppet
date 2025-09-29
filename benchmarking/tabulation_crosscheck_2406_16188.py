#!/usr/bin/env python3

import hoppet as hp
import numpy as np

def main():
    dy = 0.05    
    nloop = 4
    ymax = 17.0
    Qmin = 1.0
    Qmax = 1000.0
    dlnlnQ = dy/8.0
    order = -6

    # Start hoppet
    # VFNS variant corresponds to table 1 of the HOPPET-v2 update note
    hp.SetApproximateDGLAPN3LO(hp.n3lo_splitting_approximation_up_to_2410_08089)
    hp.StartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, hp.factscheme_MSbar)
    asQ0 = 0.35
    Q0 = np.sqrt(2.0)
    # Do the evolution. 
    hp.Evolve(asQ0, Q0, nloop, 1.0, hp.BenchmarkPDFunpol, Q0)

    # Evaluate the PDFs at some x values and print them
    xvals = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.7, 0.9]
    Q = 100.0

    print("VFNS results with the n3lo_splitting_approximation_up_to_2410_08089 splitting functions")
    print("                                  Evaluating PDFs at Q =",Q, " GeV")
    print("    x      u-ubar      d-dbar    dbar-ubar   2(ubr+dbr)     s+sbar      c+cbar      b+bbar      gluon")
    for ix in range(len(xvals)):
        pdf_array = hp.Eval(xvals[ix], Q)
        print("{:7.1E} {:11.4E} {:11.4E} {:11.4E} {:11.4E}  {:11.4E} {:11.4E} {:11.4E} {:11.4E}".format(
            xvals[ix],
            pdf_array[6 + 2] - pdf_array[6 - 2], 
            pdf_array[6 + 1] - pdf_array[6 - 1], 
            (pdf_array[6 - 1] - pdf_array[6 - 2]),
            2 * (pdf_array[6 - 1] + pdf_array[6 - 2]),
            pdf_array[6 - 3] + pdf_array[6 + 3],
            pdf_array[6 - 4] + pdf_array[6 + 4],
            pdf_array[6 - 5] + pdf_array[6 + 5],
            pdf_array[6 + 0]
        ))
    print("")
    
    # Fixed flavour variant for comparion with 2406.16188
    hp.SetFFN(4)
    hp.SetApproximateDGLAPN3LO(hp.n3lo_splitting_approximation_up_to_2310_05744)
    hp.StartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, hp.factscheme_MSbar)
    hp.Evolve(asQ0, Q0, nloop, 1.0, hp.BenchmarkPDFunpol, Q0)
    
    print("")
    print(" Table to compare with Tables 1 and 2 of 2406.16188")
    print("                                  Evaluating PDFs at Q =",Q, " GeV")
    print("    x      u-ubar      d-dbar    dbar-ubar   2(ubr+dbr)    s-sbar      s+sbar      c+cbar      gluon")
    for ix in range(len(xvals)):
        pdf_array = hp.Eval(xvals[ix], Q)
        print("{:7.1E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E}".format(
            xvals[ix],
            pdf_array[6 + 2] - pdf_array[6 - 2], 
            pdf_array[6 + 1] - pdf_array[6 - 1], 
            (pdf_array[6 - 1] - pdf_array[6 - 2]),
            2 * (pdf_array[6 - 1] + pdf_array[6 - 2]),
            pdf_array[6 + 3] - pdf_array[6 - 3],
            pdf_array[6 - 3] + pdf_array[6 + 3],
            pdf_array[6 - 4] + pdf_array[6 + 4],
            pdf_array[6 + 0]
        ))

    hp.DeleteAll()

if __name__ == "__main__":
    main()
