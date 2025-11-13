#! /usr/bin/env python3
"""
This is a small example that uses the Hoppet and LHAPDF python
interfaces. It loads an LHAPDF set and assigns it to Hoppet.
"""

import hoppet as hp
from hoppet import lhapdf
import lhapdf
import numpy as np
import argparse
import time

def main():
    # Get commandline
    parser = argparse.ArgumentParser(description="Check of an LHAPDF grid against HOPPET evolution.")
    parser.add_argument('-pdf', type=str, default="PDF4LHC21_40", help='LHAPDF set name (default PDF4LHC21_40)')
    parser.add_argument('-yorder', type=int, default=2, help='order for interpolation in y=ln1/x (default=2)')
    parser.add_argument('-lnlnQorder', type=int, default=2, help='order for interpolation in lnlnQ (default=2)')

    args = parser.parse_args()

    hlinfo = hp.lhapdf.load(args.pdf)
    p_lhapdf = lhapdf.mkPDF(args.pdf, 0)

    print("Checking that hlinfo is valid:", hlinfo is not None)

    # Overwrite the yorder and lnlnQorder interpolation orders
    hp.SetYLnlnQInterpOrders(args.yorder, args.lnlnQorder)

    # Evaluate the PDFs at some x values and print them
    xvals = [1e-5,1e-4,1e-3,1e-2,0.1,0.3,0.5,0.7,0.9]
    Q = 100.0

    print("")
    print("           Evaluating PDFs at Q =",Q, " GeV")
    print("-------------- HOPPET ------------------")
    if not hlinfo.QED:
        print("    x      u-ubar      d-dbar    2(ubr+dbr)    c+cbar       gluon")
        for ix in range(9):
            pdf_array = hp.Eval(xvals[ix], Q)
            print("{:7.1E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E}".format(
                xvals[ix],
                pdf_array[hp.iflv_u] - pdf_array[hp.iflv_ubar], 
                pdf_array[hp.iflv_d] - pdf_array[hp.iflv_dbar], 
                2 * (pdf_array[hp.iflv_dbar] + pdf_array[hp.iflv_ubar]),
                pdf_array[hp.iflv_cbar] + pdf_array[hp.iflv_c],
                pdf_array[hp.iflv_g]
            ))
        print("-------------- LHAPDF ------------------")
        for ix in range(9):
            pdf_dict = p_lhapdf.xfxQ(xvals[ix], Q)
            #print(type(pdf_dict)
            print("{:7.1E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E}".format(
                xvals[ix],
                pdf_dict[2] - pdf_dict[-2], 
                pdf_dict[1] - pdf_dict[-1], 
                2 * (pdf_dict[-1] + pdf_dict[-2]),
                pdf_dict[-4] + pdf_dict[4],
                pdf_dict[21]
            ))
    else:
        print("    x      u-ubar      d-dbar    2(ubr+dbr)    c+cbar       gluon      photon")
        for ix in range(9):
            pdf_array = hp.Eval(xvals[ix], Q)
            print("{:7.1E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E}".format(
                xvals[ix],
                pdf_array[hp.iflv_u] - pdf_array[hp.iflv_ubar], 
                pdf_array[hp.iflv_d] - pdf_array[hp.iflv_dbar], 
                2 * (pdf_array[hp.iflv_dbar] + pdf_array[hp.iflv_ubar]),
                pdf_array[hp.iflv_cbar] + pdf_array[hp.iflv_c],
                pdf_array[hp.iflv_g], pdf_array[hp.iflv_photon]
            ))
        print("-------------- LHAPDF ------------------")
        for ix in range(9):
            pdf_dict = p_lhapdf.xfxQ(xvals[ix], Q)
            #print(type(pdf_dict)
            print("{:7.1E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E} {:11.4E}".format(
                xvals[ix],
                pdf_dict[2] - pdf_dict[-2], 
                pdf_dict[1] - pdf_dict[-1], 
                2 * (pdf_dict[-1] + pdf_dict[-2]),
                pdf_dict[-4] + pdf_dict[4],
                pdf_dict[21], pdf_dict[22]
            ))
    


    #hp.WriteLHAPDFGrid("test_python",0)

    # Define grids for timing test
    npoints = 1000
    xvals_timing = np.logspace(np.log10(1e-5), np.log10(0.9), npoints)
    Qvals_timing = np.logspace(np.log10(1.0), np.log10(1000.0), npoints)

    # Timing HOPPET
    start_hoppet = time.perf_counter()
    for Q in Qvals_timing:
        for x in xvals_timing:
            pdf_array = hp.Eval(x, Q)
    end_hoppet = time.perf_counter()
    print(f"HOPPET evaluation time {(end_hoppet - start_hoppet)/npoints/npoints*1e9:.2f} ns")

    # Timing LHAPDF
    # Load the PDF from LHAPDF
    start_lhapdf = time.perf_counter()
    for Q in Qvals_timing:
        for x in xvals_timing:
            pdf_dict = p_lhapdf.xfxQ(x, Q)
            # The following call is faster but works only with a patch to LHAPDF
            # pdf_list = p_lhapdf.xfxQv(x, Q)
    end_lhapdf = time.perf_counter()
    print(f"LHAPDF evaluation time {(end_lhapdf - start_lhapdf)/npoints/npoints*1e9:.2f} ns")

    # AlphaS timings
    start_hoppet_as = time.perf_counter()
    for Q in Qvals_timing:
        as_hoppet = hp.AlphaS(Q)
    end_hoppet_as = time.perf_counter()
    print(f"HOPPET alphaS time {(end_hoppet_as - start_hoppet_as)/npoints*1e9:.2f} ns") 
    start_lhapdf_as = time.perf_counter()
    for Q in Qvals_timing:
        as_lhapdf = p_lhapdf.alphasQ(Q)
    end_lhapdf_as = time.perf_counter()
    print(f"LHAPDF alphaS time {(end_lhapdf_as - start_lhapdf_as)/npoints*1e9:.2f} ns")
    
    hp.DeleteAll()


if __name__ == "__main__":
    main()
