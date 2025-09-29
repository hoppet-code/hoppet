#! /usr/bin/env python3
"""
This is a small example that uses the Hoppet and LHAPDF python
interfaces. It loads an LHAPDF set at a low scale and evolves it up
with Hoppet.
"""

import hoppet as hp
import lhapdf
import numpy as np
import argparse
import sys
import time

def load_lhapdf_assign_hoppet(lhapdfname, Q0in, dy, nloopin = 0, Q0_just_above_mb = False, 
                                    Q0_just_above_mc = False, 
                                    exact_nnlo_nf = False, exact_nnlo_splitting = False, 
                                    n3lo_splitting = '2410', FFN = -1, assign = False,
                                    yorder=6, lnlnQorder=4):
    # Load the PDF from LHAPDF
    p_lhapdf = lhapdf.mkPDF(lhapdfname, 0)

    # Now that we have the PDF we define the interface as needed by hoppet
    def lhapdf_interface(x, Q):
        pdf = np.zeros(13)
        lhapdf = p_lhapdf.xfxQ(None, x, Q)
        # Map HOPPET indices to LHAPDF PIDs
        pid_map = [ -6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6 ]
        for i, pid in enumerate(pid_map):
            pdf[i] = lhapdf.get(pid, 0.0)
        return pdf
    
    # Get some information from the PDF like order in QCD, masses etc.
    nloop = p_lhapdf.orderQCD + 1 # LHAPDF starts at 0
    xmin = p_lhapdf.xMin
    xmax = p_lhapdf.xMax
    Qmin = np.sqrt(p_lhapdf.q2Min)
    Qmax = np.sqrt(p_lhapdf.q2Max)
    mc = p_lhapdf.quarkThreshold(4)
    mb = p_lhapdf.quarkThreshold(5)
    mt = p_lhapdf.quarkThreshold(6)
    if(p_lhapdf.hasFlavor(6) == False): mt = 2*Qmax# If no top is defined set it to a high value
    
    # Get alphas at Q0, and check that Q0 is not below lahpdf Qmin
    Q0 = max(Q0in, Qmin)
    eps = 1e-4
    if Q0_just_above_mc:
        Q0 = mc + eps
    elif Q0_just_above_mb:
        Q0 = mb + eps
    asQ0 = p_lhapdf.alphasQ(Q0)

    # Print some info to the screen
    print(f"{lhapdfname} read succesfully with the following parameters extracted: \nNumber of loops: {nloop}\nxmin: {xmin}\nxmax: {xmax}\nQmin: {Qmin}\nQmax: {Qmax}\nmc: {mc}\nmb: {mb}\nmt: {mt}")
    
    if nloopin > 0:
        nloop = nloopin
        print(f"Overriding number of loops to nloop = {nloop}")
        
    # Now we start hoppet
    print(f"Starting Hoppet with dy = {dy} and nloop = {nloop}")

    # By default we use parametrised nf thresholds and splitting
    # functions (this only applies to the NNLO part, since at N3LO we
    # are currently forced to use exact nf but approximate splitting
    # functions).
    hp.SetExactDGLAP(exact_nnlo_nf, exact_nnlo_splitting)
    print(f"Using exact NNLO nf thresholds: {exact_nnlo_nf}, exact NNLO splitting functions: {exact_nnlo_splitting}")

    # n3lo splitting function approximation
    if nloop == 4:
        if n3lo_splitting == '2310':
            hp.SetApproximateDGLAPN3LO(hp.n3lo_splitting_approximation_up_to_2310_05744)
        elif n3lo_splitting == '2404':
            hp.SetApproximateDGLAPN3LO(hp.n3lo_splitting_approximation_up_to_2404_09701)
        elif n3lo_splitting == '2410':
            hp.SetApproximateDGLAPN3LO(hp.n3lo_splitting_approximation_up_to_2410_08089) # This is the default value in hoppet at the moment
        else:
            print(f"Error: Unknown n3lo-splitting value {n3lo_splitting}")
            sys.exit(1)
        print(f"N3LO splitting function approximation: {n3lo_splitting}")

    # Right now I can't see a way to find the flavour scheme in the
    # LHAPDF interface. For now we assume it is variable unless the
    # user specifies FFN 
    if FFN > 0:
        hp.SetFFN(FFN)
        print(f"Using Fixed Flavour Number scheme with nf = {FFN}")
    else:
        hp.SetPoleMassVFN(mc,mb,mt)
        print(f"Using Pole Mass Variable Flavour Number scheme with mc = {mc}, mb = {mb}, mt = {mt}")

    hp.SetYLnlnQInterpOrders(yorder, lnlnQorder)
    print(f"Set yorder = {yorder}, lnlnQorder = {lnlnQorder}")
    hp.Start(dy, nloop)
    
    if assign:
        print(f"Assigning PDF using hoppetAssign using Q0 = {Q0} GeV with as(Q0) = {asQ0}")
        hp.SetCoupling(asQ0, Q0, nloop)
        hp.Assign(lhapdf_interface)
    else:
        print(f"Evolving PDF from Q0 = {Q0} GeV with as(Q0) = {asQ0}")
        hp.Evolve(asQ0, Q0, nloop, 1.0, lhapdf_interface, Q0)

def main():
    # Get commandline
    parser = argparse.ArgumentParser(description="Check of an LHAPDF grid against HOPPET evolution.")
    parser.add_argument('-pdf', required=True, help='LHAPDF set name (required, ex. NNPDF30_nnlo_as_0118)')
    parser.add_argument('-dy', type=float, default=0.05, help='dy for HOPPET evolution (default: 0.05)')
    parser.add_argument('-Q0', type=float, default=1.0, help='Initial Q0 value (default: Qmin from LHAPDF)')
    parser.add_argument('-yorder', type=int, default=2, help='order for interpolation in y=ln1/x (default=2)')
    parser.add_argument('-lnlnQorder', type=int, default=2, help='order for interpolation in lnlnQ (default=2)')

    args = parser.parse_args()

    load_lhapdf_assign_hoppet(args.pdf, args.Q0, args.dy, yorder=args.yorder, lnlnQorder=args.lnlnQorder, assign = True)

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
    print("")
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
    p_lhapdf = lhapdf.mkPDF(args.pdf, 0)
    start_lhapdf = time.perf_counter()
    for Q in Qvals_timing:
        for x in xvals_timing:
            pdf_dict = p_lhapdf.xfxQ(None, x, Q)
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
