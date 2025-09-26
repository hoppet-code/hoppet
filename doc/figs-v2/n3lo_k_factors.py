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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import time

def main():
    # Define flavour names for indices 0 to 12
    flavour_names = [r'$\bar{t}$', r'$\bar{b}$', r'$\bar{c}$', r'$\bar{s}$', r'$\bar{u}$', r'$\bar{d}$', r'$g$', r'$d$', r'$u$', r'$s$', r'$c$', r'$b$', r'$t$']
    # Get commandline
    parser = argparse.ArgumentParser(description="Check of an LHAPDF grid against HOPPET evolution.")
    parser.add_argument('-pdf', required=True, help='LHAPDF set name (required, ex. NNPDF30_nnlo_as_0118)')
    parser.add_argument('-dy', type=float, default=0.05, help='dy for HOPPET evolution (default: 0.05)')
    parser.add_argument('-Q0', type=float, default=1.0, help='Initial Q0 value (default: Qmin from LHAPDF)')
    parser.add_argument('-yorder', type=int, default=6, help='order for interpolation in y=ln1/x (default of -1 uses same as evolution)')
    parser.add_argument('-lnlnQorder', type=int, default=4, help='order for interpolation in lnlnQ (default=4)')

    args = parser.parse_args()

    # x points to plot
    npoints = 100
    x_log = np.logspace(np.log10(1e-5), np.log10(0.2), npoints)
    x_lin = np.linspace(1e-2, 0.9, npoints)

    # Arrays to hold results
    nnlo_log = np.zeros((npoints, 13))
    nnlo_lin = np.zeros((npoints, 13))
    n3lo_log = np.zeros((npoints, 13))
    n3lo_lin = np.zeros((npoints, 13))
    n3lo_novfn3_log = np.zeros((npoints, 13))
    n3lo_novfn3_lin = np.zeros((npoints, 13))

    nnlo_log, nnlo_lin = start_run_hoppet_fill_arrays(args.dy, 3, False, npoints, x_lin, x_log)
    n3lo_novfn3_log, n3lo_novfn3_lin = start_run_hoppet_fill_arrays(args.dy, 4, False, npoints, x_lin, x_log)
    n3lo_log, n3lo_lin = start_run_hoppet_fill_arrays(args.dy, 4, True, npoints, x_lin, x_log)

    # Arrays to hold results for PDF4LHC21_40
    pdf4lhc21_nnlo_log = np.zeros((npoints, 13))
    pdf4lhc21_nnlo_lin = np.zeros((npoints, 13))
    pdf4lhc21_n3lo_log = np.zeros((npoints, 13))
    pdf4lhc21_n3lo_lin = np.zeros((npoints, 13))
    pdf4lhc21_n3lo_novfn3_log = np.zeros((npoints, 13))
    pdf4lhc21_n3lo_novfn3_lin = np.zeros((npoints, 13))

    load_lhapdf_start_evolve_hoppet(args.pdf, args.Q0, args.dy, yorder=args.yorder, lnlnQorder=args.lnlnQorder, vfn3 = False)
    pdf4lhc21_nnlo_log, pdf4lhc21_nnlo_lin = hoppet_fill_array(npoints, x_lin, x_log)
    load_lhapdf_start_evolve_hoppet(args.pdf, args.Q0, args.dy, yorder=args.yorder, lnlnQorder=args.lnlnQorder, vfn3 = False, nloopin=4)
    pdf4lhc21_n3lo_novfn3_log, pdf4lhc21_n3lo_novfn3_lin = hoppet_fill_array(npoints, x_lin, x_log)
    load_lhapdf_start_evolve_hoppet(args.pdf, args.Q0, args.dy, yorder=args.yorder, lnlnQorder=args.lnlnQorder, vfn3 = True, nloopin=4)
    pdf4lhc21_n3lo_log, pdf4lhc21_n3lo_lin = hoppet_fill_array(npoints, x_lin, x_log)

    with PdfPages("n3lo_nnlo_ratio.pdf") as pdf:
        # Log-spaced x plot
        ratios, labels = compute_ratios_labels(n3lo_log, nnlo_log, flavour_names)
        plot_ratios(x_log, ratios, labels, 'Benchmark initial condition', 'hoppet v2.0.0', pdf=pdf, logx=True)

        # Linearly spaced x plot
        ratios_lin, labels_lin = compute_ratios_labels(n3lo_lin, nnlo_lin, flavour_names)
        plot_ratios(x_lin, ratios_lin, labels_lin, 'Benchmark initial condition', 'hoppet v2.0.0', pdf=pdf, logx=False)
        
        # Log-spaced x plot
        ratios, labels = compute_ratios_labels(n3lo_novfn3_log, nnlo_log, flavour_names)
        plot_ratios(x_log, ratios, labels, 'Benchmark initial condition (no VFN3)', 'hoppet v2.0.0', pdf=pdf, logx=True)

        # Linearly spaced x plot
        ratios_lin, labels_lin = compute_ratios_labels(n3lo_novfn3_lin, nnlo_lin, flavour_names)
        plot_ratios(x_lin, ratios_lin, labels_lin, 'Benchmark initial condition (no VFN3)', 'hoppet v2.0.0', pdf=pdf, logx=False)

        # Log-spaced x plot
        ratios, labels = compute_ratios_labels(pdf4lhc21_n3lo_log, pdf4lhc21_nnlo_log, flavour_names)
        plot_ratios(x_log, ratios, labels, 'PDF4LHC21_40 initial condition', 'hoppet v2.0.0', pdf=pdf, logx=True)

        # Linearly spaced x plot
        ratios_lin, labels_lin = compute_ratios_labels(pdf4lhc21_n3lo_lin, pdf4lhc21_nnlo_lin, flavour_names)
        plot_ratios(x_lin, ratios_lin, labels_lin, 'PDF4LHC21_40 initial condition', 'hoppet v2.0.0', pdf=pdf, logx=False)

        # Log-spaced x plot
        ratios, labels = compute_ratios_labels(pdf4lhc21_n3lo_novfn3_log, pdf4lhc21_nnlo_log, flavour_names)
        plot_ratios(x_log, ratios, labels, 'PDF4LHC21_40 initial condition (no VFN3)', 'hoppet v2.0.0', pdf=pdf, logx=True)

        # Linearly spaced x plot
        ratios_lin, labels_lin = compute_ratios_labels(pdf4lhc21_n3lo_novfn3_lin, pdf4lhc21_nnlo_lin, flavour_names)
        plot_ratios(x_lin, ratios_lin, labels_lin, 'PDF4LHC21_40 initial condition (no VFN3)', 'hoppet v2.0.0', pdf=pdf, logx=False)


    #load_lhapdf_start_evolve_hoppet(args.pdf, args.Q0, args.dy,yorder=args.yorder, lnlnQorder=args.lnlnQorder)


    hp.DeleteAll()

def start_run_hoppet_fill_arrays(dy, nloop, vfn3, npoints, x_lin, x_log):
    hp.SetN3LOnfthresholds(vfn3) # For debugging as this is much faster
    hp.Start(dy, nloop)

    Q0 = np.sqrt(2.0)
    asQ0 = 0.35
    print(f"Starting Hoppet with dy = {dy} and nloop = {nloop}")
    print(f"Evolving benchmark PDF from Q0 = {Q0} GeV with as(Q0) = {asQ0}")
    hp.Evolve(asQ0, Q0, nloop, 1.0, hp.BenchmarkPDFunpol, Q0)

    return hoppet_fill_array(npoints, x_lin, x_log)

def hoppet_fill_array(npoints, x_lin, x_log):
    # Arrays to hold results
    out_log = np.zeros((npoints, 13))
    out_lin = np.zeros((npoints, 13))

    for i, x in enumerate(x_log):
        out_log[i, :] = hp.Eval(x, 100.0)
    
    for i, x in enumerate(x_lin):
        out_lin[i, :] = hp.Eval(x, 100.0)

    return out_log, out_lin

def compute_ratios_labels(n3lo_arr, nnlo_arr, flavour_names):
    ratios = []
    labels = []
    # -5 and 5: sum (indices 1 and 11)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio_m5_p5 = np.where(
            (nnlo_arr[:, 1] + nnlo_arr[:, 11]) != 0,
            (n3lo_arr[:, 1] + n3lo_arr[:, 11]) / (nnlo_arr[:, 1] + nnlo_arr[:, 11]),
            np.nan
        )
    ratios.append(ratio_m5_p5)
    labels.append(r'$\bar{b}+b$')
    # -4 and 4: sum (indices 2 and 10)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio_m4_p4 = np.where(
            (nnlo_arr[:, 2] + nnlo_arr[:, 10]) != 0,
            (n3lo_arr[:, 2] + n3lo_arr[:, 10]) / (nnlo_arr[:, 2] + nnlo_arr[:, 10]),
            np.nan
        )
    ratios.append(ratio_m4_p4)
    labels.append(r'$\bar{c}+c$')
    # All other flavours individually (skip -6, 6, -5, 5, -4, 4)
    for i in range(13):
        if i in [0, 12, 1, 11, 2, 10]:
            continue
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.where(nnlo_arr[:, i] != 0, n3lo_arr[:, i] / nnlo_arr[:, i], np.nan)
        ratios.append(ratio)
        labels.append(flavour_names[i])
    return ratios, labels

def plot_ratios(xvals, ratios, labels, title, extra_label, pdf=None, logx=False):
    plt.figure(figsize=(10, 6))
    for ratio, label in zip(ratios, labels):
        plt.plot(xvals, ratio, label=label)
    if logx:
        plt.xscale('log')
    plt.xlabel('x')
    plt.ylabel(r'N$^3$LO / NNLO')
    plt.title(title)
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    plt.text(0.98, 0.98, extra_label, color='grey', alpha=0.6,
             fontsize=14, ha='right', va='top', transform=plt.gca().transAxes)
    plt.tight_layout()
    if pdf is not None:
        pdf.savefig()
        plt.close()
    else:
        plt.show()


def load_lhapdf_start_evolve_hoppet(lhapdfname, Q0in, dy, nloopin = 0, Q0_just_above_mb = False, 
                                    Q0_just_above_mc = False, 
                                    exact_nnlo_nf = False, exact_nnlo_splitting = False, 
                                    n3lo_splitting = '2410', FFN = -1, assign = False,
                                    yorder=6, lnlnQorder=4, vfn3 = False):
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
    hp.SetN3LOnfthresholds(vfn3)
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

if __name__ == "__main__":
    main()
