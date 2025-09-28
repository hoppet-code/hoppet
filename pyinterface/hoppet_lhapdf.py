try:
    import lhapdf
except ImportError:
    raise ImportError("The 'lhapdf' package is required for hoppet_lhapdf but is not installed.")

try:
    import numpy as np
except ImportError:
    raise ImportError("The 'numpy' package is required for hoppet_lhapdf but is not installed.")

import sys

try:
    import hoppet as hp
except ImportError:
    raise ImportError("The 'hoppet' extension is required for hoppet_lhapdf but is not installed or built.")

def load_lhapdf_assign_hoppet(lhapdfname, imem = 0):
    """
    Loads a PDF from LHAPDF, extracts relevant parameters, and assigns it to HOPPET.

    Parameters
    ----------
    lhapdfname : str
        The name of the LHAPDF set to load.
    imem : int, optional
        The member number of the PDF set (default is 0).

    This function sets up the HOPPET evolution with parameters from the LHAPDF set,
    assigns the PDF to HOPPET, and prints information about the configuration.
    """
    # Load the PDF from LHAPDF
    p_lhapdf = lhapdf.mkPDF(lhapdfname, imem)

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
    Q0 = Qmin
    asQ0 = p_lhapdf.alphasQ(Q0)

    # Print some info to the screen
    print(f"{lhapdfname} read succesfully with the following parameters extracted: \nNumber of loops: {nloop}\nxmin: {xmin}\nxmax: {xmax}\nQmin: {Qmin}\nQmax: {Qmax}\nmc: {mc}\nmb: {mb}\nmt: {mt}")
    
    # Now we start hoppet
    dy = 0.05
    print(f"Starting Hoppet with dy = {dy} and nloop = {nloop}")

    hp.SetPoleMassVFN(mc,mb,mt)
    print(f"Using Pole Mass Variable Flavour Number scheme with mc = {mc}, mb = {mb}, mt = {mt}")

    yorder = 2
    lnlnQorder = 2
    hp.SetYLnlnQInterpOrders(yorder, lnlnQorder)
    print(f"Set yorder = {yorder}, lnlnQorder = {lnlnQorder}")
    hp.Start(dy, nloop)
    
    print(f"Assigning PDF using hoppetAssign using Q0 = {Q0} GeV with as(Q0) = {asQ0}")
    hp.SetCoupling(asQ0, Q0, nloop)
    hp.Assign(lhapdf_interface)
