try:
    import lhapdf
except ImportError:
    raise ImportError("The 'lhapdf' package is required for the hoppet LHAPDF interface but is not installed.")

try:
    import numpy as np
except ImportError:
    raise ImportError("The 'numpy' package is required for the hoppet LHAPDF interface but is not installed.")

try:
    import math
except ImportError:
    raise ImportError("The 'math' package is required for the hoppet LHAPDF interface but is not installed.")

try:
    import hoppet as hp
except ImportError:
    raise ImportError("The 'hoppet' extension is required for the hoppet LHAPDF interface but is not installed or built.")

def load(lhapdfname, imem = 0):
    """
    Loads a PDF from LHAPDF, extracts relevant parameters, and assigns
    it to HOPPET.

    Parameters
    ----------
    lhapdfname : str
        The name of the LHAPDF set to load.
    imem : int, optional
        The member number of the PDF set (default is 0).

    Returns
    -------
    hoppet_lhapdf
        An instance of the hoppet_lhapdf class containing the LHAPDF
        and associated meta data.

    This function sets up the HOPPET evolution with parameters from
    the LHAPDF set, assigns the PDF to HOPPET, and prints information
    about the configuration.
    """
    # Load the PDF from LHAPDF into the hoppet_lhapdf class
    hl = hoppet_lhapdf(lhapdfname, imem)

    # Now that we have the PDF we define the interface as needed by hoppet
    def lhapdf_interface(x, Q):
        pdf = np.zeros(13)
        lhapdf = hl.pdf.xfxQ(None, x, Q)
        # Map HOPPET indices to LHAPDF PIDs
        pid_map = [ -6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6 ]
        for i, pid in enumerate(pid_map):
            pdf[i] = lhapdf.get(pid, 0.0)
        return pdf
    
    # Print some info to the screen
    print(f"{lhapdfname} read succesfully with the following parameters extracted: \nNumber of loops: {hl.nloop}\nxmin: {hl.xmin}\nxmax: {hl.xmax}\nQmin: {hl.Qmin}\nQmax: {hl.Qmax}\nmc: {hl.thresholds[4]}\nmb: {hl.thresholds[5]}\nmt: {hl.thresholds[6]}")

    # Now we start hoppet
    dy = 0.05
    ymax = float(math.ceil(np.log(1.0/hl.xmin)))
    if ymax > 15.0:
        dlnlnQ = dy/8.0
    else:
        dlnlnQ = dy/4.0

    hp.SetPoleMassVFN(hl.thresholds[4], hl.thresholds[5], hl.thresholds[6])
    print(f"Using Pole Mass Variable Flavour Number scheme with mc = {hl.thresholds[4]}, mb = {hl.thresholds[5]}, mt = {hl.thresholds[6]}")

    # Set default interpolation orders, which together with dy=0.05
    # guarantee O(10^-4) accuracy
    order = -6
    yorder = 2
    lnlnQorder = 2
    print(f"Setting interpolation orders yorder = {yorder}, lnlnQorder = {lnlnQorder}")
    hp.SetYLnlnQInterpOrders(yorder, lnlnQorder)

    print(f"Starting Hoppet with ymax = {ymax} and dy = {dy} and nloop = {hl.nloop} and dlnlnQ = {dlnlnQ} and order = {order}")
    hp.StartExtended(ymax, dy, hl.Qmin, hl.Qmax, dlnlnQ, hl.nloop, order, hp.factscheme_MSbar)

    print(f"Assigning PDF using hoppetAssign with a coupling as(Q0 = {hl.MZ}) = {hl.asMZ}")
    hp.SetCoupling(hl.asMZ, hl.MZ, hl.nloop)
    hp.Assign(lhapdf_interface)

    return hl

# Small class that contains the PDF from LHAPDF and some meta data
class hoppet_lhapdf:
    def __init__(self, pdfname, member):
        self.pdf = lhapdf.mkPDF(pdfname, member)
        self.name = pdfname
        self.Qmin = np.sqrt(self.pdf.q2Min)
        self.Qmax = np.sqrt(self.pdf.q2Max)
        self.xmin = self.pdf.xMin
        self.xmax = self.pdf.xMax
        self.nloop = self.pdf.orderQCD + 1 # LHAPDF starts at 0
        self.thresholds = {fl: self.pdf.quarkThreshold(fl) for fl in [1,2,3,4,5,6]}
        if(self.pdf.hasFlavor(6) == False): self.thresholds[6] = 2*self.Qmax # If no top is defined set it to a high value
        self.MZ = 91.188 # Z mass in GeV
        self.asMZ = self.pdf.alphasQ(self.MZ)

    def __del__(self):
        # No explicit resource management needed, but method provided for extensibility
        pass

    def __repr__(self):
        return (f"hoppet_lhapdf(name={self.name}, Qmin={self.Qmin}, Qmax={self.Qmax}, "
                f"xmin={self.xmin}, xmax={self.xmax}, nloop={self.nloop}, "
                f"thresholds={self.thresholds}, MZ={self.MZ}, asMZ={self.asMZ})")