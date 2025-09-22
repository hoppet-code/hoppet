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

    uv = N_uv * pow(x, 0.8) * pow((1 - x), 3)
    dv = N_dv * pow(x, 0.8) * pow((1 - x), 4)
    dbar = N_db * pow(x, -0.1) * pow(1 - x, 6)
    ubar = dbar * (1 - x)

    pdf = np.zeros(13) # Initialise the array with zeros

    pdf[ 0+6] = N_g * pow(x, -0.1) * pow(1 - x, 5)
    pdf[-3+6] = 0.2 * (dbar + ubar)
    pdf[ 3+6] = pdf[-3 + 6]
    pdf[ 2+6] = uv + ubar
    pdf[-2+6] = ubar
    pdf[ 1+6] = dv + dbar
    pdf[-1+6] = dbar
    
# Overkill
    pdf[ 4+6] = 0.0
    pdf[ 5+6] = 0.0
    pdf[ 6+6] = 0.0
    pdf[-4+6] = 0.0
    pdf[-5+6] = 0.0
    pdf[-6+6] = 0.0

    return pdf

def main():
    Qmax = 13000.0
    nloop_coefs = 4
    xmur = 1.0
    xmuf = 1.0
    sc_choice = hp.scale_choice_Q # Uses Q as the central scale choice
    Qmin = 1.0

    # Set heavy flavour scheme
    mc = 1.414213563  # sqrt(2.0_dp) + epsilon
    mb = 4.5
    mt = 175.0
    hp.SetPoleMassVFN(mc, mb, mt)

    # Streamlined initialization
    # including  parameters for x-grid
    order = -6 # interpolation order, not perturbative order in alphas!
    ymax  = 16.0
    dy    = 0.05  # dble_val_opt("-dy",0.1_dp)
    dlnlnQ = dy/4.0
    nloop = 3 
    minQval = min(xmuf*Qmin, Qmin)
    maxQval = max(xmuf*Qmax, Qmax)

    # initialise the grid and dglap holder, using the streamlined
    # interface for simplicity
    hp.StartExtended(ymax,dy,minQval,maxQval,dlnlnQ,nloop,order,hp.factscheme_MSbar)

    # Setup all constants and parameters needed by the structure functions
    hp.StartStrFct(nloop_coefs)
    
    asQ0 = 0.35
    Q0 = np.sqrt(2.0)
    # Do the evolution. 
    hp.Evolve(asQ0, Q0, nloop, 1.0, hera_lhc, Q0)

    hp.InitStrFct(nloop_coefs, True, xmur, xmuf)

    # write out the structure functions
    ymax = np.log(1e5) 
    Q = 100.0
    write_f1_py(Q, ymax, 10, Q, Q, True)
    write_f2_py(Q, ymax, 10, Q, Q, True)
    write_f3_py(Q, ymax, 10, Q, Q, True)

    hp.DeleteAll()

def write_f1_py(Qtest, ymax, ny, muF=None, muR=None, use_sep_orders=False):
    print(f"# Q = {Qtest:10.4f}")
    if use_sep_orders:
        print("#     x        F1Wp(LO)    F1Wm(LO)    F1Wp(NLO)   F1Wm(NLO)  F1Wp(NNLO)  F1Wm(NNLO)"
              "  F1Wp(N3LO)  F1Wm(N3LO)    F1Z(LO)    F1Z(NLO)    F1Z(NNLO)"
              "   F1Z(N3LO)    F1γ(LO)    F1γ(NLO)    F1γ(NNLO)   F1γ(N3LO)   F1γZ(LO)"
              "   F1γZ(NLO)   F1γZ(NNLO)  F1γZ(N3LO)")
    else:
        print("# x F1Wp F1Wm F1Z F1γ F1γZ")

    for iy in range(ny, 0, -1):
        ytest = iy * ymax / ny
        xval = np.exp(-ytest)
        if use_sep_orders:
            res_lo = hp.StrFctLO(xval, Qtest, muR, muF)
            print(f"{xval:12.4e} {res_lo[hp.iF1Wp]:11.4e} {res_lo[hp.iF1Wm]:11.4e}", end=' ')
            F1Z_LO = res_lo[hp.iF1Z]
            F1EM_LO = res_lo[hp.iF1EM]
            F1gZ_LO = res_lo[hp.iF1gZ]

            res_nlo = hp.StrFctNLO(xval, Qtest, muR, muF)
            print(f"{res_nlo[hp.iF1Wp]:11.4e} {res_nlo[hp.iF1Wm]:11.4e}", end=' ')
            F1Z_NLO = res_nlo[hp.iF1Z]
            F1EM_NLO = res_nlo[hp.iF1EM]
            F1gZ_NLO = res_nlo[hp.iF1gZ]

            res_nnlo = hp.StrFctNNLO(xval, Qtest, muR, muF)
            print(f"{res_nnlo[hp.iF1Wp]:11.4e} {res_nnlo[hp.iF1Wm]:11.4e}", end=' ')
            F1Z_NNLO = res_nnlo[hp.iF1Z]
            F1EM_NNLO = res_nnlo[hp.iF1EM]
            F1gZ_NNLO = res_nnlo[hp.iF1gZ]

            res_n3lo = hp.StrFctN3LO(xval, Qtest, muR, muF)
            print(f"{res_n3lo[hp.iF1Wp]:11.4e} {res_n3lo[hp.iF1Wm]:11.4e}", end=' ')
            F1Z_N3LO = res_n3lo[hp.iF1Z]
            F1EM_N3LO = res_n3lo[hp.iF1EM]
            F1gZ_N3LO = res_n3lo[hp.iF1gZ]

            print(f"{F1Z_LO:11.4e} {F1Z_NLO:11.4e} {F1Z_NNLO:11.4e} {F1Z_N3LO:11.4e} "
                  f"{F1EM_LO:11.4e} {F1EM_NLO:11.4e} {F1EM_NNLO:11.4e} {F1EM_N3LO:11.4e} "
                  f"{F1gZ_LO:11.4e} {F1gZ_NLO:11.4e} {F1gZ_NNLO:11.4e} {F1gZ_N3LO:11.4e}")
        else:
            res = hp.StrFct(xval, Qtest, muR, muF)
            print(f"{xval:12.4e} {res[hp.iF1Wp]:11.4e} {res[hp.iF1Wm]:11.4e} {res[hp.iF1Z]:11.4e} {res[hp.iF1EM]:11.4e} {res[hp.iF1gZ]:11.4e}")
    print()
    print()

def write_f2_py(Qtest, ymax, ny, muF=None, muR=None, use_sep_orders=False):
    print(f"# Q = {Qtest:10.4f}")
    if use_sep_orders:
        print("#     x        F2Wp(LO)    F2Wm(LO)    F2Wp(NLO)   F2Wm(NLO)  F2Wp(NNLO)  F2Wm(NNLO)"
              "  F2Wp(N3LO)  F2Wm(N3LO)    F2Z(LO)    F2Z(NLO)    F2Z(NNLO)"
              "   F2Z(N3LO)    F2γ(LO)    F2γ(NLO)    F2γ(NNLO)   F2γ(N3LO)   F2γZ(LO)"
              "   F2γZ(NLO)   F2γZ(NNLO)  F2γZ(N3LO)")
    else:
        print("# x F2Wp F2Wm F2Z F2γ F2γZ")

    for iy in range(ny, 0, -1):
        ytest = iy * ymax / ny
        xval = np.exp(-ytest)
        if use_sep_orders:
            res_lo = hp.StrFctLO(xval, Qtest, muR, muF)
            print(f"{xval:12.4e} {res_lo[hp.iF2Wp]:11.4e} {res_lo[hp.iF2Wm]:11.4e}", end=' ')
            F2Z_LO = res_lo[hp.iF2Z]
            F2EM_LO = res_lo[hp.iF2EM]
            F2gZ_LO = res_lo[hp.iF2gZ]

            res_nlo = hp.StrFctNLO(xval, Qtest, muR, muF)
            print(f"{res_nlo[hp.iF2Wp]:11.4e} {res_nlo[hp.iF2Wm]:11.4e}", end=' ')
            F2Z_NLO = res_nlo[hp.iF2Z]
            F2EM_NLO = res_nlo[hp.iF2EM]
            F2gZ_NLO = res_nlo[hp.iF2gZ]

            res_nnlo = hp.StrFctNNLO(xval, Qtest, muR, muF)
            print(f"{res_nnlo[hp.iF2Wp]:11.4e} {res_nnlo[hp.iF2Wm]:11.4e}", end=' ')
            F2Z_NNLO = res_nnlo[hp.iF2Z]
            F2EM_NNLO = res_nnlo[hp.iF2EM]
            F2gZ_NNLO = res_nnlo[hp.iF2gZ]

            res_n3lo = hp.StrFctN3LO(xval, Qtest, muR, muF)
            print(f"{res_n3lo[hp.iF2Wp]:11.4e} {res_n3lo[hp.iF2Wm]:11.4e}", end=' ')
            F2Z_N3LO = res_n3lo[hp.iF2Z]
            F2EM_N3LO = res_n3lo[hp.iF2EM]
            F2gZ_N3LO = res_n3lo[hp.iF2gZ]

            print(f"{F2Z_LO:11.4e} {F2Z_NLO:11.4e} {F2Z_NNLO:11.4e} {F2Z_N3LO:11.4e} "
                  f"{F2EM_LO:11.4e} {F2EM_NLO:11.4e} {F2EM_NNLO:11.4e} {F2EM_N3LO:11.4e} "
                  f"{F2gZ_LO:11.4e} {F2gZ_NLO:11.4e} {F2gZ_NNLO:11.4e} {F2gZ_N3LO:11.4e}")
        else:
            res = hp.StrFct(xval, Qtest, muR, muF)
            print(f"{xval:12.4e} {res[hp.iF2Wp]:11.4e} {res[hp.iF2Wm]:11.4e} {res[hp.iF2Z]:11.4e} {res[hp.iF2EM]:11.4e} {res[hp.iF2gZ]:11.4e}")
    print()
    print()

def write_f3_py(Qtest, ymax, ny, muF=None, muR=None, use_sep_orders=False):
    print(f"# Q = {Qtest:10.4f}")
    if use_sep_orders:
        print("#     x        F3Wp(LO)    F3Wm(LO)    F3Wp(NLO)   F3Wm(NLO)  F3Wp(NNLO)  F3Wm(NNLO)"
              "  F3Wp(N3LO)  F3Wm(N3LO)    F3Z(LO)    F3Z(NLO)    F3Z(NNLO)"
              "   F3Z(N3LO)    F3γZ(LO)"
              "   F3γZ(NLO)   F3γZ(NNLO)  F3γZ(N3LO)")
    else:
        print("# x F3Wp F3Wm F3Z F3γZ")

    for iy in range(ny, 0, -1):
        ytest = iy * ymax / ny
        xval = np.exp(-ytest)
        if use_sep_orders:
            res_lo = hp.StrFctLO(xval, Qtest, muR, muF)
            print(f"{xval:12.4e} {res_lo[hp.iF3Wp]:11.4e} {res_lo[hp.iF3Wm]:11.4e}", end=' ')
            F3Z_LO = res_lo[hp.iF3Z]
            F3gZ_LO = res_lo[hp.iF3gZ]

            res_nlo = hp.StrFctNLO(xval, Qtest, muR, muF)
            print(f"{res_nlo[hp.iF3Wp]:11.4e} {res_nlo[hp.iF3Wm]:11.4e}", end=' ')
            F3Z_NLO = res_nlo[hp.iF3Z]
            F3gZ_NLO = res_nlo[hp.iF3gZ]

            res_nnlo = hp.StrFctNNLO(xval, Qtest, muR, muF)
            print(f"{res_nnlo[hp.iF3Wp]:11.4e} {res_nnlo[hp.iF3Wm]:11.4e}", end=' ')
            F3Z_NNLO = res_nnlo[hp.iF3Z]
            F3gZ_NNLO = res_nnlo[hp.iF3gZ]
    
            res_n3lo = hp.StrFctN3LO(xval, Qtest, muR, muF)
            print(f"{res_n3lo[hp.iF3Wp]:11.4e} {res_n3lo[hp.iF3Wm]:11.4e}", end=' ')
            F3Z_N3LO = res_n3lo[hp.iF3Z]
            F3gZ_N3LO = res_n3lo[hp.iF3gZ]

            print(f"{F3Z_LO:11.4e} {F3Z_NLO:11.4e} {F3Z_NNLO:11.4e} {F3Z_N3LO:11.4e} "
                  f"{F3gZ_LO:11.4e} {F3gZ_NLO:11.4e} {F3gZ_NNLO:11.4e} {F3gZ_N3LO:11.4e}")
        else:
            res = hp.StrFct(xval, Qtest, muR, muF)
            print(f"{xval:12.4e} {res[hp.iF3Wp]:11.4e} {res[hp.iF3Wm]:11.4e} {res[hp.iF3Z]:11.4e} {res[hp.iF3gZ]:11.4e}")
    print()
    print()

if __name__ == "__main__":
    main()
