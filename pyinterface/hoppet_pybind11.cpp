#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <Python.h>
#include "hoppet.h"
#include <vector>
#include <functional>

namespace py = pybind11;

// Global variables for array management
const unsigned int pdf_len = 13; 
const unsigned int qed_pdf_len = 18;
const unsigned int str_fnc_len = 14;
const unsigned int str_fnc_flav_len = 3;
unsigned int py_pdf_len = pdf_len;
unsigned int py_str_fnc_len = str_fnc_len;

double *global_pdf = nullptr;
double *global_str_fnc = nullptr;

// Global callback storage
std::function<std::vector<double>(double, double)> pdf_callback;
std::function<std::vector<double>(double, double)> str_fnc_callback;

// Wrapper functions to bridge Python and C callbacks
static void pdf_subroutine_wrapper(const double &x, const double &Q, double *res) {
    if (pdf_callback) {
        auto result = pdf_callback(x, Q);
        for (size_t i = 0; i < py_pdf_len && i < result.size(); i++) {
            res[i] = result[i];
        }
    }
}

static void str_fnc_subroutine_wrapper(const double &x, const double &Q, double *res) {
    if (str_fnc_callback) {
        auto result = str_fnc_callback(x, Q);
        for (size_t i = 0; i < py_str_fnc_len && i < result.size(); i++) {
            res[i] = result[i];
        }
    }
}

// Initialize and free global arrays
void init_global_pdf() {
    if (global_pdf == nullptr) {
        global_pdf = new double[py_pdf_len];
    }
}

void free_global_pdf() {
    if (global_pdf != nullptr) {
        delete[] global_pdf;
        global_pdf = nullptr;
    }
}

void init_global_str_fnc() {
    if (global_str_fnc == nullptr) {
        global_str_fnc = new double[py_str_fnc_len];
    }
}

void free_global_str_fnc() {
    if (global_str_fnc != nullptr) {
        delete[] global_str_fnc;
        global_str_fnc = nullptr;
    }
}

// Convert arrays to Python lists
std::vector<double> pdf_to_vector(double *pdf) {
    std::vector<double> result(py_pdf_len);
    for (unsigned int i = 0; i < py_pdf_len; i++) {
        result[i] = pdf[i];
    }
    return result;
}

std::vector<double> str_fnc_to_vector(double *str_fnc) {
    std::vector<double> result(py_str_fnc_len);
    for (unsigned int i = 0; i < py_str_fnc_len; i++) {
        result[i] = str_fnc[i];
    }
    return result;
}

// Wrapper functions for the Python interface
void SetQED(const int & withqed, const int & qcdqed, const int & plq) {
    if (global_pdf != nullptr) {
        throw std::runtime_error("Cannot change QED settings after PDF initialization");
    }
    hoppetSetQED(withqed, qcdqed, plq);
    if (withqed) {
        py_pdf_len = qed_pdf_len;
    } else {
        py_pdf_len = pdf_len;
    }
}

void Start(const double & dy, const int & nloop) {
    init_global_pdf();
    hoppetStart(dy, nloop);
}

void StartExtended(const double & ymax, const double & dy, const double & Qmin,
                   const double & Qmax, const double & dlnlnQ, const int & nloop,
                   const int & order, const int & factscheme) {
    init_global_pdf();
    hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, factscheme);
}

void DeleteAll() {
    free_global_pdf();
    free_global_str_fnc();
    pdf_callback = nullptr;
    str_fnc_callback = nullptr;
    hoppetDeleteAll();
}

void Assign(const std::function<std::vector<double>(double, double)>& callback) {
    pdf_callback = callback;
    hoppetAssign(pdf_subroutine_wrapper);
}

void Evolve(const double & asQ0, const double & Q0alphas, const int & nloop,
            const double & muR_Q, const std::function<std::vector<double>(double, double)>& callback,
            const double & Q0pdf) {
    pdf_callback = callback;
    hoppetEvolve(asQ0, Q0alphas, nloop, muR_Q, pdf_subroutine_wrapper, Q0pdf);
}

void CachedEvolve(const std::function<std::vector<double>(double, double)>& callback) {
    pdf_callback = callback;
    hoppetCachedEvolve(pdf_subroutine_wrapper);
}

std::vector<double> Eval(const double & x, const double & Q) {
    if (global_pdf == nullptr) {
        throw std::runtime_error("Global pdf array is not initialized");
    }
    hoppetEval(x, Q, global_pdf);
    return pdf_to_vector(global_pdf);
}

std::vector<double> EvalSplit(const double & x, const double & Q, const int & iloop, const int & nf) {
    if (global_pdf == nullptr) {
        throw std::runtime_error("Global pdf array is not initialized");
    }
    hoppetEvalSplit(x, Q, iloop, nf, global_pdf);
    return pdf_to_vector(global_pdf);
}

std::vector<double> BenchmarkPDFunpol(const double & x, const double & Q) {
    if (global_pdf == nullptr) {
        throw std::runtime_error("Global pdf array is not initialized");
    }
    hoppetBenchmarkPDFunpol(x, Q, global_pdf);
    return pdf_to_vector(global_pdf);
}

void InitStrFct(const int & order_max, const int & separate_orders, const double & xR, const double & xF) {
    init_global_str_fnc();
    hoppetInitStrFct(order_max, separate_orders, xR, xF);
}

void InitStrFctFlav(const int & order_max, const int & separate_orders, const double & xR, const double & xF, const int & flavour_decomposition) {
    py_str_fnc_len = str_fnc_flav_len;
    init_global_str_fnc();
    hoppetInitStrFctFlav(order_max, separate_orders, xR, xF, flavour_decomposition);
}

std::vector<double> StrFct(const double & x, const double & Q, const double & muR_in, const double & muF_in) {
    if (global_str_fnc == nullptr) {
        throw std::runtime_error("Global structure function array is not initialized");
    }
    hoppetStrFct(x, Q, muR_in, muF_in, global_str_fnc);
    return str_fnc_to_vector(global_str_fnc);
}

std::vector<double> StrFctNoMu(const double & x, const double & Q) {
    if (global_str_fnc == nullptr) {
        throw std::runtime_error("Global structure function array is not initialized");
    }
    hoppetStrFctNoMu(x, Q, global_str_fnc);
    return str_fnc_to_vector(global_str_fnc);
}

std::vector<double> StrFctLO(const double & x, const double & Q, const double & muR_in, const double & muF_in) {
    if (global_str_fnc == nullptr) {
        throw std::runtime_error("Global structure function array is not initialized");
    }
    hoppetStrFctLO(x, Q, muR_in, muF_in, global_str_fnc);
    return str_fnc_to_vector(global_str_fnc);
}

std::vector<double> StrFctNLO(const double & x, const double & Q, const double & muR_in, const double & muF_in) {
    if (global_str_fnc == nullptr) {
        throw std::runtime_error("Global structure function array is not initialized");
    }
    hoppetStrFctNLO(x, Q, muR_in, muF_in, global_str_fnc);
    return str_fnc_to_vector(global_str_fnc);
}

std::vector<double> StrFctFlav(const double & x, const double & Q, const double & muR_in, const double & muF_in, const int & flav) {
    if (global_str_fnc == nullptr) {
        throw std::runtime_error("Global structure function array is not initialized");
    }
    hoppetStrFctFlav(x, Q, muR_in, muF_in, flav, global_str_fnc);
    return str_fnc_to_vector(global_str_fnc);
}

std::vector<double> StrFctNoMuFlav(const double & x, const double & Q, const int & flav) {
    if (global_str_fnc == nullptr) {
        throw std::runtime_error("Global structure function array is not initialized");
    }
    hoppetStrFctNoMuFlav(x, Q, flav, global_str_fnc);
    return str_fnc_to_vector(global_str_fnc);
}

std::vector<double> StrFctLOFlav(const double & x, const double & Q, const double & muR_in, const double & muF_in, const int & flav) {
    if (global_str_fnc == nullptr) {
        throw std::runtime_error("Global structure function array is not initialized");
    }
    hoppetStrFctLOFlav(x, Q, muR_in, muF_in, flav, global_str_fnc);
    return str_fnc_to_vector(global_str_fnc);
}

std::vector<double> StrFctNLOFlav(const double & x, const double & Q, const double & muR_in, const double & muF_in, const int & flav) {
    if (global_str_fnc == nullptr) {
        throw std::runtime_error("Global structure function array is not initialized");
    }
    hoppetStrFctNLOFlav(x, Q, muR_in, muF_in, flav, global_str_fnc);
    return str_fnc_to_vector(global_str_fnc);
}

std::vector<double> StrFctNNLO(const double & x, const double & Q, const double & muR_in, const double & muF_in) {
    if (global_str_fnc == nullptr) {
        throw std::runtime_error("Global structure function array is not initialized");
    }
    hoppetStrFctNNLO(x, Q, muR_in, muF_in, global_str_fnc);
    return str_fnc_to_vector(global_str_fnc);
}

std::vector<double> StrFctN3LO(const double & x, const double & Q, const double & muR_in, const double & muF_in) {
    if (global_str_fnc == nullptr) {
        throw std::runtime_error("Global structure function array is not initialized");
    }
    hoppetStrFctN3LO(x, Q, muR_in, muF_in, global_str_fnc);
    return str_fnc_to_vector(global_str_fnc);
}

PYBIND11_MODULE(hoppet_pybind11, m) {
    m.doc() = R"pbdoc(
A Higher Order Perturbative Parton Evolution Toolkit

HOPPET is a Fortran 95 package for carrying out DGLAP evolution and other 
common manipulations of parton distribution functions (PDFs).

Citation:
G.P. Salam, J. Rojo, 'A Higher Order Perturbative Parton Evolution Toolkit (HOPPET)', 
Comput. Phys. Commun. 180 (2009) 120-156, arXiv:0804.3755

and                                                       

A. Karlberg, P. Nason, G.P. Salam, G. Zanderighi & F. Dreyer (arXiv:2509.XXXXX). 

Example:

   import hoppet as hp
   import numpy as np
   
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
   
For more examples look at https://github.com/hoppet-code/hoppet/tree/master/example_py	
    )pbdoc";

    m.attr("__author__") = "Frederic Dreyer, Alexander Karlberg, Paolo Nason, Juan Rojo, Gavin Salam, Giulia Zanderighi";

    // Structure function indices constants
    m.attr("iF1Wp") = hoppet::iF1Wp;
    m.attr("iF2Wp") = hoppet::iF2Wp;
    m.attr("iF3Wp") = hoppet::iF3Wp;
    m.attr("iF1Wm") = hoppet::iF1Wm;
    m.attr("iF2Wm") = hoppet::iF2Wm;
    m.attr("iF3Wm") = hoppet::iF3Wm;
    m.attr("iF1Z") = hoppet::iF1Z;
    m.attr("iF2Z") = hoppet::iF2Z;
    m.attr("iF3Z") = hoppet::iF3Z;
    m.attr("iF1EM") = hoppet::iF1EM;
    m.attr("iF2EM") = hoppet::iF2EM;
    m.attr("iF1gZ") = hoppet::iF1gZ;
    m.attr("iF2gZ") = hoppet::iF2gZ;
    m.attr("iF3gZ") = hoppet::iF3gZ;

    // Scale choice constants
    m.attr("scale_choice_fixed") = hoppet::scale_choice_fixed;
    m.attr("scale_choice_Q") = hoppet::scale_choice_Q;
    m.attr("scale_choice_arbitrary") = hoppet::scale_choice_arbitrary;

    // Splitting function constants
    m.attr("nnlo_splitting_exact") = hoppet::nnlo_splitting_exact;
    m.attr("nnlo_splitting_param") = hoppet::nnlo_splitting_param;
    m.attr("nnlo_splitting_Nfitav") = hoppet::nnlo_splitting_Nfitav;
    m.attr("nnlo_splitting_Nfiterr1") = hoppet::nnlo_splitting_Nfiterr1;
    m.attr("nnlo_splitting_Nfiterr2") = hoppet::nnlo_splitting_Nfiterr2;

    m.attr("n3lo_splitting_exact") = hoppet::n3lo_splitting_exact;
    m.attr("n3lo_splitting_param") = hoppet::n3lo_splitting_param;
    m.attr("n3lo_splitting_Nfitav") = hoppet::n3lo_splitting_Nfitav;
    m.attr("n3lo_splitting_Nfiterr1") = hoppet::n3lo_splitting_Nfiterr1;
    m.attr("n3lo_splitting_Nfiterr2") = hoppet::n3lo_splitting_Nfiterr2;
    m.attr("n3lo_splitting_approximation_up_to_2310_05744") = hoppet::n3lo_splitting_approximation_up_to_2310_05744;
    m.attr("n3lo_splitting_approximation_up_to_2404_09701") = hoppet::n3lo_splitting_approximation_up_to_2404_09701;
    m.attr("n3lo_splitting_approximation_up_to_2410_08089") = hoppet::n3lo_splitting_approximation_up_to_2410_08089;

    // Threshold constants
    m.attr("nnlo_nfthreshold_exact") = hoppet::nnlo_nfthreshold_exact;
    m.attr("nnlo_nfthreshold_param") = hoppet::nnlo_nfthreshold_param;
    m.attr("n3lo_nfthreshold_on") = hoppet::n3lo_nfthreshold_on;
    m.attr("n3lo_nfthreshold_off") = hoppet::n3lo_nfthreshold_off;

    // Factorization scheme constants
    m.attr("factscheme_MSbar") = hoppet::factscheme_MSbar;
    m.attr("factscheme_DIS") = hoppet::factscheme_DIS;
    m.attr("factscheme_PolMSbar") = hoppet::factscheme_PolMSbar;
    m.attr("factscheme_FragMSbar") = hoppet::factscheme_FragMSbar;

    // Main interface functions
    m.def("SetQED", &SetQED, "Setup QED PDFs");
    m.def("Start", &Start, "Initialize hoppet with basic parameters");
    m.def("StartExtended", &StartExtended, "Initialize hoppet with extended parameters");
    m.def("DeleteAll", &DeleteAll, "Delete all hoppet objects");
    
    m.def("Assign", &Assign, "Assign a PDF subroutine");
    m.def("Evolve", &Evolve, "Evolve PDFs with given parameters");
    m.def("CachedEvolve", &CachedEvolve, "Evolve using cached evolution");
    
    m.def("Eval", &Eval, "Evaluate PDF at x, Q");
    m.def("EvalSplit", &EvalSplit, "Evaluate splitting function");
    m.def("BenchmarkPDFunpol", &BenchmarkPDFunpol, "Benchmark unpolarized PDF");

    // Structure function initialization and evaluation
    m.def("InitStrFct", &InitStrFct, "Initialize structure functions");
    m.def("InitStrFctFlav", &InitStrFctFlav, "Initialize flavor-decomposed structure functions");
    
    m.def("StrFct", &StrFct, "Calculate structure function");
    m.def("StrFctNoMu", &StrFctNoMu, "Calculate structure function without mu parameters");
    m.def("StrFctLO", &StrFctLO, "Calculate LO structure function");
    m.def("StrFctNLO", &StrFctNLO, "Calculate NLO structure function");
    m.def("StrFctNNLO", &StrFctNNLO, "Calculate NNLO structure function");
    m.def("StrFctN3LO", &StrFctN3LO, "Calculate N3LO structure function");
    
    m.def("StrFctFlav", &StrFctFlav, "Calculate flavor-decomposed structure function");
    m.def("StrFctNoMuFlav", &StrFctNoMuFlav, "Calculate flavor-decomposed structure function without mu");
    m.def("StrFctLOFlav", &StrFctLOFlav, "Calculate LO flavor-decomposed structure function");
    m.def("StrFctNLOFlav", &StrFctNLOFlav, "Calculate NLO flavor-decomposed structure function");

    // Direct C functions (renamed to remove hoppet prefix and underscores)
    m.def("PreEvolve", &hoppetPreEvolve, "Prepare cached evolution");
    m.def("AlphaS", &hoppetAlphaS, "Get coupling at scale Q");
    m.def("SetCoupling", &hoppetSetCoupling, "Set coupling parameters");
    m.def("EvalPID", &hoppetEvalPID, "Evaluate PDF for specific PID");
    m.def("SetFFN", &hoppetSetFFN, "Set fixed-flavor number scheme");
    m.def("SetVFN", &hoppetSetVFN, "Set variable-flavor number scheme");
    m.def("SetPoleMassVFN", &hoppetSetPoleMassVFN, "Set pole-mass VFN");
    m.def("SetMSbarMassVFN", &hoppetSetMSbarMassVFN, "Set MSbar-mass VFN");
    m.def("SetExactDGLAP", &hoppetSetExactDGLAP, "Set exact DGLAP");
    m.def("SetApproximateDGLAPN3LO", &hoppetSetApproximateDGLAPN3LO, "Set approximate N3LO DGLAP");
    m.def("SetSplittingNNLO", &hoppetSetSplittingNNLO, "Set NNLO splitting functions");
    m.def("SetSplittingN3LO", &hoppetSetSplittingN3LO, "Set N3LO splitting functions");
    m.def("SetN3LOnfthresholds", &hoppetSetN3LOnfthresholds, "Set N3LO nf thresholds");
    m.def("SetYLnlnQInterpOrders", &hoppetSetYLnlnQInterpOrders, "Set interpolation orders");
    m.def("StartStrFct", &hoppetStartStrFct, "Start structure functions");
    m.def("StartStrFctExtended", &hoppetStartStrFctExtended, "Start structure functions (extended)");
    m.def("WriteLHAPDFGrid", &hoppetWriteLHAPDFGrid, "Write LHAPDF grid");
}