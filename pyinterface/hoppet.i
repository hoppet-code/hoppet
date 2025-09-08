/* File: hoppet.i */
%module hoppet
%include "std_string.i"  // This must be included in order to handle c strings

%{
#define SWIG_FILE_WITH_INIT
#include "hoppet.h"
#include <Python.h>

// Check if the callback is a callable object and set it as a global
// variable. This routine is used in the Assign, Evolve and
// CachedEvolve functions, to handle the fact that the pdf_subrpoutine
// is expected to be a pointer. Defined as a macro to avoid
// boilerplate code.
#define CHECK_AND_SET_CALLBACK(callback) \
    if (!PyCallable_Check(callback)) { \
        PyErr_SetString(PyExc_TypeError, "Expected a callable object"); \
        return; \
    } \
    PyObject_SetAttrString(PyImport_AddModule("__main__"), "pdf_callback", callback);

// Check if the global pdf array is initialized. This routine is used
// in the Eval etc. functions to ensure that the global pdf array is
// initialized before calling the hoppetEval function. Defined as a
// macro to avoid boilerplate code.
# define CHECK_GLOBAL_PDF_INITIALIZED \
    if (global_pdf == nullptr) { \
        PyErr_SetString(PyExc_RuntimeError, "Global pdf array is not initialized"); \
        return nullptr; \
    }

const unsigned int pdf_len = 18; // 18 to have space for photon and leptons if needed
double *global_pdf = nullptr;

// Wrapper function to bridge Python and C callback
static void pdf_subroutine_wrapper(const double &x, const double &Q, double *res) {
    PyGILState_STATE gstate = PyGILState_Ensure();  // Ensure GIL for Python calls

    PyObject *callback = PyObject_GetAttrString(PyImport_AddModule("__main__"), "pdf_callback");
    if (callback && PyCallable_Check(callback)) {
        PyObject *arglist = Py_BuildValue("(dd)", x, Q);
        PyObject *result = PyObject_CallObject(callback, arglist);
        Py_DECREF(arglist);

        if (result && PySequence_Check(result)) {
            for (int i = 0; i < 13; i++) {
                PyObject *item = PySequence_GetItem(result, i);
                res[i] = PyFloat_AsDouble(item);
                Py_DECREF(item);
            }
            Py_DECREF(result);
        }
    }

    PyGILState_Release(gstate);  // Release GIL
}

// Initialize the global pdf array
void init_global_pdf() {
    if (global_pdf == nullptr) {
        global_pdf = new double[pdf_len];
    }
}

// Free the global pdf array
void free_global_pdf() {
    if (global_pdf != nullptr) {
        delete[] global_pdf;
        global_pdf = nullptr;
    }
}

// Wrapper function to convert the pdf array to a Python list
static PyObject* pdf_to_array(double *pdf) {
    PyObject *py_list = PyList_New(pdf_len);
    for (unsigned int i = 0; i < pdf_len; i++) {
        PyList_SetItem(py_list, i, PyFloat_FromDouble(pdf[i]));
    }
    return py_list;
}

// To correctly initialize the global pdf array we need to modify the
// hoppetStart and StarExtended. We also modify hoppetDeleteAll so
// that it frees the pdf.
static void Start(const double & dy, const int & nloop){
    init_global_pdf();
    hoppetStart(dy, nloop);
}

static void StartExtended(const double & ymax,   //< highest value of ln1/x user wants to access
                          const double & dy,     //< internal ln1/x grid spacing: 0.1-0.25 is a sensible range
                          const double & Qmin,   //< lower limit of Q range
                          const double & Qmax,   //< upper limit of Q range
                          const double & dlnlnQ, //< internal table spacing in lnlnQ (e.g. dy/4)
                          const int & nloop,     //< the maximum number of loops we'll want (<=3)
                          const int & order,     //< order of numerical interpolation (e.g. -6)
                          const int & factscheme){
    init_global_pdf();
    hoppetStartExtended(ymax, dy, Qmin, Qmax, dlnlnQ, nloop, order, factscheme);
}

static void DeleteAll() {
    free_global_pdf();
    hoppetDeleteAll();
}

static void Assign(PyObject *callback) {
    CHECK_AND_SET_CALLBACK(callback)
    hoppetAssign(pdf_subroutine_wrapper);
}

static void Evolve(const double & asQ0, const double & Q0alphas, const int & nloop, 
                                 const double & muR_Q, PyObject *callback, const double & Q0pdf) {
    CHECK_AND_SET_CALLBACK(callback)
    hoppetEvolve(asQ0, Q0alphas, nloop, muR_Q, pdf_subroutine_wrapper, Q0pdf);
}

static void CachedEvolve(PyObject *callback) {
    CHECK_AND_SET_CALLBACK(callback)
    hoppetCachedEvolve(pdf_subroutine_wrapper);
}

static PyObject* Eval(const double & x, const double & Q) {
    CHECK_GLOBAL_PDF_INITIALIZED
    hoppetEval(x, Q, global_pdf);
    return pdf_to_array(global_pdf);
}

static PyObject* EvalSplit(const double & x, const double & Q, const int & iloop, const int & nf) {
    CHECK_GLOBAL_PDF_INITIALIZED
    hoppetEvalSplit(x, Q, iloop, nf, global_pdf);
    return pdf_to_array(global_pdf);
}

static PyObject* StrFct(const double & x, const double & Q, const double & muR_in, const double & muF_in) {
    CHECK_GLOBAL_PDF_INITIALIZED
    hoppetStrFct(x, Q, muR_in, muF_in, global_pdf);
    return pdf_to_array(global_pdf);
}

static PyObject* StrFctNoMu(const double & x, const double & Q) {
    CHECK_GLOBAL_PDF_INITIALIZED
    hoppetStrFctNoMu(x, Q, global_pdf);
    return pdf_to_array(global_pdf);
}

static PyObject* StrFctLO(const double & x, const double & Q, const double & muR_in, const double & muF_in) {
    CHECK_GLOBAL_PDF_INITIALIZED
    hoppetStrFctLO(x, Q, muR_in, muF_in, global_pdf);
    return pdf_to_array(global_pdf);
}

static PyObject* StrFctNLO(const double & x, const double & Q, const double & muR_in, const double & muF_in) {
    CHECK_GLOBAL_PDF_INITIALIZED
    hoppetStrFctNLO(x, Q, muR_in, muF_in, global_pdf);
    return pdf_to_array(global_pdf);
}

static PyObject* StrFctFlav(const double & x, const double & Q, const double & muR_in, const double & muF_in, const int & flav) {
    CHECK_GLOBAL_PDF_INITIALIZED
    hoppetStrFctFlav(x, Q, muR_in, muF_in, flav, global_pdf);
    return pdf_to_array(global_pdf);
}

static PyObject* StrFctNoMuFlav(const double & x, const double & Q, const int & flav) {
    CHECK_GLOBAL_PDF_INITIALIZED
    hoppetStrFctNoMuFlav(x, Q, flav, global_pdf);
    return pdf_to_array(global_pdf);
}

static PyObject* StrFctLOFlav(const double & x, const double & Q, const double & muR_in, const double & muF_in, const int & flav) {
    CHECK_GLOBAL_PDF_INITIALIZED
    hoppetStrFctLOFlav(x, Q, muR_in, muF_in, flav, global_pdf);
    return pdf_to_array(global_pdf);
}

static PyObject* StrFctNLOFlav(const double & x, const double & Q, const double & muR_in, const double & muF_in, const int & flav) {
    CHECK_GLOBAL_PDF_INITIALIZED
    hoppetStrFctNLOFlav(x, Q, muR_in, muF_in, flav, global_pdf);
    return pdf_to_array(global_pdf);
}

static PyObject* StrFctNNLO(const double & x, const double & Q, const double & muR_in, const double & muF_in) {
    CHECK_GLOBAL_PDF_INITIALIZED
    hoppetStrFctNNLO(x, Q, muR_in, muF_in, global_pdf);
    return pdf_to_array(global_pdf);
}

static PyObject* StrFctN3LO(const double & x, const double & Q, const double & muR_in, const double & muF_in) {
    CHECK_GLOBAL_PDF_INITIALIZED
    hoppetStrFctN3LO(x, Q, muR_in, muF_in, global_pdf);
    return pdf_to_array(global_pdf);
}
%}

// Routines that retain their "hoppet" prefix in this list have an
// explicit functions defined above in order to correctly return a pdf
// object. The rest we rename to remove the trailing underscore,
// remove the hoppet prefix and turn into CamelCase.
%rename(PreEvolve               )      hoppetpreevolve_;     
%rename(AlphaS                  )      hoppetalphas_; 
%rename(SetCoupling             )      hoppetsetcoupling_;       
%rename(SetFFN                  )      hoppetsetffn_;       
%rename(SetVFN                  )      hoppetsetvfn_;       
%rename(SetPoleMassVFN          )      hoppetsetpolemassvfn_;       
%rename(SetMSbarMassVFN         )      hoppetsetmsbarmassvfn_;       
%rename(SetExactDGLAP           )      hoppetsetexactdglap_;
%rename(SetApproximateDGLAPN3LO )      hoppetsetapproximatedglapn3lo_;
%rename(SetSplittingNNLO        )      hoppetsetsplittingnnlo_;
%rename(SetSplittingN3LO        )      hoppetsetsplittingn3lo_;
%rename(SetN3LOnfthresholds     )      hoppetsetn3lonfthresholds_;
%rename(SetYLnlnQInterpOrders   )      hoppetsetylnlnqinterporders_;
%rename(SetQED                  )      hoppetsetqed_;
%rename(StartStrFct             )      hoppetstartstrfct_;
%rename(StartStrFctExtended     )      hoppetstartstrfctextended_;
%rename(InitStrFct              )      hoppetinitstrfct_;
%rename(InitStrFctFlav          )      hoppetinitstrfctflav_;
%rename(WriteLHAPDFGrid         )      hoppetWriteLHAPDFGrid;
%rename(hoppetWriteLHAPDFgrid   )      hoppetwritelhapdfgrid_;

// These are the functions that have an explicit interface defined
// above, so we make sure not to include the C++ versions
%ignore hoppetstart_;
%ignore hoppetstartextended_;
%ignore hoppetassign_; // The callback function is Assign
%ignore hoppetevolve_; // The callback function is Evolve
%ignore hoppetcachedevolve_; // The callback function is CachedEvolve
%ignore hoppeteval_;          
%ignore hoppetevalsplit_;
%ignore hoppetdeleteall_;
%ignore hoppetstrfct_;
%ignore hoppetstrfctnomu_;
%ignore hoppetstrfctlo_;
%ignore hoppetstrfctnlo_;
%ignore hoppetstrfctflav_;
%ignore hoppetstrfctnomuflav_;
%ignore hoppetstrfctloflav_;
%ignore hoppetstrfctnloflav_;
%ignore hoppetstrfctnnlo_;
%ignore hoppetstrfctn3lo_;
%ignore hoppetwritelhapdfwithlen_;


%include "hoppet.h"

%inline %{
    void init_global_pdf();
    void free_global_pdf();
    PyObject* pdf_to_array(double *pdf);
    void Start(const double & dy, const int & nloop);
    void StartExtended(const double & ymax, const double & dy, const double & Qmin, 
                        const double & Qmax, const double & dlnlnQ, const int & nloop, 
                        const int & order, const int & factscheme);
    void DeleteAll();
    void Assign(PyObject *callback);
    void Evolve(const double & asQ0, const double & Q0alphas, const int & nloop, 
                              const double & muR_Q, PyObject *callback, const double & Q0pdf);
    void CachedEvolve(PyObject *callback);
    PyObject* Eval(const double & x, const double & Q);
    PyObject* EvalSplit(const double & x, const double & Q, const int & iloop, const int & nf);
    PyObject* StrFct(const double & x, const double & Q, const double & muR_in, const double & muF_in);
    PyObject* StrFctNoMu(const double & x, const double & Q);
    PyObject* StrFctLO(const double & x, const double & Q, const double & muR_in, const double & muF_in);
    PyObject* StrFctNLO(const double & x, const double & Q, const double & muR_in, const double & muF_in);
    PyObject* StrFctFlav(const double & x, const double & Q, const double & muR_in, const double & muF_in, const int & flav);
    PyObject* StrFctNoMuFlav(const double & x, const double & Q, const int & flav);
    PyObject* StrFctLOFlav(const double & x, const double & Q, const double & muR_in, const double & muF_in, const int & flav);
    PyObject* StrFctNLOFlav(const double & x, const double & Q, const double & muR_in, const double & muF_in, const int & flav);
    PyObject* StrFctNNLO(const double & x, const double & Q, const double & muR_in, const double & muF_in);
    PyObject* StrFctN3LO(const double & x, const double & Q, const double & muR_in, const double & muF_in);
%}
