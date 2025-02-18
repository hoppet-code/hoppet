/* File: hoppet_v1.i */
%module hoppet_v1

%{
#define SWIG_FILE_WITH_INIT
#include "hoppet_v1.h"
#include <iostream>
#include <Python.h>

extern void Assign(void (* pdf_subroutine)(const double &, 
                                                 const double &, double *));

extern void Evolve(const double &,
                         const double &,
                         const int    &,
                         const double &,
                         void (* pdf_subroutine)(const double &, 
                                                 const double &, double *),
                         const double &);

extern void CachedEvolve(void (* pdf_subroutine)(const double &, 
                                                       const double &, double *));                         

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

static void Assign(PyObject *callback) {
    if (!PyCallable_Check(callback)) {
        PyErr_SetString(PyExc_TypeError, "Expected a callable object");
        return;
    }
    PyObject_SetAttrString(PyImport_AddModule("__main__"), "pdf_callback", callback);
    hoppetAssign(pdf_subroutine_wrapper);
}

static void Evolve(const double & asQ0, const double & Q0alphas, const int & nloop, 
                                 const double & muR_Q, PyObject *callback, const double & Q0pdf) {
    if (!PyCallable_Check(callback)) {
        PyErr_SetString(PyExc_TypeError, "Expected a callable object");
        return;
    }
    PyObject_SetAttrString(PyImport_AddModule("__main__"), "pdf_callback", callback);
    hoppetEvolve(asQ0, Q0alphas, nloop, muR_Q, pdf_subroutine_wrapper, Q0pdf);
}

static void CachedEvolve(PyObject *callback) {
    if (!PyCallable_Check(callback)) {
        PyErr_SetString(PyExc_TypeError, "Expected a callable object");
        return;
    }
    PyObject_SetAttrString(PyImport_AddModule("__main__"), "pdf_callback", callback);
    hoppetCachedEvolve(pdf_subroutine_wrapper);
}

// This function simply creates a pointer to a new array of doubles as neede by hoppetEval etc.
double *new_pdf(){
    return new double[13];
}

// Wrapper function to convert the pdf array to a Python list
static PyObject* pdf_to_array(double *pdf) {
    PyObject *py_list = PyList_New(13);
    for (int i = 0; i < 13; i++) {
        PyList_SetItem(py_list, i, PyFloat_FromDouble(pdf[i]));
    }
    return py_list;
}

%}

%rename(Start          )          hoppetstart_;
%rename(StartExtended  )          hoppetstartextended_;
%rename(hoppetAssign         )          hoppetassign_; // The callback function is Assign
%rename(hoppetEvolve         )          hoppetevolve_; // The callback function is Evolve
%rename(PreEvolve      )          hoppetpreevolve_;     
%rename(hoppetCachedEvolve   )          hoppetcachedevolve_; // The callback function is CachedEvolve
%rename(AlphaS         )          hoppetalphas_; 
%rename(SetFFN         )          hoppetsetffn_;       
%rename(SetVFN         )          hoppetsetvfn_;       
%rename(SetPoleMassVFN )          hoppetsetpolemassvfn_;       
%rename(SetMSbarMassVFN)          hoppetsetmsbarmassvfn_;       
%rename(SetExactDGLAP  )          hoppetsetexactdglap_;
%rename(Eval           )          hoppeteval_;          
%rename(EvalSplit      )          hoppetevalsplit_;
%rename(SetQED         )          hoppetsetqed_;
%rename(DeleteAll      )          hoppetdeleteall_;

%rename(StartStrFct        )      hoppetstartstrfct_;
%rename(StartStrFctExtended)      hoppetstartstrfctextended_;
%rename(InitStrFct         )      hoppetinitstrfct_;
%rename(InitStrFctFlav     )      hoppetinitstrfctflav_;
%rename(StrFct             )      hoppetstrfct_;
%rename(StrFctNoMu         )      hoppetstrfctnomu_;
%rename(StrFctLO           )      hoppetstrfctlo_;
%rename(StrFctNLO          )      hoppetstrfctnlo_;
%rename(StrFctFlav         )      hoppetstrfctflav_;
%rename(StrFctNoMuFlav     )      hoppetstrfctnomuflav_;
%rename(StrFctLOFlav       )      hoppetstrfctloflav_;
%rename(StrFctNLOFlav      )      hoppetstrfctnloflav_;
%rename(StrFctNNLO         )      hoppetstrfctnnlo_;
%rename(StrFctN3LO         )      hoppetstrfctn3lo_;

%include "hoppet_v1.h"

%inline %{
    void Assign(PyObject *callback);
    void Evolve(const double & asQ0, const double & Q0alphas, const int & nloop, 
                              const double & muR_Q, PyObject *callback, const double & Q0pdf);
    void CachedEvolve(PyObject *callback);
    double *new_pdf();
    PyObject* pdf_to_array(double *pdf);
%}