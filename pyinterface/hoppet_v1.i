/* File: hoppet_v1.i */
%module hoppet_v1

%{
#define SWIG_FILE_WITH_INIT
#include "hoppet_v1.h"
#include <iostream>
#include <Python.h>

extern void hoppetAssign(void (* pdf_subroutine)(const double &, 
                                                 const double &, double *));

extern void hoppetEvolve(const double &,
                         const double &,
                         const int    &,
                         const double &,
                         void (* pdf_subroutine)(const double &, 
                                                 const double &, double *),
                         const double &);

extern void hoppetEval(const double &,
                       const double &,
                       double *);
                         
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

// Wrapper function to convert the pdf array to a Python list
static PyObject* pdf_wrapper(double *pdf) {
    PyObject *py_list = PyList_New(13);
    for (int i = 0; i < 13; i++) {
        PyList_SetItem(py_list, i, PyFloat_FromDouble(pdf[i]));
    }
    return py_list;
}

static void hoppetAssign_wrapper(PyObject *callback) {
    if (!PyCallable_Check(callback)) {
        PyErr_SetString(PyExc_TypeError, "Expected a callable object");
        return;
    }
    PyObject_SetAttrString(PyImport_AddModule("__main__"), "pdf_callback", callback);
    hoppetAssign(pdf_subroutine_wrapper);
}

static void hoppetEvolve_wrapper(const double & asQ0, const double & Q0alphas, const int & nloop, 
                                 const double & muR_Q, PyObject *callback, const double & Q0pdf) {
    if (!PyCallable_Check(callback)) {
        PyErr_SetString(PyExc_TypeError, "Expected a callable object");
        return;
    }
    PyObject_SetAttrString(PyImport_AddModule("__main__"), "pdf_callback", callback);
    hoppetEvolve(asQ0, Q0alphas, nloop, muR_Q, pdf_subroutine_wrapper, Q0pdf);
}

static PyObject* hoppetEval_wrapper(const double &x, const double &Q, double *pdf) {
    hoppetEval(x, Q, pdf);
    return pdf_wrapper(pdf);
}

double *new_pdf(){
    return new double[13];
}

%}

%rename(hoppetStart          )          hoppetstart_;
%rename(hoppetStartExtended  )          hoppetstartextended_;
%rename(hoppetAssign         )          hoppetassign_;
%rename(hoppetEvolve         )          hoppetevolve_;
%rename(hoppetPreEvolve      )          hoppetpreevolve_;     
%rename(hoppetCachedEvolve   )          hoppetcachedevolve_;
%rename(hoppetAlphaS         )          hoppetalphas_; 
%rename(hoppetSetFFN         )          hoppetsetffn_;       
%rename(hoppetSetVFN         )          hoppetsetvfn_;       
%rename(hoppetSetPoleMassVFN )          hoppetsetpolemassvfn_;       
%rename(hoppetSetMSbarMassVFN)          hoppetsetmsbarmassvfn_;       
%rename(hoppetSetExactDGLAP  )          hoppetsetexactdglap_;
%rename(hoppetEval           )          hoppeteval_;          
%rename(hoppetEvalSplit      )          hoppetevalsplit_;
%rename(hoppetSetQED         )          hoppetsetqed_;
%rename(hoppetDeleteAll      )          hoppetdeleteall_;

%rename(hoppetStartStrFct        )      hoppetstartstrfct_;
%rename(hoppetStartStrFctExtended)      hoppetstartstrfctextended_;
%rename(hoppetInitStrFct         )      hoppetinitstrfct_;
%rename(hoppetInitStrFctFlav     )      hoppetinitstrfctflav_;
%rename(hoppetStrFct             )      hoppetstrfct_;
%rename(hoppetStrFctNoMu         )      hoppetstrfctnomu_;
%rename(hoppetStrFctLO           )      hoppetstrfctlo_;
%rename(hoppetStrFctNLO          )      hoppetstrfctnlo_;
%rename(hoppetStrFctFlav         )      hoppetstrfctflav_;
%rename(hoppetStrFctNoMuFlav     )      hoppetstrfctnomuflav_;
%rename(hoppetStrFctLOFlav       )      hoppetstrfctloflav_;
%rename(hoppetStrFctNLOFlav      )      hoppetstrfctnloflav_;
%rename(hoppetStrFctNNLO         )      hoppetstrfctnnlo_;
%rename(hoppetStrFctN3LO         )      hoppetstrfctn3lo_;

%include "hoppet_v1.h"

%inline %{
    void hoppetAssign_wrapper(PyObject *callback);
    void hoppetEvolve_wrapper(const double & asQ0, const double & Q0alphas, const int & nloop, 
                              const double & muR_Q, PyObject *callback, const double & Q0pdf);
    PyObject* hoppetEval_wrapper(const double &x, const double &Q, double *pdf);
    double *new_pdf();
%}