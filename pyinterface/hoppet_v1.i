/* File: hoppet_v1.i */
%module hoppet_v1

%{
#define SWIG_FILE_WITH_INIT
#include "hoppet_v1.h"
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