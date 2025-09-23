// -*- C++ -*-
// C++ include file for hoppet's vanilla interface
#ifndef __HOPPET__
#define __HOPPET__

// define nicer forms of standard f77 naming
#define hoppetStart                    hoppetstart_
#define hoppetStartExtended            hoppetstartextended_
#define hoppetAssign                   hoppetassign_
#define hoppetSetCoupling              hoppetsetcoupling_
#define hoppetEvolve                   hoppetevolve_        
#define hoppetPreEvolve                hoppetpreevolve_     
#define hoppetCachedEvolve             hoppetcachedevolve_
#define hoppetAlphaS                   hoppetalphas_ 
#define hoppetSetFFN                   hoppetsetffn_       
#define hoppetSetVFN                   hoppetsetvfn_       
#define hoppetSetPoleMassVFN           hoppetsetpolemassvfn_       
#define hoppetSetMSbarMassVFN          hoppetsetmsbarmassvfn_       
#define hoppetSetExactDGLAP            hoppetsetexactdglap_
#define hoppetSetApproximateDGLAPN3LO  hoppetsetapproximatedglapn3lo_
#define hoppetSetSplittingNNLO         hoppetsetsplittingnnlo_
#define hoppetSetSplittingN3LO         hoppetsetsplittingn3lo_
#define hoppetSetN3LOnfthresholds      hoppetsetn3lonfthresholds_
#define hoppetSetYLnlnQInterpOrders    hoppetsetylnlnqinterporders_
#define hoppetEval                     hoppeteval_          
#define hoppetEvalPID                  hoppetevalpid_          
#define hoppetEvalSplit                hoppetevalsplit_
#define hoppetSetQED                   hoppetsetqed_
#define hoppetDeleteAll                hoppetdeleteall_
#define hoppetBenchmarkPDFunpol        hoppetbenchmarkpdfunpol_

/// The fortran subroutines pertaining to the below structure function 
/// interfaces can be found at the end of structure_functions.f90
#define hoppetStartStrFct              hoppetstartstrfct_
#define hoppetStartStrFctExtended      hoppetstartstrfctextended_
#define hoppetInitStrFct               hoppetinitstrfct_
#define hoppetInitStrFctFlav           hoppetinitstrfctflav_
#define hoppetStrFct                   hoppetstrfct_
#define hoppetStrFctNoMu               hoppetstrfctnomu_
#define hoppetStrFctLO                 hoppetstrfctlo_
#define hoppetStrFctNLO                hoppetstrfctnlo_
#define hoppetStrFctFlav               hoppetstrfctflav_
#define hoppetStrFctNoMuFlav           hoppetstrfctnomuflav_
#define hoppetStrFctLOFlav             hoppetstrfctloflav_
#define hoppetStrFctNLOFlav            hoppetstrfctnloflav_
#define hoppetStrFctNNLO               hoppetstrfctnnlo_
#define hoppetStrFctN3LO               hoppetstrfctn3lo_

#include <string>

namespace hoppet {
  /// indices for the different structure functions
  const int iF1Wp = 1+6; //< F1 W+ : D + Ubar                                                       
  const int iF2Wp = 2+6; //< F2 W+ : D + Ubar                                                      
  const int iF3Wp = 3+6; //< F3 W+ : D + Ubar                                                      
  const int iF1Wm =-1+6; //< F1 W- : Dbar + U                                                      
  const int iF2Wm =-2+6; //< F2 W- : Dbar + U                                                      
  const int iF3Wm =-3+6; //< F3 W- : Dbar + U                                                      
  const int iF1Z  = 4+6; //< F1 Z  : (D + Dbar) * v_i^2a_i^2_down + (U + Ubar) * v_i^2a_i^2_up     
  const int iF2Z  = 5+6; //< F2 Z  : (D + Dbar) * v_i^2a_i^2_down + (U + Ubar) * v_i^2a_i^2_up     
  const int iF3Z  = 6+6; //< F3 Z  : (D + Dbar) * 2v_ia_i_down + (U + Ubar) * 2v_ia_i_up           
  const int iF1EM =-4+6; //< F1 γ  : (D + Dbar) * e2_down + (U + Ubar) * e2_up                     
  const int iF2EM =-5+6; //< F2 γ  : (D + Dbar) * e2_down + (U + Ubar) * e2_up                     
  const int iF1gZ = 0+6; //< F1 γZ : (D + Dbar) * e_down * 2v_i_down + (U + Ubar) * e_up * 2v_i_up 
  const int iF2gZ =-6+6; //< F2 γZ : (D + Dbar) * e_down * 2v_i_down + (U + Ubar) * e_up * 2v_i_up
  const int iF3gZ = 7+6; //< F3 γZ : (D + Dbar) * e_down * 2a_i_down + (U + Ubar) * e_up * 2a_i_up

  const int scale_choice_fixed     = 0; ///< muR,muF scales
					///predetermined in the
					///hoppetStartStrFct call
  const int scale_choice_Q         = 1; ///< muR,muF scales equal to
					///xR,xQ * Q (xR,xQ to be set
					///in hoppetStartStrFct)
  const int scale_choice_arbitrary = 2; ///< muR,muF scales can be
					///chosen freely in the
					///hoppetStrFctLO (etc.) and
					///hoppetStrFct calls
  const int nnlo_splitting_exact = -2;
  const int nnlo_splitting_param = -1;
  // these three should keep their numerical values because of a
  // correspondence with vogt imod values
  const int nnlo_splitting_Nfitav   =  0;
  const int nnlo_splitting_Nfiterr1 =  1;
  const int nnlo_splitting_Nfiterr2 =  2;

  const int n3lo_splitting_exact = -2;
  const int n3lo_splitting_param = -1;
  // these three should keep their numerical values because of a
  // correspondence with vogt imod values
  const int n3lo_splitting_Nfitav   =  0;
  const int n3lo_splitting_Nfiterr1 =  1;
  const int n3lo_splitting_Nfiterr2 =  2;
  // As of 2024-04-16 there are several approximations available of
  // the n3lo splitting functions. To maintain some backwards
  // compatibility we have a switch below that allows the user to pick
  // which set of approximations to use. The have progressively more
  // moments.
  const int n3lo_splitting_approximation_up_to_2310_05744 = 100;
  //< Uses non-singlet of 1610.07477+1707.08315, pure-singlet (qq) of
  //2302.07593, qg of 2307.04158 and gq and gg of 2310.05744
  const int n3lo_splitting_approximation_up_to_2404_09701
            = 101; //< Replaces gq with that of 2404.09701
  const int n3lo_splitting_approximation_up_to_2410_08089
            = 102; //< Additionally replaces gg with that of 2410.08089
  
  const int nnlo_nfthreshold_exact = -12;
  const int nnlo_nfthreshold_param = -11;

  const int n3lo_nfthreshold_on  = 1;
  const int n3lo_nfthreshold_off = 0;

  const int factscheme_MSbar    = 1;
  const int factscheme_DIS      = 2;
  const int factscheme_PolMSbar = 3;
  const int factscheme_FragMSbar = 4;

}

extern "C" {


  // The following does not work
  // void hoppetSetQED(const bool & withqed, const bool & qcdqed, const bool & plq);
  // There must be a problem in passing boolean from c++ to fortran. Passing them
  // as integer works, but is it safe? (TBD)
  void hoppetSetQED(const int & withqed, const int & qcdqed, const int & plq);

  /// initialise the underlying grid, splitting functions and pdf-table
  /// objects, using the dy grid spacing and splitting functions up to
  /// nloop loops; all other parameters are set automatically
  void hoppetStart(const double & dy, const int & nloop);

  /// an extended interface for starting hoppet
  void hoppetStartExtended(
       const double & ymax,   //< highest value of ln1/x user wants to access
       const double & dy,     //< internal ln1/x grid spacing: 0.1-0.25 is a sensible range
       const double & Qmin,   //< lower limit of Q range
       const double & Qmax,   //< upper limit of Q range
       const double & dlnlnQ, //< internal table spacing in lnlnQ (e.g. dy/4)
       const int & nloop,     //< the maximum number of loops we'll want (<=3)
       const int & order,     //< order of numerical interpolation (e.g. -6)
       const int & factscheme //< one of the factschemes defined above
       );


  /// Set things up to be a fixed-flavour number scheme with the given
  /// fixed_nf number of flavours
  void hoppetSetFFN(const int & fixed_nf);


  /// Set things up to be a variable-flavour number scheme with the given
  /// quark (pole) masses. Now deprecated; use hoppetSetPoleMassVFN instead
  void  hoppetSetVFN(const double &mc, const double & mb, const double & mt);

  /// Set things up to be a variable-flavour number scheme with the
  /// given quark (pole) masses. Thresholds are crossed at the pole
  /// masses, both for the coupling and the PDF evolution.
  void  hoppetSetPoleMassVFN(const double &mc, const double & mb, const double & mt);

  /// Set things up to be a variable-flavour number scheme with the given
  /// quark (MSbar) masses. Thresholds are crossed at the MSbar
  /// masses, both for the coupling and the PDF evolution.
  void  hoppetSetMSbarMassVFN(const double &mc, const double & mb, const double & mt);

  
  /// Arrange for the use of exact NNLO splitting and mass-threshold
  /// functions.
  void hoppetSetExactDGLAP(const int & exact_nfthreshold, const int & exact_splitting);

  /// Arrange for the use of various approximate N3LO splitting functions.
  void hoppetSetApproximateDGLAPN3LO(const int & splitting_variant);

  /// Arrange for the use of various NNLO splitting functions.
  void hoppetSetSplittingNNLO(const int & splitting_variant);

  /// Arrange for the use of various N3LO splitting functions.
  void hoppetSetSplittingN3LO(const int & splitting_variant);

  ///  Arrange for the use of N3LO mass thresholds or not.
  void hoppetSetN3LOnfthresholds(const int & variant);

  ///  Override the default interpolation order in y and lnlnQ.
  void hoppetSetYLnlnQInterpOrders(const int & yorder, const int & lnlnQorder);

  /// Given a pdf_subroutine with the interface shown below, initialise
  /// our internal pdf table.
  void hoppetAssign(void (* pdf_subroutine)(const double & x, 
                                            const double & Q, double * res) );

//=====================================================================
// Set up the strong coupling such that alphas(Q)=alphas_Q, with the
// given number of loops (nloop).
//
// The user should have set the quark masses or requested a FFN scheme
// prior to calling this function.
//
// This function is provided mainly for use in conjunction with hoppetAssign.
// In particular, it has the side effect of modifying the structure of the
// PDF tables to make sure they know about the mass thresholds.
//
// If QED has been requested, a QED coupling will also be set up
// (its value is not currently configurable from this interface).
//
// If you call hoppetEvolve (below), there is no need to separately call
// hoppetSetCoupling.
  void hoppetSetCoupling(const double & asQ0,
                         const double & Q0alphas,
                         const int    & nloop);


  /// Given a pdf_subroutine with the interface shown below, fill the 
  /// table by evolving the PDF from scale Q0pdf, with alphas provided 
  /// at scale Q0alphas
  void hoppetEvolve(const double & asQ0,
                    const double & Q0alphas,
                    const int    & nloop,
                    const double & muR_Q,
                    void (* pdf_subroutine)(const double & x, 
                                            const double & Q, double * res),
                    const double & Q0pdf);


  /// Prepare a cached evolution
  void hoppetPreEvolve(const double & asQ0, 
                       const double & Q0alphas, 
                       const int    & nloop, 
                       const double & muR_Q, 
                       const double & Q0pdf);


  /// Carry out a cached evolution based on the initial condition
  /// that can be obtained from pdf_subroutine at the scale Q0pdf set in
  /// hoppetPreEvolve
  void hoppetCachedEvolve(void (*pdf_subroutine)(const double & x, 
                                     const double & Q, double * res));

  /// Return the coupling at scale Q
  double hoppetAlphaS(const double & Q);

  /// Return in f[0..12+5] the value of the internally stored pdf at
  /// the given x,Q, with the usual LHApdf meanings for the indices
  /// -6:6 + photon and lepton.
  void hoppetEval(const double & x,
                  const double & Q,
                  double * f);


  /// Return in f the value of the internally stored pdf at the given
  /// x,Q, and pid with the usual LHApdf meanings for the indices
  /// -6:6.
  double hoppetEvalPID(const double & x,
		       const double & Q,
		       const int & pid);


  /// Return in f[0..12] the value of 
  ///
  ///    [P(iloop,nf) \otimes pdf] (x,Q)
  ///
  /// where P(iloop,nf) is the iloop-splitting function for the given
  /// value of nf, and pdf is our internally stored pdf.
  ///
  /// The normalisation is such that the nloop dglap evolution equation is
  ///
  ///     dpdf/dlnQ^2 = sum_{iloop=1}^nloop 
  ///                        (alphas/(2*pi))^iloop * P(iloop,nf) \otimes pdf
  ///
  /// Note that each time nf changes relative to a previous call for the
  /// same iloop, the convolution has to be repeated for the whole
  /// table. So for efficient results when requiring multiple nf values,
  /// calls with the same nf value should be grouped together.
  ///
  /// In particular, for repeated calls with the same value of nf, the
  /// convolutions are carried out only on the first call (i.e. once for
  /// each value of iloop). Multiple calls with different values for
  /// iloop can be carried out without problems.
  ///
  void hoppetEvalSplit(const double & x,
                       const double & Q,
                       const int    & iloop,
                       const int    & nf,
                       double * f);

  //---------------------------------------------------------------------
  /// Delete all hoppet objects associated with the streamlined interface,
  /// including grid definitions, PDF tables, couplings, etc.
  ///
  /// NB: this does not delete anything associated with the structure function
  /// part of the interface
  void hoppetDeleteAll();
  
  ///----------------------------------------------------------------------
  /// Setup of constants and parameters needed for structure functions
  void hoppetStartStrFct(const int & order_max);
			 
  
  ///----------------------------------------------------------------------
  /// Setup of constants and parameters needed for structure functions
  void hoppetStartStrFctExtended(const int & order_max,
				 const int & nflav,
				 const int & scale_choice,
				 const double & constant_mu,
				 const int & param_coefs,
				 const double & wmass,
				 const double & zmass);
  
  /// Initialize the structure functions up to specified order
  /// this requires the PDF to have been set up beforehand, and filled in tables(0)
  void hoppetInitStrFct(const int & order_max,
			const int & separate_orders,
			const double & xR,
			const double & xF);

  /// Initialize the structure functions up to specified order
  /// this requires the PDF to have been set up beforehand, and filled in tables(0)
  void hoppetInitStrFctFlav(const int & order_max,
			    const int & separate_orders,
			    const double & xR,
			    const double & xF,
			    const int & flavour_decomposition);


  /// calculate the structure function at x, Q, muR, muF.
  /// This is the sum over all orders. 
  /// The result is placed in F, which must be an array of size at least 14
  /// with the indices as defined above (iF1Wp, etc.)
  void hoppetStrFct(const double & x,
		    const double & Q,
		    const double & muR_in,
		    const double & muF_in,
		    double * F);

  /// calculate the structure function at x, Q, with muR and muF as
  /// requested in hoppetStartStrFct. This is the sum over all orders
  void hoppetStrFctNoMu(const double & x,
		    const double & Q,
 		    double * F);

  /// F_LO
  /// calculate the leading order structure function at x, muF
  ///
  void hoppetStrFctLO(const double & x,
		       const double & Q,
		       const double & muR_in,
		       const double & muF_in,
		       double * F);
  /// F_NLO
  /// calculate the next-to-leading order structure function at x, muF
  ///
  void hoppetStrFctNLO(const double & x,
			const double & Q,
			const double & muR_in,
			const double & muF_in,
			double * F);
  /// calculate the structure function at x, Q, muR, muF.
  /// This is the sum over all orders. 
  /// The result is placed in F, which must be an array of size at least 14
  /// with the indices as defined above (iF1Wp, etc.)
  void hoppetStrFctFlav(const double & x,
			const double & Q,
			const double & muR_in,
			const double & muF_in,
			const double & iflav,
			double * F);

  /// calculate the structure function at x, Q, with muR and muF as
  /// requested in hoppetStartStrFct. This is the sum over all orders
  void hoppetStrFctNoMuFlav(const double & x,
			    const double & Q,
			    const double & iflav,
			    double * F);

  /// F_LO
  /// calculate the leading order structure function at x, muF
  ///
  void hoppetStrFctLOFlav(const double & x,
			  const double & Q,
			  const double & muR_in,
			  const double & muF_in,
			  const double & iflav,
			  double * F);
  /// F_NLO
  /// calculate the next-to-leading order structure function at x, muF
  ///
  void hoppetStrFctNLOFlav(const double & x,
			   const double & Q,
			   const double & muR_in,
			   const double & muF_in,
			   const double & iflav,
			   double * F);
  /// F_NNLO
  /// calculate the next-to-next-to-leading order structure function at x, muF
  ///
  void hoppetStrFctNNLO(const double & x,
			const double & Q,
			const double & muR_in,
			const double & muF_in,
			double * F);
  /// F_N3LO
  /// calculate the next-to-next-to-next-to-leading order structure function at x, muF
  ///
  void hoppetStrFctN3LO(const double & x,
			const double & Q,
			const double & muR_in,
			const double & muF_in,
			double * F);
  
  void hoppetwritelhapdfwithlen_(int & basename_len,
				 const char* basename,
				 const int & pdf_index);
  
  /// Write out the contents of tables(0) (assumed to be the PDF) in the
  /// LHAPDF format
  void hoppetWriteLHAPDFGrid(const std::string& basename, const int& pdf_index) {
    // Convert std::string to C-style string and get length
    const char* basename_cstr = basename.c_str();
    int basename_len = basename.length();
    
    // Call Fortran subroutine (pass integers by reference)
    hoppetwritelhapdfwithlen_(basename_len, basename_cstr, pdf_index);
  }

  void hoppetBenchmarkPDFunpol(const double & x,
			       const double &Q,
			       double * xpdf);

}

#endif // __HOPPET__
