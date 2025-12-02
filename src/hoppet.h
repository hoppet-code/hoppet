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
#define hoppetAlphaQED                 hoppetalphaqed_ 
#define hoppetSetFFN                   hoppetsetffn_       
#define hoppetSetVFN                   hoppetsetvfn_       
#define hoppetSetPoleMassVFN           hoppetsetpolemassvfn_       
#define hoppetSetMSbarMassVFN          hoppetsetmsbarmassvfn_       
#define hoppetSetExactDGLAP            hoppetSetExactDGLAP_c
#define hoppetSetApproximateDGLAPN3LO  hoppetsetapproximatedglapn3lo_
#define hoppetSetSplittingNNLO         hoppetsetsplittingnnlo_
#define hoppetSetSplittingN3LO         hoppetsetsplittingn3lo_
#define hoppetSetN3LOnfthresholds      hoppetsetn3lonfthresholds_
#define hoppetSetYLnlnQInterpOrders    hoppetsetylnlnqinterporders_
#define hoppetEval                     hoppeteval_          
#define hoppetEvalFortranIFlv          hoppetevaliflv_          
#define hoppetEvalSplit                hoppetevalsplit_
#define hoppetSetQED                   hoppetSetQED_c
#define hoppetDeleteAll                hoppetdeleteall_
#define hoppetBenchmarkPDFunpol        hoppetbenchmarkpdfunpol_

/// The fortran subroutines pertaining to the below structure function 
/// interfaces can be found at the end of structure_functions.f90
#define hoppetStartStrFct              hoppetstartstrfct_
#define hoppetStartStrFctExtended      hoppetStartStrFctExtended_c
#define hoppetInitStrFct               hoppetInitStrFct_c
#define hoppetInitStrFctFlav           hoppetInitStrFctFlav_c
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
#include <iostream>
#include <stdbool.h>

namespace hoppet {
  /// @defgroup HOPPET_CC_CONSTS Constants used when integrating splitting functions
  /// @{
  /// @brief Constants indicating which part of a splitting function to return.
  ///
  /// When converting a splitting function to a grid_conv object, the function
  /// being integrated must return one specific part of the splitting function.
  /// These constants label the different parts
  constexpr int cc_REAL=1;      ///< the sum of regular + plus parts
  constexpr int cc_VIRT=2;      ///< minus the "plus"
  constexpr int cc_REALVIRT=3;  ///< just the regular part
  constexpr int cc_DELTA=4;     ///< the delta function part
  /// @}

  /// @defgroup HOPPET_IFLV Flavour indices
  /// @{
  /// These are the indices used to address the different parton flavours
  /// in HOPPET. Note that they are shifted by +6 compared to the Fortran
  /// numbering, i.e. they start from zero
  constexpr int iflv_g = 0+6;
  constexpr int iflv_d = 1+6, iflv_u = 2+6, iflv_s = 3+6, iflv_c = 4+6;
  constexpr int iflv_b = 5+6, iflv_t = 6+6;
  constexpr int iflv_dbar = -1+6;
  constexpr int iflv_ubar = -2+6;
  constexpr int iflv_sbar = -3+6;
  constexpr int iflv_cbar = -4+6;
  constexpr int iflv_bbar = -5+6;
  constexpr int iflv_tbar = -6+6;
  constexpr int iflv_photon = 8+6;
  constexpr int iflv_electron = 9+6; ///< this is the sum of e- and e+
  constexpr int iflv_muon = 10+6;    ///< this is the sum of mu- and mu+
  constexpr int iflv_tau = 11+6;     ///< this is the sum of tau- and tau+
  constexpr int iflv_min = iflv_tbar;     ///< lowest actual flavour index for QCD-only evolution
  constexpr int iflv_max = iflv_t;        ///< highest actual flavour index for QCD-only evolution
  constexpr int iflv_info = iflv_max + 1; ///< index for storing extra representation info
  constexpr int ncompmin = iflv_min;  ///< lowest stored index for QCD-only evolution
  constexpr int ncompmax = iflv_info; ///< highest stored index for QCD-only evolution
  /// the value of iflv_min in fortran, useful for calculating offsets in some cases
  constexpr int iflv_min_fortran = -6;
  /// @}

  /// @defgroup HOPPET_STRFCT_CONSTS Structure function indices
  /// @{
  /// @brief Indices of the different structure functions
  const int iF1Wp = 1+6; ///< F1 W+ : D + Ubar                                                       
  const int iF2Wp = 2+6; ///< F2 W+ : D + Ubar                                                      
  const int iF3Wp = 3+6; ///< F3 W+ : D + Ubar                                                      
  const int iF1Wm =-1+6; ///< F1 W- : Dbar + U                                                      
  const int iF2Wm =-2+6; ///< F2 W- : Dbar + U                                                      
  const int iF3Wm =-3+6; ///< F3 W- : Dbar + U                                                      
  const int iF1Z  = 4+6; ///< F1 Z  : (D + Dbar) * v_i^2a_i^2_down + (U + Ubar) * v_i^2a_i^2_up     
  const int iF2Z  = 5+6; ///< F2 Z  : (D + Dbar) * v_i^2a_i^2_down + (U + Ubar) * v_i^2a_i^2_up     
  const int iF3Z  = 6+6; ///< F3 Z  : (D + Dbar) * 2v_ia_i_down + (U + Ubar) * 2v_ia_i_up           
  const int iF1EM =-4+6; ///< F1 γ  : (D + Dbar) * e2_down + (U + Ubar) * e2_up                     
  const int iF2EM =-5+6; ///< F2 γ  : (D + Dbar) * e2_down + (U + Ubar) * e2_up                     
  const int iF1gZ = 0+6; ///< F1 γZ : (D + Dbar) * e_down * 2v_i_down + (U + Ubar) * e_up * 2v_i_up 
  const int iF2gZ =-6+6; ///< F2 γZ : (D + Dbar) * e_down * 2v_i_down + (U + Ubar) * e_up * 2v_i_up
  const int iF3gZ = 7+6; ///< F3 γZ : (D + Dbar) * e_down * 2a_i_down + (U + Ubar) * e_up * 2a_i_up  
  /// @}

  /// @defgroup HOPPET_SFSCALE_CONSTS Constants that control scale choices in structure functions
  /// @{
  /// @brief Constants that control the choice of \f$\mu_R\f$ and \f$\mu_F\f$ scales in structure functions
  
  /// muR,muF scales predetermined in the hoppetStartStrFct call
  constexpr int scale_choice_fixed     = 0; 
  /// muR,muF scales equal to xR,xQ * Q (xR,xQ to be set in hoppetStartStrFct)
  constexpr int scale_choice_Q         = 1; 
  /// muR,muF scales can be chosen freely in the hoppetStrFctLO (etc.)
  /// and hoppetStrFct calls
  constexpr int scale_choice_arbitrary = 2; 
  /// @}

  constexpr int nnlo_splitting_exact = -2;
  constexpr int nnlo_splitting_param = -1;

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

  constexpr int n3lo_nfthreshold_libOME        = 1;
  constexpr int n3lo_nfthreshold_exact_fortran = 2;
  constexpr int n3lo_nfthreshold_off = 0;
  constexpr int n3lo_nfthreshold_on = n3lo_nfthreshold_libOME;

  const int factscheme_MSbar    = 1;
  const int factscheme_DIS      = 2;
  const int factscheme_PolMSbar = 3;
  const int factscheme_FragMSbar = 4;

}

extern "C" {

  /// hoppet's global variable indicating which part of the
  /// splitting function is being requested when converting to 
  /// a grid_conv object
  extern int hoppet_global_cc_piece;

  /// Set the configuration of QED evolution (off by default)
  ///
  /// @param withqed If true, include QED evolution alongside QCD
  /// @param qcdqed If true, include mixed QCD+QED terms , \f$O(\alpha_s \alpha)\f$
  /// @param plq If true, in the \f$O(\alpha^2)\f$ \f$P_{\ell q}\f$ QED term
  ///
  /// This is to be called before hoppetStart
  void hoppetSetQED(const bool & withqed, const bool & qcdqed, const bool & plq);

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
  ///
  /// @param fixed_nf Value of fixed number of flavours (e.g. 5)
  void hoppetSetFFN(const int & fixed_nf);


  /// Set things up to be a variable-flavour number scheme with the
  /// given quark (pole) masses. Now deprecated; use
  /// hoppetSetPoleMassVFN instead, which is what is being called
  /// undder the hood.
  ///
  /// @param mc Charm mass
  /// @param mb Bottom mass
  /// @param mt Top mass
  void  hoppetSetVFN(const double &mc, const double & mb, const double & mt);

  /// Set things up to be a variable-flavour number scheme with the
  /// given quark (pole) masses. Thresholds are crossed at the pole
  /// masses, both for the coupling and the PDF evolution.
  ///
  /// @param mc Charm mass
  /// @param mb Bottom mass
  /// @param mt Top mass
  void  hoppetSetPoleMassVFN(const double &mc, const double & mb, const double & mt);

  /// Set things up to be a variable-flavour number scheme with the given
  /// quark (MSbar) masses. Thresholds are crossed at the MSbar
  /// masses, both for the coupling and the PDF evolution.
  ///
  /// @param mc Charm mass
  /// @param mb Bottom mass
  /// @param mt Top mass
  void  hoppetSetMSbarMassVFN(const double &mc, const double & mb, const double & mt);

  
  /// Arrange for the use of exact NNLO splitting and mass-threshold
  /// functions.
  ///
  /// @param exact_nfthreshold If True use the exact NNLO mass
  /// threshold functions. If False use the faster parametrisations
  /// @param exact_splitting If True use the exact NNLO splitting
  /// functions. If False use the faster parametrisations
  ///
  /// To be called before hoppetStart
  void hoppetSetExactDGLAP(const bool & exact_nfthreshold, const bool & exact_splitting);

  /// Arrange for the use of various approximate N3LO splitting functions.
  ///
  /// @param splitting_variant One of the hoppet.n3lo_splitting_approximation_* options.
  ///
  /// To be called before hoppetStart
  void hoppetSetApproximateDGLAPN3LO(const int & splitting_variant);

  /// Arrange for the use of various NNLO splitting functions.
  ///
  /// @param splitting_variant One of the hoppet.nnlo_splitting_* options.
  ///
  /// To be called before hoppetStart
  void hoppetSetSplittingNNLO(const int & splitting_variant);

  /// Change the variant for the N3LO splitting functions.
  ///
  /// @param splitting_variant One of the hoppet.n3lo_splitting_* options. 
  ///
  /// To be called before hoppetStart
  void hoppetSetSplittingN3LO(const int & splitting_variant);

  ///  Arrange for the use of N3LO mass thresholds or not.
  ///
  /// @param variant One of the hoppet.n3lo_nfthresholds_* options.
  ///
  /// To be called before hoppetStart
  void hoppetSetN3LOnfthresholds(const int & variant);

  ///  Override the default interpolation order in y and lnlnQ.
  ///
  /// @param yorder The interpolation order in y (default: 5)
  /// @param lnlnQorder The interpolation order in lnlnQ (default: 4)
  ///
  /// 2 corresponds to quadratic, 3 to cubic etc.
  void hoppetSetYLnlnQInterpOrders(const int & yorder, const int & lnlnQorder);

  /// Given a pdf_subroutine with the interface shown below, initialise
  /// our internal pdf table.
  void hoppetAssign(void (* pdf_subroutine)(const double & x, 
                                            const double & Q, double * res) );

  /// Set up the strong coupling such that alphas(Q)=alphas_Q, with the
  /// given number of loops (nloop).
  ///
  /// The user should have set the quark masses or requested a FFN scheme
  /// prior to calling this function.
  ///
  /// This function is provided mainly for use in conjunction with
  /// hoppetAssign (C++) :func:`Assign` (Python).  In particular, it
  /// has the side effect of modifying the structure of the PDF tables
  /// to make sure they know about the mass thresholds.
  ///
  /// If QED has been requested, a QED coupling will also be set up
  /// (its value is not currently configurable from this interface).
  ///
  /// If you call hoppetEvolve (C++) :func:`Evolve` (Python), there is
  /// no need to separately call this routine.
  ///
  /// @param asQ0 alphas at the scale Q0
  /// @param Q0alphas the scale Q0
  /// @param nloop Perturbative order (1: LO, 2: NLO, 3: NNLO, 4: N3LO)
  void hoppetSetCoupling(const double & asQ0,
                         const double & Q0alphas,
                         const int    & nloop);


  /// Given a pdf_subroutine with the interface shown below, fill the 
  /// table by evolving the PDF from scale Q0pdf, with alphas provided 
  /// at scale Q0alphas
  ///
  /// @param asQ0           value of $\alpha_s$ at scale Q0alphas
  /// @param Q0alphas       Scale at which $\alpha_s$ has been supplied [GeV]
  /// @param nloop          Number of loops for evolution (1=LO, 2=NLO, 3=NNLO, 4=N3LO)
  /// @param muR_Q          Ratio of renormalisation scheme to factorisation scale during evolution
  /// @param pdf_subroutine Function pointer to the user-defined PDF function
  /// @param Q0pdf          Initial scale for PDF [GeV]
  ///
  void hoppetEvolve(const double & asQ0,
                    const double & Q0alphas,
                    const int    & nloop,
                    const double & muR_Q,
                    void (* pdf_subroutine)(const double & x, 
                                            const double & Q, double * res),
                    const double & Q0pdf);


  /// Prepare a cached evolution.  Once this has been called, one can
  /// use hoppetCachedEvolve (C++) or :func:`CachedEvolve` (Python) to
  /// carry out the cached evolution of a specific initial condition.
  ///
  /// @param asQ0     Strong coupling at initial scale Q0alphas
  /// @param Q0alphas Initial scale for alpha_s [GeV]
  /// @param nloop    Number of loops for evolution (1=LO, 2=NLO, 3=NNLO, 4=N3LO)
  /// @param muR_Q    Ratio of renormalisation scheme to factorisation scale during evolution
  /// @param Q0pdf    Initial scale for PDF [GeV]
  ///
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

  /// Return the strong coupling at scale Q
  ///
  /// @param Q Scale in GeV
  /// @return The strong coupling
  double hoppetAlphaS(const double & Q);

  /// Return the QED coupling at scale Q
  ///
  /// @param Q Scale in GeV
  /// @return The QED coupling
  double hoppetAlphaQED(const double & Q);

  /// Return in f[0..N] the value of the internally stored pdf at
  /// the given x,Q. The indices of f correspond to the iflv
  /// constants defined above, e.g. f[iflv_g] is the gluon, etc.
  ///
  /// The array f should be of size at least:
  ///
  /// - 13 for QCD only
  /// - 18 for QCD+photon+leptons (the lepton indices give the sum of lepton and antilepton)
  ///
  void hoppetEval(const double & x,
                  const double & Q,
                  double * f);

  /// the interface to get a single flavour, which takes iflv
  /// starting from -6. This is the direct Fortran interface.
  /// The C++/Python interface is hoppetEvalIFlv below, which takes the
  /// hoppet iflv constants, which have been shifted by +6
  /// to make them non-negative and match the indices to be used
  /// with hoppetEval
  double hoppetEvalFortranIFlv(const double & x,
           const double & Q,
           const int & iflv_starts_minus6);

  /// Return xf(x,Q) for the flavour indicated by iflv, which should
  /// be one of the hoppet iflv_* constants (iflv_g, iflv_d,
  /// iflv_ubar, etc.)
  ///
  /// @param x Longitudinal momentum fraction (0 < x < 1)
  /// @param Q Scale in GeV
  /// @param iflv One of the hoppet.iflv_g, hoppet.iflv_d, etc.
  ///
  /// @return xf(x,Q) for the flavour indicated by iflv
  inline double hoppetEvalIFlv(const double & x,
		       const double & Q,
		       const int & iflv) {return hoppetEvalFortranIFlv(x,Q,iflv-6);}


  /// Return in f[0..12] the value of 
  ///
  ///    [P(iloop,nf) \f$\otimes\f$ pdf] (x,Q)
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
  
  /// Minimal setup of structure functions
  ///
  /// @param order_max highest order in QCD to compute (1: LO, 2: NLO, 3: NNLO, 4: N3LO)
  ///
  void hoppetStartStrFct(const int & order_max);
			 
  
  /// Setup of constants and parameters needed for structure functions
  ///
  /// @param     order_max      highest order in QCD to compute (1: LO, 2: NLO, 3: NNLO, 4: N3LO)
  /// @param     nflav          integer number of flavours (if negative use variable flavour)
  /// @param     scale_choice   (0: fixed scale, 1: use Q, 2: use arbitrary scale)
  /// @param     constant_mu    if scale_choice = scale_choice_fixed (= 0) then this is the fixed scale
  /// @param     param_coefs    if .true. use parametrised coefficients functions
  /// @param     wmass          Mass of the W boson
  /// @param     zmass          Mass of the z boson
  ///
  void hoppetStartStrFctExtended(const int & order_max,
				 const int & nflav,
				 const int & scale_choice,
				 const double & constant_mu,
				 const bool & param_coefs,
				 const double & wmass,
				 const double & zmass);
  
  /// Initialize the structure functions up to specified order
  /// this requires the PDF to have been set up beforehand, and filled in tables(0)
  void hoppetInitStrFct(const int & order_max,
			const bool & separate_orders,
			const double & xR,
			const double & xF);

  /// Initialize the structure functions up to specified order
  /// this requires the PDF to have been set up beforehand, and filled in tables(0)
  void hoppetInitStrFctFlav(const int & order_max,
			    const bool & separate_orders,
			    const double & xR,
			    const double & xF);


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
			const int & iflav,
			double * F);

  /// calculate the structure function at x, Q, with muR and muF as
  /// requested in hoppetStartStrFct. This is the sum over all orders
  void hoppetStrFctNoMuFlav(const double & x,
			    const double & Q,
			    const int & iflav,
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
  

#ifdef __cplusplus
extern "C" {
#endif

// Proper C prototype for direct Fortran call (see bind(C) in Fortran)
  void hoppetwritelhapdfwithlen_c(const int& basename_len, const char* basename, const int& pdf_index);

#ifdef __cplusplus
}
#endif
  
  /// Write out the contents of tables(0) (assumed to be the PDF) in the
  /// LHAPDF format
  ///
  /// @param basename The basename of the output
  /// @param pdf_index The intended index in the LHAPDF meaning. If index is 0 both the PDF and the .info file is saved to disk
  inline void hoppetWriteLHAPDFGrid(const std::string& basename, const int& pdf_index) {
    // Fortran expects a char array, its length, and the index (all as pointers)
    const char* basename_cstr = basename.c_str();
    int basename_len = basename.length();

    hoppetwritelhapdfwithlen_c(basename_len, basename_cstr, pdf_index);

  }

  void hoppetBenchmarkPDFunpol(const double & x,
			       const double &Q,
			       double * xpdf);


  /// Take a pointer to an array of C characters and if the version
  /// string length is < maxlen-1, fill the C-array with a null
  /// terminated string with the version, and return 0. If the version
  /// string length is >= maxlen-1, then do not fill the C-array,
  /// and return the required length (including the null terminator).
  int hoppetVersionC(char *cchar, int maxlen);

}

/// @brief  Return the version of hoppet as a std::string
/// @return The version string
inline std::string hoppetVersion() {
  const int maxlen = 100; // should be more than enough
  char cchar[maxlen];
  int errcode_or_maxlen = hoppetVersionC(cchar, maxlen);
  if (errcode_or_maxlen == 0) {
    return std::string(cchar);
  } else {
    // we need a longer buffer
    char * cchar2 = new char[errcode_or_maxlen];
    errcode_or_maxlen = hoppetVersionC(cchar2, errcode_or_maxlen);
    std::string result(cchar2);
    delete[] cchar2;
    return result;
  }
}


#endif // __HOPPET__
