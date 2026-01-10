#include "inc/hoppet_config.inc"

!! module with helpers for probes for setting derived splitting matrices etc.
module hoppet_probes
  use types; use consts_dp
  implicit none

  !! normalisations for each of the different flavour probes, designed
  !! so that when we change representation from evolution to human there
  !! are no accidental cancellations; the gluon one is used separately
  !! from the others, so it can retain a simple normalisation
  real(dp), parameter :: probe_norm(-2:2) = (/sqrt(5.0_dp), pi, one, sqrt(2.0_dp), one/)
end module hoppet_probes

!======================================================================
!  
!! This module provide a type for matrices of splitting functions, i.e. split_mat
!! objects. It also provides associated 
!!
!! - subroutines for manipulating the splitting matrices (similarly 
!!   to the subroutine for grid_conv objects)
!!
!! - subroutines for initialising the splitting matrices with the LO
!!   NLO and NNLO splitting functions
!!
!! - subroutines for initialising Polarised LO and NLO splitting matrices. 
!!
!! - tools for getting "derived" splitting matrices.
!!
!! Note that it also contains some older routines for coefficient functions,
!! but these are to be considered deprecated since v1.2 of hoppet, since
!! all structure associated with coefficient functions can (and usually should)
!! be held in a split_mat object.
!======================================================================
module hoppet_split_mat
  use types; use consts_dp; use splitting_functions
#ifdef HOPPET_ENABLE_N3LO_FORTRAN_MTM  
  use mass_thresholds_n3lo
#endif
  !use qcd
  use dglap_choices
  use warnings_and_errors
  use pdf_representation
  use convolution
  use assertions
  implicit none
  private

  !-----------------------------------------------------------------
  ! holds all the components of the splitting functions
  type split_mat
     !private
     !-- These are the singlet evolution pieces
     !   qg is defined including a 2nf factor
     type(grid_conv)           :: singlet(iflv_g:iflv_sigma,iflv_g:iflv_sigma)
     type(grid_conv), pointer  :: gg, qq, gq, qg
     !-- These are the non-singlet pieces     ! 
     type(grid_conv) :: NS_plus   ! multiplies [(q_i + qbar_i) - (q_j + qbar_j)]
     type(grid_conv) :: NS_minus  ! multiplies [(q_i - qbar_i) - (q_j - qbar_j)]
     type(grid_conv) :: NS_V      ! multiplies [\sum_i (q_i - qbar_i)]
     ! the nf value for which this splitting matrix is defined
     integer            :: nf_int
  end type split_mat
  public             :: split_mat


  !-----------------------------------------------------------------
  ! holds the components of a coefficient function
  !
  ! NB: as of v1.2 this should be considered obsolete and will be
  ! removed in future versions: split_mat objects should be used instead.
  type coeff_mat
     !private
     !-- need a grid def. Sometimes no g or q will be defined
     !   but still want the coefficient function associated with a grid
     !   def
     !   delta is the magnitude of a delta function piece
     !   quark charge contributions are included by the routines here.
     !   for gluons the implicitness is \sum_{q,qbar} e_q^2
     type(grid_def)     :: grid
     type(grid_conv)    :: g, q
     real(dp)           :: delta
     logical            :: HO
  end type coeff_mat
  public             :: coeff_mat


  !---------- avoid need for main program to access convolution
  !  give it both possible names. Maybe not necessary?
  public :: grid_def
  interface cobj_InitGridDef
     ! perhaps illegal, so do it explicitly?
     !module procedure InitGridDef
     module procedure conv_InitGridDef_single, conv_InitGridDef_multi
  end interface
  public :: cobj_InitGridDef
  public :: InitGridDef

  public :: cobj_InitSplitLinks
  public :: InitSplitMatLO, InitSplitMatNLO
  public :: InitSplitMatNNLO, InitSplitMatN3LO
  public :: InitSplitMatTimeNLO
  public :: InitSplitMatPolLO, InitSplitMatPolNLO

  interface InitSplitMat
     module procedure InitSplitMat_sm_fact, InitSplitMat_zero
  end interface
  public :: InitSplitMat!, Delete_sm
  !public :: AddWithCoeff_sm, Multiply_sm
  !public :: cobj_PConv, cobj_PConv_1d
  
  !-------- things for splitting functions --------------------
  interface cobj_InitCoeff
     module procedure cobj_InitCoeffLO, cobj_InitCoeffHO, cobj_InitCoeff_cf
  end interface
  public :: cobj_InitCoeff, cobj_InitCoeffLO, cobj_InitCoeffHO
  !public :: cobj_CConv
  !public :: cobj_AddCoeff, cobj_DelCoeff
  public :: GetDerivedSplitMatProbes, AllocSplitMat, SetDerivedSplitMat


  interface operator(*)
     module procedure cobj_PConv, cobj_PConv_1d, cobj_CConv
  end interface
  public :: operator(*)
  interface operator(.conv.)
     module procedure cobj_PConv, cobj_PConv_1d, cobj_CConv
  end interface
  public  :: operator(.conv.)
  public  :: cobj_Eval2LConv
  !-- avoid access through here
  private :: cobj_PConv, cobj_PConv_1d, cobj_CConv

  interface SetToZero
     module procedure SetToZero_sm
  end interface
  public :: SetToZero

  interface Multiply
     module procedure Multiply_sm
  end interface
  public :: Multiply

  interface AddWithCoeff
     module procedure AddWithCoeff_sm
  end interface
  public :: AddWithCoeff

  interface SetToConvolution
     module procedure SetToConvolution_sm
  end interface
  public :: SetToConvolution

  interface SetToCommutator
     module procedure SetToCommutator_sm
  end interface
  public :: SetToCommutator

  interface Delete
     module procedure Delete_sm, cobj_DelCoeff
  end interface
  public :: Delete


contains

  !======================================================================
  !! make sure all required internal links are set up for the splitting 
  !! matrix 
  subroutine cobj_InitSplitLinks(P)
    type(split_mat), target,  intent(inout) :: P
    
    if ((iflv_sigma - iflv_g) /= 1) call wae_error(&
         &'cobj_InitSplitLinks:', 'local gluon and sigma ids not consistent')
    
    !-- NB singlet matrix is not written with usual convention of 
    !   (sigma g), but rather (g sigma)
    P%gg => P%singlet(iflv_g,iflv_g)
    P%gq => P%singlet(iflv_g,iflv_sigma)
    P%qg => P%singlet(iflv_sigma,iflv_g)
    P%qq => P%singlet(iflv_sigma,iflv_sigma)
  end subroutine cobj_InitSplitLinks
  

  !======================================================================
  !! Initialise a LO unpolarised splitting matrix, with the nf value that
  !! is current from the qcd module. 
  subroutine InitSplitMatLO(grid, P)
    use qcd
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: P

    !-- info to describe the splitting function
    !P%loops  = 1
    P%nf_int = nf_int

    call cobj_InitSplitLinks(P)

    call InitGridConv(grid, P%gg, sf_Pgg)
    call InitGridConv(grid, P%qq, sf_Pqq)
    call InitGridConv(grid, P%gq, sf_Pgq)
    call InitGridConv(grid, P%qg, sf_Pqg)

    !-- now fix up pieces so that they can be directly used as a matrix
    call Multiply(P%qg, 2*nf)

    !-- PqqV +- PqqbarV
    call InitGridConv(P%NS_plus,  P%qq)
    call InitGridConv(P%NS_minus, P%qq)

    !-- PNSminus + nf * (PqqS - PqqbarS)
    call InitGridConv(P%NS_V, P%NS_minus)
  end subroutine InitSplitMatLO


  !======================================================================
  !! Initialise a NLO unpolarised splitting matrix, with the nf value that
  !! is current from the qcd module. (MSbar scheme)
  subroutine InitSplitMatNLO(grid, P, factscheme)
    use qcd
    type(grid_def),    intent(in)    :: grid
    type(split_mat),   intent(inout) :: P
    integer, optional, intent(in)    :: factscheme
    integer :: factscheme_local
    !-----------------------------------------
    type(grid_conv) :: P1qqV, P1qqbarV, P1qqS
    !-- needed for DIS schemes
    type(grid_conv) :: Cq,Cg
    type(grid_conv) :: CqPqq, CqPqg, CqPgq, CqPgg
    type(grid_conv) :: CgPgg, CgPqg, CgPgq, CgPqq


    factscheme_local = default_or_opt(factscheme_default, factscheme)
    if (factscheme_local /= factscheme_MSbar) then
       ! NB: do not support DIS scheme here since it involves 
       !     determination of Pqq etc at LO (already done once, so do 
       !     not want to repeat) -- rather do this stuff in
       !     dglap_holders, where Pqq(LO) will in any case be available.
       !     (Is this "chickening out"?)
       write(0,*) 'InitSplitMatNLO: unsupported fact scheme', factscheme
       call wae_error('InitSplitMatNLO: stopping')
    end if
    
    !-- info to describe the splitting function
    !P%loops = 2
    P%nf_int = nf_int
    
    call cobj_InitSplitLinks(P)

    !-- these are the building blocks
    call InitGridConv(grid, P1qqV, sf_P1qqV)
    call InitGridConv(grid, P1qqbarV, sf_P1qqbarV)
    call InitGridConv(grid, P1qqS, sf_P1qqS)

    !-- PqqV + PqqbarV
    call InitGridConv(P%NS_plus, P1qqV)
    call AddWithCoeff(P%NS_plus, P1qqbarV, one)
    !-- PqqV - PqqbarV
    call InitGridConv(P%NS_minus, P1qqV)
    call AddWithCoeff(P%NS_minus, P1qqbarV, -one)

    !-- PNSminus + nf * (PqqS - PqqbarS) 
    !   [NB at NLO, PqqS = PqqbarS] 
    call InitGridConv(P%NS_V, P%NS_minus)
    
    !-- Pqq in matrix:  PNS_plus + nf*(PqqS + PqqbarS)
    !   [NB at NLO, PqqS = PqqbarS] 
    call InitGridConv(P%qq, P%NS_plus)
    call AddWithCoeff(P%qq, P1qqS, two*nf)

    !-- rest of singlet matrix
    call InitGridConv(grid, P%gq, sf_P1gq)
    call InitGridConv(grid, P%gg, sf_P1gg)
    call InitGridConv(grid, P%qg, sf_P1qg)
    !-- recall that the way it is defined it needs a factor 2nf
    call Multiply(P%qg, two*nf)

    !-- tidy up 
    call Delete(P1qqV)
    call Delete(P1qqbarV)
    call Delete(P1qqS)

  end subroutine InitSplitMatNLO



  !======================================================================
  !! Initialise a NNLO unpolarised splitting matrix, with the nf value that
  !! is current from the qcd module. (MSbar scheme)
  subroutine InitSplitMatNNLO(grid, P, factscheme)
    use qcd; use convolution_communicator
    type(grid_def),    intent(in)    :: grid
    type(split_mat),   intent(inout) :: P
    integer, optional, intent(in) :: factscheme
    integer :: factscheme_local
    !-----------------------------------------
    type(grid_conv) :: P2NSS
    real(dp) :: dummy

    !call wae_error('InitSplitMatNNLO: NNLO not yet implemented')
    factscheme_local = default_or_opt(factscheme_default, factscheme)
    if (factscheme_local /= factscheme_MSbar) then
       write(0,*) 'InitSplitMatNNLO: unsupported fact scheme', factscheme
       call wae_error('InitSplitMatNNLO: stopping')
    end if
    
    !-- info to describe the splitting function
    !P%loops = 3
    P%nf_int = nf_int

    call cobj_InitSplitLinks(P)

    call InitGridConv(grid, P%NS_plus, sf_P2NSPlus)
    call InitGridConv(grid, P%NS_minus, sf_P2NSMinus)
    
    !-- if understanding of convention is right then P_V = P_- + P_S
    call InitGridConv(P%NS_V, P%NS_minus)
    call InitGridConv(grid, P2NSS, sf_P2NSS)
    call AddWithCoeff(P%NS_V, P2NSS)    
    call Delete(P2NSS)

    !-- now the singlet functions
    call InitGridConv(grid, P%qg, sf_P2qg2nf)
    call InitGridConv(grid, P%gg, sf_P2gg)
    call InitGridConv(grid, P%gq, sf_P2gq)
    !-- qq = "pure-singlet" + P+
    call InitGridConv(grid, P%qq, sf_P2PS)
    call AddWithCoeff(P%qq, P%NS_plus)
  end subroutine InitSplitMatNNLO


  !======================================================================
  !! Initialise a N3LO unpolarised splitting matrix, with the nf value that
  !! is current from the qcd module. (MSbar scheme)
  subroutine InitSplitMatN3LO(grid, P, factscheme)
    use qcd; use convolution_communicator
    type(grid_def),    intent(in)    :: grid
    type(split_mat),   intent(inout) :: P
    integer, optional, intent(in) :: factscheme
    integer :: factscheme_local, nf_store
    !-----------------------------------------
    type(grid_conv) :: P3NSS
    real(dp) :: dummy
    integer, save :: nwarn_nf_lo = 5, nwarn_nf_hi = 5

    !call wae_error('InitSplitMatN3LO: N3LO not yet implemented')
    factscheme_local = default_or_opt(factscheme_default, factscheme)
    if (factscheme_local /= factscheme_MSbar) then
       write(0,*) 'InitSplitMatN3LO: unsupported fact scheme', factscheme
       call wae_error('InitSplitMatN3LO: stopping')
    end if
    
    !-- info to describe the splitting function
    !P%loops = 3
    ! The splitting matrices are currently hard coded for nf=3,4,5
    nf_store = nf_int
    if (nf_int < 3) then
      nf_int = 3
      call wae_warn(nwarn_nf_lo, "InitSplitMatN3LO: nf_int < 3, setting to 3; nf_int was", intval = nf_store)
    else if (nf_int > 5 .and. n3lo_splitting_approximation .lt.n3lo_splitting_approximation_up_to_2512) then
      nf_int = 5
      call wae_warn(nwarn_nf_hi, "InitSplitMatN3LO: nf_int > 5, setting to 5; nf_int was", intval = nf_store)
    end if
    !nf_int = min(max(nf_int,3),5)
    P%nf_int = nf_store

    call cobj_InitSplitLinks(P)

    call InitGridConv(grid, P%NS_plus, sf_P3NSPlus)
    call InitGridConv(grid, P%NS_minus, sf_P3NSMinus)
    
    !-- if understanding of convention is right then P_V = P_- + P_S
    call InitGridConv(P%NS_V, P%NS_minus)
    call InitGridConv(grid, P3NSS, sf_P3NSS)
    call AddWithCoeff(P%NS_V, P3NSS)    
    call Delete(P3NSS)

    !-- now the singlet functions
    call InitGridConv(grid, P%qg, sf_P3qg2nf)
    call InitGridConv(grid, P%gg, sf_P3gg)
    call InitGridConv(grid, P%gq, sf_P3gq)
    !-- qq = "pure-singlet" + P+
    call InitGridConv(grid, P%qq, sf_P3PS)
    call AddWithCoeff(P%qq, P%NS_plus)

    ! Restore nf_int
    nf_int = nf_store
  end subroutine InitSplitMatN3LO

  !======================================================================
  !! Initialise a NLO unpolarised splitting matrix, with the nf value that
  !! is current from the qcd module. (MSbar scheme)
  subroutine InitSplitMatTimeNLO(grid, P)
    use qcd
    type(grid_def),    intent(in)    :: grid
    type(split_mat),   intent(inout) :: P
    !-----------------------------------------
    type(grid_conv) :: P1qqV, P1qqbarV

    !-- info to describe the splitting function
    !P%loops = 2
    P%nf_int = nf_int
    
    call cobj_InitSplitLinks(P)

    !-- these are the building blocks
    call InitGridConv(grid, P1qqV, sf_P1qqV)
    ! add in the timelike correction
    call AddWithCoeff(P1qqV, sf_TmSP1qqNS)
    call InitGridConv(grid, P1qqbarV, sf_P1qqbarV)

    !-- PqqV + PqqbarV
    call InitGridConv(P%NS_plus, P1qqV)
    call AddWithCoeff(P%NS_plus, P1qqbarV, one)
    !-- PqqV - PqqbarV
    call InitGridConv(P%NS_minus, P1qqV)
    call AddWithCoeff(P%NS_minus, P1qqbarV, -one)

    !-- PNSminus + nf * (PqqS - PqqbarS) 
    !   [NB at NLO, PqqS = PqqbarS] 
    call InitGridConv(P%NS_V, P%NS_minus)
    
    !--  singlet matrix
    ! Note that for the singlet matrix here, we are using the ESW notation
    ! which swaps qg/gq names in the matrix, i.e.
    !
    ! - the matrix qg entry will multiply a gluon fragementation function
    !   so it corresponds to q->g splitting which is similar to the initial-state
    !   Pgq; the q entry is the sum over all quark fragmentation functions; since
    !   each quark can branch to a gluon we need a 2nf factor
    !
    ! - the matrix gq entry will multiply a quark fragmentation function, so
    !   it corresponds to a g->q splitting, which is similar to the initial-state
    !   Pqg; the ESW definition includes a factor 2nf which should not be there
    !   because the singlet already includes the sum over all flavours.
    !
    call InitGridConv(grid, P%qq, sf_TP1qq)
    call InitGridConv(grid, P%gg, sf_TP1gg)
    call InitGridConv(grid, P%gq, sf_TP1qg) ! note swap of qg/gq for time-like case
    call InitGridConv(grid, P%qg, sf_TP1gq) ! note swap of qg/gq for time-like case
    !-- recall that the way it is defined it needs a factor 2nf
    call Multiply(P%qg, two*nf)
    call Multiply(P%gq, one/(two*nf))

    !-- tidy up 
    call Delete(P1qqV)
    call Delete(P1qqbarV)

  end subroutine InitSplitMatTimeNLO



  !======================================================================
  !! Initialise a LO Polarised splitting matrix, with the nf value that
  !! is current from the qcd module.
  subroutine InitSplitMatPolLO(grid, P)
    use qcd
    type(grid_def),    intent(in)    :: grid
    type(split_mat),   intent(inout) :: P

    !-- info to describe the splitting function
    !P%loops  = 1
    P%nf_int = nf_int

    call cobj_InitSplitLinks(P)

    call InitGridConv(grid, P%gg, sf_DPgg)
    call InitGridConv(grid, P%qq, sf_DPqq)
    call InitGridConv(grid, P%gq, sf_DPgq)
    call InitGridConv(grid, P%qg, sf_DPqg)

    !-- now fix up pieces so that they can be directly used as a matrix
    call Multiply(P%qg, 2*nf)

    !-- PqqV +- PqqbarV
    call InitGridConv(P%NS_plus,  P%qq)
    call InitGridConv(P%NS_minus, P%qq)

    !-- PNSminus + nf * (PqqS - PqqbarS)
    call InitGridConv(P%NS_V, P%NS_minus)
  end subroutine InitSplitMatPolLO

  !======================================================================
  !! Initialise a NLO Polarised splitting matrix, with the nf value that
  !! is current from the qcd module.
  subroutine InitSplitMatPolNLO(grid, P, factscheme)
    use qcd
    type(grid_def),    intent(in)    :: grid
    type(split_mat),   intent(inout) :: P
    integer, optional, intent(in)    :: factscheme
    integer :: factscheme_local
    !-----------------------------------------
    type(grid_conv) :: DP1qqV, DP1qqbarV, DP1qqS

    factscheme_local = default_or_opt(factscheme_Poldefault, factscheme)
    if (factscheme_local /= factscheme_PolMSbar) then
       ! NB: do not support DIS scheme here since it involves 
       !     determination of Pqq etc at LO (already done once, so do 
       !     not want to repeat) -- rather do this stuff in
       !     dglap_holders, where Pqq(LO) will in any case be available.
       !     (Is this "chickening out"?)
       write(0,*) 'InitSplitMatPolNLO: unsupported fact scheme', factscheme
       call wae_error('InitSplitMatPolNLO: stopping')
    end if
    
    !-- info to describe the splitting function
    !P%loops = 2
    P%nf_int = nf_int
    
    call cobj_InitSplitLinks(P)

    !-- these are the building blocks
    call InitGridConv(grid, DP1qqV, sf_DP1qqV)
    call InitGridConv(grid, DP1qqbarV, sf_DP1qqbarV)
    call InitGridConv(grid, DP1qqS, sf_DP1qqS)

    !-- PqqV + PqqbarV
    call InitGridConv(P%NS_plus, DP1qqV)
    call AddWithCoeff(P%NS_plus, DP1qqbarV, one)
    !-- PqqV - PqqbarV
    call InitGridConv(P%NS_minus, DP1qqV)
    call AddWithCoeff(P%NS_minus, DP1qqbarV, -one)

    !-- PNSminus + nf * (PqqS - PqqbarS) 
    !   [NB at NLO, PqqS = PqqbarS] 
    call InitGridConv(P%NS_V, P%NS_minus)
    
    !-- Pqq in matrix:  PNS_plus + nf*(PqqS + PqqbarS)
    !   [NB at NLO, PqqS = PqqbarS] 
    call InitGridConv(P%qq, P%NS_plus)
    call AddWithCoeff(P%qq, DP1qqS, two*nf)

    !-- rest of singlet matrix
    call InitGridConv(grid, P%gq, sf_DP1gq)
    call InitGridConv(grid, P%gg, sf_DP1gg)
    call InitGridConv(grid, P%qg, sf_DP1qg)
    !-- recall that the way it is defined it needs a factor 2nf
    call Multiply(P%qg, two*nf)

    !-- tidy up 
    call Delete(DP1qqV)
    call Delete(DP1qqbarV)
    call Delete(DP1qqS)

  end subroutine InitSplitMatPolNLO


  !----------------------------------------------------------------------
  !! initialize a splitting function matrix with another one (potentially 
  !! multiplied by some factor). Memory allocation (if needed) is automatically
  !! handled by the subsiduary routines.
  subroutine InitSplitMat_sm_fact(P, Pin, factor)
    type(split_mat),  intent(inout)        :: P
    type(split_mat),  intent(in)           :: Pin
    real(dp),         intent(in), optional :: factor
    !if (nf_d /= (nf_int+1)/2) write(0,*) 'WARNING: non-standard nf_d'

    !P%loops  = Pin%loops
    P%nf_int = Pin%nf_int
    call cobj_InitSplitLinks(P)

    call InitGridConv(P%gg,    Pin%gg,     factor)
    call InitGridConv(P%qq,    Pin%qq,     factor)
    call InitGridConv(P%gq,    Pin%gq,     factor)
    call InitGridConv(P%qg,    Pin%qg,     factor)

    !-- from here on, one could exploit info in loops so as to 
    !   reduce number of multiplications... But I am too lazy in this first
    !   go
    call InitGridConv(P%NS_plus,  Pin%NS_plus,  factor)
    call InitGridConv(P%NS_minus, Pin%NS_minus, factor)
    call InitGridConv(P%NS_V,     Pin%NS_V,     factor)
    
  end subroutine InitSplitMat_sm_fact

  !! Initialise the splitting matrix P to zero, with a structure
  !! corresponding to the value of nf as supplied
  subroutine InitSplitMat_zero(grid, P, nf)
    use convolution, only: InitGridConv, grid_def
    type(grid_def),   intent(in)    :: grid
    type(split_mat),  intent(inout) :: P
    integer,          intent(in)    :: nf
    P%nf_int = nf
    call cobj_InitSplitLinks(P)

    call InitGridConv(grid, P%gg)
    call InitGridConv(grid, P%qq)
    call InitGridConv(grid, P%gq)
    call InitGridConv(grid, P%qg)

    call InitGridConv(grid, P%NS_plus )
    call InitGridConv(grid, P%NS_minus)
    call InitGridConv(grid, P%NS_V    )
    
  end subroutine InitSplitMat_zero


  !----------------------------------------------------------------------
  !! init a splitting function set with another one (potentially 
  !! multiplied by some factor)
  subroutine AddWithCoeff_sm(P, Padd, factor)
    type(split_mat),     intent(inout) :: P
    type(split_mat),     intent(in)    :: Padd
    real(dp),       intent(in), optional :: factor
    !if (nf_d /= (nf_int+1)/2) write(0,*) 'WARNING: non-standard nf_d'

    !P%loops = max(P%loops, Padd%loops)
    P%nf_int = assert_eq(P%nf_int, Padd%nf_int, &
         &'AddWithCoeff_sm: nf must be the same')
    call AddWithCoeff(P%gg,       Padd%gg,       factor)
    call AddWithCoeff(P%qq,       Padd%qq,       factor)
    call AddWithCoeff(P%gq,       Padd%gq,       factor)
    call AddWithCoeff(P%qg,       Padd%qg,       factor)
    call AddWithCoeff(P%NS_plus,  Padd%NS_plus,  factor)
    call AddWithCoeff(P%NS_minus, Padd%NS_minus, factor)
    call AddWithCoeff(P%NS_V,     Padd%NS_V,     factor)

  end subroutine AddWithCoeff_sm


  !----------------------------------------------------------------------
  !! Multiply the splitting matrix P by factor.
  subroutine Multiply_sm(P, factor)
    type(split_mat),  intent(inout) :: P
    real(dp),         intent(in)    :: factor
    !if (nf_d /= (nf_int+1)/2) write(0,*) 'WARNING: non-standard nf_d'

    call Multiply(P%gg,       factor)
    call Multiply(P%qq,       factor)
    call Multiply(P%gq,       factor)
    call Multiply(P%qg,       factor)
    call Multiply(P%NS_plus,  factor)
    call Multiply(P%NS_minus, factor)
    call Multiply(P%NS_V,     factor)

  end subroutine Multiply_sm


  !======================================================================
  !! Set the splitting matrix P to zero
  subroutine SetToZero_sm(P)
    type(split_mat),  intent(inout) :: P

    call SetToZero(P%gg)
    call SetToZero(P%qq)
    call SetToZero(P%gq)
    call SetToZero(P%qg)
    call SetToZero(P%NS_plus)
    call SetToZero(P%NS_minus)
    call SetToZero(P%NS_V)

  end subroutine SetToZero_sm


  !======================================================================
  !! Set the splitting matrix PA to be equal to PB .conv. PC
  subroutine SetToConvolution_sm(PA, PB, PC)
    type(split_mat),  intent(inout) :: PA
    type(split_mat),  intent(in)    :: PB, PC

    PA%nf_int = assert_eq(PB%nf_int, PC%nf_int, 'SetToConvolution_sm')
    call cobj_InitSplitLinks(PA)

    call SetToConvolution(PA%singlet , PB%singlet , PC%singlet )
    call SetToConvolution(PA%NS_plus , PB%NS_plus , PC%NS_plus )
    call SetToConvolution(PA%NS_minus, PB%NS_minus, PC%NS_minus)
    call SetToConvolution(PA%NS_V    , PB%NS_V    , PC%NS_V    )

  end subroutine SetToConvolution_sm


  !----------------------------------------------------------------------
  !! Set the splitting matrix PA to be equal to [PB, PC]
  subroutine SetToCommutator_sm(PA, PB, PC)
    type(split_mat),  intent(inout) :: PA
    type(split_mat),  intent(in)    :: PB, PC

    PA%nf_int = assert_eq(PB%nf_int, PC%nf_int, 'SetToCommutator_sm')
    call cobj_InitSplitLinks(PA)

    call SetToCommutator(PA%singlet, PB%singlet, PC%singlet)
    ! Recall that following calls also set things to zero
    call InitGridConv(PB%NS_plus%grid, PA%NS_plus )
    call InitGridConv(PB%NS_plus%grid, PA%NS_minus)
    call InitGridConv(PB%NS_plus%grid, PA%NS_V    )

  end subroutine SetToCommutator_sm


  !----------------------------------------------------------------------
  !! Clear a splitting function (dealloc memory)
  subroutine Delete_sm(P)
    type(split_mat),     intent(inout) :: P

    call Delete(P%gg); nullify(P%gg)
    call Delete(P%qq); nullify(P%qq)
    call Delete(P%gq); nullify(P%gq)
    call Delete(P%qg); nullify(P%qg)
    call Delete(P%NS_plus)
    call Delete(P%NS_minus)
    call Delete(P%NS_V)
  end subroutine Delete_sm

  !----------------------------------------------------------------------
  !! Return the convolution of the splitting matrix P with q_in.
  !!
  !! If q_in is in the human representation is is returned in the
  !! human representation (though internal conversions are carried
  !! out). If it is in an evolution representation then it is
  !! _assumed_ to be in the representation that corresponds to the nf
  !! value of P (and no consistency checking is performed) and no
  !! change of representation is carried out.
  !!
  !! This function also has a 1d overloaded version
  function cobj_PConv(P,q_in) result(Pxq)
    type(split_mat), intent(in)         :: P
    real(dp),        intent(in), target :: q_in(0:,ncompmin:)
    real(dp)                     :: Pxq(0:ubound(q_in,dim=1),ncompmin:ncompmax)
    !-- for when we have to change rep
    integer                       :: pdfr
    !type(pdf_rep)                :: prep
    real(dp), allocatable, target :: q_ev(:,:)
    real(dp), pointer             :: q(:,:)
    integer :: i

    if (ncomponents < P%nf_int) then
       call wae_error('cobj_Pconv:',&
            &'ncomponents in representation is < nf in P.')
    end if

    pdfr = GetPdfRep(q_in)
    if (pdfr /= pdfr_Human) then
       ! check we have the right number of flavours
       if (GetPdfRep(q_in) /= p%nf_int) then
          call wae_error('cobj_Pconv:',&
               &'n-flavours of q_in representation /= n-flavours for splitting function = ', &
               &intval=GetPdfRep(q_in), dbleval=one*p%nf_int)
       end if
       q => q_in
    else
       allocate(q_ev(0:ubound(q_in,dim=1),ncompmin:ncompmax))
       !prep = DefaultEvlnRep(p%nf_int)
       call CopyHumanPdfToEvln(p%nf_int, q_in, q_ev)
       q => q_ev
    end if
    
    PxQ(:, iflv_V)     = P%NS_V * q(:,iflv_V)
    !PxQ(:, iflv_g)     = P%gg * q(:,iflv_g) + P%gq * q(:,iflv_sigma)
    !PxQ(:, iflv_sigma) = P%qg * q(:,iflv_g) + P%qq * q(:,iflv_sigma)
    PxQ(:,iflv_g:iflv_sigma) = P%singlet .conv. q(:,iflv_g:iflv_sigma)

    !-- we need nf-1 NS (+ and -) pieces
    do i = 2, P%nf_int
       PxQ(:,+i) = P%NS_plus  * q(:,+i)
       PxQ(:,-i) = P%NS_minus * q(:,-i)
    end do

    !-- everything else should be set to zero
    do i = P%nf_int + 1, ncomponents
       PxQ(:,+i) = zero
       PxQ(:,-i) = zero
    end do

    !call LabelPdfAsRep(Pxq,pdfr_Evln)
    call LabelPdfAsRep(Pxq,P%nf_int)
    if (pdfr == pdfr_Human) then
       q_ev = Pxq ! avoid an extra temporary... (real cleverness might
                  ! have been to avoid this copy...?)
       !call CopyEvlnPdfToHuman(prep, q_ev, Pxq) ! does labelling as well
       call CopyEvlnPdfToHuman(p%nf_int, q_ev, Pxq) ! does labelling as well
       deallocate(q_ev)
    end if
    
  end function cobj_PConv

  !-----------------------------------------------------------------
  !! In some situations the 1d overloaded version is useful.
  !! Do just for PConv, but think also about doing it for others
  !! later on...
  function cobj_PConv_1d(P,q_in) result(Pxq)
    type(split_mat), intent(in)         :: P
    real(dp),        intent(in), target :: q_in(0:,ncompmin:,:)
    real(dp)                            :: &
         &Pxq(0:ubound(q_in,dim=1),ncompmin:ncompmax,size(q_in,dim=3))
    integer :: i
    do i = 1, size(q_in,dim=3)
       Pxq(:,:,i) = cobj_PConv(P,q_in(:,:,i))
    end do
  end function cobj_PConv_1d
  


  !------------------------------------------------------------
  !! Allocate space for a set of splitting functions, labelling
  !! also with the appropriate nf and nloops...
  subroutine AllocSplitMat(grid,P,nf_in)!, nloops)
    type(grid_def),    intent(in)   :: grid
    type(split_mat),   intent(out)  :: P
    integer,           intent(in)   :: nf_in
    !integer, optional, intent(in)   :: nloops
    
    P%nf_int = nf_in
    !P%loops  = default_or_opt(0,nloops)
    call cobj_InitSplitLinks(P)
    call AllocGridConv(grid,P%gg      )
    call AllocGridConv(grid,P%qq      )
    call AllocGridConv(grid,P%gq      )
    call AllocGridConv(grid,P%qg      )
    call AllocGridConv(grid,P%NS_plus )
    call AllocGridConv(grid,P%NS_minus)
    call AllocGridConv(grid,P%NS_V    )
  end subroutine AllocSplitMat

  
  !-----------------------------------------------------------------
  !! Returns the set of probes needed to establish a matrix of
  !! "derived" effective splitting functions.
  subroutine GetDerivedSplitMatProbes(grid, nf_in, probes)
    use hoppet_probes, only: probe_norm
    type(grid_def), intent(in) :: grid
    integer,        intent(in) :: nf_in
    real(dp),       pointer    :: probes(:,:,:)
    !-----------
    real(dp), pointer :: probes_1d(:,:)
    integer :: nprobes, nprobes_1d, iprobe

    call GetDerivedProbes(grid,probes_1d)

    nprobes_1d = ubound(probes_1d,2)

    ! need to double the number of probes because we need to deduce
    ! the 2x2 matrix for singlet evolution (so once we'll do it with
    ! the singlet quark set to zero; the other time we'll do it with the
    ! singlet gluon set to zero).
    nprobes    = 2*nprobes_1d

    allocate(probes(0:ubound(probes_1d,1),ncompmin:ncompmax, 1:nprobes))

    probes = zero
    do iprobe = 1, nprobes
       ! make sure representation is correct...
       !call LabelPdfAsRep(probes(:,:,iprobe),pdfr_Evln)
       call LabelPdfAsRep(probes(:,:,iprobe),nf_in)
    end do

    probes(:,    iflv_V, 1:nprobes_1d) = probe_norm(    iflv_V) * probes_1d ! NS_V
    probes(:,         2, 1:nprobes_1d) = probe_norm(         2) * probes_1d ! NS+
    probes(:,        -2, 1:nprobes_1d) = probe_norm(        -2) * probes_1d ! NS-
    probes(:,iflv_sigma, 1:nprobes_1d) = probe_norm(iflv_sigma) * probes_1d ! Quark column of singlet

    ! gluon column of singlet
    probes(:,iflv_g, nprobes_1d+1:nprobes) = probe_norm(iflv_g) * probes_1d       

    ! normally this would be part of housekeeping job of convolution
    ! module (SetDerivedConv), but since probes_1d is lost from 
    ! view at this point, we have to do it ourselves...
    deallocate(probes_1d)
  end subroutine GetDerivedSplitMatProbes

  
  !---------------------------------------------------------------------
  !! Given an allocated split_mat and the results of operating on the probes
  !! determine the resulting "derived" splitting function.
  subroutine SetDerivedSplitMat(P,probes)
    use hoppet_probes, only: probe_norm
    type(split_mat),  intent(inout) :: P
    real(dp),         pointer       :: probes(:,:,:)
    !-----------------------------------------
    integer :: nprobes_1d,nprobes, il, ih

    nprobes   = size(probes,dim=3)
    nprobes_1d = nprobes/2

    il = 1; ih = nprobes_1d
    call SetDerivedConv_nodealloc(P%NS_V,     probes(:,    iflv_V,il:ih)/probe_norm(    iflv_V))
    call SetDerivedConv_nodealloc(P%NS_plus,  probes(:,         2,il:ih)/probe_norm(         2))
    call SetDerivedConv_nodealloc(P%NS_minus, probes(:,        -2,il:ih)/probe_norm(        -2))
    call SetDerivedConv_nodealloc(P%gq,       probes(:,    iflv_g,il:ih)/probe_norm(iflv_sigma))
    call SetDerivedConv_nodealloc(P%qq,       probes(:,iflv_sigma,il:ih)/probe_norm(iflv_sigma))
    
    il = nprobes_1d+1; ih = nprobes
    call SetDerivedConv_nodealloc(P%gg,       probes(:,iflv_g    ,il:ih)/probe_norm(iflv_g))
    call SetDerivedConv_nodealloc(P%qg,       probes(:,iflv_sigma,il:ih)/probe_norm(iflv_g))

    deallocate(probes)
  end subroutine SetDerivedSplitMat
  

  !======================================================================
  ! From here onwards we have things to do with coefficient functions
  subroutine cobj_InitCoeffLO(grid,C, factor)
    type(grid_def), intent(in)    :: grid
    type(coeff_mat),     intent(inout) :: C
    real(dp),       intent(in), optional :: factor

    C%grid = grid
    C%HO = .false.
    if (present(factor)) then
       C%delta = factor
    else
       C%delta = one
    end if
  end subroutine cobj_InitCoeffLO
  
  !-- HO. coeff function is not quite so simple. Note that we do simply 
  !   call it NLO, because it could also be relevant at O(as^2).
  !   NOTE UNUSUAL ORDER...
  subroutine cobj_InitCoeffHO(grid, C, coeff_g, coeff_q)
    use convolution_communicator
    type(grid_def), intent(in)    :: grid
    type(coeff_mat),     intent(inout) :: C
    interface
       function coeff_g(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: coeff_g
       end function coeff_g
    end interface
    interface
       function coeff_q(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: coeff_q
       end function coeff_q
    end interface
    real(dp) :: sanity1, sanity2

    C%grid = grid
    C%HO = .true.
    C%delta = zero
    call InitGridConv(grid,C%g, coeff_g)
    call InitGridConv(grid,C%q, coeff_q)
    !-- a sanity check: quite often the arguments for the
    !   quark and gluon pieces get exchanged.
    !   Cannot always ensure correctness (e.g. FL), but can 
    !   detect some cases of misuse
    cc_piece = cc_VIRT
    sanity1 = coeff_g(one)
    cc_piece = cc_DELTA
    sanity2 = coeff_g(zero)
    if (sanity1 /= zero .or. sanity2 /= zero) then
       write(0,*) 'WARNING in cobj_InitCoeffHO **********************'
       write(0,*) 'gluon coefficient function has virtual corrections'
       write(0,*) 'this could be a sign that quark and gluon cf fns have&
            & been exchanged'
       write(0,*) '**************************************************'
    end if
  end subroutine cobj_InitCoeffHO

  !-- initialise C as being factor*Cin ----------------------
  subroutine cobj_InitCoeff_cf(C, Cin, factor)
    type(coeff_mat),     intent(inout) :: C
    type(coeff_mat),     intent(in)    :: Cin
    real(dp),       intent(in), optional :: factor

    C%grid = Cin%grid
    C%HO = Cin%HO
    C%delta = Cin%delta
    if (present(factor)) C%delta = C%delta * factor
    if (C%HO) then
       call InitGridConv(C%q, Cin%q, factor)
       call InitGridConv(C%g, Cin%g, factor)
    end if
  end subroutine cobj_InitCoeff_cf

  
  !-- initialise C as being factor*Cin ----------------------
  subroutine cobj_AddCoeff(C, Cin, factor)
    type(coeff_mat),     intent(inout) :: C
    type(coeff_mat),     intent(in)    :: Cin
    real(dp),       intent(in), optional :: factor

    call ValidateGD(C%grid,Cin%grid,'cobj_AddCoeff')
    if (present(factor)) then 
       C%delta = C%delta  + Cin%delta * factor
    else
       C%delta = C%delta  + Cin%delta
    end if
    
    if (Cin%HO) then
       if (C%HO) then
          call AddWithCoeff(C%q, Cin%q, factor)
          call AddWithCoeff(C%g, Cin%g, factor)
       else
          call InitGridConv(C%q, Cin%q, factor)
          call InitGridConv(C%g, Cin%g, factor)
       end if
    end if
    C%HO = C%HO .or. Cin%HO
  end subroutine cobj_AddCoeff

  !----------------------------------------------------------------------
  ! Clear a coefficient function (dealloc memory)
  subroutine cobj_DelCoeff(C)
    type(coeff_mat),     intent(inout) :: C

    if (C%HO) then
       call Delete(C%g)
       call Delete(C%q)
    end if
  end subroutine cobj_DelCoeff



  !-------------------------------------------------------------------
  function cobj_CConv(C,q) result(Cxq)
    type(coeff_mat), intent(in) :: C
    real(dp),   intent(in) :: q(0:,0:)
    real(dp)               :: Cxq(0:ubound(q,dim=1))
    !--------------------------------------------------
    real(dp) :: ud(0:ubound(q,dim=1))
    real(dp), parameter :: uq2 = (4.0_dp/9.0_dp)
    real(dp), parameter :: dq2 = (1.0_dp/9.0_dp)

    call wae_error('cobj_CConv:', &
         &'coefficient functions not yet suppoprted with new representation')
    Cxq = zero
    !REPLACE ud = q(:,iflv_u)*uq2 + q(:,iflv_d)*dq2
    !REPLACE if (C%HO) then
    !REPLACE    !-- CHECK: factor of two in front of nf_u and nf_d is there both 
    !REPLACE    !   in ESW and the ADS paper.
    !REPLACE    Cxq = (two*(nf_u*uq2 + nf_d*dq2)) * (C%g .conv. q(:,iflv_g))&
    !REPLACE         & + (C%q .conv. ud)
    !REPLACE    if (C%delta /= zero) Cxq = Cxq + C%delta * ud
    !REPLACE else
    !REPLACE    Cxq = C%delta * ud
    !REPLACE end if
  end function cobj_CConv

  !-------------------------------------------------------------------
  ! Calculates [(C2 - y^2/(1+(1-y)^2) CL) .conv. pdf] .atx. xbj
  ! put it here to have access to .conv. for coefficient functions
  function cobj_Eval2LConv(C2, CL, pdf, xbj, ybj) result(res)
    use types; use consts_dp; use convolution
    type(coeff_mat), intent(in)  :: C2, CL
    real(dp),   intent(in)  :: pdf(0:,0:), xbj, ybj
    real(dp)                :: res
    !----------------------------------------------
    type(gdval) :: gdx
    type(coeff_mat)  :: C2L

    call cobj_InitCoeff(C2L, C2)
    call cobj_AddCoeff (C2L, CL, - ybj**2/(1+(1-ybj)**2) )
    gdx = xbj .with. C2%grid
    res = EvalGridQuant(gdx%grid, (C2L .conv. pdf), -log(gdx%val))
    !res = ((C2L.conv.pdf) .atx. gdx ) 
    call cobj_DelCoeff (C2L)
  end function cobj_Eval2LConv
    
end module hoppet_split_mat
!end module dglap_objects_hidden


!-----------------------------------------------------------------
!! module that defines mass threshold matrices
!! and associated manipulations
module hoppet_mass_threshold_mat
  use types; use consts_dp
  use convolution
  use hoppet_split_mat
  implicit none

  private

  !-----------------------------------------------------------------
  !! A type to hold the general structure of a mass threshold matrix,
  !! consisting of a standard splitting matrix acting on the light flavours,
  !! plus three extra pieces to handle the heavy quark contributions.
  !!
  !! Note: a subset of the structure is discussed in 
  !! See logbooks/2025-12-flavour-struct for some discussion
  !! of how this type will work Eqs.(32-35) of 1403.6356, and 1406.4654
  !! plus new terms for Q-Qbar heavy-quark contributions in 2512.13508
  type :: mass_threshold_mat
    type(split_mat) :: P_light !!< the part that acts on the light flavours
    type(grid_conv) :: PShq !!< (h+hbar) from singlet(nflight)
    type(grid_conv) :: PShg !!< (h+hbar) from gluon  (nflight)
    type(grid_conv) :: NShV !!< (h-hbar) from valence(nflight) i.e. sum_i (q_i-qbar_i)
  end type mass_threshold_mat

  public :: mass_threshold_mat

  interface InitMTM
    module procedure InitMTM_zero, InitMTM_from_split_mat, InitMTM_from_new_mtm
  end interface
  public :: InitMTM

  interface operator(*)
    module procedure ConvMTMNew
  end interface
  interface operator(.conv.)
    module procedure ConvMTMNew
  end interface
  interface Delete
    module procedure DelNewMTM
  end interface
  interface Multiply
    module procedure Multiply_mtm
  end interface
  interface AddWithCoeff
    module procedure AddWithCoeff_mtm_mtm, AddWithCoeff_mtm_splitmat
  end interface
  public :: operator(*), operator(.conv.), Delete, Multiply, AddWithCoeff

  public :: GetDerivedMTMProbes, SetDerivedMTM

  interface SetToConvolution
    module procedure SetToConvolution_mtm_sm, SetToConvolution_sm_mtm
  end interface
  public :: SetToConvolution

contains

  !! Initialise an mtm to zero, with the given value of nf_light
  subroutine InitMTM_zero(grid, mtm, nf_light)
    type(grid_def),           intent(in)    :: grid      !< the grid to use
    type(mass_threshold_mat), intent(inout) :: mtm       !< the mass threshold matrix to initialise
    integer,                  intent(in)    :: nf_light  !< number of light flavours
    !---------------------------------------------
    call InitSplitMat(grid, mtm%P_light, nf=nf_light)

    call InitGridConv(grid, mtm%PShg)
    call InitGridConv(grid, mtm%PShq)
    call InitGridConv(grid, mtm%NShV)
  end subroutine InitMTM_zero

  !! Initialise a mass threshold object from nf-1 -> nf from a splitting
  !! matrix with nf flavours.
  subroutine InitMTM_from_split_mat(mtm, P)
    type(mass_threshold_mat), intent(out) :: mtm
    type(split_mat),              intent(in)  :: P
    !---------------------------------------------
    integer :: nf_light_int
    real(dp) :: nf_light_dp, nf_heavy_dp
    type(grid_conv) :: ps_plus, ps_minus

    nf_light_int = P%nf_int - 1
    nf_light_dp  = real(nf_light_int, dp)
    nf_heavy_dp  = real(nf_light_int + 1, dp)

    call cobj_InitSplitLinks(mtm%P_light)
    mtm%P_light%nf_int = nf_light_int

    ! derivation of how to do this is given in logbooks/2025-12-flavour-struct
    ! section 2

    ! set up things that carry over directly
    call InitGridConv(mtm%P_light%NS_plus , P%NS_plus )
    call InitGridConv(mtm%P_light%NS_minus, P%NS_minus)
    call InitGridConv(mtm%P_light%gq      , P%gq      )
    call InitGridConv(mtm%P_light%gg      , P%gg      )

    ! set up terms that just need an overall factor
    call InitGridConv(mtm%P_light%qg      , P%qg      )
    call Multiply    (mtm%P_light%qg, nf_light_dp / nf_heavy_dp)

    call InitGridConv(mtm%PShg, P%qg)
    call Multiply    (mtm%PShg, one / nf_heavy_dp)

    ! now get the "pure singlet" (plus and minus -- i.e. pure singlet and pure valence)
    call InitGridConv(ps_plus , P%qq)
    call AddWithCoeff(ps_plus , P%NS_plus , -one)

    call InitGridConv(ps_minus, P%NS_V)
    call AddWithCoeff(ps_minus, P%NS_minus, -one)

    ! and use those to set up the remaining terms
    call InitGridConv(mtm%P_light%qq  , ps_plus)
    call Multiply    (mtm%P_light%qq  , nf_light_dp / nf_heavy_dp)
    call AddWithCoeff(mtm%P_light%qq  , P%NS_plus)
  
    call InitGridConv(mtm%P_light%NS_V, ps_minus)
    call Multiply    (mtm%P_light%NS_V, nf_light_dp / nf_heavy_dp)
    call AddWithCoeff(mtm%P_light%NS_V, P%NS_minus)

    call InitGridConv(mtm%PShq, ps_plus)
    call Multiply    (mtm%PShq, one / nf_heavy_dp)

    call InitGridConv(mtm%NShV, ps_minus)
    call Multiply    (mtm%NShV, one / nf_heavy_dp)

    call Delete(ps_plus )
    call Delete(ps_minus)

  end subroutine InitMTM_from_split_mat

  !----------------------------------------------
  !! Initialise an mtm from another one
  subroutine InitMTM_from_new_mtm(mtm, mtm_in)
    type(mass_threshold_mat), intent(inout) :: mtm
    type(mass_threshold_mat), intent(in)  :: mtm_in

    call InitSplitMat(mtm%P_light, mtm_in%P_light)

    call InitGridConv(mtm%PShg, mtm_in%PShg)
    call InitGridConv(mtm%PShq, mtm_in%PShq)
    call InitGridConv(mtm%NShV, mtm_in%NShV)
  end subroutine InitMTM_from_new_mtm

  !----------------------------------------------
  !! return the convolution of mtm with q
  function ConvMTMNew(mtm,q) result(Pxq)
    use pdf_representation
    type(mass_threshold_mat), intent(in) :: mtm
    real(dp),   intent(in) :: q(0:,ncompmin:)
    real(dp)               :: Pxq(0:ubound(q,dim=1),ncompmin:ubound(q,dim=2))
    real(dp) :: singlet(0:ubound(q,dim=1)), valence(0:ubound(q,dim=1)), qbar_sum(0:ubound(q,dim=1))
    real(dp) :: dq_from_singlet(0:ubound(q,dim=1)) !!< addition to each flavour from singlet
    real(dp) :: h_plus(0:ubound(q,dim=1)), h_minus(0:ubound(q,dim=1))
    integer :: i, nf_light, nf_heavy

    if (GetPdfRep(q) /= pdfr_Human) call wae_error('ConvMTMNew',&
         &'q is not in Human representation')

    PxQ = mtm%P_light * q

    nf_light = mtm%P_light%nf_int

    qbar_sum = sum(q(:,-nf_light:-1),dim=2)
    singlet  = sum(q(:,1:nf_light),dim=2)   ! not yet the singlet, just the sum of light quarks
    valence  = singlet - qbar_sum
    singlet  = singlet + qbar_sum

    h_plus  = mtm%PShg * q(:,0) + mtm%PShq * singlet    
    h_minus = mtm%NShV * valence
    Pxq(:,  nf_light+1 ) = half * (h_plus + h_minus)
    Pxq(:,-(nf_light+1)) = half * (h_plus - h_minus)
  end function ConvMTMNew

  !----------------------------------------------
  !! mutliply the mtm by fact
  subroutine Multiply_mtm(mtm, fact)
    type(mass_threshold_mat), intent(inout) :: mtm
    real(dp),                     intent(in)    :: fact
    call Multiply(mtm%P_light, fact)
    call Multiply(mtm%PShg, fact)
    call Multiply(mtm%PShq, fact)
    call Multiply(mtm%NShV, fact)
  end subroutine Multiply_mtm

  !----------------------------------------------
  !! apply mtm += factor * mtm_to_add (factor is optional)
  subroutine AddWithCoeff_mtm_mtm(mtm, mtm_to_add, factor)
    type(mass_threshold_mat), intent(inout) :: mtm
    type(mass_threshold_mat), intent(in)    :: mtm_to_add
    real(dp),                     intent(in), optional :: factor

    call AddWithCoeff(mtm%P_light, mtm_to_add%P_light, factor)
    call AddWithCoeff(mtm%PShg, mtm_to_add%PShg, factor)
    call AddWithCoeff(mtm%PShq, mtm_to_add%PShq, factor)
    call AddWithCoeff(mtm%NShV, mtm_to_add%NShV, factor)
  end subroutine AddWithCoeff_mtm_mtm

  !----------------------------------------------
  !! apply mtm += factor * sm_to_add (factor is optional)
  !!
  !! only works if sm_to_add has either nf_int = mtm%P_light%nf_int
  !! or nf_int = mtm%P_light%nf_int + 1
  subroutine AddWithCoeff_mtm_splitmat(mtm, sm_to_add, factor)
    use warnings_and_errors
    use hoppet_to_string, only : to_string
    type(mass_threshold_mat), intent(inout) :: mtm
    type(split_mat),              intent(in)    :: sm_to_add
    real(dp),                     intent(in), optional :: factor
    type(mass_threshold_mat) :: mtm_to_add

    if (sm_to_add%nf_int == mtm%P_light%nf_int) then
      call AddWithCoeff(mtm%P_light, sm_to_add, factor)
    else if (sm_to_add%nf_int == mtm%P_light%nf_int+1) then
      call InitMTM(mtm_to_add, sm_to_add)
      call AddWithCoeff(mtm, mtm_to_add, factor)
      call Delete(mtm_to_add)
    else
      call wae_error('AddWithCoeff_mtm_splitmat: nf_int mismatch between mtm (nf_light='&
                     // to_string(mtm%P_light%nf_int) //&
                      ') and split_mat to add (nf_int=' // to_string(sm_to_add%nf_int) // ')')
    end if
  end subroutine AddWithCoeff_mtm_splitmat

  !-------------------------------------------------------
  !! Set MTM_A = MTM_B .conv. P_C, on condition that
  !! (MTM_B%P_light%nf_int == P_C%nf_int)
  subroutine SetToConvolution_mtm_sm(MTM_A, MTM_B, P_C)
    use warnings_and_errors
    use hoppet_to_string
    type(mass_threshold_mat), intent(inout) :: MTM_A
    type(mass_threshold_mat), intent(in)    :: MTM_B
    type(split_mat),              intent(in)    :: P_C
    !---------------------------------------------
    type(grid_conv) :: qg_gq, gq_qg
    integer :: nf_light_int

    if (MTM_B%P_light%nf_int /= P_C%nf_int) &
         call wae_error('SetToConvolution_mtm_sm: nf_int mismatch between MTM_B (P_light%nf_int='&
                       // to_string(MTM_B%P_light%nf_int) // ') and P_C (nf_int=' &
                       // to_string(P_C%nf_int) // ')')

    call SetToConvolution(MTM_A%P_light, MTM_B%P_light, P_C)
    call SetToConvolution(MTM_A%NShV, MTM_B%NShV, P_C%NS_V)

    ! A%hq = B%hq * C%qq + B%hg * C%gq 
    call SetToConvolution(MTM_A%PShq, MTM_B%PShq, P_C%qq)
    call SetToConvolution(qg_gq, MTM_B%PShg, P_C%gq)
    call AddWithCoeff    (MTM_A%PShq, qg_gq)

    ! A%hg = B%hg * C%gg + B%hq * C%qg 
    call SetToConvolution(MTM_A%PShg, MTM_B%PShg, P_C%gg)
    call SetToConvolution(gq_qg, MTM_B%PShq, P_C%qg)
    call AddWithCoeff    (MTM_A%PShg, gq_qg)

    call Delete(qg_gq)
    call Delete(gq_qg)
  end subroutine SetToConvolution_mtm_sm


  !-------------------------------------------------------
  !! Set MTM_A = P_B .conv. MTM_C, on condition that
  !! (P_B%P_light%nf_int == MTM_C%P_light%nf_int+1)
  subroutine SetToConvolution_sm_mtm(MTM_A, P_B, MTM_C)
    use assertions
    use hoppet_to_string
    use convolution
    type(mass_threshold_mat), intent(inout) :: MTM_A
    type(split_mat),              intent(in)    :: P_B
    type(mass_threshold_mat), intent(in)    :: MTM_C
    !---------------------------------------------
    type(grid_conv) :: tmp1, tmp2, P_B_PS_plus, P_B_PS_minus
    integer  :: nf_light_int, nf_heavy_int
    real(dp) :: nf_light_dp , nf_heavy_dp

    nf_light_int = assert_eq(P_B%nf_int-1, MTM_C%P_light%nf_int,&
                             "SetToConvolution_sm_mtm: nf_int mismatch between P_B and MTM_C")
    nf_heavy_int = nf_light_int + 1
    nf_light_dp = real(nf_light_int, dp)
    nf_heavy_dp = real(nf_heavy_int, dp)

    call cobj_InitSplitLinks(MTM_A%P_light)
    MTM_A%P_light%nf_int = nf_light_int

    ! set up the pure singlet terms
    call InitGridConv(P_B_PS_plus , P_B%qq) 
    call AddWithCoeff(P_B_PS_plus , P_B%NS_plus , -one)
    call InitGridConv(P_B_PS_minus, P_B%NS_V) 
    call AddWithCoeff(P_B_PS_minus, P_B%NS_minus , -one)

    ! Get two non-singlet components
    ! (here and below, see section 3 of 202-12-flavour-struct logbook)
    call SetToConvolution(MTM_A%P_light%NS_plus , P_B%NS_plus , MTM_C%P_light%NS_plus )
    call SetToConvolution(MTM_A%P_light%NS_minus, P_B%NS_minus, MTM_C%P_light%NS_minus)

    ! now get the singlet-plus
    call SetToConvolution(MTM_A%P_light%qq, P_B%NS_plus, MTM_C%P_light%qq) ! (P_+ * T_S+)
    call InitGridConv(tmp1, MTM_C%P_light%qq)   ! (T_S+ 
    call AddWithCoeff(tmp1, MTM_C%PShq)         !  + T_h+S+)
    call SetToConvolution(tmp2, P_B_PS_plus, tmp1) ! (P_S+ * (T_S+ + T_h+S+))
    call AddWithCoeff(MTM_A%P_light%qq, tmp2, nf_light_dp / nf_heavy_dp) ! += nl/nh * (P_S+ * (T_S+ + T_h+S+))
    call SetToConvolution(tmp2, P_B%qg, MTM_C%P_light%gq) ! (P_qg * T_gq)
    call AddWithCoeff(MTM_A%P_light%qq, tmp2, nf_light_dp / nf_heavy_dp) ! += nl/nh * (P_qg * T_gq)


    ! and then the singlet-minus (i.e. valence)
    call SetToConvolution(MTM_A%P_light%NS_V, P_B%NS_minus, MTM_C%P_light%NS_V) ! (P_- * T_S-)
    call InitGridConv(tmp1, MTM_C%P_light%NS_V)   ! (T_S- 
    call AddWithCoeff(tmp1, MTM_C%NShV)         !  + T_h-S-)
    call SetToConvolution(tmp2, P_B_PS_minus, tmp1) ! (P_S- * (T_S- + T_h-S-))
    call AddWithCoeff(MTM_A%P_light%NS_V, tmp2, nf_light_dp / nf_heavy_dp) ! nl/nh * (P_S- * (T_S- + T_h-S-))

    ! the glue-glue piece
    call SetToConvolution(MTM_A%P_light%gg, P_B%gg, MTM_C%P_light%gg)
    call InitGridConv(tmp1, MTM_C%P_light%qg)
    call AddWithCoeff(tmp1, MTM_C%PShg)
    call SetToConvolution(tmp2, P_B%gq, tmp1)
    call AddWithCoeff(MTM_A%P_light%gg, tmp2)

    ! the gluon-quark
    call SetToConvolution(MTM_A%P_light%gq, P_B%gg, MTM_C%P_light%gq)
    call InitGridConv(tmp1, MTM_C%P_light%qq)
    call AddWithCoeff(tmp1, MTM_C%PShq)
    call SetToConvolution(tmp2, P_B%gq, tmp1)
    call AddWithCoeff(MTM_A%P_light%gq, tmp2)

    ! PShq piece
    ! 1/nh (P_PS+ * T_{S+})
    call SetToConvolution(MTM_A%PShq, P_B_PS_plus, MTM_C%P_light%qq)
    call Multiply(MTM_A%PShq, one / nf_heavy_dp)
    ! += (P_+ + 1/nh P_{PS+}) * T_{h+S_+}
    call InitGridConv(tmp1, P_B%NS_plus)
    call AddWithCoeff(tmp1, P_B_PS_plus, one / nf_heavy_dp)
    call SetToConvolution(tmp2, tmp1, MTM_C%PShq)
    call AddWithCoeff(MTM_A%PShq, tmp2)
    ! += 1/nh (P_qg * T_gq)
    call SetToConvolution(tmp2, P_B%qg, MTM_C%P_light%gq) ! (P_qg * T_gq)
    call AddWithCoeff(MTM_A%PShq, tmp2, one / nf_heavy_dp) ! += 1/nh * (P_qg * T_gq)

    ! NShV piece
    ! 1/nh (P_PS- * T_{S-})
    call SetToConvolution(MTM_A%NShV, P_B_PS_minus, MTM_C%P_light%NS_V)
    call Multiply(MTM_A%NShV, one / nf_heavy_dp)
    ! += (P_+ + 1/nh P_{PS-}) * T_{h-S_-}
    call InitGridConv(tmp1, P_B%NS_minus)
    call AddWithCoeff(tmp1, P_B_PS_minus, one / nf_heavy_dp)
    call SetToConvolution(tmp2, tmp1, MTM_C%NShV)
    call AddWithCoeff(MTM_A%NShV, tmp2)

    ! PShg piece
    ! (P_+ + 1/nh P_{PS+})*T_hg
    call InitGridConv(tmp1,P_B%NS_plus)
    call AddWithCoeff(tmp1,P_B_PS_plus, one/nf_heavy_dp)
    call SetToConvolution(MTM_A%PShg, tmp1, MTM_C%PShg)
    ! += 1/nh * (P_{S_+g}*T_{gg} + P_{PS+} * T_{S+g}
    call SetToConvolution(tmp1, P_B%qg, MTM_C%P_light%gg)
    call SetToConvolution(tmp2, P_B_PS_plus, MTM_C%P_light%qg)
    call AddWithCoeff(tmp1,tmp2)
    call AddWithCoeff(MTM_A%PShg, tmp1, one/nf_heavy_dp)
    
    ! qg piece
    ! T_{qg} = P_+ * T_{S_+g}
    call SetToConvolution(MTM_A%P_light%qg, P_B%NS_plus, MTM_C%P_light%qg)
    ! tmp2 = P_{PS+} * (T_{S_+g} + T_{h+g})
    call InitGridConv(tmp1, MTM_C%PShg)
    call AddWithCoeff(tmp1, MTM_C%P_light%qg)
    call SetToConvolution(tmp2, P_B_PS_plus, tmp1)
    ! tmp2 += P_{S_+g} * T_{gg}
    call SetToConvolution(tmp1, P_B%qg, MTM_C%P_light%gg)
    call AddWithCoeff(tmp2,tmp1)
    ! T_{qg} += nl/nh * tmp2
    call AddWithCoeff(MTM_A%P_light%qg, tmp2, nf_light_dp/nf_heavy_dp)

    call Delete(tmp1)
    call Delete(tmp2)
    call Delete(P_B_PS_plus)
    call Delete(P_B_PS_minus)
  end subroutine SetToConvolution_sm_mtm


  !! allocate and fill in probes to be used for deriving an 
  !! MTM object from a set of convolutions
  subroutine GetDerivedMTMProbes(grid, nf_light, probes)
    use pdf_representation
    type(grid_def), intent(in) :: grid
    integer,        intent(in) :: nf_light
    real(dp),       pointer    :: probes(:,:,:)
    real(dp), allocatable :: tmp(:,:)
    integer :: i

    call GetDerivedSplitMatProbes(grid, nf_light, probes)
    ! MTMs don't yet work in evolution representation, so convert
    ! the probes to human
    allocate(tmp(size(probes,dim=1), size(probes,dim=2)))
    do i = 1, size(probes,dim=3)
      call CopyEvlnPdfToHuman(nf_light, probes(:,:,i), tmp)
      probes(:,:,i) = tmp
    end do
  end subroutine GetDerivedMTMProbes

  !! Set the derived MTM from the probes, assuming
  !! mtm has been already allocated for the correct nf_light
  !!
  !! The probes are automatically deallocated
  subroutine SetDerivedMTM(mtm, probes)
    use pdf_representation
    use hoppet_probes, only : probe_norm
    type(mass_threshold_mat), intent(inout) :: mtm
    real(dp),   pointer,          intent(in) :: probes(:,:,:)
    real(dp), allocatable :: tmp(:,:)
    integer :: nprobes_1d,nprobes, il, ih, nh
    integer :: i

    nh = mtm%P_light%nf_int + 1

    ! convert probes back to evolution representation as needed by 
    ! SetDerivedSplitMat below
    allocate(tmp(size(probes,dim=1), size(probes,dim=2)))
    do i = 1, size(probes,dim=3)
      tmp = probes(:,:,i)
      call CopyHumanPdfToEvln(mtm%P_light%nf_int, tmp, probes(:,:,i))
    end do

    ! indices for probes here, and below, assume the standard structure
    ! used in GetDerivedSplitMatProbes
    nprobes   = size(probes,dim=3)
    nprobes_1d = nprobes/2

    il = 1; ih = nprobes_1d
    call SetDerivedConv_nodealloc(mtm%PShq, (probes(:,+nh,il:ih)+probes(:,-nh,il:ih))/probe_norm((iflv_sigma)))
    call SetDerivedConv_nodealloc(mtm%NShV, (probes(:,+nh,il:ih)-probes(:,-nh,il:ih))/probe_norm((iflv_V    )))
    il = nprobes_1d+1; ih = nprobes
    call SetDerivedConv_nodealloc(mtm%PShg, (probes(:,+nh,il:ih)+probes(:,-nh,il:ih))/probe_norm((iflv_g    )))

    call SetDerivedSplitMat(mtm%P_light, probes)
  end subroutine SetDerivedMTM
  !-------------------------------------------------------
  subroutine DelNewMTM(mtm)
    type(mass_threshold_mat), intent(inout) :: mtm
    call Delete(mtm%P_light)
    call Delete(mtm%PShq)
    call Delete(mtm%PShg)
    call Delete(mtm%NShV)
  end subroutine DelNewMTM

end module hoppet_mass_threshold_mat


!! Module with routines to initialise mass threshold matrix at NNLO/N3LO
module hoppet_init_mtm
  use types; use consts_dp; use splitting_functions
#ifdef HOPPET_ENABLE_N3LO_FORTRAN_MTM  
  use mass_thresholds_n3lo
#endif
  use qcd
  use hoppet_split_mat, only: cobj_InitSplitLinks
  use hoppet_mass_threshold_mat
  use dglap_choices
  use warnings_and_errors
  use pdf_representation
  use convolution
  use assertions

  implicit none

  private

  public :: InitMTMNNLO, InitMTMN3LO
  public :: InitMTMLibOME, InitMTMN3LOExactFortran

contains


  subroutine InitMTMNNLO(grid,MTM)
    use qcd
    type(grid_def),           intent(in)  :: grid
    type(mass_threshold_mat), intent(out) :: MTM
    integer :: nf_light
    !logical, parameter :: vogt_A2PShg = .false.
    !logical, parameter :: vogt_A2PShg = .true.

    ! our convention is that nf_int is the number of flavours including the
    ! heavy quark
    call cobj_InitSplitLinks(MTM%P_light)
    MTM%P_light%nf_int = nf_int - 1

    call InitGridConv(grid, MTM%PSHq, sf_A2PShq)
    select case (nnlo_nfthreshold_variant)
    case(nnlo_nfthreshold_param)
       call InitGridConv(grid, MTM%PSHg, sf_A2PShg_vogt)
    case(nnlo_nfthreshold_exact)
       call InitGridConv(grid, MTM%PSHg, sf_A2PShg)
    case default
       call wae_error('InitMTMNNLO', 'Unknown nnlo_threshold_variant',&
            &intval=nnlo_nfthreshold_variant)
    end select
   
    call InitGridConv(grid, MTM%P_light%NS_plus, sf_A2NSqq_H)
    call InitGridConv(MTM%P_light%NS_minus, MTM%P_light%NS_plus)  ! the plus and minus non-singlet pieces are the same
    call InitGridConv(MTM%P_light%NS_V,     MTM%P_light%NS_minus) ! the valence and minus non-singlet pieces are the same
    call InitGridConv(grid, MTM%P_light%gg, sf_A2Sgg_H)
    call InitGridConv(grid, MTM%P_light%gq, sf_A2Sgq_H)

    !! ! things specific to MSbar (CCN32-57)
    !! MTM%P_light%gg_extra_MSbar_delta = 8.0_dp/three * CF * TR
    !! ! here we need an extra -4CF * P_{q+qbar,g} = -8CF*P_{qg}
    !! call InitGridConv(grid, MTM%PShg_MSbar, sf_Pqg)
    !! call Multiply(MTM%PShg_MSbar, -8.0_dp*CF)
    !! call AddWithCoeff(MTM%PShg_MSbar, MTM%PSHg)

    ! pieces that are zero at NNLO
    call InitGridConv(MTM%P_light%qq, MTM%P_light%NS_plus)  ! qq = NS_plus + PSqq_H, but PSqq_H = 0 at NNLO

    call InitGridConv(grid, MTM%P_light%qg)
    call InitGridConv(grid, MTM%NShV)

    ! just store info that it is NNLO. For now it is obvious, but
    ! one day when we have NNNLO it may be useful, to indicate
    ! which structures exist and which do not. 
    ! (Mind you inexistent structures have yet to be implemented here...)
    !! MTM%loops = 3
    !! !-- no default value
    !! MTM%nf_int = nf_int
    !! !-- by default 
    !! MTM%masses_are_MSbar = .false.
  end subroutine InitMTMNNLO

  !! Initialise an N3LO MTM object based on the choice in 
  !! the global variable n3lo_nfthreshold
  subroutine InitMTMN3LO(grid,MTM_N3LO)
    type(grid_def),           intent(in)  :: grid
    type(mass_threshold_mat), intent(out) :: MTM_N3LO
    integer, save :: nwarn_n3lo_nfthreshold = 1

    if(n3lo_nfthreshold .eq. n3lo_nfthreshold_libOME_2510) then
      call InitMTMLibOME(grid,MTM_N3LO,nloop=4, LM=zero, include_NShV=.false.)

    else if(n3lo_nfthreshold .eq. n3lo_nfthreshold_libOME_2512) then
      call InitMTMLibOME(grid,MTM_N3LO,nloop=4, LM=zero, include_NShV=.true.)

    else if(n3lo_nfthreshold .eq. n3lo_nfthreshold_exact_fortran) then
      call InitMTMN3LOExactFortran(grid,MTM_N3LO)

    else if (n3lo_nfthreshold .eq. n3lo_nfthreshold_off) then
      call InitMTMNNLO(grid,MTM_N3LO) !! dummy initialisation to NNLO, just to get started
      !! MTM_N3LO%loops = 4
      call Multiply(MTM_N3LO, zero)

      call wae_warn(nwarn_n3lo_nfthreshold, 'InitMTMN3LO: n3lo_nfthreshold = off not recommended for nloop=4')

    else
      call wae_error('InitMTMN3LO', 'Unknown n3lo_threshold_variant',&
           &intval=n3lo_nfthreshold)
    end if

  end subroutine InitMTMN3LO


  !! Initialise an MTM object using the libome routines of arXiv:2510.02175
  !!
  !! \param grid    The grid to use
  !! \param MTM     The MTM object to initialise
  !! \param nloop   The perturbative (hoppet) nloops; libOME order = nloop-1 = # of powers of alpha_s
  !! \param LM      log(m^2/mu^2) where m is the heavy quark mass and mu the factorisation scale
  !!
  !! The number of light flavours is the global nf_int-1
  !!
  !! The expected relative accuracy on the underlying libOME is at least
  !! as good as (2048 * epsilon(1.0_dp))~4.5e-13,  which will usually be
  !! significantly better than the standard global integration precision
  !! (convolution's DefaultConvolutionEps(), which defaults to 1e-7)
  subroutine InitMTMLibOME(grid, MTM, nloop, LM, include_NShV)
    use iso_c_binding, only: c_double, c_int
    use qcd, only: nf_int
    use hoppet_libome_fortran
    type(grid_def),           intent(in)  :: grid
    type(mass_threshold_mat), intent(out) :: MTM
    integer,                  intent(in)  :: nloop
    real(dp),  optional,      intent(in)  :: LM
    logical,   optional,      intent(in)  :: include_NShV    
    !-----------
    integer, save :: nwarn_NShV_zero = 2
    real(c_double) :: LM_c, nf_light_d
    integer(c_int) :: order_c
    logical, save :: first_time = .true.

    if (first_time) then
      first_time = .false.
      write(6,'(a)') 'Initialisation of Mass Threshold Matrices with code from libOME'
      write(6,'(a)') '*     J. Ablinger, A. Behring, J. Blümlein, A. De Freitas,'
      write(6,'(a)') '*     A. von Manteuffel, C. Schneider, and K. Schönwald,'
      write(6,'(a)') '*     "The Single-Mass Variable Flavor Number Scheme at Three-Loop Order",'
      write(6,'(a)') '*     arXiv:2510.02175 (DESY 24-037), 2512.13508, https://gitlab.com/libome/libome'
    end if

    ! our convention is that nf_int is the number of flavours including the
    ! heavy quark
    call cobj_InitSplitLinks(MTM%P_light)
    MTM%P_light%nf_int = nf_int - 1

    nf_light_d = real(nf_int-1,c_double)
    if (present(LM)) then
      LM_c = real(LM,c_double)
    else
      LM_c = zero
    end if
    order_c = nloop - 1

    call InitGridConv(grid, MTM%PShq            , conv_OME(AQqPS_ptr     , order=order_c, nf_light=nf_light_d, LM=LM_c))
    call InitGridConv(grid, MTM%PShg            , conv_OME(AQg_ptr       , order=order_c, nf_light=nf_light_d, LM=LM_c))
    call InitGridConv(grid, MTM%P_light%NS_plus , conv_OME(AqqQNSEven_ptr, order=order_c, nf_light=nf_light_d, LM=LM_c))
    call InitGridConv(grid, MTM%P_light%NS_minus, conv_OME(AqqQNSOdd_ptr , order=order_c, nf_light=nf_light_d, LM=LM_c))
    call InitGridConv(grid, MTM%P_light%gg      , conv_OME(AggQ_ptr      , order=order_c, nf_light=nf_light_d, LM=LM_c))
    call InitGridConv(grid, MTM%P_light%gq      , conv_OME(AgqQ_ptr      , order=order_c, nf_light=nf_light_d, LM=LM_c))
    call InitGridConv(grid, MTM%P_light%qg      , conv_OME(AqgQ_ptr      , order=order_c, nf_light=nf_light_d, LM=LM_c))
    call InitGridConv(grid, MTM%P_light%qq      , conv_OME(AqqQPS_ptr    , order=order_c, nf_light=nf_light_d, LM=LM_c))
    call AddWithCoeff(MTM%P_light%qq, MTM%P_light%NS_plus)  ! qq = NS_plus + PSqq_H

    if (default_or_opt(.true., include_NShV)) then
      call InitGridConv(grid, MTM%NShV  , conv_OME(AQqPSs_ptr   , order=order_c, nf_light=nf_light_d, LM=LM_c))
    else
      if (nloop >= 4) call wae_warn(nwarn_NShV_zero, &
                                    'InitMTMLibOME: initialising N3LO or higher with NShV component set to zero')
      call InitGridConv(grid, MTM%NShV)
    end if
    !call InitGridConv(grid, MTM%NShV)
    !call InitGridConv(grid, MTM%NShV  , conv_OME(AQqPSs_ptr     , order=order_c, nf_light=nf_light_d, LM=LM_c))

    ! the valence and minus non-singlet pieces are the same up to the orders (alphas^3) available in libOME
    call InitGridConv(MTM%P_light%NS_V,     MTM%P_light%NS_minus) 

    ! set the MSbar piece to zero -- MSbar thresholds not yet supported, but still needs
    ! to be there to avoid errors in the code.
    !call InitGridConv(grid, MTM%PShg_MSbar)
    !MTM%masses_are_MSbar = .false.
    !MTM%P_light%gg_extra_MSbar_delta = zero
    !MTM%loops  = nloop
    !MTM%nf_int = nf_int

  end subroutine InitMTMLibOME


  subroutine InitMTMN3LOExactFortran(grid,MTM_N3LO)
    type(grid_def),           intent(in)  :: grid
    type(mass_threshold_mat), intent(out) :: MTM_N3LO
    logical, save :: first_time = .true.
    integer, save :: nwarn_NShV = 1
    !! 103 uses two alphas points to separate out aalphas**3 from alphas**2, 
    !! 203 uses just one extreme one and is faster while giving consistent results
    integer, parameter :: ord3 = 203 
    
#ifndef HOPPET_ENABLE_N3LO_FORTRAN_MTM  
    call wae_error('InitMTMN3LO', 'n3lo_nfthreshold = exact_fortran requested but code not compiled with HOPPET_ENABLE_N3LO_FORTRAN_MTM')
#else
    if (first_time) then
        write(6,'(a)') 'Initialisation of Mass Threshold Matrices at N3LO. The paper'
        write(6,'(a)') '*     J. Ablinger, A. Behring, J. Blümlein, A. De Freitas,'
        write(6,'(a)') '*     A. von Manteuffel, C. Schneider, and K. Schönwald,'
        write(6,'(a)') '*     "The Single-Mass Variable Flavor Number Scheme at Three-Loop Order",'
        write(6,'(a)') '*     arXiv:2510.02175 (DESY 24-037) and 2512.13508'
        write(6,'(a)') '* shall be cited at any use. Corresponding code from J. Bluemlein, March 19 2024 and '
        write(6,'(a)') '* K. Schönwald, March 15 2024, including results from arXiv:2207.00027, arXiv:2403.00513'
        write(6,'(a)') '*'
        write(6,'(a)') "The initialisation could take up to a minute, depending on the machine."
        first_time = .false.
    end if

    ! our convention is that nf_int is the number of flavours including the
    ! heavy quark
    call cobj_InitSplitLinks(MTM_N3LO%P_light)
    MTM_N3LO%P_light%nf_int = nf_int - 1
    
    write(6,'(a,i1,a,i1)') 'Initialising N3LO mass (exact Fortran) threshold matrices for nf = ', nf_int-1,' -> ', nf_int
    call InitGridConv(grid, MTM_N3LO%PShq  , JB( PS, null(), null(), order=ord3))
    call AddWithCoeff(      MTM_N3LO%PShq  , JB(APS, null(), null(), order=  3))
    
    call InitGridConv(grid, MTM_N3LO%PShg  , JB( QG, null(), null(), order=ord3))
    call AddWithCoeff(      MTM_N3LO%PShg  , JB(AQGtotal, null(), null(), order=  3)) ! The full thing

    call InitGridConv(grid, MTM_N3LO%P_light%NS_plus, JB( NSREG,  NSPLU,  NSDEL, order=ord3))
    call AddWithCoeff(      MTM_N3LO%P_light%NS_plus, JB(ANSREG, ANSPLU, ANSDEL, order=  3))

    
    call InitGridConv(grid, MTM_N3LO%P_light%NS_minus, JB(OREG_znfasLL, OPLUS_znfasLL, ODEL_znfasLL, ord3))

    !
    call InitGridConv(grid, MTM_N3LO%P_light%gg, JB( GGREG,  GGPLU,  GGDEL, order=ord3))
    call AddWithCoeff(      MTM_N3LO%P_light%gg, JB(AGGREG, AGGPLU, AGGDEL, order=  3))
    !
    call InitGridConv(grid, MTM_N3LO%P_light%gq  , JB( GQ, null(), null(), order=ord3))
    call AddWithCoeff(      MTM_N3LO%P_light%gq  , JB(AGQ, null(), null(), order=  3))
    !
    call InitGridConv(grid, MTM_N3LO%P_light%qg  , JB(QGL, null(), null(), order=ord3))
    ! there is no AQGL, all is included in QGL

    !! there is no APSL, all is included in PSL
    call InitGridConv(grid, MTM_N3LO%P_light%qq, JB(PSL, null(), null(), order=ord3))
    call AddWithCoeff(MTM_N3LO%P_light%qq, MTM_N3LO%P_light%NS_plus)  ! qq = NS_plus + PSqq_H
    
    ! the NShV piece is not provided in the original fortran code, so we set it to zero here
    !call wae_warn(nwarn_NShV, 'InitMTMN3LOExactFortran: the N3LO exact Fortran MTM does not provide the NShV piece; setting it to zero')
    call InitGridConv(grid, MTM_N3LO%NShV, JB(AQqPSs3_znfasLL, null(), null(), order=3))

    ! the valence and minus non-singlet pieces are the same
    call InitGridConv(MTM_N3LO%P_light%NS_V, MTM_N3LO%P_light%NS_minus) 

    ! set the MSbar piece to zero -- MSbar thresholds not yet supported, but still needs
    ! to be there to avoid errors in the code.
    !call InitGridConv(grid, MTM_N3LO%PShg_MSbar)
    !MTM_N3LO%masses_are_MSbar = .false.
    !MTM_N3LO%P_light%gg_extra_MSbar_delta = zero
    !MTM_N3LO%loops = 4
    !MTM_N3LO%nf_int = nf_int

#endif
  end subroutine InitMTMN3LOExactFortran

end module hoppet_init_mtm
  
!! Module that provides a single point of usage for split_mat,
!! mass_threshold_mat and related routines
module dglap_objects
  use hoppet_split_mat
  use hoppet_mass_threshold_mat   
  use hoppet_init_mtm

end module dglap_objects