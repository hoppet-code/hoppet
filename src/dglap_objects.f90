!======================================================================
! 
!! This module provide a type for matrices of splitting functions as
!! well as for coefficient functions and mass-threshold matrices. The
!! two latter ones are liable to evolve.
!!
!! The module also provides:
!!
!! - subroutines for manipulating the splitting matrices (similarly 
!!   to the subroutine for grid_conv objects)
!!
!! - subroutines for initialising the splitting matrices with the LO
!!   NLO and NNLO splitting functions
!!
!! - subroutines for initialising Polarised LO and NLO splitting matrices. 
!!
!! - mass thresholds, coefficient functions, 
!!
!! - tools for getting "derived" splitting matrices.
!!
!======================================================================
module dglap_objects
  use types; use consts_dp; use splitting_functions
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
     !-- These are the non-singlet pieces
     type(grid_conv)    :: NS_plus, NS_minus, NS_V
     !-- LO -> loops=1, etc...
     !integer            :: loops  ! [not needed anymore?]
     integer            :: nf_int
  end type split_mat
  public             :: split_mat

  !---------------------------------------------------------------
  ! for going from nf to nf+1. Currently contains pieces
  ! required for O(as^2), but not the full general structure.
  type mass_threshold_mat
     type(grid_conv) :: PShq, PShg
     type(grid_conv) :: NSqq_H, Sgg_H, Sgq_H
     ! pieces related to the case of thresholds at MSbar masses
     type(grid_conv) :: PShg_MSbar
     real(dp)        :: Sgg_H_extra_MSbar_delta
     ! LOOPS == 1+POWER OF AS2PI, NF_INT = nf including heavy quark.
     ! This is potentially adaptable.
     integer         :: loops, nf_int
     logical         :: masses_are_MSbar
  end type mass_threshold_mat
  public :: mass_threshold_mat

  !-----------------------------------------------------------------
  ! holds the components of a coefficient function
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
  public :: InitSplitMatNNLO
  public :: InitSplitMatPolLO, InitSplitMatPolNLO

  public :: InitSplitMat!, Delete_sm
  !public :: AddWithCoeff_sm, Multiply_sm
  !public :: cobj_PConv, cobj_PConv_1d

  public :: InitMTMNNLO, SetNfMTM !, cobj_ConvMTM  , cobj_DelMTM
  public :: SetMassSchemeMTM
  
  !-------- things for splitting functions --------------------
  interface cobj_InitCoeff
     module procedure cobj_InitCoeffLO, cobj_InitCoeffHO, cobj_InitCoeff_cf
  end interface
  public :: cobj_InitCoeff, cobj_InitCoeffLO, cobj_InitCoeffHO
  !public :: cobj_CConv
  !public :: cobj_AddCoeff, cobj_DelCoeff
  public :: GetDerivedSplitMatProbes, AllocSplitMat, SetDerivedSplitMat


  interface operator(*)
     module procedure cobj_PConv, cobj_PConv_1d, cobj_CConv, cobj_ConvMTM
  end interface
  public :: operator(*)
  interface operator(.conv.)
     module procedure cobj_PConv, cobj_PConv_1d, cobj_CConv, cobj_ConvMTM
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
     module procedure Delete_sm, cobj_DelMTM,  cobj_DelCoeff
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

    ! NO LONGER NECESSARY
!!$    !-- dummy to initialize Vogt routines (needed in exact cases 
!!$    !   for the A3 piece to be set up). Do it better later on if it works?
!!$    cc_piece = cc_real
!!$    dummy = sf_P2NSMinus(0.5_dp)
!!$    dummy = sf_P2gg(0.5_dp)

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
  !! multiplied by some factor). Memory allocation should be automatically
  !! handled by the subsiduary routines.
  subroutine InitSplitMat(P, Pin, factor)
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
    
  end subroutine InitSplitMat


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

    probes(:,iflv_V, 1:nprobes_1d) = probes_1d       ! NS_V
    probes(:,   2, 1:nprobes_1d) = probes_1d       ! NS+
    probes(:,  -2, 1:nprobes_1d) = probes_1d       ! NS-
    probes(:,iflv_sigma, 1:nprobes_1d) = probes_1d   ! Quark column of singlet

    ! gluon column of singlet
    probes(:,iflv_g, nprobes_1d+1:nprobes) = probes_1d       

    ! normally this would be part of housekeeping job of convolution
    ! module (SetDerivedConv), but since probes_1d is lost from 
    ! view at this point, we have to do it ourselves...
    deallocate(probes_1d)
  end subroutine GetDerivedSplitMatProbes

  
  !---------------------------------------------------------------------
  !! Given an allocated split_mat and the results of operating on the probes
  !! determine the resulting "derived" splitting function.
  subroutine SetDerivedSplitMat(P,probes)
    type(split_mat),  intent(inout) :: P
    real(dp),         pointer       :: probes(:,:,:)
    !-----------------------------------------
    integer :: nprobes_1d,nprobes, il, ih

    nprobes   = size(probes,dim=3)
    nprobes_1d = nprobes/2

    il = 1; ih = nprobes_1d
    call SetDerivedConv_nodealloc(P%NS_V,     probes(:,iflv_V,il:ih))
    call SetDerivedConv_nodealloc(P%NS_plus,  probes(:,2,il:ih))
    call SetDerivedConv_nodealloc(P%NS_minus, probes(:,-2,il:ih))
    call SetDerivedConv_nodealloc(P%gq,       probes(:,iflv_g,il:ih))
    call SetDerivedConv_nodealloc(P%qq,       probes(:,iflv_sigma,il:ih))
    
    il = nprobes_1d+1; ih = nprobes
    call SetDerivedConv_nodealloc(P%gg,       probes(:,iflv_g,il:ih))
    call SetDerivedConv_nodealloc(P%qg,       probes(:,iflv_sigma,il:ih))

    deallocate(probes)
  end subroutine SetDerivedSplitMat
  
  !======================================================================
  !         MASS THRESHOLDS
  !======================================================================
  ! Here we keep all things to do with mass thresholds

  subroutine InitMTMNNLO(grid,MTM)
    use qcd
    type(grid_def),           intent(in)  :: grid
    type(mass_threshold_mat), intent(out) :: MTM
    !logical, parameter :: vogt_A2PShg = .false.
    !logical, parameter :: vogt_A2PShg = .true.
    integer, save :: warn_param = 3

    call InitGridConv(grid, MTM%PSHq, sf_A2PShq)
    select case (nnlo_nfthreshold_variant)
    case(nnlo_nfthreshold_param)
       call InitGridConv(grid, MTM%PSHg, sf_A2PShg_vogt)
       call wae_warn(warn_param,'InitMTMNNLO:&
            & using parametrisation  (less accuracte) for A2PShg')
    case(nnlo_nfthreshold_exact)
       call InitGridConv(grid, MTM%PSHg, sf_A2PShg)
    case default
       call wae_error('InitMTMNNLO', 'Unknown nnlo_threshold_variant',&
            &intval=nnlo_nfthreshold_variant)
    end select
   
    call InitGridConv(grid, MTM%NSqq_H, sf_A2NSqq_H)
    call InitGridConv(grid, MTM%Sgg_H, sf_A2Sgg_H)
    call InitGridConv(grid, MTM%Sgq_H, sf_A2Sgq_H)

    ! things specific to MSbar (CCN32-57)
    MTM%Sgg_H_extra_MSbar_delta = 8.0_dp/three * CF * TR
    ! here we need an extra -4CF * P_{q+qbar,g} = -8CF*P_{qg}
    call InitGridConv(grid, MTM%PShg_MSbar, sf_Pqg)
    call Multiply(MTM%PShg_MSbar, -8.0_dp*CF)
    call AddWithCoeff(MTM%PShg_MSbar, MTM%PSHg)


    ! just store info that it is NNLO. For now it is obvious, but
    ! one day when we have NNNLO it may be useful, to indicate
    ! which structures exist and which do not. 
    ! (Mind you inexistent structures have yet to be implemented here...)
    MTM%loops = 3
    !-- no default value
    MTM%nf_int = 0
    !-- by default 
    MTM%masses_are_MSbar = .false.
  end subroutine InitMTMNNLO

  !---------------------------------------------------------------------
  ! want to be able to set nf, defined as number of flavours
  ! after heavy matching
  subroutine SetNfMTM(MTM,nf_lcl)
    type(mass_threshold_mat), intent(inout) :: MTM
    integer,                intent(in)    :: nf_lcl
    if (MTM%loops > 3) then
       call wae_Error('SetNfMTM: MTM had loops > 3; nf is probably fixed')
    end if
    MTM%nf_int = nf_lcl
  end subroutine SetNfMTM
  
  !---------------------------------------------------------------------
  ! set the quark mass scheme for the MTM matching
  subroutine SetMassSchemeMTM(MTM, masses_are_MSbar)
    type(mass_threshold_mat), intent(inout) :: MTM
    logical,                  intent(in)    :: masses_are_MSbar
    MTM%masses_are_MSbar = masses_are_MSbar
  end subroutine SetMassSchemeMTM
  

  !----------------------------------------------------------------------
  ! Returns the amount to be added to go from nf-1 to nf flavours.
  ! Will try some tests to make sure that nf is consistent?
  !
  ! PDFs are assumed to be in "HUMAN" REPRESENTATION.
  function cobj_ConvMTM(MTM,q) result(Pxq)
    type(mass_threshold_mat), intent(in) :: MTM
    real(dp),   intent(in) :: q(0:,ncompmin:)
    real(dp)               :: Pxq(0:ubound(q,dim=1),ncompmin:ncompmax)
    real(dp) :: singlet(0:ubound(q,dim=1))
    integer :: i, nf_light, nf_heavy

    !-- general sanity checks
    if (MTM%loops <=0 .or. MTM%nf_int <=0) call wae_error('cobj_ConvMTM:',&
         &'Mass threshold matrix is undefined')

    if (GetPdfRep(q) /= pdfr_Human) call wae_error('cobj_ConvMTM',&
         &'q is not in Human representation')

    nf_heavy = MTM%nf_int
    nf_light = nf_heavy - 1
    !write(0,*) 'Doing a MT convolution with nf_heavy =',nf_heavy

    if (ncomponents < nf_heavy) call wae_error('cobj_ConvMTM:',&
            &'ncomponents in representation is < nf in MTM.')
    if (iflv_g /= 0) call wae_Error('cobj_ConvMTM:','iflv_g/=0')
    
    !-- sanity check on nf
    !if (any(q(:,-nf_heavy)/=zero) .or. any(q(:,nf_heavy)/=zero)) &
    !     & call wae_error('cobj_ConvMTM:',&
    !     &'Distribution already has non-zero components at nf_heavy')
    
    singlet = sum(q(:,-nf_light:-1),dim=2) + sum(q(:,1:nf_light),dim=2)

    if (MTM%masses_are_MSbar) then
      Pxq(:,nf_heavy) = half*(&
           &(MTM%PShq .conv. singlet) + (MTM%PShg_MSbar .conv. q(:,iflv_g)) )
    else 
      Pxq(:,nf_heavy) = half*(&
           &(MTM%PShq .conv. singlet) + (MTM%PShg .conv. q(:,iflv_g)) )
    end if
    
    Pxq(:,-nf_heavy) = Pxq(:,nf_heavy)
    Pxq(:,iflv_g)     = (MTM%Sgq_H.conv. singlet) + (MTM%Sgg_H.conv. q(:,iflv_g))

    ! the extra delta-function piece in the MSbar scheme
    if (MTM%masses_are_MSbar) then
      Pxq(:,iflv_g) = Pxq(:,iflv_g) + MTM%Sgg_H_extra_MSbar_delta*q(:,iflv_g)
    end if
    
    
    do i = -ncomponents, ncomponents
       if (abs(i) > nf_heavy) then
          Pxq(:,i) = zero
       else if (i == iflv_g .or. abs(i) == nf_heavy) then
          cycle
       else
          Pxq(:,i) = MTM%NSqq_H .conv. q(:,i)
       end if
    end do

    call LabelPdfAsRep(Pxq,pdfr_Human)
  end function cobj_ConvMTM
  

  !---- usual cleaning up --------------------------------------------
  subroutine cobj_DelMTM(MTM)
    type(mass_threshold_mat), intent(inout) :: MTM
    call Delete(MTM%PSHq)
    call Delete(MTM%PSHg)
    call Delete(MTM%NSqq_H)
    call Delete(MTM%Sgg_H)
    call Delete(MTM%Sgq_H)
    MTM%loops = -1
    MTM%nf_int = -1
  end subroutine cobj_DelMTM


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
    
end module dglap_objects
!end module dglap_objects_hidden
