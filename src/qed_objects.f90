module qed_objects
  use types; use consts_dp
  use convolution; use qcd
  use pdf_representation
  use qed_coupling_module
  use warnings_and_errors
  implicit none
  private

  ! a leading-order splitting matrix (multiplies alpha/2pi)
  type qed_split_mat_lo
     type(grid_conv) :: Pqq_01, Pqy_01, Pyq_01, Pyy_01
     integer         :: nu, nd, nl, nf
  end type qed_split_mat_lo

  ! a NLO splitting matrix (multiplies (alpha alpha_s)/(2pi)^2)
  type qed_split_mat_nlo
     type(grid_conv) :: Pqy_11, Pyy_11, Pgy_11
     type(grid_conv) :: Pqg_11, Pyg_11, Pgg_11
     type(grid_conv) :: PqqV_11, PqqbarV_11, Pgq_11, Pyq_11
     integer         :: nu, nd, nl, nf
  end type qed_split_mat_nlo

  ! a container for both LO and NLO splitting matrices
  type qed_split_mat
     type(qed_split_mat_lo)  :: lo
     type(qed_split_mat_nlo) :: nlo
   contains
     procedure :: SetNf => QEDSplitMatSetNf
  end type qed_split_mat
  public :: qed_split_mat
  public :: InitQEDSplitMat
  public :: QEDSplitMatSetNf
  
  interface operator(*)
     module procedure :: conv_qed_lo, conv_qed_nlo
  end interface operator(*)
  public :: operator(*)

  ! 8 is photon, 9,10,11 are leptons (each is lepton + anti-lepton)
  integer, parameter, public :: ncompmaxLeptons = 11
  integer, parameter, public :: ncompmaxPhoton  = 8

  public :: AllocPDFWithPhoton, AllocPDFWithLeptons

  integer, parameter, public :: iflv_photon   =  8
  ! NB: convention here is that electron, muon and tau contain the sum
  !     of the particle and anti-particle distributions
  integer, parameter, public :: iflv_electron =  9
  integer, parameter, public :: iflv_muon     = 10
  integer, parameter, public :: iflv_tau      = 11
  
contains

  !------------------------------------------
  subroutine AllocPDFWithPhoton(grid, pdf)
    type(grid_def), intent(in) :: grid
    real(dp),       pointer    :: pdf(:,:)
    
    call AllocGridQuant(grid,  pdf, ncompmin, ncompmaxPhoton)
  end subroutine AllocPDFWithPhoton

  !------------------------------------------
  subroutine AllocPDFWithLeptons(grid, pdf)
    type(grid_def), intent(in) :: grid
    real(dp),       pointer    :: pdf(:,:)
    
    call AllocGridQuant(grid,  pdf, ncompmin, ncompmaxLeptons)
  end subroutine AllocPDFWithLeptons
  
    
  
  !------------------------------------------
  subroutine InitQEDSplitMat(grid, qed_split)
    use qed_splitting_functions
    type(grid_def),      intent(in)    :: grid
    type(qed_split_mat), intent(inout) :: qed_split
    integer :: nl, nd, nu
    real(dp) :: fy_over_fg

    call InitGridConv(grid, qed_split%lo%Pqq_01, sf_Pqq_01)
    call InitGridConv(grid, qed_split%lo%Pyq_01, sf_Pyq_01)
    call InitGridConv(grid, qed_split%lo%Pqy_01, sf_Pqy_01)
    call InitGridConv(grid, qed_split%lo%Pyy_01, sf_Pyy_01)

    call InitGridConv(grid, qed_split%nlo%Pqg_11, sf_Pqg_11)
    call InitGridConv(grid, qed_split%nlo%Pyg_11, sf_Pyg_11)
    call InitGridConv(grid, qed_split%nlo%Pgg_11, sf_Pgg_11)

    ! get the fy splitting functions from the fg ones 
    fy_over_fg = CF * CA / TR 
    call InitGridConv(qed_split%nlo%Pqy_11, qed_split%nlo%Pqg_11, fy_over_fg)
    call InitGridConv(qed_split%nlo%Pgy_11, qed_split%nlo%Pyg_11, fy_over_fg)
    call InitGridConv(qed_split%nlo%Pyy_11, qed_split%nlo%Pgg_11, fy_over_fg)

    call InitGridConv(grid, qed_split%nlo%PqqV_11, sf_PqqV_11)
    call InitGridConv(grid, qed_split%nlo%PqqbarV_11, sf_PqqbarV_11)
    call InitGridConv(grid, qed_split%nlo%Pgq_11, sf_Pgq_11)
    call InitGridConv(qed_split%nlo%Pyq_11, qed_split%nlo%Pgq_11) ! identical...

    ! set some vaguely sensible values (the splitting functions
    ! themselves don't depend on nf at the order we are dealing with)
    nl = 3
    nd = (nf_int+1) / 2
    nu = (nf_int  ) / 2
    call qed_split%SetNf(nl, nd, nu)
  end subroutine InitQEDSplitMat

  !----------------------------------------------------------------------
  subroutine QEDSplitMatSetNf(qed_sm, nl, nd, nu)
    class(qed_split_mat) :: qed_sm
    integer, intent(in) :: nl, nu, nd

    qed_sm%lo%nl = nl
    qed_sm%lo%nd = nd
    qed_sm%lo%nu = nu
    qed_sm%lo%nf = nd+nu
    
    qed_sm%nlo%nl = nl
    qed_sm%nlo%nd = nd
    qed_sm%nlo%nu = nu
    qed_sm%nlo%nf = nd+nu
  end subroutine QEDSplitMatSetNf
  

  !----------------------------------------------------------------------
  function conv_qed_lo(qed_lo, gq) result(gout)
    type(qed_split_mat_lo), intent(in) :: qed_lo
    real(dp),               intent(in) :: gq(0:, ncompmin:)
    real(dp)                           :: gout(0:ubound(gq,dim=1), ncompmin:ubound(gq,dim=2))
    !---------------------------------------
    real(dp) :: flvsum(0:ubound(gq,dim=1)), flvout(0:ubound(gq,dim=1))
    integer  :: i
    ! the charge, colour, etc. factor when branching from a flavour;
    ! we allow for the leptons here, even if not using them, to
    ! allow for photon branching to them...
    real(dp) :: chg2_fromflv(ncompmin:ncompmaxLeptons)
    ! the charge, colour, etc. factor when branching to a flavour
    real(dp) :: chg2_toflv   (ncompmin:ncompmaxLeptons)
    

    ! for now we force this...
    if (GetPdfRep(gq(:,:ncompmax)) /= pdfr_Human) call wae_error('conv_qed_lo',&
         &'gq is not in "Human" format')
    if (ubound(gq,dim=2) < ncompmaxPhoton) call wae_error('conv_qed_lo',&
         &'gq appears not to have photon component')
    
    ! set up the squared electric charges of each of the components
    chg2_fromflv = zero
    chg2_fromflv(-1:(-2*qed_lo%nd+1):-2) = e_dn2
    chg2_fromflv( 1:( 2*qed_lo%nd-1):+2) = e_dn2
    chg2_fromflv(-2:(-2*qed_lo%nu ):-2) = e_up2
    chg2_fromflv( 2:( 2*qed_lo%nu ):+2) = e_up2
    chg2_fromflv(9:9+qed_lo%nl-1) = one

    ! then set up a charge * multiplicity for branching to a given component
    chg2_toflv = zero
    ! include a factor of CA for the quarks
    chg2_toflv(-6:6) = CA * chg2_fromflv(-6:6)
    ! include a factor of 2 for the leptons, because we sum leptons and
    ! anti-leptons
    chg2_toflv(9:11) = chg2_fromflv(9:11) * two

    !write(6,*) nf_int
    !write(6,*) chg2_fromflv
    !write(6,*) chg2_toflv
    
    ! now get the full "photon"-emission power; if the PDF doesn't
    ! include leptons, then the branching from them is not included
    ! (even if nl /= 0 in the splitting matrix)
    flvsum = zero
    do i = ncompmin, min(ubound(gq,dim=2), ncompmaxLeptons)
       if (chg2_fromflv(i) /= zero) flvsum = flvsum + chg2_fromflv(i) * gq(:,i)
    end do
    
    ! now set the result
    gout = zero
    ! first the evolution of the photon from quarks and leptons
    gout(:,8) = qed_lo%Pyq_01 * flvsum
    ! Then the derivative of the photons associated with their "decay" to fermions.
    ! Include a factor of 1/2 because we sum explicitly over quarks and anti-quarks
    ! whereas the original normalisation from eqs.(21&22) of 1512.00612 involves
    ! only a sum over flavours, not anti-flavours
    gout(:,8) = gout(:,8) + (half * sum(chg2_toflv)) * (qed_lo%Pyy_01 * gq(:,8))

    ! now the evolution of the quarks and leptons:
    ! first calculate the generic splitting from a photon to a fermion
    ! without any charge or colour factors
    flvout    = qed_lo%Pqy_01 * gq(:,8)
    ! then add it in with appropriate charges
    do i = ncompmin, min(ubound(gq,dim=2), ncompmaxLeptons)
       if (chg2_toflv(i) /= zero) then
          gout(:,i) = chg2_toflv(i) * flvout
          gout(:,i) = gout(:,i) + chg2_fromflv(i) * (qed_lo%Pqq_01 * gq(:,i))
       end if
    end do
    
  end function conv_qed_lo
         
  
  !----------------------------------------------------------------------
  function conv_qed_nlo(qed_nlo, gq) result(gout)
    type(qed_split_mat_nlo), intent(in) :: qed_nlo
    real(dp),                intent(in) :: gq(0:, ncompmin:)
    real(dp)                            :: gout(0:ubound(gq,dim=1), ncompmin:ubound(gq,dim=2))
    !---------------------------------------
    real(dp) :: flvsum(0:ubound(gq,dim=1)), flvout(0:ubound(gq,dim=1))
    integer  :: i
    ! the charge, colour, etc. factor when branching from a flavour;
    ! we allow for the leptons here, even if not using them, to
    ! allow for photon branching to them...
    real(dp) :: chg2_fromflv(ncompmin:ncompmaxPhoton)
    ! the sum over charges (just flavours)
    real(dp) :: chg2sum

    ! for now we force this...
    if (GetPdfRep(gq(:,:ncompmax)) /= pdfr_Human) call wae_error('conv_qed_nlo',&
         &'gq is not in "Human" format')
    if (ubound(gq,dim=2) < ncompmaxPhoton) call wae_error('conv_qed_nlo',&
         &'gq appears not to have photon component')

    ! set up the squared electric charges of each of the components
    ! only use quarks here, because leptons don't contribute at O(alpha alpha_s)
    chg2_fromflv = zero
    chg2_fromflv(-1:(-2*qed_nlo%nd+1):-2) = e_dn2
    chg2_fromflv( 1:( 2*qed_nlo%nd-1):+2) = e_dn2
    chg2_fromflv(-2:(-2*qed_nlo%nu ):-2) = e_up2
    chg2_fromflv( 2:( 2*qed_nlo%nu ):+2) = e_up2
    chg2sum = half * sum(chg2_fromflv) 


    !write(6,*) nf_int
    !write(6,*) chg2_fromflv
    !write(6,*) chg2_toflv
    
    ! now get the full "photon"-emission power; if the PDF doesn't
    ! include leptons, then the branching from them is not included
    ! (even if nl /= 0 in the splitting matrix)
    flvsum = zero
    do i = ncompmin, min(ubound(gq,dim=2), ncompmaxPhoton)
       if (chg2_fromflv(i) /= zero) flvsum = flvsum + chg2_fromflv(i) * gq(:,i)
    end do
    
    ! now set the result
    gout = zero
    ! first the evolution of the gluon and photon from quarks
    gout(:,0) = qed_nlo%Pgq_11 * flvsum
    gout(:,8) = qed_nlo%Pyq_11 * flvsum
    ! Then the derivative of the photons associated with their "decay" to fermions.
    ! Include a factor of 1/2 because we sum explicitly over quarks and anti-quarks
    ! whereas the original normalisation from eqs.(21&22) of 1512.00612 involves
    ! only a sum over flavours, not anti-flavours
    gout(:,8) = gout(:,8) + chg2sum * (qed_nlo%Pyy_11 * gq(:,8) + qed_nlo%Pyg_11 * gq(:,0))
    gout(:,0) = gout(:,0) + chg2sum * (qed_nlo%Pgg_11 * gq(:,0) + qed_nlo%Pgy_11 * gq(:,8))

    ! now the evolution of the quarks:
    ! first calculate the generic splitting from a photon to a fermion
    ! without any charge or colour factors
    flvout    = (qed_nlo%Pqy_11 * gq(:,8) + qed_nlo%Pqg_11 * gq(:,0))
    ! then add it in with appropriate charges
    do i = ncompmin, min(ubound(gq,dim=2), ncompmaxPhoton)
       if (chg2_fromflv(i) /= zero) then
          gout(:,i) = chg2_fromflv(i) * flvout
          gout(:,i) = gout(:,i) + chg2_fromflv(i) * (&
               &qed_nlo%PqqV_11 * gq(:,i) + qed_nlo%PqqbarV_11 * gq(:,-i))
       end if
    end do
    
  end function conv_qed_nlo
  
end module qed_objects
