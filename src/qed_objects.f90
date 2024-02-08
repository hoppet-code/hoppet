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

  ! a NNLO splitting matrix (multiplies ( alpha/(2pi) )^2) 
  ! contains only Plq splitting!
  type qed_split_mat_nnlo 
     type(grid_conv) :: Plq_02
     integer         :: nu, nd, nl, nf
  end type qed_split_mat_nnlo

  ! a container for both LO and NLO splitting matrices
  type qed_split_mat
     type(qed_split_mat_lo)  :: lo
     type(qed_split_mat_nlo) :: nlo
     type(qed_split_mat_nnlo) :: nnlo
   ! only with F2003
   !contains
   !  procedure :: SetNf => QEDSplitMatSetNf
  end type qed_split_mat
  public :: qed_split_mat
  public :: InitQEDSplitMat
  public :: QEDSplitMatSetNf

  
  interface operator(*)
     module procedure conv_qed_lo, conv_qed_nlo, conv_qed_nnlo
  end interface operator(*)
  public :: operator(*)
  interface operator(.conv.)
     module procedure conv_qed_lo, conv_qed_nlo, conv_qed_nnlo
  end interface operator(.conv.)
  public :: operator(.conv.)

  interface Copy
     module procedure Copy_qed_split_mat_lo, Copy_qed_split_mat_nlo, Copy_qed_split_mat_nnlo, Copy_qed_split_mat
  end interface Copy
  public :: Copy
  
  interface Delete
     module procedure Delete_qed_split_mat_lo, Delete_qed_split_mat_nlo, Delete_qed_split_mat_nnlo, Delete_qed_split_mat
  end interface
  public :: Delete


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

  logical, public :: OffLeptonChargesFrom=.false.
  logical, public :: OffLeptonChargesTo=.false.
  logical, public :: OffNlGammaGamma=.false.
  
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

    call InitGridConv(grid, qed_split%nnlo%Plq_02, sf_Plq_02)

    ! set some vaguely sensible values (the splitting functions
    ! themselves don't depend on nf at the order we are dealing with)
    nl = 3
    nd = (nf_int+1) / 2
    nu = (nf_int  ) / 2
    call QEDSplitMatSetNf(qed_split, nl, nd, nu)
  end subroutine InitQEDSplitMat

  !----------------------------------------------------------------------
  subroutine QEDSplitMatSetNf(qed_sm, nl, nd, nu)
    type(qed_split_mat) :: qed_sm
    integer, intent(in) :: nl, nu, nd

    qed_sm%lo%nl = nl
    qed_sm%lo%nd = nd
    qed_sm%lo%nu = nu
    qed_sm%lo%nf = nd+nu
    
    qed_sm%nlo%nl = nl
    qed_sm%nlo%nd = nd
    qed_sm%nlo%nu = nu
    qed_sm%nlo%nf = nd+nu

    qed_sm%nnlo%nl = nl
    qed_sm%nnlo%nd = nd
    qed_sm%nnlo%nu = nu
    qed_sm%nnlo%nf = nd+nu    
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
    real(dp) :: chg2_toflv   (ncompmin:ncompmaxLeptons), sum_chg2_tofl
    

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
    if(OffLeptonChargesFrom) then
       chg2_fromflv(9:9+qed_lo%nl-1) = 0
    else
       chg2_fromflv(9:9+qed_lo%nl-1) = one       
    endif

    ! then set up a charge * multiplicity for branching to a given component
    chg2_toflv = zero
    ! include a factor of CA=NC for the quarks
    chg2_toflv(-6:6) = CA * chg2_fromflv(-6:6)

    ! include a factor of 2 for the leptons, because we sum leptons and
    ! anti-leptons
    chg2_toflv(9:9+qed_lo%nl-1) = two

    if(OffNlGammaGamma) then
       sum_chg2_tofl=sum(chg2_toflv(-6:6))
    else
       sum_chg2_tofl=sum(chg2_toflv)
    endif

    if(OffLeptonChargesTo) then
       chg2_toflv(9:9+qed_lo%nl-1) = 0
    endif

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
    gout(:,8) = gout(:,8) + (half * sum_chg2_tofl) * (qed_lo%Pyy_01 * gq(:,8))

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

  !----------------------------------------------------------------------
  function conv_qed_nnlo(qed_nnlo, gq) result(gout)
    type(qed_split_mat_nnlo), intent(in) :: qed_nnlo
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
    ! the charge, colour, etc. factor when branching to a flavour
    real(dp) :: chg2_toflv   (ncompmin:ncompmaxLeptons), sum_chg2_tofl

    ! for now we force this...
    if (GetPdfRep(gq(:,:ncompmax)) /= pdfr_Human) call wae_error('conv_qed_nnlo',&
         &'gq is not in "Human" format')
    if (ubound(gq,dim=2) < ncompmaxPhoton) call wae_error('conv_qed_nnlo',&
         &'gq appears not to have photon component')

    ! set up the squared electric charges of each of the components
    ! only use quarks here, because leptons don't contribute at O(alpha alpha_s)
    chg2_fromflv = zero
    chg2_fromflv(-1:(-2*qed_nnlo%nd+1):-2) = e_dn2
    chg2_fromflv( 1:( 2*qed_nnlo%nd-1):+2) = e_dn2
    chg2_fromflv(-2:(-2*qed_nnlo%nu ):-2) = e_up2
    chg2_fromflv( 2:( 2*qed_nnlo%nu ):+2) = e_up2

    ! then set up a charge * multiplicity for branching to a given component
    chg2_toflv = zero
    ! turn on only leptons
    ! include a factor of 2 for the leptons, because we sum leptons and
    ! anti-leptons
    chg2_toflv(9:9+qed_nnlo%nl-1) = two

    ! now get the full "photon"-emission power
    flvsum = zero
    do i = ncompmin, min(ubound(gq,dim=2), ncompmaxPhoton)
       if (chg2_fromflv(i) /= zero) flvsum = flvsum + chg2_fromflv(i) * gq(:,i)
    end do
    
    ! now set the result
    gout = zero
    ! we add only the splitting from a quark to a lepton        
    do i = ncompmin, min(ubound(gq,dim=2), ncompmaxLeptons)
       if (chg2_toflv(i) /= zero) then
          gout(:,i) = qed_nnlo%Plq_02 * flvsum * chg2_toflv(i)
       end if
    end do
        
  end function conv_qed_nnlo

  !----------------------------------------------------------------------
  subroutine Copy_qed_split_mat_lo(in,out)
    type(qed_split_mat_lo), intent(inout) :: in
    type(qed_split_mat_lo), intent(out)   :: out
    out%nu = in%nu
    out%nd = in%nd
    out%nl = in%nl
    out%nf = in%nf
    call InitGridConv(out%Pqq_01, in%Pqq_01)
    call InitGridConv(out%Pqy_01, in%Pqy_01)
    call InitGridConv(out%Pyq_01, in%Pyq_01)
    call InitGridConv(out%Pyy_01, in%Pyy_01)    
  end subroutine Copy_qed_split_mat_lo
  
  !----------------------------------------------------------------------
  subroutine Copy_qed_split_mat_nlo(in,out)
    type(qed_split_mat_nlo), intent(inout) :: in
    type(qed_split_mat_nlo), intent(out)   :: out
    out%nu = in%nu
    out%nd = in%nd
    out%nl = in%nl
    out%nf = in%nf

    call InitGridConv(out%Pqy_11    , in%Pqy_11    )
    call InitGridConv(out%Pyy_11    , in%Pyy_11    )
    call InitGridConv(out%Pgy_11    , in%Pgy_11    )
    call InitGridConv(out%Pqg_11    , in%Pqg_11    )
    call InitGridConv(out%Pyg_11    , in%Pyg_11    )
    call InitGridConv(out%Pgg_11    , in%Pgg_11    )
    call InitGridConv(out%PqqV_11   , in%PqqV_11   )
    call InitGridConv(out%PqqbarV_11, in%PqqbarV_11)
    call InitGridConv(out%Pgq_11    , in%Pgq_11    )
    call InitGridConv(out%Pyq_11    , in%Pyq_11    )
  end subroutine Copy_qed_split_mat_nlo

  !----------------------------------------------------------------------
  subroutine Copy_qed_split_mat_nnlo(in,out)
    type(qed_split_mat_nnlo), intent(inout) :: in
    type(qed_split_mat_nnlo), intent(out)   :: out
    out%nu = in%nu
    out%nd = in%nd
    out%nl = in%nl
    out%nf = in%nf
    call InitGridConv(out%Plq_02, in%Plq_02)
  end subroutine Copy_qed_split_mat_nnlo
  
  !----------------------------------------------------------------------
  subroutine Copy_qed_split_mat(in,out)
    type(qed_split_mat), intent(inout) :: in
    type(qed_split_mat), intent(out)   :: out
    call Copy(in%lo , out%lo )
    call Copy(in%nlo, out%nlo)
    call Copy(in%nnlo, out%nnlo)
  end subroutine Copy_qed_split_mat

  
  !----------------------------------------------------------------------
  subroutine Delete_qed_split_mat_lo(in)
   type(qed_split_mat_lo), intent(inout) :: in
   call Delete(in%Pqq_01)
   call Delete(in%Pqy_01)
   call Delete(in%Pyq_01)
   call Delete(in%Pyy_01)    
 end subroutine Delete_qed_split_mat_lo
 
 !----------------------------------------------------------------------
 subroutine Delete_qed_split_mat_nlo(in)
   type(qed_split_mat_nlo), intent(inout) :: in

   call Delete(in%Pqy_11    )
   call Delete(in%Pyy_11    )
   call Delete(in%Pgy_11    )
   call Delete(in%Pqg_11    )
   call Delete(in%Pyg_11    )
   call Delete(in%Pgg_11    )
   call Delete(in%PqqV_11   )
   call Delete(in%PqqbarV_11)
   call Delete(in%Pgq_11    )
   call Delete(in%Pyq_11    )
 end subroutine Delete_qed_split_mat_nlo

 !----------------------------------------------------------------------
 subroutine Delete_qed_split_mat_nnlo(in)
   type(qed_split_mat_nnlo), intent(inout) :: in
   call Delete(in%Plq_02)
 end subroutine Delete_qed_split_mat_nnlo
 
 !----------------------------------------------------------------------
 subroutine Delete_qed_split_mat(in)
   type(qed_split_mat), intent(inout) :: in
   call Delete(in%lo)
   call Delete(in%nlo)
   call Delete(in%nnlo)
 end subroutine Delete_qed_split_mat

end module qed_objects
