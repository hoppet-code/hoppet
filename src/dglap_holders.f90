!======================================================================
!! Module containing the definition of the dglap_holder type, together
!! with various related subroutines.
!!
!! The dglap_holder type is intended to contain all splitting
!! functions and coefficient functions related to PDF evolution and
!! F2/FL convolutions. As of (5/5/2006 -- actually for a long time
!! now), the coefficient function part is not fully implemented,
!! whereas all the splitting functions are.
!!
module dglap_holders
  use types; use consts_dp
  use dglap_objects
  use pdf_representation
  use assertions; use warnings_and_errors
  implicit none
  private

  !-- 
  !! Everything needed to calculate a cross section, given structure
  !! functions and alpha_s.
  !!
  ! nomenclature is: objects are leading order, unless there is a suffix n, 
  ! in which case they are NLO.
  type dglap_holder
     type(grid_def) :: grid
     !*********** for interpolated splitting functions ****************
     !* type(split_mat), pointer :: interpP(:, :)  ![alphas,nf]
     !* real(dp),        pointer :: alphas_values(:,:) ![alphas,nf]
     !* type(split_max), pointer :: currentP => null()
     !* logical                  :: using_interp
     !*********** for interpolated splitting functions ****************

     !-----------------------  nloop, nf  --------------
     type(split_mat), pointer :: allP(:,   :)  
     
     type(split_mat), pointer :: P_LO, P_NLO, P_NNLO, P_N3LO ! pointers to allP(:,nf)
     type(coeff_mat), pointer :: allC(:,:)
     type(coeff_mat), pointer :: C2, C2_1, CL_1

     !----------------------------------  nloop,nf
     type(mass_threshold_mat), pointer :: allMTM(:,:) => null()
     type(mass_threshold_mat), pointer :: MTM2 ! legacy, deprecated, pointer to MTM(3,nf)
     type(mass_threshold_mat), pointer :: MTM_NNLO, MTM_N3LO ! will be pointers to MTM(3,nf) and MTM(4,nf)
     logical    :: MTM2_exists=.false.
     integer    :: factscheme, nloop
     !--------------------------------  nf  ------------
     !type(pdf_rep), pointer :: all_prep(:)
     !type(pdf_rep), pointer :: prep
     integer                :: nf
  !contains
  !   procedure :: SetNf => SetNfDglapHolder
  end type dglap_holder
  
  !-- this is potentially makeshift?
  integer, parameter :: nC = 3

  public :: dglap_holder, InitDglapHolder, SetNfDglapHolder

  interface Delete
     module procedure holder_Delete
  end interface
  public :: Delete

contains

  !*********** for interpolated splitting functions ****************
  !!!! MEMBER FUNCTIONS !!!!!
  !!! InitDglapHolder_from_table(...)
  !!!
  !!! SetCurrentP(alphas [,nf])  -> sets currentP
  !!!
  !*********** for interpolated splitting functions ****************
  
  !-------------------------------------------------------
  ! Sets up eveything needed to calculate cross sections
  ! i.e. splitting functions and coefficient functions
  subroutine InitDglapHolder(grid, dh, factscheme, nloop, nflo, nfhi)
    use coefficient_functions; use convolution
    use qcd; use dglap_choices
    type(grid_def),     intent(in)    :: grid
    type(dglap_holder), intent(inout), target :: dh
    integer, optional,  intent(in)    :: factscheme
    integer, optional,  intent(in)    :: nloop
    integer, optional,  intent(in)    :: nflo, nfhi
    !-- holds temporary results
    type(grid_conv) :: dconv
    !-- holds all possible combinations of coefficient and splitting functions
    !   needed for DIS schemes
    type(grid_conv) :: Cq,Cg
    type(grid_conv) :: CqPqq, CqPqg, CqPgq, CqPgg
    type(grid_conv) :: CgPqq, CgPqg, CgPgq, CgPgg
    !-- more compactly written version of DIS scheme
    type(grid_conv) :: CC(iflv_g:iflv_sigma,iflv_g:iflv_sigma)
    type(grid_conv) :: tmp2d(iflv_g:iflv_sigma,iflv_g:iflv_sigma)
    logical :: newDIS = .true.
    !logical :: newDIS = .false.
    integer :: nfstore, nflcl
    integer, save :: nfragwarn = 5
    logical :: need_MTM

    dh%factscheme = default_or_opt(factscheme_default, factscheme)
    dh%nloop     = default_or_opt(2, nloop)
    if (dh%nloop > 4 .or. dh%nloop < 1) then
       call wae_error('InitDglapHolder: nloop must be between 1 and 4')
    end if
    dh%grid = grid


    if (present(nflo) .and. present(nfhi)) then
       allocate(dh%allP(dh%nloop,nflo:nfhi))
       !allocate(dh%all_prep(nflo:nfhi))
       allocate(dh%allC(nC,nflo:nfhi))
    else
       !-- otherwise, use whatever nf is currently set
       allocate(dh%allP(dh%nloop,nf_int:nf_int))
       !allocate(dh%all_prep(nf_int:nf_int))
       allocate(dh%allC(nC,nf_int:nf_int))
    end if

    ! mass steps not always supported. Work out if they are
    need_MTM = dh%nloop >= 3 .and. dh%factscheme == factscheme_MSbar .and. mass_steps_on &
            &.and. lbound(dh%allP,dim=2) /= ubound(dh%allP,dim=2)
    if (need_MTM) allocate(dh%allMTM(3:nloop,nflo:nfhi))

    !-- want to reset it at end
    nfstore = nf_int
    
    do nflcl = lbound(dh%allP,dim=2), ubound(dh%allP,dim=2)
       !-- this sets up all the pointers to make it look like a fixed-nf
       !   dh!
       call SetNfDglapHolder(dh, nflcl)

       !----- this is the remainder, essentially unchanged ----
       select case (dh%factscheme)
       case (factscheme_MSbar)
          call InitSplitMatLO (grid, dh%P_LO)
          if (dh%nloop >= 2) call InitSplitMatNLO(grid, dh%P_NLO, dh%factscheme)
          if (dh%nloop >= 3) then
             call InitSplitMatNNLO(grid, dh%P_NNLO, dh%factscheme)
             if (need_MTM) then
               if (nflcl == lbound(dh%allP,dim=2)+1) then
                  call InitMTMNNLO(grid,dh%MTM_NNLO)
               else if (nflcl > lbound(dh%allP,dim=2) + 1) then
                  ! try to be efficient by recycling MTMs from lower nf
                  ! since all that changes is nf_int variable inside MTM_NNLO
                  call InitMTM(dh%MTM_NNLO, dh%allMTM(3, lbound(dh%allP,dim=2) + 1))
                  call SetNfMTM(dh%MTM_NNLO, nflcl)
               end if
             end if
          end if
          if (dh%nloop >= 4) then
             call InitSplitMatN3LO(grid, dh%P_N3LO, dh%factscheme)
             ! we do not need this for the lower nf, because it takes
             ! us from nf-1->nf
             if (nflcl /= lbound(dh%allP,dim=2)) call InitMTMN3LO(grid,dh%MTM_N3LO)
             !-- do this once, and only if really needed
!             if (lbound(dh%allP,dim=2) /= ubound(dh%allP,dim=2) &
!                  &.and. mass_steps_on &
!                  &.and. nflcl == lbound(dh%allP,dim=2)) then
!                dh%MTM2_exists = .true.
!             end if
          end if

          call cobj_InitCoeff(grid, dh%C2)
          call cobj_InitCoeff(grid, dh%C2_1, cf_CgF2MSbar, cf_CqF2MSbar)
          call cobj_InitCoeff(grid, dh%CL_1,  cf_CgFL, cf_CqFL)
       case (factscheme_DIS)
          call InitSplitMatLO (grid, dh%P_LO)
          if (dh%nloop >= 2) &
               &call InitSplitMatNLO(grid, dh%P_NLO, factscheme_MSbar)
          if (dh%nloop >= 3) write(0,*) &
               &'DIS factorisation scheme not supported for 3 loops or more'
          call cobj_InitCoeff(grid, dh%C2)
          call cobj_InitCoeff(dh%C2_1, dh%C2, zero)
          call cobj_InitCoeff(grid, dh%CL_1,  cf_CgFL, cf_CqFL)
          
          !-- now convert MSbar splitting functions into DIS scheme -------
          ! See CCN21-6 (and also CCN17-61)
          
          if (newDIS) then
             !
             ! NB THIS VERSION OF THE DIS SCHEME IS NOT YET IN AGREEMENT
             !    WITH THE OLDER VERSION...
             !    ACTUALLY: THIS IS PROBABLY NO LONGER TRUE?
             !-- create the matrix C for use in 
             !   P_matrix -> P_matrix + [C,P]
             !           /  Cq   Cg \
             !     C ==  |          |   where Cg includes 2nf factor
             !           \ -Cq  -Cg /
             call InitGridConv(grid,CC(iflv_sigma,iflv_sigma),cf_CqF2MSbar)
             call InitGridConv(grid,CC(iflv_sigma,iflv_g),    cf_CgF2MSbar)
             call Multiply(CC(iflv_sigma,iflv_g), two*nf)
             call InitGridConv(CC(iflv_g,:),CC(iflv_sigma,:),-one)
             !-- now get a temporary to hold the commutator
             call AllocGridConv(grid,tmp2d)
             !-- work out the commutator: tmp2d=[C,P]
             !write(0,*) 'lb',lbound(dh%P_LO%singlet), ubound(dh%P_LO%singlet)
             !write(0,*) 'dh%P_LO%singlet(0,0)%nsub',dh%P_LO%singlet(0,0)%grid%nsub
             !write(0,*) 'dh%P_LO%singlet(0,1)%nsub',dh%P_LO%singlet(0,1)%grid%nsub
             !-----------------------------------------------
             ! putting the explicit singlet bounds somehow eliminates
             ! an ifc memory error, which was associated with 
             ! gcb in SetToCommutator obtaining ubound=(/220,1/)
             ! No understanding of origin of error and locally
             ! singlet has right bounds
             !
             ! result however seems to be wrong.
             !
             ! Message: be wary of intel here
             !call SetToCommutator(tmp2d,CC,dh%P_LO%singlet(iflv_g:iflv_sigma,iflv_g:iflv_sigma))
             call SetToCommutator(tmp2d,CC,dh%P_LO%singlet)
             !-- add it to P1
             call AddWithCoeff(dh%P_NLO%singlet,tmp2d)
             !-- add the beta function pieces as well
             call AddWithCoeff(dh%P_NLO%singlet,CC, -twopi_beta0)
             call AddWithCoeff(dh%P_NLO%NS_plus, CC(iflv_sigma,iflv_sigma),&
                  & -twopi_beta0)
             !-- quark number conservation remains OK because
             !   Cq has the om=0 moment equal to zero.
             call AddWithCoeff(dh%P_NLO%NS_minus, CC(iflv_sigma,iflv_sigma), &
                  &-twopi_beta0)
             call AddWithCoeff(dh%P_NLO%NS_V, CC(iflv_sigma,iflv_sigma), &
                  &-twopi_beta0)
             
             !-- clean up
             call Delete(CC)
             call Delete(tmp2d)
          else
             call InitGridConv(grid, Cq, cf_CqF2MSbar)
             call InitGridConv(grid, Cg, cf_CgF2MSbar)
             call Multiply(Cg, two*nf)
             !   where possible put the smoother distribution to the right
             !   (only makes a difference when pdfs are non zero at x=1?)
             call conv_ConvConv(CqPqq, Cq, dh%P_LO%qq)
             call conv_ConvConv(CqPqg, Cq, dh%P_LO%qg)
             call conv_ConvConv(CqPgq, Cq, dh%P_LO%gq)
             call conv_ConvConv(CqPgg, Cq, dh%P_LO%gg)
             call conv_ConvConv(CgPqq, dh%P_LO%qq, Cg)
             call conv_ConvConv(CgPqg, dh%P_LO%qg, Cg)
             call conv_ConvConv(CgPgq, dh%P_LO%gq, Cg)
             call conv_ConvConv(CgPgg, dh%P_LO%gg, Cg)
             
             !
             !   First deal with P_matrix -> P_matrix + [C,P]
             !           /  Cq   Cg \
             !     C ==  |          |   where Cg includes 2nf factor
             !           \ -Cq  -Cg /
             !
             !   Pqq -> Pqq + (Cg Pgq + Cq Pqg)
             call AddWithCoeff(dh%P_NLO%qq, CgPgq, one)
             call AddWithCoeff(dh%P_NLO%qq, CqPqg, one)
             !   Pqg -> Pqg + (Cq Pqg + Cg Pgg - Cg Pqq + Cg Pqg)
             call AddWithCoeff(dh%P_NLO%qg, CqPqg, one)
             call AddWithCoeff(dh%P_NLO%qg, CgPgg, one)
             call AddWithCoeff(dh%P_NLO%qg, CgPqq,-one)
             call AddWithCoeff(dh%P_NLO%qg, CgPqg, one)
             !   Pgq -> Pqg + (-Cq Pqq - Cg Pgq - Cq Pgq + Cq Pgg)
             call AddWithCoeff(dh%P_NLO%gq, CqPqq,-one)
             call AddWithCoeff(dh%P_NLO%gq, CgPgq,-one)
             call AddWithCoeff(dh%P_NLO%gq, CqPgq,-one)
             call AddWithCoeff(dh%P_NLO%gq, CqPgg, one)
             !   Pgg -> Pgg + (-Cq Pqg - Cg Pgq)
             call AddWithCoeff(dh%P_NLO%gg, CqPqg,-one)
             call AddWithCoeff(dh%P_NLO%gg, CgPgq,-one)
             !
             !   Now deal with P_matrix -> P_matrix - beta0 * C
             !                 P_+      -> P_+      - beta0 * C_q
             call AddWithCoeff(dh%P_NLO%qq, Cq, -twopi_beta0)
             call AddWithCoeff(dh%P_NLO%qg, Cg, -twopi_beta0)
             call AddWithCoeff(dh%P_NLO%gq, Cq, +twopi_beta0)
             call AddWithCoeff(dh%P_NLO%gg, Cg, +twopi_beta0)
             call AddWithCoeff(dh%P_NLO%NS_plus, Cq, -twopi_beta0)
             !-- quark number conservation remains OK because
             !   Cq has the om=0 moment equal to zero.
             call AddWithCoeff(dh%P_NLO%NS_minus, Cq, -twopi_beta0)
             call AddWithCoeff(dh%P_NLO%NS_V, Cq, -twopi_beta0)
             !
             !   tidy up
             !write(0,*) 'Hey:',dh%P_LO%singlet(0,0)%subgc(1)%conv(0:3,1)
             !write(0,*) 'Hey:',dh%P_LO%singlet(0,1)%subgc(1)%conv(0:3,1)
             !write(0,*) 'Hey:',dh%P_LO%singlet(1,0)%subgc(1)%conv(0:3,1)
             !write(0,*) 'Hey:',dh%P_LO%singlet(1,1)%subgc(1)%conv(0:3,1)
             !write(0,*) 'Hey:',               Cq%subgc(1)%conv(0:3,1)
             !write(0,*) 'Hey:',               Cg%subgc(1)%conv(0:3,1)
             
             call Delete(Cq)
             call Delete(Cg)
             call Delete(CqPqq)
             call Delete(CqPqg)
             call Delete(CqPgq)
             call Delete(CqPgg)
             call Delete(CgPqq)
             call Delete(CgPqg)
             call Delete(CgPgq)
             call Delete(CgPgg)
             write(0,*) 'result:',dh%P_NLO%singlet(iflv_sigma,iflv_sigma)%conv(10,1)
          end if
          
       case (factscheme_PolMSbar)
          write(0,*) "SETTING UP POLARIZED EVOLUTION"
          call InitSplitMatPolLO (grid, dh%P_LO)
          if (dh%nloop >= 2) call InitSplitMatPolNLO(grid, &
               &dh%P_NLO, dh%factscheme)
          if (dh%nloop >= 3) call wae_error('InitDglapHolder',&
               &'nloop >= 3 not supported for polarized case')

       case (factscheme_FragMSbar)
          if (dh%nloop >= 3) call wae_error('InitDglapHolder',&
               &'nloop >= 3 not supported for fragmentation case')

          ! recall that equations for timelike evolution (fragmentation
          ! functions) are different from those for the spacelike case,
          ! as given for example in Owens, Phys.Lett 76B (1978) 80 and
          ! so read
          !
          ! d/dlnQ^2 Dq = as/2pi * (      Pqq Dq + Pgq Dg)
          ! d/dlnQ^2 Dg = as/2pi * (sum_q Pqg Dq + Pgg Dg)
          !
          ! those the actual Pij are identical to before; so we've transposed
          ! the matrix and for the singlet formulation also reshuffled
          ! some 2nf factors
          !
          ! start off with default space-like matrix
          call InitSplitMatLO (grid, dh%P_LO)
          ! get 2nf factors correct
          call Multiply(dh%P_LO%qg, one / (2*nflcl))
          call Multiply(dh%P_LO%gq, one * 2*nflcl)
          ! exchange qg and gq entries
          dconv = dh%P_LO%qg       ! tmp
          dh%P_LO%qg = dh%P_LO%gq
          dh%P_LO%gq = dconv

          if (dh%nloop >= 2) then
             call InitSplitMatTimeNLO (grid, dh%P_NLO)
             call wae_warn(nfragwarn, &
                  & 'DANGER: nloop=2 fragmentation flavour thresholds not implemented')
          end if
          
       case default
          write(0,*) 'factorisation scheme ',dh%factscheme,&
               &' is not currently supported' 
          stop
       end select
       
       !-- used for converting to and from the evolution representation
       !dh%prep = DefaultEvlnRep(nf_int)
!!$       dh%prep%nf    = nf_int
!!$       !dh%prep%ibase = nf_int
!!$       dh%prep%ibase = 1
    end do

    

    !-- be "clean", put nf back to where it was ----------
    if (nfstore <= ubound(dh%allP,2) .and. nfstore >= lbound(dh%allP,2)) then
       ! set it in the dglap holder as well if valid
       call SetNfDglapHolder(dh,nfstore)
    else
       ! set it just in the global qcd module, because it is not
       ! valid here...
       call qcd_SetNf(nfstore)
    end if
    
  end subroutine InitDglapHolder


  !--------------------------------------------------------------
  ! set up all pointers so that it looks like an old holder in the
  ! good-old fixed-nf days; it also sets the global nf.
  ! (NB: perhaps that is not too good?)
  !
  ! If quark_masses_are_MSbar is present, the the quark mass 
  ! scheme is also set for the mass-threshold-matrix
  subroutine SetNfDglapHolder(dh, nflcl, quark_masses_are_MSbar)
    use qcd
    !class(dglap_holder), intent(inout) :: dh
    type(dglap_holder), intent(inout) :: dh
    integer, intent(in) :: nflcl
    logical, optional, intent(in) :: quark_masses_are_MSbar

    if (nflcl < lbound(dh%allP,dim=2) .or. nflcl > ubound(dh%allP,dim=2)) then
      call wae_Error('SetNfDglapHolder: tried to set unsupported nf; requested nf was',intval=nflcl)
    end if 
    
    !-- want general nf to be consistent. Not really all that nice a 
    !   way of doing things, but cannot think of a better solution 
    !   given the current structure...
    call qcd_SetNf(nflcl)

    !-- set up links so that the remainder of the routine
    !   can stay unchanged
    dh%P_LO  => dh%allP(1,nflcl)
    if (dh%nloop >= 2) then
       dh%P_NLO => dh%allP(2,nflcl)
    else
       nullify(dh%P_NLO)
    end if
    
    if (dh%nloop >= 3) then
       dh%P_NNLO => dh%allP(3,nflcl)
       if (associated(dh%allMTM) .and. nflcl >= lbound(dh%allMTM,2)) then
         dh%MTM_NNLO => dh%allMTM(3,nflcl)
         dh%MTM2     => dh%allMTM(3,nflcl)
         if (present(quark_masses_are_MSbar)) call SetMassSchemeMTM(dh%MTM_NNLO, quark_masses_are_MSbar)
         dh%MTM2_exists = .true.
       else
         nullify(dh%MTM_NNLO)
         nullify(dh%MTM2    )
         dh%MTM2_exists = .false.
       end if
    else
      nullify(dh%P_NNLO)
      nullify(dh%MTM_NNLO)
      nullify(dh%MTM2)
      dh%MTM2_exists = .false.
   end if
    
    if (dh%nloop >= 4) then
       dh%P_N3LO => dh%allP(4,nflcl)
       if (associated(dh%allMTM) .and. nflcl >= lbound(dh%allMTM,2)) then
         dh%MTM_N3LO => dh%allMTM(4,nflcl)
         if (present(quark_masses_are_MSbar)) call SetMassSchemeMTM(dh%MTM_N3LO, quark_masses_are_MSbar)
       else
         nullify(dh%MTM_N3LO)
       end if
    else
       nullify(dh%P_N3LO)
       nullify(dh%MTM_N3LO)
    end if
    
    dh%C2 => dh%allC(1,nflcl)
    dh%C2_1 => dh%allC(2,nflcl)
    dh%CL_1 => dh%allC(3,nflcl)
    !dh%prep => dh%all_prep(nflcl)
    dh%nf = nflcl
  end subroutine SetNfDglapHolder

  !======================================================================
  !! Attempt to free up all memory associated with this holder
  subroutine holder_Delete(dh)
    type(dglap_holder), intent(inout) :: dh
    !--------------------------------------
    integer :: nflcl, iloop, i

    do nflcl = lbound(dh%allP,2), ubound(dh%allP,2)
       do i = 1, size(dh%allC(:,nflcl))
          call Delete(dh%allC(i,nflcl))
       end do
       do iloop = 1, dh%nloop
          call Delete(dh%allP(iloop,nflcl))
       end do
    end do
    
    if (associated(dh%allMTM)) then
      do nflcl = lbound(dh%allMTM,2), ubound(dh%allMTM,2)
         do iloop = lbound(dh%allMTM,1), ubound(dh%allMTM,1)
            call Delete(dh%allMTM(iloop,nflcl))
         end do
      end do
      deallocate(dh%allMTM)
    end if 

    deallocate(dh%allP)
    deallocate(dh%allC)
    !deallocate(dh%all_prep)
  end subroutine holder_Delete
  
  
end module dglap_holders
