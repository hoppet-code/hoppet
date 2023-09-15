
module coefficient_functions_holder_internal
  use types; 
  use convolution_communicator
  use dglap_objects
  use coefficient_functions 
  implicit none
 
  type coef_holder
    integer        :: nloop
    integer        :: nf, nflo, nfhi
    real(dp)       :: C2LO,   CLLO,   C3LO
    type(split_mat), pointer :: C2NLO,  CLNLO,  C3NLO
    type(split_mat), pointer :: C2NNLO, CLNNLO, C3NNLO
    type(split_mat), pointer :: C2N3LO, CLN3LO, C3N3LO
    type(split_mat), pointer :: C2N3LO_fl11, CLN3LO_fl11
    type(split_mat), pointer :: C_NLO(:,:), C_NNLO(:,:), C_N3LO(:,:)
  end type coef_holder

end module coefficient_functions_holder_internal

! Module that contains the coefficient functions up to N3LO
!
! Useful references:
!
! - http://www.liv.ac.uk/~avogt/coeff.html
!   (this is the source of the xc*.f 2- and 3-loop coefficient functions)
!
! - hep-ph/9907472 and hep-ph/0006154, which contain the
!   parametrisations and some amount of definition of
!   what it is we're doing.
!
! - Original van Neerven and Zijlstra papers
module coefficient_functions_holder
  use types; use warnings_and_errors
  use consts_dp; use convolution_communicator
  use qcd
  use convolution
  use dglap_objects
  use coefficient_functions 
  use coefficient_functions_holder_internal
  implicit none
  !private
  
  public :: coef_holder, InitCoefHolder, SetNfCoefHolder
  
  real(dp), parameter  :: tiny = 1D-10
  
  !! to facilitate conditional compilation of exact coefficient functions
  !! (which add about 20s to the compilation time), we define an
  !! external function that handles the exact coefficient functions.
  !! If support is not compiled in, it should simply give an error.
  interface
     subroutine InitCoefHolderExact(grid, ch)
        use types; use coefficient_functions_holder_internal
        implicit none
        type(grid_def), intent(in) :: grid
        type(coef_holder), intent(inout), target :: ch
     end subroutine InitCoefHolderExact
  end interface

contains

  !----------------------------------------------------------------------
  ! Set up all the coefficient functions
  subroutine InitCoefHolder(grid, ch, nloop, use_exact_cf, nflo, nfhi)
    type(grid_def),    intent(in) :: grid
    type(coef_holder), intent(inout), target :: ch
    integer, optional, intent(in) :: nloop
    logical, optional, intent(in) :: use_exact_cf
    integer, optional, intent(in) :: nflo, nfhi
    !------------------------------
    integer :: nfstore, nf_lcl
    logical :: use_exact

    ! default settings
    ch%nloop  = 4
    use_exact = .false.
    ch%nf     = nf_int
    ch%nflo   = nf_int
    ch%nfhi   = nf_int

    ! if specified, use inputs
    if (present(nloop)) ch%nloop = nloop
    if (present(use_exact_cf)) use_exact = use_exact_cf
    if (present(nflo)) ch%nflo = nflo
    if (present(nfhi)) ch%nfhi = nfhi
    
    if (ch%nloop > 4 .or. ch%nloop < 1) then
       call wae_error('InitCoefHolder: nloop must be between 1 and 4')
    end if
    
    ! now allocate the arrays
    if (nloop.ge.2) then
       allocate(ch%C_NLO(3,ch%nflo:ch%nfhi))
    endif
    if (nloop.ge.3) then
       allocate(ch%C_NNLO(3,ch%nflo:ch%nfhi))
    endif
    if (nloop.ge.4) then
       allocate(ch%C_N3LO(5,ch%nflo:ch%nfhi))
    endif

    ! first initialise the LO coefficient "functions" (just a number, since delta-fn)
    ! they don't depend on nf, so set them up now
    ch%C2LO = one
    ch%CLLO = zero
    ch%C3LO = one

    ! store nf value so that we can reset it afterwards
    nfstore = nf_int
    
    ! loop over nf values and set up coefficient functions
    do nf_lcl = ch%nflo, ch%nfhi
       call SetNfCoefHolder(ch, nf_lcl)
       if (ch%nloop.ge.2) then
          call InitC2NLO(grid, ch%C2NLO)
          call InitCLNLO(grid, ch%CLNLO)
          call InitC3NLO(grid, ch%C3NLO)
       end if

       if (use_exact_cf) then
         call InitCoefHolderExact(grid, ch)
       else
         if (nloop.ge.3) then
            ! either with the exact functions
            !if (use_exact_cf) then
            !   call InitC2NNLO_e(grid, ch%C2NNLO)
            !   call InitCLNNLO_e(grid, ch%CLNNLO)
            !   call InitC3NNLO_e(grid, ch%C3NNLO)
            !else
               ! or the parametrised version
               call InitC2NNLO(grid, ch%C2NNLO)
               call InitCLNNLO(grid, ch%CLNNLO)
               call InitC3NNLO(grid, ch%C3NNLO) 
            !endif
         endif
         if (nloop.ge.4) then
            !if (use_exact_cf) then
            !   call InitC2N3LO_e(grid, ch%C2N3LO)
            !   call InitCLN3LO_e(grid, ch%CLN3LO)
            !   call InitC3N3LO_e(grid, ch%C3N3LO)
            !   ! including the fl11 terms for Z/photon exchanges
            !   call InitC2N3LO_fl11_e(grid, ch%C2N3LO_fl11)
            !   call InitCLN3LO_fl11_e(grid, ch%CLN3LO_fl11)
            !else
               call InitC2N3LO(grid, ch%C2N3LO)
               call InitCLN3LO(grid, ch%CLN3LO)
               call InitC3N3LO(grid, ch%C3N3LO)
               ! including the fl11 terms for Z/photon exchanges
               call InitC2N3LO_fl11(grid, ch%C2N3LO_fl11)
               call InitCLN3LO_fl11(grid, ch%CLN3LO_fl11)
            !endif
         endif
      end if
    end do

    ! reset nf value to initial value
    ! and set all the pointers accordingly
    call SetNfCoefHolder(ch, nfstore)
    
  end subroutine InitCoefHolder

  !----------------------------------------------------------------------
  ! Set the internal nf to specified value
  subroutine SetNfCoefHolder(ch, nflcl)
    type(coef_holder), intent(inout), target :: ch
    integer, intent(in) :: nflcl

    if ((nflcl.lt.ch%nflo) .or. (nflcl.gt.ch%nfhi)) then
       call wae_Error('SetNfCoefHolder: tried to set unsupported nf; request val was',intval=nflcl)
    end if
    
    ch%nf   = nflcl
    call qcd_SetNf(nflcl)

    ! now set the pointers
    if (ch%nloop.ge.2) then
       ch%C2NLO => ch%C_NLO(1,nflcl)
       ch%CLNLO => ch%C_NLO(2,nflcl)
       ch%C3NLO => ch%C_NLO(3,nflcl)
    endif
    if (ch%nloop.ge.3) then
       ch%C2NNLO => ch%C_NNLO(1,nflcl)
       ch%CLNNLO => ch%C_NNLO(2,nflcl)
       ch%C3NNLO => ch%C_NNLO(3,nflcl)
    endif
    if (ch%nloop.ge.4) then
       ch%C2N3LO => ch%C_N3LO(1,nflcl)
       ch%CLN3LO => ch%C_N3LO(2,nflcl)
       ch%C3N3LO => ch%C_N3LO(3,nflcl)
       ch%C2N3LO_fl11 => ch%C_N3LO(4,nflcl)
       ch%CLN3LO_fl11 => ch%C_N3LO(5,nflcl)
    endif
    
  end subroutine SetNfCoefHolder

  
  !======================================================================
  ! NLO coefficient functions
  !======================================================================

  !======================================================================
  ! NLO C2 coefficient function as a splitting matrix
  subroutine InitC2NLO(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%qq, cf_CqF2MSbar)
    call InitGridConv(grid, C%qg, cf_CgF2MSbar)
    call InitGridConv(grid, C%gg)  ! zero
    call InitGridConv(grid, C%gq)  ! zero

    !-- now fix up pieces so that they can be directly used as a matrix
    call Multiply(C%qg, 2*nf)

    ! all non-singlet pieces are equal to singlet ones
    call InitGridConv(C%NS_plus,  C%qq)
    call InitGridConv(C%NS_minus, C%qq)
    call InitGridConv(C%NS_V, C%NS_minus)
  end subroutine InitC2NLO
  
  !======================================================================
  ! NLO CL coefficient function as a splitting matrix
  subroutine InitCLNLO(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%qq, cf_CqFL)
    call InitGridConv(grid, C%qg, cf_CgFL)
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

    !-- now fix up pieces so that they can be directly used as a matrix
    call Multiply(C%qg, 2*nf)

    ! all non-singlet pieces are equal to singlet ones
    call InitGridConv(C%NS_plus,  C%qq)
    call InitGridConv(C%NS_minus, C%qq)
    call InitGridConv(C%NS_V, C%NS_minus)
  end subroutine InitCLNLO
  
  !======================================================================
  ! NLO C3 coefficient function as a splitting matrix
  subroutine InitC3NLO(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%qq, cf_CqF3)
    call InitGridConv(grid, C%qg)  ! zero
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

    !-- now fix up pieces so that they can be directly used as a matrix
    call Multiply(C%qg, 2*nf)

    ! all non-singlet pieces are equal to singlet ones
    call InitGridConv(C%NS_plus,  C%qq)
    call InitGridConv(C%NS_minus, C%qq)
    call InitGridConv(C%NS_V, C%NS_minus)

    ! The following line helps tests to see if it makes a difference
    ! if we set the singlet qq piece to zero.
    !
    ! Conclusion:
    ! - it makes no difference for Z
    ! - it makes a difference for W with an odd number of (light) flavours
    ! - but W with an odd number of flavours isn't really a consistent thing
    !
    ! Since we're going to be running with an even # of flavours for W, we
    ! should be good whatever we put here.
    !
    !call InitGridConv(grid, C%qq)  ! zero
  end subroutine InitC3NLO
  
  !======================================================================
  ! NNLO parametrised coefficient functions
  !======================================================================

  !======================================================================
  ! NNLO C2 coefficient function as a splitting matrix
  subroutine InitC2NNLO(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%NS_plus, cfNNLO_F2NS_plus)
    call InitGridConv(grid, C%NS_minus, cfNNLO_F2NS_minus)
    ! valence shouldn't come in to F2, so it shouldn't matter
    ! what we put here (explicitly verified for F2Z and F2Wp*F2Wm)
    call InitGridConv(C%NS_V, C%NS_minus) 
    !call InitGridConv(grid, C%NS_V) 

    ! singlet qq is NS+ + pure_singlet
    call InitGridConv(C%qq, C%NS_plus)    ! start with NS+
    call AddWithCoeff(C%qq, cfNNLO_F2PS)  ! add in pure singlet piece
    call InitGridConv(grid, C%qg, cfNNLO_F2G)   ! the part driven by the gluon
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

  end subroutine InitC2NNLO

  !======================================================================
  ! NNLO CL coefficient function as a splitting matrix
  subroutine InitCLNNLO(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%NS_plus, cfNNLO_FLNS_plus)
    call InitGridConv(grid, C%NS_minus, cfNNLO_FLNS_minus)
    ! valence shouldn't come in to F2, so it shouldn't matter
    ! what we put here (explicitly verified for F2Z and F2Wp*F2Wm)
    call InitGridConv(C%NS_V, C%NS_minus) 
    !call InitGridConv(grid, C%NS_V) 

    ! singlet qq is NS+ + pure_singlet
    call InitGridConv(C%qq, C%NS_plus)    ! start with NS+
    call AddWithCoeff(C%qq, cfNNLO_FLPS)  ! add in pure singlet piece
    call InitGridConv(grid, C%qg, cfNNLO_FLG)   ! the part driven by the gluon
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

  end subroutine InitCLNNLO

  !======================================================================
  ! NNLO C3 coefficient function as a splitting matrix
  subroutine InitC3NNLO(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    ! singlet piece should have no impact on Z case and for the W case
    ! has no impact for even numbers of flavours (odd numbers of light
    ! flavours don't make sense with W, except potentially inside
    ! certain loops). So we simply set it to zero
    call InitGridConv(grid, C%qq)  ! zero
    call InitGridConv(grid, C%qg)  ! zero
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

    !-- now fix up pieces so that they can be directly used as a matrix
    call Multiply(C%qg, 2*nf)

    ! all non-singlet pieces are equal to singlet ones
    call InitGridConv(grid, C%NS_plus, cfNNLO_F3NS_plus)
    call InitGridConv(grid, C%NS_minus, cfNNLO_F3NS_minus)
    ! NS_minus and NS_V are identical according to sentence before
    ! Eq.(3.15) of 1109.3717v2
    call InitGridConv(C%NS_V, C%NS_minus) 

  end subroutine InitC3NNLO
  
  
  !======================================================================
  ! N3LO parametrised coefficient functions
  !======================================================================

  !======================================================================
  ! N3LO C2 coefficient function as a splitting matrix
  subroutine InitC2N3LO(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%NS_plus, cfN3LO_F2NS_plus)
    call InitGridConv(grid, C%NS_minus, cfN3LO_F2NS_minus)
    ! valence shouldn't come in to F2, so it shouldn't matter
    ! what we put here (explicitly verified for F2Z and F2Wp*F2Wm)
    call InitGridConv(C%NS_V, C%NS_minus) 
    !call InitGridConv(grid, C%NS_V) 

    ! singlet qq is NS+ + pure_singlet
    call InitGridConv(C%qq, C%NS_plus)    ! start with NS+
    call AddWithCoeff(C%qq, cfN3LO_F2PS)  ! add in pure singlet piece
    call InitGridConv(grid, C%qg, cfN3LO_F2G)   ! the part driven by the gluon
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

  end subroutine InitC2N3LO

  !======================================================================
  ! N3LO C2 fl11 coefficient function as a splitting matrix
  ! Only the part that is different for Z boson exchange
  subroutine InitC2N3LO_fl11(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%NS_plus, cfN3LO_F2NS_plus_fl11)
    call InitGridConv(grid, C%NS_minus, cfN3LO_F2NS_minus_fl11)
    ! valence shouldn't come in to F2, so it shouldn't matter
    ! what we put here (explicitly verified for F2Z and F2Wp*F2Wm)
    call InitGridConv(C%NS_V, C%NS_minus) 
    !call InitGridConv(grid, C%NS_V) 

    ! singlet qq is NS+ + pure_singlet
    call InitGridConv(C%qq, C%NS_plus)    ! start with NS+
    call AddWithCoeff(C%qq, cfN3LO_F2PS_fl11)  ! add in pure singlet piece
    call InitGridConv(grid, C%qg, cfN3LO_F2G_fl11)   ! the part driven by the gluon
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

  end subroutine InitC2N3LO_fl11

  !======================================================================
  ! N3LO CL coefficient function as a splitting matrix
  subroutine InitCLN3LO(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%NS_plus, cfN3LO_FLNS_plus)
    call InitGridConv(grid, C%NS_minus, cfN3LO_FLNS_minus)
    ! valence shouldn't come in to F2, so it shouldn't matter
    ! what we put here (explicitly verified for F2Z and F2Wp*F2Wm)
    call InitGridConv(C%NS_V, C%NS_minus) 
    !call InitGridConv(grid, C%NS_V) 

    ! singlet qq is NS+ + pure_singlet
    call InitGridConv(C%qq, C%NS_plus)    ! start with NS+
    call AddWithCoeff(C%qq, cfN3LO_FLPS)  ! add in pure singlet piece
    call InitGridConv(grid, C%qg, cfN3LO_FLG)   ! the part driven by the gluon
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

  end subroutine InitCLN3LO


  !======================================================================
  ! N3LO CL coefficient function as a splitting matrix
  subroutine InitCLN3LO_fl11(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%NS_plus, cfN3LO_FLNS_plus_fl11)
    call InitGridConv(grid, C%NS_minus, cfN3LO_FLNS_minus_fl11)
    ! valence shouldn't come in to F2, so it shouldn't matter
    ! what we put here (explicitly verified for F2Z and F2Wp*F2Wm)
    call InitGridConv(C%NS_V, C%NS_minus) 
    !call InitGridConv(grid, C%NS_V) 

    ! singlet qq is NS+ + pure_singlet
    call InitGridConv(C%qq, C%NS_plus)    ! start with NS+
    call AddWithCoeff(C%qq, cfN3LO_FLPS_fl11)  ! add in pure singlet piece
    call InitGridConv(grid, C%qg, cfN3LO_FLG_fl11)   ! the part driven by the gluon
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

  end subroutine InitCLN3LO_fl11

  !======================================================================
  ! N3LO C3 coefficient function as a splitting matrix
  subroutine InitC3N3LO(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    ! singlet piece should have no impact on Z case and for the W case
    ! has no impact for even numbers of flavours (odd numbers of light
    ! flavours don't make sense with W, except potentially inside
    ! certain loops). So we simply set it to zero
    call InitGridConv(grid, C%qq)  ! zero
    call InitGridConv(grid, C%qg)  ! zero
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

    !-- now fix up pieces so that they can be directly used as a matrix
    call Multiply(C%qg, 2*nf)

    ! all non-singlet pieces are equal to singlet ones
    call InitGridConv(grid, C%NS_plus, cfN3LO_F3NS_plus)
    call InitGridConv(grid, C%NS_minus, cfN3LO_F3NS_minus)
    ! NS_minus and NS_V are not identical at N3LO
    call InitGridConv(grid, C%NS_V, cfN3LO_F3NS_val) 

  end subroutine InitC3N3LO
  

  !======================================================================
  ! this is the part that needs to be added to Hoppet's
  ! cf_CqF2MSbar in order to get the F3 coefficient function
  !
  ! It is taken from Eq.(155) of hep-ph/9912355, keeping in mind that
  ! the result there multiplies as/4pi (we multiply as/2pi, so include
  ! an extra factor of 1/2).
  function cf_CqF3minusF2(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lnx, ln1mx

    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = res - CF * (1+x)
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function cf_CqF3minusF2


  function cf_CqF3 (y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    
    res = cf_CqF3minusF2(y) + cf_CqF2MSbar(y)
    
  end function cf_CqF3


  !======================================================================
  !======================================================================
  ! NNLO PIECES (parametrised version)
  !======================================================================
  !======================================================================

  !----------------------------------------------------------------------
  ! This is the coefficient function for the sum F3(nu)+F3(nubar),
  ! corresponding to C3NSP-C3NSN in W. van Neerven's program.
  function cfNNLO_F3NS_minus(y) result(res)
    use xc3ns2p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = C3NM2A(x, nf_int) + C3NS2B(x, nf_int)
    case(cc_REALVIRT)
       res = C3NM2A(x, nf_int)
    case(cc_VIRT)
       res = - C3NS2B(x, nf_int)
    case(cc_DELTA)
       res = C3NM2C(zero, nf_int)
    end select
    
    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_F3NS_minus

  !----------------------------------------------------------------------
  ! This is the coefficient function for the sum F3(nu)-F3(nubar),
  ! corresponding to C3NSP-C3NSN in W. van Neerven's program.
  function cfNNLO_F3NS_plus(y) result(res)
    use xc3ns2p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = C3NP2A(x, nf_int) + C3NS2B(x, nf_int)
    case(cc_REALVIRT)
       res = C3NP2A(x, nf_int)
    case(cc_VIRT)
       res = - C3NS2B(x, nf_int)
    case(cc_DELTA)
       res = C3NP2C(zero, nf_int)
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_F3NS_plus
  

  !----------------------------------------------------------------------
  ! This is the "plus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the electromagnetic F2, and, in
  ! addition to the singlet, that involves (u+ubar) - (d+dbar) type
  ! contributions, which is precisely what a NS+ coefficient function
  ! multiplies.
  function cfNNLO_F2NS_plus(y) result(res)
    use xc2ns2p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = C2NN2A(x, nf_int) + C2NS2B(x, nf_int)
    case(cc_REALVIRT)
       res = C2NN2A(x, nf_int)
    case(cc_VIRT)
       res = - C2NS2B(x, nf_int)
    case(cc_DELTA)
       res = C2NN2C(zero, nf_int)
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_F2NS_plus

  !----------------------------------------------------------------------
  ! This is the "minus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the charged-current F2, and, in
  ! addition to the singlet, that involves (u-ubar) - (d-dbar) type
  ! contributions, which is precisely what a NS- coefficient function
  ! multiplies.
  function cfNNLO_F2NS_minus(y) result(res)
    use xc2ns2p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = C2NC2A(x, nf_int) + C2NS2B(x, nf_int)
    case(cc_REALVIRT)
       res = C2NC2A(x, nf_int)
    case(cc_VIRT)
       res = - C2NS2B(x, nf_int)
    case(cc_DELTA)
       res = C2NC2C(zero, nf_int)
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_F2NS_minus

  !----------------------------------------------------------------------
  ! Pure singlet piece of F2 NNLO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfNNLO_F2PS(y) result(res)
    use xc2sg2p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = C2S2A(x, nf_int)
    case(cc_REALVIRT)
       res = C2S2A(x, nf_int)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_F2PS

  !----------------------------------------------------------------------
  ! Pure singlet piece of F2 NNLO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfNNLO_F2G(y) result(res)
    use xc2sg2p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = C2G2A(x, nf_int)
    case(cc_REALVIRT)
       res = C2G2A(x, nf_int)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = C2G2C(zero, nf_int)
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_F2G

  !----------------------------------------------------------------------
  ! This is the "plus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the electromagnetic FL, and, in
  ! addition to the singlet, that involves (u+ubar) - (d+dbar) type
  ! contributions, which is precisely what a NS+ coefficient function
  ! multiplies.
  function cfNNLO_FLNS_plus(y) result(res)
    use xclns2p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = CLNN2A(x, nf_int)
    case(cc_REALVIRT)
       res = CLNN2A(x, nf_int)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = CLNN2C(zero)
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_FLNS_plus

  !----------------------------------------------------------------------
  ! This is the "minus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the charged-current FL, and, in
  ! addition to the singlet, that involves (u-ubar) - (d-dbar) type
  ! contributions, which is precisely what a NS- coefficient function
  ! multiplies.
  function cfNNLO_FLNS_minus(y) result(res)
    use xclns2p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = CLNC2A(x, nf_int)
    case(cc_REALVIRT)
       res = CLNC2A(x, nf_int)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = CLNC2C(zero)
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_FLNS_minus

  !----------------------------------------------------------------------
  ! Pure singlet piece of FL NNLO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfNNLO_FLPS(y) result(res)
    use xclsg2p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = CLS2A(x, nf_int)
    case(cc_REALVIRT)
       res = CLS2A(x, nf_int)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_FLPS

  !----------------------------------------------------------------------
  ! Pure singlet piece of FL NNLO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfNNLO_FLG(y) result(res)
    use xclsg2p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = CLG2A(x, nf_int)
    case(cc_REALVIRT)
       res = CLG2A(x, nf_int)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_FLG

  

  
  !======================================================================
  !======================================================================
  ! N3LO PIECES (parametrised)
  !======================================================================
  !======================================================================

  !----------------------------------------------------------------------
  ! This is the coefficient function for the sum F3(nu)+F3(nubar),
  ! corresponding to C3NSP-C3NSN in W. van Neerven's program.
  function cfN3LO_F3NS_val(y) result(res)
    use xc3ns3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    
    select case(cc_piece)
    case(cc_REAL)
       res = C3NM3A(x, -y, nf_int, 1) + C3NS3B(x, -y, nf_int)
    case(cc_REALVIRT)
       res = C3NM3A(x, -y, nf_int, 1)
    case(cc_VIRT)
       res = - C3NS3B(x, -y, nf_int)
    case(cc_DELTA)
       !res = C3NS3C(zero, nf_int)
       res = C3NM3C(zero, nf_int) ! FD: Why is this C3NM3C instead of C3NS3C as in the NNLO case ???
                                  !     and where is the diff piece between C3NM3C and C3NP3C ?
    end select
    
    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F3NS_val

  !----------------------------------------------------------------------
  ! This is the coefficient function for the sum F3(nu)+F3(nubar),
  ! corresponding to C3NSP-C3NSN in W. van Neerven's program.
  function cfN3LO_F3NS_minus(y) result(res)
    use xc3ns3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    
    select case(cc_piece)
    case(cc_REAL)
       res = C3NM3A(x, -y, nf_int, 0) + C3NS3B(x, -y, nf_int)
    case(cc_REALVIRT)
       res = C3NM3A(x, -y, nf_int, 0)
    case(cc_VIRT)
       res = - C3NS3B(x, -y, nf_int)
    case(cc_DELTA)
       !res = C3NS3C(zero, nf_int)
       res = C3NM3C(zero, nf_int) ! FD: Why is this C3NM3C instead of C3NS3C as in the NNLO case ???
                                  !     and where is the diff piece between C3NM3C and C3NP3C ?
    end select
    
    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F3NS_minus

  !----------------------------------------------------------------------
  ! This is the coefficient function for the sum F3(nu)-F3(nubar),
  ! corresponding to C3NSP-C3NSN in W. van Neerven's program.
  function cfN3LO_F3NS_plus(y) result(res)
    use xc3ns3p
    use xcdiff3pnew
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    
    select case(cc_piece)
    case(cc_REAL)
       res = C3NM3A(x, -y, nf_int, 0) + C3Q3DFP(x, -y, nf_int) + C3NS3B(x, -y, nf_int)
    case(cc_REALVIRT)
       res = C3NM3A(x, -y, nf_int, 0) + C3Q3DFP(x, -y, nf_int)
    case(cc_VIRT)
       res = - C3NS3B(x, -y, nf_int)
    case(cc_DELTA)
       !res = C3NS3C(zero, nf_int)
       res = C3NM3C(zero, nf_int) + c3q3dfPC (zero, nf_int) 
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F3NS_plus
  

  !----------------------------------------------------------------------
  ! This is the "plus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the electromagnetic F2, and, in
  ! addition to the singlet, that involves (u+ubar) - (d+dbar) type
  ! contributions, which is precisely what a NS+ coefficient function
  ! multiplies.
  function cfN3LO_F2NS_plus(y) result(res)
    use xc2ns3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    select case(cc_piece)
    case(cc_REAL)
       res = C2NP3A(x, -y, nf_int, 1) + C2NS3B(x, -y, nf_int)
    case(cc_REALVIRT)
       res = C2NP3A(x, -y, nf_int, 1)
    case(cc_VIRT)
       res = - C2NS3B(x, -y, nf_int)
    case(cc_DELTA)
       res = C2NP3C(zero, nf_int, 1) ! FD: Why is this C2NP3C instead of C2NS3C as in the NNLO case ???       
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2NS_plus
  
  !----------------------------------------------------------------------
  ! This is the "plus" non-singlet coefficient function, but only the fl11 part
  ! 
  function cfN3LO_F2NS_plus_fl11(y) result(res)
    use xc2ns3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    select case(cc_piece)
    case(cc_REAL)
       res = C2NP3A(x, -y, nf_int, 0)
    case(cc_REALVIRT)
       res = C2NP3A(x, -y, nf_int, 0)
    case(cc_VIRT)
       res = res
    case(cc_DELTA)
       res = C2NP3C(zero, nf_int, 0) ! FD: Why is this C2NP3C instead of C2NS3C as in the NNLO case ???       
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2NS_plus_fl11

  !----------------------------------------------------------------------
  ! This is the "minus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the charged-current F2, and, in
  ! addition to the singlet, that involves (u-ubar) - (d-dbar) type
  ! contributions, which is precisely what a NS- coefficient function
  ! multiplies.
  function cfN3LO_F2NS_minus(y) result(res)
    use xc2ns3p
    use xcdiff3pnew
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    
    select case(cc_piece)
    case(cc_REAL)
       res = C2NP3A(x, -y, nf_int, 1) - C2Q3DFP(x, -y, nf_int) + C2NS3B(x, -y, nf_int)
    case(cc_REALVIRT)
       res = C2NP3A(x, -y, nf_int, 1) - C2Q3DFP(x, -y, nf_int)
    case(cc_VIRT)
       res = - C2NS3B(x, -y, nf_int)
    case(cc_DELTA)
       res = C2NP3C(zero, nf_int, 1) - c2q3dfPC (zero, nf_int) ! FD: Why is this C2NP3C instead of C2NS3C as in the NNLO case ???
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2NS_minus


  !----------------------------------------------------------------------
  ! This is the "minus" non-singlet coefficient function but only the fl11 part
  ! 
  function cfN3LO_F2NS_minus_fl11(y) result(res)
    use xc2ns3p
    ! use xcdiff3pnew
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    
    select case(cc_piece)
    case(cc_REAL)
       res = C2NP3A(x, -y, nf_int, 0)
    case(cc_REALVIRT)
       res = C2NP3A(x, -y, nf_int, 0)
    case(cc_VIRT)
       res = res
    case(cc_DELTA)
       res = C2NP3C(zero, nf_int, 0) ! FD: Why is this C2NP3C instead of C2NS3C as in the NNLO case ???
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2NS_minus_fl11

  !----------------------------------------------------------------------
  ! Pure singlet piece of F2 N3LO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfN3LO_F2PS(y) result(res)
    use xc2sg3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = C2S3A(x, -y, nf_int, 1)
    case(cc_REALVIRT)
       res = C2S3A(x, -y, nf_int, 1)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2PS
  
  !----------------------------------------------------------------------
  ! Pure singlet piece of F2 N3LO coefficient function.
  ! This is just the fl11 piece
  function cfN3LO_F2PS_fl11(y) result(res)
    use xc2sg3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = C2S3A(x, -y, nf_int, 0)
    case(cc_REALVIRT)
       res = C2S3A(x, -y, nf_int, 0)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = C2S3C(zero, nf_int) !  This piece has an fl11 contribution only
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2PS_fl11

  !----------------------------------------------------------------------
  ! Pure singlet piece of F2 N3LO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfN3LO_F2G(y) result(res)
    use xc2sg3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = C2G3A(x, -y, nf_int, 1)
    case(cc_REALVIRT)
       res = C2G3A(x, -y, nf_int, 1)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = C2G3C(zero, nf_int)
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2G

  !----------------------------------------------------------------------
  ! Pure singlet piece of F2 N3LO coefficient function.
  ! this is just the fl11 piece
  function cfN3LO_F2G_fl11(y) result(res)
    use xc2sg3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = C2G3A(x, -y, nf_int, 0)
    case(cc_REALVIRT)
       res = C2G3A(x, -y, nf_int, 0)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2G_fl11

  !----------------------------------------------------------------------
  ! This is the "plus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the electromagnetic FL, and, in
  ! addition to the singlet, that involves (u+ubar) - (d+dbar) type
  ! contributions, which is precisely what a NS+ coefficient function
  ! multiplies.
  function cfN3LO_FLNS_plus(y) result(res)
    use xclns3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = CLNP3A(x, -y, nf_int, 1)
    case(cc_REALVIRT)
       res = CLNP3A(x, -y, nf_int, 1)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = CLNP3C(zero, nf_int)
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLNS_plus
  
  !----------------------------------------------------------------------
  ! This is the "plus" non-singlet coefficient function.
  ! this is the fl11 part.
  function cfN3LO_FLNS_plus_fl11(y) result(res)
    use xclns3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = CLNP3A(x, -y, nf_int, 0)
    case(cc_REALVIRT)
       res = CLNP3A(x, -y, nf_int, 0)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLNS_plus_fl11

  !----------------------------------------------------------------------
  ! This is the "minus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the charged-current FL, and, in
  ! addition to the singlet, that involves (u-ubar) - (d-dbar) type
  ! contributions, which is precisely what a NS- coefficient function
  ! multiplies.
  function cfN3LO_FLNS_minus(y) result(res)
    use xclns3p
    use xcdiff3pnew
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    select case(cc_piece)
    case(cc_REAL)
       res = CLNP3A(x, -y, nf_int, 1) - CLQ3DFP(x, -y, nf_int)
    case(cc_REALVIRT)
       res = CLNP3A(x, -y, nf_int, 1) - CLQ3DFP(x, -y, nf_int)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = CLNP3C(zero, nf_int)
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLNS_minus
  
  !----------------------------------------------------------------------
  ! This is the "minus" non-singlet coefficient function
  ! this is the fl11 part
  function cfN3LO_FLNS_minus_fl11(y) result(res)
    use xclns3p
    use xcdiff3pnew
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    select case(cc_piece)
    case(cc_REAL)
       res = CLNP3A(x, -y, nf_int, 0) 
    case(cc_REALVIRT)
       res = CLNP3A(x, -y, nf_int, 0) 
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLNS_minus_fl11

  !----------------------------------------------------------------------
  ! Pure singlet piece of FL N3LO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfN3LO_FLPS(y) result(res)
    use xclsg3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = CLS3A(x, -y, nf_int, 1)
    case(cc_REALVIRT)
       res = CLS3A(x, -y, nf_int, 1)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLPS

  !----------------------------------------------------------------------
  ! Pure singlet piece of FL N3LO coefficient function.
  ! this is the fl11 part
  function cfN3LO_FLPS_fl11(y) result(res)
    use xclsg3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = CLS3A(x, -y, nf_int, 0)
    case(cc_REALVIRT)
       res = CLS3A(x, -y, nf_int, 0)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLPS_fl11

  !----------------------------------------------------------------------
  ! Pure singlet piece of FL N3LO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfN3LO_FLG(y) result(res)
    use xclsg3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = CLG3A(x, -y, nf_int, 1)
    case(cc_REALVIRT)
       res = CLG3A(x, -y, nf_int, 1)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLG

  !----------------------------------------------------------------------
  ! Pure singlet piece of FL N3LO coefficient function.
  ! this is the fl11 part
  function cfN3LO_FLG_fl11(y) result(res)
    use xclsg3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = CLG3A(x, -y, nf_int, 0)
    case(cc_REALVIRT)
       res = CLG3A(x, -y, nf_int, 0)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLG_fl11


end module coefficient_functions_holder
