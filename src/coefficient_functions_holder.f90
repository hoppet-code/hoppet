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
  use types; use consts_dp; use convolution_communicator
  use qcd
  use convolution_communicator
  use dglap_objects
  use coefficient_functions 
  implicit none
  public :: InitC2NLO, InitCLNLO, InitC3NLO, InitC2NNLO, InitCLNNLO, InitC3NNLO
  public :: InitC2NNLO_e, InitCLNNLO_e, InitC3NNLO_e, InitC2N3LO, InitCLN3LO, InitC3N3LO
  
contains
  
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
  ! NNLO exact coefficient functions
  !======================================================================

  !======================================================================
  ! attempt to write the NNLO C2 coefficient function as a splitting
  ! matrix
  subroutine InitC2NNLO_e(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%NS_plus, cfNNLO_F2NS_plus_e)
    call InitGridConv(grid, C%NS_minus, cfNNLO_F2NS_minus_e)
    ! valence shouldn't come in to F2, so it shouldn't matter
    ! what we put here (explicitly verified for F2Z and F2Wp*F2Wm)
    call InitGridConv(C%NS_V, C%NS_minus) 
    !call InitGridConv(grid, C%NS_V) 

    ! singlet qq is NS+ + pure_singlet
    call InitGridConv(C%qq, C%NS_plus)    ! start with NS+
    call AddWithCoeff(C%qq, cfNNLO_F2PS_e)  ! add in pure singlet piece
    call InitGridConv(grid, C%qg, cfNNLO_F2G_e)   ! the part driven by the gluon
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

  end subroutine InitC2NNLO_e

  !======================================================================
  ! attempt to write the NNLO CL coefficient function as a splitting
  ! matrix
  subroutine InitCLNNLO_e(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%NS_plus, cfNNLO_FLNS_plus_e)
    call InitGridConv(grid, C%NS_minus, cfNNLO_FLNS_minus_e)
    ! valence shouldn't come in to F2, so it shouldn't matter
    ! what we put here (explicitly verified for F2Z and F2Wp*F2Wm)
    call InitGridConv(C%NS_V, C%NS_minus) 
    !call InitGridConv(grid, C%NS_V) 

    ! singlet qq is NS+ + pure_singlet
    call InitGridConv(C%qq, C%NS_plus)    ! start with NS+
    call AddWithCoeff(C%qq, cfNNLO_FLPS_e)  ! add in pure singlet piece
    call InitGridConv(grid, C%qg, cfNNLO_FLG_e)   ! the part driven by the gluon
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

  end subroutine InitCLNNLO_e

  !======================================================================
  ! attempt to write the NNLO C3 coefficient function as a splitting
  ! matrix
  subroutine InitC3NNLO_e(grid, C)
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
    call InitGridConv(grid, C%NS_plus, cfNNLO_F3NS_plus_e)
    call InitGridConv(grid, C%NS_minus, cfNNLO_F3NS_minus_e)
    ! NS_minus and NS_V are identical according to sentence before
    ! Eq.(3.15) of 1109.3717v2
    call InitGridConv(C%NS_V, C%NS_minus) 

  end subroutine InitC3NNLO_e
  
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
    call AddWithCoeff(C%qq, cfN3LO_F2PS)  ! add in pure singlet piece
    call InitGridConv(grid, C%qg, cfN3LO_F2G)   ! the part driven by the gluon
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
  ! NNLO PIECES (exact version)
  !======================================================================
  !======================================================================


  !----------------------------------------------------------------------
  ! This is the coefficient function for the sum F3(nu)+F3(nubar),
  ! corresponding to C3NSP-C3NSN in W. van Neerven's program.
  function cfNNLO_F3NS_minus_e(y) result(res)
    use xc3ns2e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    call set_c3soft(nf_int)
    
    select case(cc_piece)
    case(cc_REAL)
       res = X3NM2A(x, nf_int) + X3NS2B(x, nf_int)
    case(cc_REALVIRT)
       res = X3NM2A(x, nf_int)
    case(cc_VIRT)
       res = - X3NS2B(x, nf_int)
    case(cc_DELTA)
       res = X3NS2C(zero, nf_int)
    end select
    
    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_F3NS_minus_e

  !----------------------------------------------------------------------
  ! This is the coefficient function for the sum F3(nu)-F3(nubar),
  ! corresponding to C3NSP-C3NSN in W. van Neerven's program.
  function cfNNLO_F3NS_plus_e(y) result(res)
    use xc3ns2e
    use xcdiff2e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    call set_c3soft(nf_int)
    
    select case(cc_piece)
    case(cc_REAL)
       res = X3NM2A(x, nf_int) + XC3DFF2(x) + X3NS2B(x, nf_int)
    case(cc_REALVIRT)
       res = X3NM2A(x, nf_int) + XC3DFF2(x)
    case(cc_VIRT)
       res = - X3NS2B(x, nf_int)
    case(cc_DELTA)
       res = X3NS2C(zero, nf_int)
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_F3NS_plus_e
  

  !----------------------------------------------------------------------
  ! This is the "plus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the electromagnetic F2, and, in
  ! addition to the singlet, that involves (u+ubar) - (d+dbar) type
  ! contributions, which is precisely what a NS+ coefficient function
  ! multiplies.
  function cfNNLO_F2NS_plus_e(y) result(res)
    use xc2ns2e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    call set_c2soft(nf_int)
    select case(cc_piece)
    case(cc_REAL)
       res = X2NP2A(x, nf_int) + X2NS2B(x, nf_int)
    case(cc_REALVIRT)
       res = X2NP2A(x, nf_int)
    case(cc_VIRT)
       res = - X2NS2B(x, nf_int)
    case(cc_DELTA)
       res = X2NS2C(zero, nf_int)
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_F2NS_plus_e

  !----------------------------------------------------------------------
  ! This is the "minus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the charged-current F2, and, in
  ! addition to the singlet, that involves (u-ubar) - (d-dbar) type
  ! contributions, which is precisely what a NS- coefficient function
  ! multiplies.
  function cfNNLO_F2NS_minus_e(y) result(res)
    use xc2ns2e
    use xcdiff2e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    call set_c2soft(nf_int)
    
    select case(cc_piece)
    case(cc_REAL)
       res = X2NP2A(x, nf_int) - XC2DFF2(x) + X2NS2B(x, nf_int)
    case(cc_REALVIRT)
       res = X2NP2A(x, nf_int) - XC2DFF2(x)
    case(cc_VIRT)
       res = - X2NS2B(x, nf_int)
    case(cc_DELTA)
       res = X2NS2C(zero, nf_int)
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_F2NS_minus_e

  !----------------------------------------------------------------------
  ! Pure singlet piece of F2 NNLO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfNNLO_F2PS_e(y) result(res)
    use xc2sg2e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = X2S2A(x, nf_int)
    case(cc_REALVIRT)
       res = X2S2A(x, nf_int)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_F2PS_e

  !----------------------------------------------------------------------
  ! Pure singlet piece of F2 NNLO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfNNLO_F2G_e(y) result(res)
    use xc2sg2e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = X2G2A(x, nf_int)
    case(cc_REALVIRT)
       res = X2G2A(x, nf_int)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_F2G_e

  !----------------------------------------------------------------------
  ! This is the "plus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the electromagnetic FL, and, in
  ! addition to the singlet, that involves (u+ubar) - (d+dbar) type
  ! contributions, which is precisely what a NS+ coefficient function
  ! multiplies.
  function cfNNLO_FLNS_plus_e(y) result(res)
    use xclns2e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = XLNP2A(x, nf_int)
    case(cc_REALVIRT)
       res = XLNP2A(x, nf_int)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_FLNS_plus_e

  !----------------------------------------------------------------------
  ! This is the "minus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the charged-current FL, and, in
  ! addition to the singlet, that involves (u-ubar) - (d-dbar) type
  ! contributions, which is precisely what a NS- coefficient function
  ! multiplies.
  function cfNNLO_FLNS_minus_e(y) result(res)
    use xclns2e
    use xcdiff2e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = XLNP2A(x, nf_int) - XCLDFF2(x)
    case(cc_REALVIRT)
       res = XLNP2A(x, nf_int) - XCLDFF2(x)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_FLNS_minus_e

  !----------------------------------------------------------------------
  ! Pure singlet piece of FL NNLO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfNNLO_FLPS_e(y) result(res)
    use xclsg2e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = XLS2A(x, nf_int)
    case(cc_REALVIRT)
       res = XLS2A(x, nf_int)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_FLPS_e

  !----------------------------------------------------------------------
  ! Pure singlet piece of FL NNLO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfNNLO_FLG_e(y) result(res)
    use xclsg2e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = XLG2A(x, nf_int)
    case(cc_REALVIRT)
       res = XLG2A(x, nf_int)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.25_dp ! since our convention is to multiply (as/2pi)^2, theirs is to multiply (as/4pi)^2
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfNNLO_FLG_e

  
  !======================================================================
  !======================================================================
  ! N3LO PIECES
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
    use xcdiff3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    
    select case(cc_piece)
    case(cc_REAL)
       res = C3NM3A(x, -y, nf_int, 0) + C3Q3DF(x, -y, nf_int, 0) + C3NS3B(x, -y, nf_int)
    case(cc_REALVIRT)
       res = C3NM3A(x, -y, nf_int, 0) + C3Q3DF(x, -y, nf_int, 0)
    case(cc_VIRT)
       res = - C3NS3B(x, -y, nf_int)
    case(cc_DELTA)
       !res = C3NS3C(zero, nf_int)
       res = C3NM3C(zero, nf_int) ! FD: Why is this C3NM3C instead of C3NS3C as in the NNLO case ???
                                  !     and where is the diff piece between C3NM3C and C3NP3C ?
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
    use xcdiff3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    
    select case(cc_piece)
    case(cc_REAL)
       res = C2NP3A(x, -y, nf_int, 1) - C2Q3DF(x, -y, nf_int, 0) + C2NS3B(x, -y, nf_int)
    case(cc_REALVIRT)
       res = C2NP3A(x, -y, nf_int, 1) - C2Q3DF(x, -y, nf_int, 0)
    case(cc_VIRT)
       res = - C2NS3B(x, -y, nf_int)
    case(cc_DELTA)
       res = C2NP3C(zero, nf_int, 1) ! FD: Why is this C2NP3C instead of C2NS3C as in the NNLO case ???
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2NS_minus


  !----------------------------------------------------------------------
  ! This is the "minus" non-singlet coefficient function but only the fl11 part
  ! 
  function cfN3LO_F2NS_minus_fl11(y) result(res)
    use xc2ns3p
    ! use xcdiff3p
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
       res = C2S3C(zero, nf_int)
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
       res = C2G3C(zero, nf_int)
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
       res = CLNP3C(zero, nf_int)
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
    use xcdiff3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    select case(cc_piece)
    case(cc_REAL)
       res = CLNP3A(x, -y, nf_int, 1) - CLQ3DF(x, -y, nf_int, 0)
    case(cc_REALVIRT)
       res = CLNP3A(x, -y, nf_int, 1) - CLQ3DF(x, -y, nf_int, 0)
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
    use xcdiff3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    select case(cc_piece)
    case(cc_REAL)
       res = CLNP3A(x, -y, nf_int, 0) - CLQ3DF(x, -y, nf_int, 0)
    case(cc_REALVIRT)
       res = CLNP3A(x, -y, nf_int, 0) - CLQ3DF(x, -y, nf_int, 0)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = CLNP3C(zero, nf_int)
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
