
module coefficient_functions_holder_exact
  use types; use warnings_and_errors
  use consts_dp; use convolution_communicator
  use qcd
  use convolution
  use dglap_objects
  use coefficient_functions_holder 
  use coefficient_functions_holder_internal
  implicit none

contains

  ! when we have the exact coefficient functions, this is the actual
  ! implementation of their initialisation. It assumes that the 
  ! the coef_holder, cf, already have nloop and the current nf set.
  !
  ! NB: to aid with conditional compilation, this should only be 
  !     be called from InitCoefHolderExact
  subroutine InitCoefHolderExact_impl(grid, ch)
    type(grid_def), intent(in) :: grid
    type(coef_holder), intent(inout), target :: ch
    !------------------------------
    write(6,*) 'Initialising exact coefficient functions for nf =', nf_int

    if (ch%nloop.ge.3) then
       write(6,*) 'Initialising C2NNLO'
       call InitC2NNLO_e(grid, ch%C2NNLO)
       write(6,*) 'Initialising CLNNLO'
       call InitCLNNLO_e(grid, ch%CLNNLO)
       write(6,*) 'Initialising C3NNLO'
       call InitC3NNLO_e(grid, ch%C3NNLO)
    end if 
        
    if (ch%nloop.ge.4) then
       write(6,*) 'Initialising C2N3LO'
       call InitC2N3LO_e(grid, ch%C2N3LO)
       write(6,*) 'Initialising CLN3LO'
       call InitCLN3LO_e(grid, ch%CLN3LO)
       write(6,*) 'Initialising C3N3LO'
       call InitC3N3LO_e(grid, ch%C3N3LO)
       ! including the fl11 terms for Z/photon exchanges
       write(6,*) 'Initialising C2N3LO_fl11'
       call InitC2N3LO_fl11_e(grid, ch%C2N3LO_fl11)
       write(6,*) 'Initialising CLN3LO_fl11'
       call InitCLN3LO_fl11_e(grid, ch%CLN3LO_fl11)
    end if
    
  end subroutine InitCoefHolderExact_impl

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
  ! N3LO exact coefficient functions
  !======================================================================

  !======================================================================
  ! N3LO C2 coefficient function as a splitting matrix
  subroutine InitC2N3LO_e(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%NS_plus, cfN3LO_F2NS_plus_e)
    call InitGridConv(grid, C%NS_minus, cfN3LO_F2NS_minus_e)
    ! valence shouldn't come in to F2, so it shouldn't matter
    ! what we put here (explicitly verified for F2Z and F2Wp*F2Wm)
    call InitGridConv(C%NS_V, C%NS_minus) 
    !call InitGridConv(grid, C%NS_V) 

    ! singlet qq is NS+ + pure_singlet
    call InitGridConv(C%qq, C%NS_plus)    ! start with NS+
    call AddWithCoeff(C%qq, cfN3LO_F2PS_e)  ! add in pure singlet piece
    call InitGridConv(grid, C%qg, cfN3LO_F2G_e)   ! the part driven by the gluon
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

  end subroutine InitC2N3LO_e

  !======================================================================
  ! N3LO C2 fl11 coefficient function as a splitting matrix
  ! Only the part that is different for Z boson exchange
  subroutine InitC2N3LO_fl11_e(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%NS_plus, cfN3LO_F2NS_plus_fl11_e)
    call InitGridConv(grid, C%NS_minus, cfN3LO_F2NS_minus_fl11_e)
    ! valence shouldn't come in to F2, so it shouldn't matter
    ! what we put here (explicitly verified for F2Z and F2Wp*F2Wm)
    call InitGridConv(C%NS_V, C%NS_minus) 
    !call InitGridConv(grid, C%NS_V) 

    ! singlet qq is NS+ + pure_singlet
    call InitGridConv(C%qq, C%NS_plus)    ! start with NS+
    call AddWithCoeff(C%qq, cfN3LO_F2PS_fl11_e)  ! add in pure singlet piece
    call InitGridConv(grid, C%qg, cfN3LO_F2G_fl11_e)   ! the part driven by the gluon
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

  end subroutine InitC2N3LO_fl11_e

  !======================================================================
  ! N3LO CL coefficient function as a splitting matrix
  subroutine InitCLN3LO_e(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%NS_plus, cfN3LO_FLNS_plus_e)
    call InitGridConv(grid, C%NS_minus, cfN3LO_FLNS_minus_e)
    ! valence shouldn't come in to F2, so it shouldn't matter
    ! what we put here (explicitly verified for F2Z and F2Wp*F2Wm)
    call InitGridConv(C%NS_V, C%NS_minus) 
    !call InitGridConv(grid, C%NS_V) 

    ! singlet qq is NS+ + pure_singlet
    call InitGridConv(C%qq, C%NS_plus)    ! start with NS+
    call AddWithCoeff(C%qq, cfN3LO_FLPS_e)  ! add in pure singlet piece
    call InitGridConv(grid, C%qg, cfN3LO_FLG_e)   ! the part driven by the gluon
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

  end subroutine InitCLN3LO_e


  !======================================================================
  ! N3LO CL coefficient function as a splitting matrix
  subroutine InitCLN3LO_fl11_e(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    C%nf_int = nf_int
    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%NS_plus, cfN3LO_FLNS_plus_fl11_e)
    call InitGridConv(grid, C%NS_minus, cfN3LO_FLNS_minus_fl11_e)
    ! valence shouldn't come in to F2, so it shouldn't matter
    ! what we put here (explicitly verified for F2Z and F2Wp*F2Wm)
    call InitGridConv(C%NS_V, C%NS_minus) 
    !call InitGridConv(grid, C%NS_V) 

    ! singlet qq is NS+ + pure_singlet
    call InitGridConv(C%qq, C%NS_plus)    ! start with NS+
    call AddWithCoeff(C%qq, cfN3LO_FLPS_fl11_e)  ! add in pure singlet piece
    call InitGridConv(grid, C%qg, cfN3LO_FLG_fl11_e)   ! the part driven by the gluon
    call InitGridConv(grid, C%gq)  ! zero
    call InitGridConv(grid, C%gg)  ! zero

  end subroutine InitCLN3LO_fl11_e

  !======================================================================
  ! N3LO C3 coefficient function as a splitting matrix
  subroutine InitC3N3LO_e(grid, C)
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
    call InitGridConv(grid, C%NS_plus, cfN3LO_F3NS_plus_e)
    !call SetDefaultConvolutionEps(1d-4)
    call InitGridConv(grid, C%NS_minus, cfN3LO_F3NS_minus_e)
    !call SetDefaultConvolutionEps(1d-7)
    ! NS_minus and NS_V are not identical at N3LO
    call InitGridConv(grid, C%NS_V, cfN3LO_F3NS_val_e) 

  end subroutine InitC3N3LO_e



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
  ! N3LO PIECES (exact)
  !======================================================================
  !======================================================================

  !----------------------------------------------------------------------
  ! This is the coefficient function for the sum F3(nu)+F3(nubar),
  ! corresponding to C3NSP-C3NSN in W. van Neerven's program.
  function cfN3LO_F3NS_val_e(y) result(res)
    use xc3ns3p
    use xc3ns3e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    CALL SET_C3SOFT_N3LO(nf_int)

    select case(cc_piece)
    case(cc_REAL)
       res = X3NM3A(x, -y, nf_int, 1) + X3NS3B(x, -y, nf_int)
    case(cc_REALVIRT)
       res = X3NM3A(x, -y, nf_int, 1)
    case(cc_VIRT)
       res = - X3NS3B(x, -y, nf_int)
    case(cc_DELTA)
       res = X3NS3C(zero, -y, nf_int) 
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F3NS_val_e

  !----------------------------------------------------------------------
  ! This is the coefficient function for the sum F3(nu)+F3(nubar),
  ! corresponding to C3NSP-C3NSN in W. van Neerven's program.
  function cfN3LO_F3NS_minus_e(y) result(res)
    use xc3ns3e
    use xc3ns3p
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    CALL SET_C3SOFT_N3LO(nf_int)

    select case(cc_piece)
    case(cc_REAL)
       res = X3NM3A(x, -y, nf_int, 0) + X3NS3B(x, -y, nf_int)
    case(cc_REALVIRT)
       res = X3NM3A(x, -y, nf_int, 0)
    case(cc_VIRT)
       res = - X3NS3B(x, -y, nf_int)
    case(cc_DELTA)
       res = X3NS3C(zero, -y, nf_int) 
    end select
    
    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F3NS_minus_e

  !----------------------------------------------------------------------
  ! This is the coefficient function for the sum F3(nu)-F3(nubar),
  ! corresponding to C3NSP-C3NSN in W. van Neerven's program.
  function cfN3LO_F3NS_plus_e(y) result(res)
    use xc3ns3p
    use xc3ns3e
    use xcdiff3pnew
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    CALL SET_C3SOFT_N3LO(nf_int)

    select case(cc_piece)
    case(cc_REAL)
       res = X3NM3A(x, -y, nf_int, 0) + C3Q3DFP(x, -y, nf_int) + X3NS3B(x, -y, nf_int)
    case(cc_REALVIRT)
       res = X3NM3A(x, -y, nf_int, 0) + C3Q3DFP(x, -y, nf_int)
    case(cc_VIRT)
       res = - X3NS3B(x, -y, nf_int)
    case(cc_DELTA)
       res = X3NS3C(zero, -y, nf_int) + c3q3dfPC (zero, nf_int)
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F3NS_plus_e
  

  !----------------------------------------------------------------------
  ! This is the "plus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the electromagnetic F2, and, in
  ! addition to the singlet, that involves (u+ubar) - (d+dbar) type
  ! contributions, which is precisely what a NS+ coefficient function
  ! multiplies.
  function cfN3LO_F2NS_plus_e(y) result(res)
    use xc2ns3p
    use xc2ns3e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    CALL SET_C2SOFT_N3LO(nf_int)
    select case(cc_piece)
    case(cc_REAL)
       res = X2NP3A(x, -y, nf_int, 1) + X2NS3B(x, -y, nf_int)
    case(cc_REALVIRT)
       res = X2NP3A(x, -y, nf_int, 1)
    case(cc_VIRT)
       res = - X2NS3B(x, -y, nf_int)
    case(cc_DELTA)
       res = X2NP3C(zero, nf_int, 1)        
    end select
    
    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2NS_plus_e
  
  !----------------------------------------------------------------------
  ! This is the "plus" non-singlet coefficient function, but only the fl11 part
  ! 
  function cfN3LO_F2NS_plus_fl11_e(y) result(res)
    use xc2ns3e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    call SET_C2SOFT_N3LO(nf_int)
    select case(cc_piece)
    case(cc_REAL)
       res = X2NP3A(x, -y, nf_int, 0)
    case(cc_REALVIRT)
       res = X2NP3A(x, -y, nf_int, 0)
    case(cc_VIRT)
       res = res
    case(cc_DELTA)
       res = X2NP3C(zero, nf_int, 0) 
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2NS_plus_fl11_e

  !----------------------------------------------------------------------
  ! This is the "minus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the charged-current F2, and, in
  ! addition to the singlet, that involves (u-ubar) - (d-dbar) type
  ! contributions, which is precisely what a NS- coefficient function
  ! multiplies.
  function cfN3LO_F2NS_minus_e(y) result(res)
    use xc2ns3p
    use xc2ns3e
    use xcdiff3pnew
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    CALL SET_C2SOFT_N3LO(nf_int)
    select case(cc_piece)
    case(cc_REAL)
       res = X2NP3A(x, -y, nf_int, 1) - C2Q3DFP(x, -y, nf_int) + X2NS3B(x, -y, nf_int)
    case(cc_REALVIRT)
       res = X2NP3A(x, -y, nf_int, 1) - C2Q3DFP(x, -y, nf_int)
    case(cc_VIRT)
       res = - X2NS3B(x, -y, nf_int)
    case(cc_DELTA)
       res = X2NP3C(zero, nf_int, 1) - c2q3dfPC (zero, nf_int) 
    end select
 
    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2NS_minus_e


  !----------------------------------------------------------------------
  ! This is the "minus" non-singlet coefficient function but only the fl11 part
  ! 
  function cfN3LO_F2NS_minus_fl11_e(y) result(res)
    use xc2ns3e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    call SET_C2SOFT_N3LO(nf_int)
    select case(cc_piece)
    case(cc_REAL)
       res = X2NP3A(x, -y, nf_int, 0)
    case(cc_REALVIRT)
       res = X2NP3A(x, -y, nf_int, 0)
    case(cc_VIRT)
       res = res
    case(cc_DELTA)
       res = X2NP3C(zero, nf_int, 0) 
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2NS_minus_fl11_e

  !----------------------------------------------------------------------
  ! Pure singlet piece of F2 N3LO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfN3LO_F2PS_e(y) result(res)
    use xc2ps3e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = X2S3A(x, nf_int, 1)
    case(cc_REALVIRT)
       res = X2S3A(x, nf_int, 1)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2PS_e
  
  !----------------------------------------------------------------------
  ! Pure singlet piece of F2 N3LO coefficient function.
  ! This is just the fl11 piece
  function cfN3LO_F2PS_fl11_e(y) result(res)
    use xc2ps3e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = X2S3A(x, nf_int, 0)
    case(cc_REALVIRT)
       res = X2S3A(x, nf_int, 0)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = X2S3C (x, nf_int)
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2PS_fl11_e

  !----------------------------------------------------------------------
  ! Pure singlet piece of F2 N3LO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfN3LO_F2G_e(y) result(res)
    use xc2gl3e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = X2G3A(x, nf_int, 1)
    case(cc_REALVIRT)
       res = X2G3A(x, nf_int, 1)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2G_e

  !----------------------------------------------------------------------
  ! Pure singlet piece of F2 N3LO coefficient function.
  ! this is just the fl11 piece
  function cfN3LO_F2G_fl11_e(y) result(res)
    use xc2gl3e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = X2G3A(x, nf_int, 0)
    case(cc_REALVIRT)
       res = X2G3A(x, nf_int, 0)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_F2G_fl11_e

  !----------------------------------------------------------------------
  ! This is the "plus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the electromagnetic FL, and, in
  ! addition to the singlet, that involves (u+ubar) - (d+dbar) type
  ! contributions, which is precisely what a NS+ coefficient function
  ! multiplies.
  function cfN3LO_FLNS_plus_e(y) result(res)
    use xclns3e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    select case(cc_piece)
    case(cc_REAL)
       res = XLNP3A(x, nf_int, 1)
    case(cc_REALVIRT)
       res = XLNP3A(x, nf_int, 1)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLNS_plus_e
  
  !----------------------------------------------------------------------
  ! This is the "plus" non-singlet coefficient function.
  ! this is the fl11 part.
  function cfN3LO_FLNS_plus_fl11_e(y) result(res)
    use xclns3e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = XLNP3A(x, nf_int, 0)
    case(cc_REALVIRT)
       res = XLNP3A(x, nf_int, 0)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLNS_plus_fl11_e

  !----------------------------------------------------------------------
  ! This is the "minus" non-singlet coefficient function (Vogt's
  ! comments say it applies to the charged-current FL, and, in
  ! addition to the singlet, that involves (u-ubar) - (d-dbar) type
  ! contributions, which is precisely what a NS- coefficient function
  ! multiplies.
  function cfN3LO_FLNS_minus_e(y) result(res)
    use xclns3e
    use xcdiff3pnew
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    select case(cc_piece)
    case(cc_REAL)
       res = XLNP3A(x, nf_int, 1) - CLQ3DFP(x, -y, nf_int)
    case(cc_REALVIRT)
       res = XLNP3A(x, nf_int, 1) - CLQ3DFP(x, -y, nf_int)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero!
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLNS_minus_e
  
  !----------------------------------------------------------------------
  ! This is the "minus" non-singlet coefficient function
  ! this is the fl11 part
  function cfN3LO_FLNS_minus_fl11_e(y) result(res)
    use xclns3p
    use xclns3e
    use xcdiff3pnew
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero
    select case(cc_piece)
    case(cc_REAL)
       res = XLNP3A(x, nf_int, 0) 
    case(cc_REALVIRT)
       res = XLNP3A(x, nf_int, 0) 
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLNS_minus_fl11_e

  !----------------------------------------------------------------------
  ! Pure singlet piece of FL N3LO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfN3LO_FLPS_e(y) result(res)
    use xclps3e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = XLS3A(x, nf_int, 1)
    case(cc_REALVIRT)
       res = XLS3A(x, nf_int, 1)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLPS_e

  !----------------------------------------------------------------------
  ! Pure singlet piece of FL N3LO coefficient function.
  ! this is the fl11 part
  function cfN3LO_FLPS_fl11_e(y) result(res)
    use xclps3e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = XLS3A(x, nf_int, 0)
    case(cc_REALVIRT)
       res = XLS3A(x, nf_int, 0)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLPS_fl11_e

  !----------------------------------------------------------------------
  ! Pure singlet piece of FL N3LO coefficient function. To get the full
  ! singlet result, this is to be added to NS_plus (cf. relation just after
  ! Eq.(3.18) of 1109.3717v2
  function cfN3LO_FLG_e(y) result(res)
    use xclgl3e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = XLG3A(x, nf_int, 1)
    case(cc_REALVIRT)
       res = XLG3A(x, nf_int, 1)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLG_e

  !----------------------------------------------------------------------
  ! Pure singlet piece of FL N3LO coefficient function.
  ! this is the fl11 part
  function cfN3LO_FLG_fl11_e(y) result(res)
    use xclgl3e
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL)
       res = XLG3A(x, nf_int, 0)
    case(cc_REALVIRT)
       res = XLG3A(x, nf_int, 0)
    case(cc_VIRT)
       res = zero
    case(cc_DELTA)
       res = zero
    end select

    res = res * 0.125_dp ! since our convention is to multiply (as/2pi)^3, theirs is to multiply (as/4pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function cfN3LO_FLG_fl11_e

    
end module coefficient_functions_holder_exact


!=============================================================================
! This just calls the function in coefficient_functions_holder_exact
subroutine InitCoefHolderExact(grid, ch)
  use types; use coefficient_functions_holder_exact, dummy => InitCoefHolderExact
  implicit none
  type(grid_def), intent(in) :: grid
  type(coef_holder), intent(inout), target :: ch
  call InitCoefHolderExact_impl(grid, ch)
end subroutine InitCoefHolderExact
