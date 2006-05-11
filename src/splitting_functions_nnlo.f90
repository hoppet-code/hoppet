

!======================================================================
! interface to exact 3-loop splitting functions, as programmed by Vogt
! (see xpij2e.f and xpns2e.f for details).
!
! There are several things to remember: 
!    . Vogt y == my x
!    . Vogt is normalised to powers of (as/4pi)^3 
!    . colour-factor dependence IS available (taken from qcd module)
!      but only for CF, CA -- TR is fixed to be 1/2
!    . systematic separation in A,B,C functions needs to be
!      taken into account; C function needs to be called with x=0?x
!    . check inclusion of 2nf in qg piece
! 
!----------------------------------------------------------------------
module splitting_functions_nnlo_e
  use types; use consts_dp; use convolution_communicator
  use qcd; use warnings_and_errors
  use xpij2e; use xpns2e
  implicit none
  private

  public :: sf_P2gg, sf_P2qg2nf, sf_P2PS, sf_P2gq
  public :: sf_P2NSPlus, sf_P2NSMinus, sf_P2NSS
  ! these are of help elsewhere in identifying (run-time) the sets of 
  ! splitting functions that have been used?
  public :: name_xpij2, name_xpns2

contains
  function sf_P2gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = X2GGA(x, nf_int) + X2GGB(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - X2GGB(x, nf_int)
    case(cc_DELTA)
       res = X2GGC(zero, nf_int)
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2gg

  ! ..This is the (regular) pure-singlet splitting functions P_ps^(2).    
  !    P_qq^(2) is obtained by adding the non-singlet quantity P_NS^(2)+. 
  function sf_P2PS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = X2PSA(x, nf_int) + zero
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2PS

  !--------------------------------------------------------------------
  ! this already includes the factor of 2nf
  function sf_P2qg2nf(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = X2QGA(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2qg2nf


  function sf_P2gq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = X2GQA(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2gq

  !------------------------- non-singlet pieces ------------------------
  function sf_P2NSPlus(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = X2NSPA(x, nf_int) + X2NSB(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - X2NSB(x, nf_int)
    case(cc_DELTA)
       res = X2NSC(zero, nf_int)
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2NSPlus


  function sf_P2NSMinus(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = X2NSMA(x, nf_int) + X2NSB(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - X2NSB(x, nf_int)
    case(cc_DELTA)
       res = X2NSC(zero, nf_int)
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2NSMinus

  !-- according to comments in Vogt code P_S = P_V - P_-
  function sf_P2NSS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = X2NSSA(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2NSS

  !----------------------------------------------------------------
  ! The Vogt expressions are valid only for standard colour factors
  ! This routine makes sure that the appropriate conditions hold
  subroutine sf_VogtValidate
    if (tr == half) return
    call wae_Error('sf_VogtValidate: &
         &colour factors must be set to default values', &
         &'in order to use the Vogt splitting function parameterisations')
  end subroutine sf_VogtValidate

end module splitting_functions_nnlo_e


!======================================================================
! interface to approx 3-loop splitting functions, as parameterised by Vogt
! (see xpij2p.f and xpns2p.f for details).
!
! There are several things to remember: 
!    . Vogt y == my x
!    . Vogt is normalised to powers of (as/4pi)^3 
!    . colour-factor dependence is not available
!    . systematic separation in A,B,C functions needs to be
!      taken into account; C function needs to be called with x=0?x
!    . check inclusion of 2nf in qg piece
! 
!----------------------------------------------------------------------
module splitting_functions_nnlo_p
  use types; use consts_dp; use convolution_communicator
  use qcd; use warnings_and_errors
  use xpij2p; use xpns2p
  implicit none
  private

  public :: sf_P2gg, sf_P2qg2nf, sf_P2PS, sf_P2gq
  public :: sf_P2NSPlus, sf_P2NSMinus, sf_P2NSS
  ! these are of help elsewhere in identifying (run-time) the sets of 
  ! splitting functions that have been used?
  public :: name_xpij2, name_xpns2

contains
  function sf_P2gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P2GGA(x, nf_int) + P2GGB(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - P2GGB(x, nf_int)
    case(cc_DELTA)
       res = P2GGC(zero, nf_int)
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2gg

  ! ..This is the (regular) pure-singlet splitting functions P_ps^(2).    
  !    P_qq^(2) is obtained by adding the non-singlet quantity P_NS^(2)+. 
  function sf_P2PS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P2PSA(x, nf_int) + zero
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2PS

  !--------------------------------------------------------------------
  ! this already includes the factor of 2nf
  function sf_P2qg2nf(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P2QGA(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2qg2nf


  function sf_P2gq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P2GQA(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2gq

  !------------------------- non-singlet pieces ------------------------
  function sf_P2NSPlus(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P2NSPA(x, nf_int) + P2NSB(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - P2NSB(x, nf_int)
    case(cc_DELTA)
       res = P2NSPC(zero, nf_int)
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2NSPlus


  function sf_P2NSMinus(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P2NSMA(x, nf_int) + P2NSB(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - P2NSB(x, nf_int)
    case(cc_DELTA)
       res = P2NSMC(zero, nf_int)
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2NSMinus

  !-- according to comments in Vogt code P_S = P_V - P_-
  function sf_P2NSS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P2NSSA(x, nf_int)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2NSS

  !----------------------------------------------------------------
  ! The Vogt expressions are valid only for standard colour factors
  ! This routine makes sure that the appropriate conditions hold
  subroutine sf_VogtValidate
    if (ca == three .and. tr == half .and. cf == four/three) return
    call wae_Error('sf_VogtValidate: &
         &colour factors must be set to default values', &
         &'in order to use the Vogt splitting function parameterisations')
  end subroutine sf_VogtValidate

end module splitting_functions_nnlo_p



!======================================================================
! interface to guessed 3-loop splitting functions, as parameterised by Vogt
! (see xpij2n.f and xpns2n.f for details).
!
! There are several things to remember: 
!    . Vogt y == my x
!    . Vogt is normalised to powers of (as/4pi)^3 
!    . colour-factor dependence is not available
!    . systematic separation in A,B,C functions needs to be
!      taken into account; C function needs to be called with x=0?x
!    . check inclusion of 2nf in qg piece
! 
!----------------------------------------------------------------------
module splitting_functions_nnlo_n
  use types; use consts_dp; use convolution_communicator
  use qcd; use dglap_choices; use warnings_and_errors
  use xpij2n; use xpns2n
  implicit none
  private

  public :: sf_P2gg, sf_P2qg2nf, sf_P2PS, sf_P2gq
  public :: sf_P2NSPlus, sf_P2NSMinus, sf_P2NSS
  ! these are of help elsewhere in identifying (run-time) the sets of 
  ! splitting functions that have been used?
  public :: name_xpij2, name_xpns2

contains
  function sf_P2gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P2GGA(x, nf_int, nnlo_splitting_variant) + P2GGB(x, nf_int, nnlo_splitting_variant)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - P2GGB(x, nf_int, nnlo_splitting_variant)
    case(cc_DELTA)
       res = P2GGC(zero, nf_int, nnlo_splitting_variant)
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2gg

  ! ..This is the (regular) pure-singlet splitting functions P_ps^(2).    
  !    P_qq^(2) is obtained by adding the non-singlet quantity P_NS^(2)+. 
  function sf_P2PS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P2PSA(x, nf_int, nnlo_splitting_variant) + zero
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2PS

  !--------------------------------------------------------------------
  ! this already includes the factor of 2nf
  function sf_P2qg2nf(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P2QGA(x, nf_int, nnlo_splitting_variant)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2qg2nf


  function sf_P2gq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P2GQA(x, nf_int, nnlo_splitting_variant)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2gq

  !------------------------- non-singlet pieces ------------------------
  function sf_P2NSPlus(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P2NSPA(x, nf_int, nnlo_splitting_variant) + P2NSPB(x, nf_int, nnlo_splitting_variant)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - P2NSPB(x, nf_int, nnlo_splitting_variant)
    case(cc_DELTA)
       res = P2NSPC(zero, nf_int, nnlo_splitting_variant)
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2NSPlus


  function sf_P2NSMinus(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P2NSMA(x, nf_int, nnlo_splitting_variant) + P2NSMB(x, nf_int, nnlo_splitting_variant)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - P2NSMB(x, nf_int, nnlo_splitting_variant)
    case(cc_DELTA)
       res = P2NSMC(zero, nf_int, nnlo_splitting_variant)
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2NSMinus

  !-- according to comments in Vogt code P_S = P_V - P_-
  function sf_P2NSS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P2NSSA(x, nf_int, nnlo_splitting_variant)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 8.0_dp   ! compensate (as/4pi)^3 -> (as/2pi)^3
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P2NSS

  !----------------------------------------------------------------
  ! The Vogt expressions are valid only for standard colour factors
  ! This routine makes sure that the appropriate conditions hold
  subroutine sf_VogtValidate
    if (ca == three .and. tr == half .and. cf == four/three) return
    call wae_Error('sf_VogtValidate: &
         &colour factors must be set to default values', &
         &'in order to use the Vogt splitting function parameterisations')
  end subroutine sf_VogtValidate

end module splitting_functions_nnlo_n


!======================================================================
! Modules for extracting selecting the relevant splitting functions, 
! either "n" version, "p" versions. May be extended to
! exact versions at some later stage...
!======================================================================
module sf_nnlo_n_rename
  use splitting_functions_nnlo_n,&
       &sf_n_P2gg      => sf_P2gg     ,&
       &sf_n_P2PS      => sf_P2PS     ,&
       &sf_n_P2qg2nf   => sf_P2qg2nf  ,&
       &sf_n_P2gq      => sf_P2gq     ,&
       &sf_n_P2NSPlus  => sf_P2NSPlus ,&
       &sf_n_P2NSMinus => sf_P2NSMinus,&
       &sf_n_P2NSS     => sf_P2NSS
  implicit none
  private
  public :: sf_n_P2gg, sf_n_P2qg2nf, sf_n_P2PS, sf_n_P2gq
  public :: sf_n_P2NSPlus, sf_n_P2NSMinus, sf_n_P2NSS
end module sf_nnlo_n_rename


module sf_nnlo_p_rename
  use splitting_functions_nnlo_p,&
       &sf_p_P2gg      => sf_P2gg     ,&
       &sf_p_P2PS      => sf_P2PS     ,&
       &sf_p_P2qg2nf   => sf_P2qg2nf  ,&
       &sf_p_P2gq      => sf_P2gq     ,&
       &sf_p_P2NSPlus  => sf_P2NSPlus ,&
       &sf_p_P2NSMinus => sf_P2NSMinus,&
       &sf_p_P2NSS     => sf_P2NSS
  implicit none
  private
  public :: sf_p_P2gg, sf_p_P2qg2nf, sf_p_P2PS, sf_p_P2gq
  public :: sf_p_P2NSPlus, sf_p_P2NSMinus, sf_p_P2NSS
end module sf_nnlo_p_rename


module sf_nnlo_e_rename
  use splitting_functions_nnlo_e,&
       &sf_e_P2gg      => sf_P2gg     ,&
       &sf_e_P2PS      => sf_P2PS     ,&
       &sf_e_P2qg2nf   => sf_P2qg2nf  ,&
       &sf_e_P2gq      => sf_P2gq     ,&
       &sf_e_P2NSPlus  => sf_P2NSPlus ,&
       &sf_e_P2NSMinus => sf_P2NSMinus,&
       &sf_e_P2NSS     => sf_P2NSS
  implicit none
  private
  public :: sf_e_P2gg, sf_e_P2qg2nf, sf_e_P2PS, sf_e_P2gq
  public :: sf_e_P2NSPlus, sf_e_P2NSMinus, sf_e_P2NSS
end module sf_nnlo_e_rename


module splitting_functions_nnlo
  use types; use qcd; use dglap_choices
  use warnings_and_errors
  use sf_nnlo_e_rename
  use sf_nnlo_p_rename
  use sf_nnlo_n_rename
  public :: sf_P2gg, sf_P2qg2nf, sf_P2PS, sf_P2gq
  public :: sf_P2NSPlus, sf_P2NSMinus, sf_P2NSS
contains
  
  
  real(dp) function sf_P2gg     (y) result(res)
    real(dp), intent(in) :: y
    select case(nnlo_splitting_variant)
    case(nnlo_splitting_exact)
       res = sf_e_P2gg(y)
    case(nnlo_splitting_param)
       res = sf_p_P2gg(y)
    case(nnlo_splitting_nfitav, nnlo_splitting_nfiterr1, nnlo_splitting_nfiterr2)
       res = sf_n_P2gg(y)
    case default
       call wae_error('splitting_functions_nnlo','unrecognized imod',&
            &intval=nnlo_splitting_variant)
    end select
  end function sf_P2gg
  
  real(dp) function sf_P2PS     (y) result(res)
    real(dp), intent(in) :: y
    select case(nnlo_splitting_variant)
    case(nnlo_splitting_exact)
       res = sf_e_P2PS(y)
    case(nnlo_splitting_param)
       res = sf_p_P2PS(y)
    case(nnlo_splitting_nfitav, nnlo_splitting_nfiterr1, nnlo_splitting_nfiterr2)
       res = sf_n_P2PS(y)
    case default
       call wae_error('splitting_functions_nnlo','unrecognized imod',&
            &intval=nnlo_splitting_variant)
    end select
  end function sf_P2PS
  
  real(dp) function sf_P2qg2nf  (y) result(res)
    real(dp), intent(in) :: y
    select case(nnlo_splitting_variant)
    case(nnlo_splitting_exact)
       res = sf_e_P2qg2nf(y)
    case(nnlo_splitting_param)
       res = sf_p_P2qg2nf(y)
    case(nnlo_splitting_nfitav, nnlo_splitting_nfiterr1, nnlo_splitting_nfiterr2)
       res = sf_n_P2qg2nf(y)
    case default
       call wae_error('splitting_functions_nnlo','unrecognized imod',&
            &intval=nnlo_splitting_variant)
    end select
  end function sf_P2qg2nf
  
  real(dp) function sf_P2gq     (y) result(res)
    real(dp), intent(in) :: y
    select case(nnlo_splitting_variant)
    case(nnlo_splitting_exact)
       res = sf_e_P2gq(y)
    case(nnlo_splitting_param)
       res = sf_p_P2gq(y)
    case(nnlo_splitting_nfitav, nnlo_splitting_nfiterr1, nnlo_splitting_nfiterr2)
       res = sf_n_P2gq(y)
    case default
       call wae_error('splitting_functions_nnlo','unrecognized imod',&
            &intval=nnlo_splitting_variant)
    end select
  end function sf_P2gq
  
  real(dp) function sf_P2NSPlus (y) result(res)
    real(dp), intent(in) :: y
    select case(nnlo_splitting_variant)
    case(nnlo_splitting_exact)
       res = sf_e_P2NSPlus(y)
    case(nnlo_splitting_param)
       res = sf_p_P2NSPlus(y)
    case(nnlo_splitting_nfitav, nnlo_splitting_nfiterr1, nnlo_splitting_nfiterr2)
       res = sf_n_P2NSPlus(y)
    case default
       call wae_error('splitting_functions_nnlo','unrecognized imod',&
            &intval=nnlo_splitting_variant)
    end select
  end function sf_P2NSPlus
  
  real(dp) function sf_P2NSMinus(y) result(res)
    real(dp), intent(in) :: y
    select case(nnlo_splitting_variant)
    case(nnlo_splitting_exact)
       res = sf_e_P2NSMinus(y)
    case(nnlo_splitting_param)
       res = sf_p_P2NSMinus(y)
    case(nnlo_splitting_nfitav, nnlo_splitting_nfiterr1, nnlo_splitting_nfiterr2)
       res = sf_n_P2NSMinus(y)
    case default
       call wae_error('splitting_functions_nnlo','unrecognized imod',&
            &intval=nnlo_splitting_variant)
    end select
  end function sf_P2NSMinus

  real(dp) function sf_P2NSS(y) result(res)
    real(dp), intent(in) :: y
    select case(nnlo_splitting_variant)
    case(nnlo_splitting_exact)
       res = sf_e_P2NSS(y)
    case(nnlo_splitting_param)
       res = sf_p_P2NSS(y)
    case(nnlo_splitting_nfitav, nnlo_splitting_nfiterr1, nnlo_splitting_nfiterr2)
       res = sf_n_P2NSS(y)
    case default
       call wae_error('splitting_functions_nnlo','unrecognized imod',&
            &intval=nnlo_splitting_variant)
    end select
  end function sf_P2NSS
end module splitting_functions_nnlo
