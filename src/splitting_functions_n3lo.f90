!======================================================================
! This file contains the 4-loop splitting functions as computed by
! Vogt et al. Currently (Jan 2024) only the "guessed" versions are
! available but we keep both the exact and the parametrised versions
! here to enable their implementation later on (although for now they
! return an error). The references for the guessed splitting functions
! are given in 2302.07593, 2307.04158, and 2310.05744.
!
! The non-singlet splitting functions are taken from 1707.08315, with
! some pieces exact, some parametrised and some guessed.

!======================================================================
! interface to exact 3-loop splitting functions, as programmed by Vogt
! (see xpij2e.f and xpns2e.f for details).
!
! There are several things to remember: 
!    . Vogt y == my x
!    . Vogt is normalised to powers of (as/4pi)^4 
!    . colour-factor dependence IS available (taken from qcd module)
!      but only for CF, CA -- TR is fixed to be 1/2
!    . systematic separation in A,B,C functions needs to be
!      taken into account; C function needs to be called with x=0?x
!    . check inclusion of 2nf in qg piece
! 
!----------------------------------------------------------------------
module splitting_functions_n3lo_e
  use types; use consts_dp; use convolution_communicator
  use qcd; use warnings_and_errors
  use xpij2e; use xpns2e
  use dglap_choices
  implicit none
  private

  public :: sf_P3gg, sf_P3qg2nf, sf_P3PS, sf_P3gq
  public :: sf_P3NSPlus, sf_P3NSMinus, sf_P3NSS
  ! these are of help elsewhere in identifying (run-time) the sets of 
  ! splitting functions that have been used?
  public :: name_xpij2, name_xpns2

contains
  function sf_P3gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call wae_Error('Exact N3LO splitting functions not yet available.')

    call sf_VogtValidate
    x = exp(-y)
    res = zero

!    select case(cc_piece)
!    case(cc_REAL,cc_REALVIRT)
!       res = X2GGA(x, nf_int) + X2GGB(x, nf_int)
!    end select
!    select case(cc_piece)
!    case(cc_VIRT,cc_REALVIRT)
!       res = res - X2GGB(x, nf_int)
!    case(cc_DELTA)
!       res = X2GGC(zero, nf_int)
!    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3gg

  ! ..This is the (regular) pure-singlet splitting functions P_ps^(2).    
  !    P_qq^(2) is obtained by adding the non-singlet quantity P_NS^(2)+. 
  function sf_P3PS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call wae_Error('Exact N3LO splitting functions not yet available.')

    call sf_VogtValidate
    x = exp(-y)
    res = zero

!    select case(cc_piece)
!    case(cc_REAL,cc_REALVIRT)
!       res = X2PSA(x, nf_int) + zero
!    end select
!    select case(cc_piece)
!    case(cc_VIRT,cc_REALVIRT)
!       res = res - zero
!    case(cc_DELTA)
!       res = zero
!    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3PS

  !--------------------------------------------------------------------
  ! this already includes the factor of 2nf
  function sf_P3qg2nf(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call wae_Error('Exact N3LO splitting functions not yet available.')

    call sf_VogtValidate
    x = exp(-y)
    res = zero

!    select case(cc_piece)
!    case(cc_REAL,cc_REALVIRT)
!       res = X2QGA(x, nf_int)
!    end select
!    select case(cc_piece)
!    case(cc_VIRT,cc_REALVIRT)
!       res = res - zero
!    case(cc_DELTA)
!       res = zero
!    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3qg2nf


  function sf_P3gq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call wae_Error('Exact N3LO splitting functions not yet available.')

    call sf_VogtValidate
    x = exp(-y)
    res = zero

!    select case(cc_piece)
!    case(cc_REAL,cc_REALVIRT)
!       res = X2GQA(x, nf_int)
!    end select
!    select case(cc_piece)
!    case(cc_VIRT,cc_REALVIRT)
!       res = res - zero
!    case(cc_DELTA)
!       res = zero
!    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3gq

  !------------------------- non-singlet pieces ------------------------
  function sf_P3NSPlus(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call wae_Error('Exact N3LO splitting functions not yet available.')

    call sf_VogtValidate
    x = exp(-y)
    res = zero

!    select case(cc_piece)
!    case(cc_REAL,cc_REALVIRT)
!       res = X2NSPA(x, nf_int) + X2NSB(x, nf_int)
!    end select
!    select case(cc_piece)
!    case(cc_VIRT,cc_REALVIRT)
!       res = res - X2NSB(x, nf_int)
!    case(cc_DELTA)
!       res = X2NSC(zero, nf_int)
!    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3NSPlus


  function sf_P3NSMinus(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call wae_Error('Exact N3LO splitting functions not yet available.')

    call sf_VogtValidate
    x = exp(-y)
    res = zero

!    select case(cc_piece)
!    case(cc_REAL,cc_REALVIRT)
!       res = X2NSMA(x, nf_int) + X2NSB(x, nf_int)
!    end select
!    select case(cc_piece)
!    case(cc_VIRT,cc_REALVIRT)
!       res = res - X2NSB(x, nf_int)
!    case(cc_DELTA)
!       res = X2NSC(zero, nf_int)
!    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3NSMinus

  !-- according to comments in Vogt code P_S = P_V - P_-
  function sf_P3NSS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call wae_Error('Exact N3LO splitting functions not yet available.')

    call sf_VogtValidate
    x = exp(-y)
    res = zero

!    select case(cc_piece)
!    case(cc_REAL,cc_REALVIRT)
!       res = X2NSSA(x, nf_int)
!    end select
!    select case(cc_piece)
!    case(cc_VIRT,cc_REALVIRT)
!       res = res - zero
!    case(cc_DELTA)
!       res = zero
!    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3NSS

  !----------------------------------------------------------------
  ! The Vogt expressions are valid only for standard colour factors
  ! This routine makes sure that the appropriate conditions hold
  subroutine sf_VogtValidate
    if (tr == half) return
    call wae_Error('sf_VogtValidate: &
         &colour factors must be set to default values', &
         &'in order to use the Vogt splitting function parameterisations')
  end subroutine sf_VogtValidate

end module splitting_functions_n3lo_e


module splitting_functions_n3lo_p
  use types; use consts_dp; use convolution_communicator
  use qcd; use warnings_and_errors
  use xpij2p; use xpns2p
  implicit none
  private

  public :: sf_P3gg, sf_P3qg2nf, sf_P3PS, sf_P3gq
  public :: sf_P3NSPlus, sf_P3NSMinus, sf_P3NSS
  ! these are of help elsewhere in identifying (run-time) the sets of 
  ! splitting functions that have been used?
  public :: name_xpij2, name_xpns2

contains
  function sf_P3gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call wae_Error('Parametrised N3LO splitting functions not yet available.')

    call sf_VogtValidate
    x = exp(-y)
    res = zero

!    select case(cc_piece)
!    case(cc_REAL,cc_REALVIRT)
!       res = P3GGA(x, nf_int) + P3GGB(x, nf_int)
!    end select
!    select case(cc_piece)
!    case(cc_VIRT,cc_REALVIRT)
!       res = res - P3GGB(x, nf_int)
!    case(cc_DELTA)
!       res = P3GGC(zero, nf_int)
!    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3gg

  ! ..This is the (regular) pure-singlet splitting functions P_ps^(2).    
  !    P_qq^(2) is obtained by adding the non-singlet quantity P_NS^(2)+. 
  function sf_P3PS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call wae_Error('Parametrised N3LO splitting functions not yet available.')

    call sf_VogtValidate
    x = exp(-y)
    res = zero

!    select case(cc_piece)
!    case(cc_REAL,cc_REALVIRT)
!       res = P3PSA(x, nf_int) + zero
!    end select
!    select case(cc_piece)
!    case(cc_VIRT,cc_REALVIRT)
!       res = res - zero
!    case(cc_DELTA)
!       res = zero
!    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3PS

  !--------------------------------------------------------------------
  ! this already includes the factor of 2nf
  function sf_P3qg2nf(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call wae_Error('Parametrised N3LO splitting functions not yet available.')

    call sf_VogtValidate
    x = exp(-y)
    res = zero

!    select case(cc_piece)
!    case(cc_REAL,cc_REALVIRT)
!       res = P3QGA(x, nf_int)
!    end select
!    select case(cc_piece)
!    case(cc_VIRT,cc_REALVIRT)
!       res = res - zero
!    case(cc_DELTA)
!       res = zero
!    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3qg2nf


  function sf_P3gq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call wae_Error('Parametrised N3LO splitting functions not yet available.')

    call sf_VogtValidate
    x = exp(-y)
    res = zero

!    select case(cc_piece)
!    case(cc_REAL,cc_REALVIRT)
!       res = P3GQA(x, nf_int)
!    end select
!    select case(cc_piece)
!    case(cc_VIRT,cc_REALVIRT)
!       res = res - zero
!    case(cc_DELTA)
!       res = zero
!    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3gq

  !------------------------- non-singlet pieces ------------------------
  function sf_P3NSPlus(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call wae_Error('Parametrised N3LO splitting functions not yet available.')

    call sf_VogtValidate
    x = exp(-y)
    res = zero

!    select case(cc_piece)
!    case(cc_REAL,cc_REALVIRT)
!       res = P3NSPA(x, nf_int) + P3NSB(x, nf_int)
!    end select
!    select case(cc_piece)
!    case(cc_VIRT,cc_REALVIRT)
!       res = res - P3NSB(x, nf_int)
!    case(cc_DELTA)
!       res = P3NSPC(zero, nf_int)
!    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3NSPlus


  function sf_P3NSMinus(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call wae_Error('Parametrised N3LO splitting functions not yet available.')

    call sf_VogtValidate
    x = exp(-y)
    res = zero

!    select case(cc_piece)
!    case(cc_REAL,cc_REALVIRT)
!       res = P3NSMA(x, nf_int) + P3NSB(x, nf_int)
!    end select
!    select case(cc_piece)
!    case(cc_VIRT,cc_REALVIRT)
!       res = res - P3NSB(x, nf_int)
!    case(cc_DELTA)
!       res = P3NSMC(zero, nf_int)
!    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3NSMinus

  !-- according to comments in Vogt code P_S = P_V - P_-
  function sf_P3NSS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call wae_Error('Parametrised N3LO splitting functions not yet available.')

!    call sf_VogtValidate
!    x = exp(-y)
!    res = zero
!
!    select case(cc_piece)
!    case(cc_REAL,cc_REALVIRT)
!       res = P3NSSA(x, nf_int)
!    end select
!    select case(cc_piece)
!    case(cc_VIRT,cc_REALVIRT)
!       res = res - zero
!    case(cc_DELTA)
!       res = zero
!    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3NSS

  !----------------------------------------------------------------
  ! The Vogt expressions are valid only for standard colour factors
  ! This routine makes sure that the appropriate conditions hold
  subroutine sf_VogtValidate
    if (ca == three .and. tr == half .and. cf == four/three) return
    call wae_Error('sf_VogtValidate: &
         &colour factors must be set to default values', &
         &'in order to use the Vogt splitting function parameterisations')
  end subroutine sf_VogtValidate

end module splitting_functions_n3lo_p



!======================================================================
! interface to guessed 4-loop splitting functions, as parameterised by Vogt
! (see xpgg3a.f90 xpgq3a.f90 xpps3a.f90 xpqg3a.f90 for details).
!
! There are several things to remember: 
!    . Vogt y == my x
!    . Vogt is normalised to powers of (as/4pi)^4 
!    . colour-factor dependence is not available
!    . systematic separation in A,B,C functions needs to be
!      taken into account; C function needs to be called with x=0?x
!    . check inclusion of 2nf in qg piece
! 
!----------------------------------------------------------------------
module splitting_functions_n3lo_n
  use types; use consts_dp; use convolution_communicator
  use qcd; use dglap_choices; use warnings_and_errors
  use xpgg3a; use xpgg3a_2410; use xpgq3a; use xpgq3a_2404; use&
       & xpqg3a; use xpps3a; use xpns3s; use xpns3p; use xpns3m
  implicit none
  private

  public :: sf_P3gg, sf_P3qg2nf, sf_P3PS, sf_P3gq
  public :: sf_P3NSPlus, sf_P3NSMinus, sf_P3NSS
  ! these are of help elsewhere in identifying (run-time) the sets of 
  ! splitting functions that have been used?
  public :: name_xpgg3, name_xpgq3, name_xpqg3, name_xpps3,&
       & name_xpns3s, name_xpns3m, name_xpns3p

contains
  function sf_P3gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    if(n3lo_splitting_approximation .eq. n3lo_splitting_approximation_up_to_2410_08089) then
       select case(cc_piece)
       case(cc_REAL,cc_REALVIRT)
          res = P3GGA_2410(x, nf_int, n3lo_splitting_variant) + P3GGB_2410(x, nf_int)
       end select
       select case(cc_piece)
       case(cc_VIRT,cc_REALVIRT)
          res = res - P3GGB_2410(x, nf_int)
       case(cc_DELTA)
          res = P3GGC_2410(zero, nf_int, n3lo_splitting_variant)
       end select
    else
       select case(cc_piece)
       case(cc_REAL,cc_REALVIRT)
          res = P3GGA(x, nf_int, n3lo_splitting_variant) + P3GGB(x, nf_int)
       end select
       select case(cc_piece)
       case(cc_VIRT,cc_REALVIRT)
          res = res - P3GGB(x, nf_int)
       case(cc_DELTA)
          res = P3GGC(zero, nf_int)
       end select
    endif
    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3gg

  ! ..This is the (regular) pure-singlet splitting functions P_ps^(2).    
  !    P_qq^(2) is obtained by adding the non-singlet quantity P_NS^(2)+. 
  function sf_P3PS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P3PSA(x, nf_int, n3lo_splitting_variant) + zero
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3PS

  !--------------------------------------------------------------------
  ! this already includes the factor of 2nf
  function sf_P3qg2nf(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P3QGA(x, nf_int, n3lo_splitting_variant)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3qg2nf


  function sf_P3gq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       if(n3lo_splitting_approximation .eq. n3lo_splitting_approximation_up_to_2310_05744) then
          res = P3GQA(x, nf_int, n3lo_splitting_variant)
       else
          res = P3GQA_2404(x, nf_int, n3lo_splitting_variant)
       endif
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3gq

  !------------------------- non-singlet pieces ------------------------
  function sf_P3NSPlus(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P3NSPA(x, nf_int, n3lo_splitting_variant) + P3NSPB(x, nf_int, n3lo_splitting_variant)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - P3NSPB(x, nf_int, n3lo_splitting_variant)
    case(cc_DELTA)
       res = P3NSPC(zero, nf_int, n3lo_splitting_variant)
    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3NSPlus


  function sf_P3NSMinus(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P3NSMA(x, nf_int, n3lo_splitting_variant) + P3NSMB(x, nf_int, n3lo_splitting_variant)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - P3NSMB(x, nf_int, n3lo_splitting_variant)
    case(cc_DELTA)
       res = P3NSMC(zero, nf_int, n3lo_splitting_variant)
    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3NSMinus

  !-- according to comments in Vogt code P_S = P_V - P_-
  function sf_P3NSS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x

    call sf_VogtValidate
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = P3NSSA(x, nf_int, n3lo_splitting_variant)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    res = res / 16.0_dp   ! compensate (as/4pi)^4 -> (as/2pi)^4
    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P3NSS

  !----------------------------------------------------------------
  ! The Vogt expressions are valid only for standard colour factors
  ! This routine makes sure that the appropriate conditions hold
  subroutine sf_VogtValidate
    if (ca == three .and. tr == half .and. cf == four/three) return
    call wae_Error('sf_VogtValidate: &
         &colour factors must be set to default values', &
         &'in order to use the Vogt splitting function parameterisations')
  end subroutine sf_VogtValidate

end module splitting_functions_n3lo_n


!======================================================================
! Modules for extracting selecting the relevant splitting functions, 
! either "n" version, "p" versions. May be extended to
! exact versions at some later stage...
!======================================================================
module sf_n3lo_n_rename
  use splitting_functions_n3lo_n,&
       &sf_n_P3gg      => sf_P3gg     ,&
       &sf_n_P3PS      => sf_P3PS     ,&
       &sf_n_P3qg2nf   => sf_P3qg2nf  ,&
       &sf_n_P3gq      => sf_P3gq     ,&
       &sf_n_P3NSPlus  => sf_P3NSPlus ,&
       &sf_n_P3NSMinus => sf_P3NSMinus,&
       &sf_n_P3NSS     => sf_P3NSS
  implicit none
  private
  public :: sf_n_P3gg, sf_n_P3qg2nf, sf_n_P3PS, sf_n_P3gq
  public :: sf_n_P3NSPlus, sf_n_P3NSMinus, sf_n_P3NSS
end module sf_n3lo_n_rename


module sf_n3lo_p_rename
  use splitting_functions_n3lo_p,&
       &sf_p_P3gg      => sf_P3gg     ,&
       &sf_p_P3PS      => sf_P3PS     ,&
       &sf_p_P3qg2nf   => sf_P3qg2nf  ,&
       &sf_p_P3gq      => sf_P3gq     ,&
       &sf_p_P3NSPlus  => sf_P3NSPlus ,&
       &sf_p_P3NSMinus => sf_P3NSMinus,&
       &sf_p_P3NSS     => sf_P3NSS
  implicit none
  private
  public :: sf_p_P3gg, sf_p_P3qg2nf, sf_p_P3PS, sf_p_P3gq
  public :: sf_p_P3NSPlus, sf_p_P3NSMinus, sf_p_P3NSS
end module sf_n3lo_p_rename


module sf_n3lo_e_rename
  use splitting_functions_n3lo_e,&
       &sf_e_P3gg      => sf_P3gg     ,&
       &sf_e_P3PS      => sf_P3PS     ,&
       &sf_e_P3qg2nf   => sf_P3qg2nf  ,&
       &sf_e_P3gq      => sf_P3gq     ,&
       &sf_e_P3NSPlus  => sf_P3NSPlus ,&
       &sf_e_P3NSMinus => sf_P3NSMinus,&
       &sf_e_P3NSS     => sf_P3NSS
  implicit none
  private
  public :: sf_e_P3gg, sf_e_P3qg2nf, sf_e_P3PS, sf_e_P3gq
  public :: sf_e_P3NSPlus, sf_e_P3NSMinus, sf_e_P3NSS
end module sf_n3lo_e_rename


module splitting_functions_n3lo
  use types; use qcd; use dglap_choices
  use warnings_and_errors
  use sf_n3lo_e_rename
  use sf_n3lo_p_rename
  use sf_n3lo_n_rename
  public :: sf_P3gg, sf_P3qg2nf, sf_P3PS, sf_P3gq
  public :: sf_P3NSPlus, sf_P3NSMinus, sf_P3NSS
contains
  
  
  real(dp) function sf_P3gg     (y) result(res)
    real(dp), intent(in) :: y
    select case(n3lo_splitting_variant)
    case(n3lo_splitting_exact)
       res = sf_e_P3gg(y)
    case(n3lo_splitting_param)
       res = sf_p_P3gg(y)
    case(n3lo_splitting_nfitav, n3lo_splitting_nfiterr1, n3lo_splitting_nfiterr2)
       res = sf_n_P3gg(y)
    case default
       call wae_error('splitting_functions_n3lo','unrecognized imod',&
            &intval=n3lo_splitting_variant)
    end select
  end function sf_P3gg
  
  real(dp) function sf_P3PS     (y) result(res)
    real(dp), intent(in) :: y
    select case(n3lo_splitting_variant)
    case(n3lo_splitting_exact)
       res = sf_e_P3PS(y)
    case(n3lo_splitting_param)
       res = sf_p_P3PS(y)
    case(n3lo_splitting_nfitav, n3lo_splitting_nfiterr1, n3lo_splitting_nfiterr2)
       res = sf_n_P3PS(y)
    case default
       call wae_error('splitting_functions_n3lo','unrecognized imod',&
            &intval=n3lo_splitting_variant)
    end select
  end function sf_P3PS
  
  real(dp) function sf_P3qg2nf  (y) result(res)
    real(dp), intent(in) :: y
    select case(n3lo_splitting_variant)
    case(n3lo_splitting_exact)
       res = sf_e_P3qg2nf(y)
    case(n3lo_splitting_param)
       res = sf_p_P3qg2nf(y)
    case(n3lo_splitting_nfitav, n3lo_splitting_nfiterr1, n3lo_splitting_nfiterr2)
       res = sf_n_P3qg2nf(y)
    case default
       call wae_error('splitting_functions_n3lo','unrecognized imod',&
            &intval=n3lo_splitting_variant)
    end select
  end function sf_P3qg2nf
  
  real(dp) function sf_P3gq     (y) result(res)
    real(dp), intent(in) :: y
    select case(n3lo_splitting_variant)
    case(n3lo_splitting_exact)
       res = sf_e_P3gq(y)
    case(n3lo_splitting_param)
       res = sf_p_P3gq(y)
    case(n3lo_splitting_nfitav, n3lo_splitting_nfiterr1, n3lo_splitting_nfiterr2)
       res = sf_n_P3gq(y)
    case default
       call wae_error('splitting_functions_n3lo','unrecognized imod',&
            &intval=n3lo_splitting_variant)
    end select
  end function sf_P3gq
  
  real(dp) function sf_P3NSPlus (y) result(res)
    real(dp), intent(in) :: y
    select case(n3lo_splitting_variant)
    case(n3lo_splitting_exact)
       res = sf_e_P3NSPlus(y)
    case(n3lo_splitting_param)
       res = sf_p_P3NSPlus(y)
    case(n3lo_splitting_nfitav, n3lo_splitting_nfiterr1, n3lo_splitting_nfiterr2)
       res = sf_n_P3NSPlus(y)
    case default
       call wae_error('splitting_functions_n3lo','unrecognized imod',&
            &intval=n3lo_splitting_variant)
    end select
  end function sf_P3NSPlus
  
  real(dp) function sf_P3NSMinus(y) result(res)
    real(dp), intent(in) :: y
    select case(n3lo_splitting_variant)
    case(n3lo_splitting_exact)
       res = sf_e_P3NSMinus(y)
    case(n3lo_splitting_param)
       res = sf_p_P3NSMinus(y)
    case(n3lo_splitting_nfitav, n3lo_splitting_nfiterr1, n3lo_splitting_nfiterr2)
       res = sf_n_P3NSMinus(y)
    case default
       call wae_error('splitting_functions_n3lo','unrecognized imod',&
            &intval=n3lo_splitting_variant)
    end select
  end function sf_P3NSMinus

  real(dp) function sf_P3NSS(y) result(res)
    real(dp), intent(in) :: y
    select case(n3lo_splitting_variant)
    case(n3lo_splitting_exact)
       res = sf_e_P3NSS(y)
    case(n3lo_splitting_param)
       res = sf_p_P3NSS(y)
    case(n3lo_splitting_nfitav, n3lo_splitting_nfiterr1, n3lo_splitting_nfiterr2)
       res = sf_n_P3NSS(y)
    case default
       call wae_error('splitting_functions_n3lo','unrecognized imod',&
            &intval=n3lo_splitting_variant)
    end select
  end function sf_P3NSS
end module splitting_functions_n3lo
