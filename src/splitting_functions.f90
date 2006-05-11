!======================================================================
! Consider the following possibility: want special types for
!      set of LO splitting functions 
!      set of coefficient functions

!----------------------------------------------------------------------
! NB: compared to all normal splitting functions we have multiplied 
! by x, to allow for integration measure. NB: this applies also 
! to the plus functions, which do not care about integration measures?
!
! \int dx f(x)_+  g(x) = \int dy [exp(-y)f(exp(-y))]_+ g(exp(-y))
!
!
module splitting_functions
  use types; use consts_dp; use convolution_communicator
  use coefficient_functions; use qcd; use warnings_and_errors
  !! early estimate of splitting function based on fit to moments
  !!use splitting_functions_nnlo_n
  !! parametrization of exact splitting functions
  !!use splitting_functions_nnlo_p
  ! this one gives a set that depends on imod.
  use splitting_functions_nnlo
  implicit none
  private

  public :: sf_Pgg, sf_Pqq, sf_Pgq, sf_Pqg
  public :: sf_P1qqV, sf_P1gg, sf_P1qqbarV, sf_P1qqS, sf_P1qg, sf_P1gq
  public :: sf_P1qqBryan, sf_P1qgBryan, sf_P1minus
  public :: sf_P1fromg, sf_P1fromq
  public :: sf_P1qqV_DIS, sf_P1qg_DIS

  public :: sf_A2PShq, sf_A2PShg, sf_A2PShg_vogt
  public :: sf_A2NSqq_H, sf_A2Sgg_H, sf_A2Sgq_H

  public :: sf_P2gg, sf_P2qg2nf, sf_P2PS, sf_P2gq
  public :: sf_P2NSPlus, sf_P2NSMinus, sf_P2NSS

  public :: sf_DPqq, sf_DPqg, sf_DPgq, sf_DPgg 
  public :: sf_DP1qqV, sf_DP1qqbarV, sf_DP1qqS
  public :: sf_DP1qg, sf_DP1gq, sf_DP1gg

!!$  ! these are of help elsewhere in identifying (run-time) the sets of 
!!$  ! splitting functions that have been used?
!!$  public :: name_xpij2, name_xpns2

  real(dp),  parameter :: four_thirds = four/three

contains
  !----------------------------------------------------------------------
  function sf_Pgg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = two*ca*(x/(one-x) + (one-x)/x + x*(1-x))
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - two*ca*one/(one-x)
    case(cc_DELTA)
       res = (11.0_dp*ca - four*nf*tr)/6.0_dp
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_Pgg

  !----------------------------------------------------------------------
  function sf_Pqq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = cf*(one+x**2)/(one-x)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - cf*two/(one-x)
    case(cc_DELTA)
       res = cf*three*half
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_Pqq

  !----------------------------------------------------------------------
  function sf_Pgq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = cf*(one + (one-x)**2)/x
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_Pgq

  !----------------------------------------------------------------------
  function sf_Pqg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = tr*(x**2 + (one-x)**2)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_Pqg
  


  !======================================================================
  ! From here onwards, the NLO splitting functions
  function sf_P1qqV(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lnx, ln1mx, pqq
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnx = log(x); ln1mx = log(one - x)
       pqq = two/(one-x) - one - x
!!$       res =    cf**2*( -(two*lnx*ln1mx + 1.5_dp*lnx)*pqq&
!!$            & -(1.5_dp + 3.5_dp*x)*lnx - half*(one+x)*lnx**2 - 5*(1-x))&
!!$            & + cf*ca*( (half*lnx**2 + 11/6.0_dp * lnx&
!!$            &          + 67/18.0_dp - pi**2/6.0_dp) * pqq + (1+x)*lnx&
!!$            &        + 20/three*(1-x))&
!!$            & + cf*tf*(-(two/three*lnx + 10/9.0_dp)*pqq - four/three*(1-x))
res=    CF*Tf*((-1.1111111111111112_dp - (2*lnx)/3._dp)*pqq -          &
     &     (4*(1 - x))/3._dp) +                                        &
     &  CA*CF*((3.7222222222222223_dp + (11*lnx)/6._dp + lnx**2/2._dp -& 
     &        Pi**2/6._dp)*pqq + (20*(1 - x))/3._dp + lnx*(1 + x)) +   &
     &  CF**2*(((-3*lnx)/2._dp - 2*ln1mx*lnx)*pqq - 5*(1 - x) -        &
     &     (lnx**2*(1 + x))/2._dp - lnx*(1.5_dp + (7*x)/2._dp))
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       pqq = -two/(1-x)
       res = res + CA*CF*(3.7222222222222223_dp - Pi**2/6._dp)*pqq - (10*CF&
           & *pqq*Tf)/9._dp
    case(cc_DELTA)
       res = -(CF*(0.16666666666666666_dp + (2*Pi**2)/9._dp)*Tf) + CA&
           & *CF*(0.7083333333333334_dp + (11*Pi**2)/18._dp - 3*zeta3) +&
           & CF**2*(0.375_dp - Pi**2/2._dp + 6*zeta3)
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function sf_P1qqV


  !======================================================================
  ! DIS version of the above.
  ! Currently the formulae are wrong. Use the complete formulae from 
  ! module dglap_holders (where double convolutions are used to work them out).
  function sf_P1qqV_DIS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    
    write(0,*) 'sf_P1qqV_DIS: DIS scheme splitting functions &
         &currently not supported' 
    stop
    res = sf_P1qqV(y) - cf_CqF2MSbar(y)*(11*CA - 2*nf)/6.0_dp
  end function sf_P1qqV_DIS

  

  !----------------------------------------------------------------------
  function sf_P1gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lnx, ln1mx, pgg, S2x, pggmx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnx = log(x); ln1mx = log(one - x)
       pgg   = (one/(one-x) + one/x -two + x*(one-x))
       pggmx = (one/(one+x) - one/x -two - x*(one+x))
       S2x   = sf_S2(x)
       res = CF*Tf*(-16 + 4/(3._dp*x) + 8*x + (20*x**2)/3._dp - lnx**2*(2&
           & + 2*x)- lnx*(6 + 10*x)) + CA*Tf*(2 - (20*pgg)/9._dp - 2*x -&
           & (4*lnx*(1 + x))/3._dp + (26*(-(1/x) + x**2))/9._dp) + CA**2&
           & *(pgg*(7.444444444444445_dp - 4*ln1mx*lnx + lnx**2 - Pi**2&
           & /3._dp) + 2*pggmx*S2x + (27*(1 - x))/2._dp + 4*lnx**2*(1 + x)&
           & + (67*(-(1/x) + x**2))/9._dp - lnx*(8.333333333333334_dp -&
           & (11*x)/3._dp + (44*x**2)/3._dp))
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       pgg = -one/(1-x)
       res = res + CA**2*pgg*(7.444444444444445_dp - Pi**2/3._dp) - (20*CA&
           & *pgg*Tf)/9._dp
    case(cc_DELTA)
       res = (-4*CA*Tf)/3._dp - CF*Tf + CA**2*(2.6666666666666665_dp&
           & + 3*zeta3)
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function sf_P1gg


  !----------------------------------------------------------------------
  function sf_P1qqbarV(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lnx, pqqmx, S2x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnx = log(x)
       pqqmx = two/(one+x) - one + x
       S2x   = sf_S2(x)
       res = CF*(-CA/2._dp + CF)*(2*pqqmx*S2x + 4*(1 - x) + 2*lnx*(1 + x))
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res + zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P1qqbarV

  !----------------------------------------------------------------------
  function sf_P1qqS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lnx, ln1mx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnx = log(x); ln1mx = log(one - x)
       res = (CF*TR*(20 - 9*(2 - lnx + lnx**2)*x - 9*(-6 - 5*lnx + lnx&
           & **2)*x**2 + 8*(-7 + 3*lnx)*x**3))/(9._dp*x) 
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res + zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P1qqS


  !----------------------------------------------------------------------
  function sf_P1qg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lnx, ln1mx, pqg, pqgmx, S2x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnx = log(x); ln1mx = log(one - x)
       pqg   = x**2 + (one-x)**2
       pqgmx = x**2 + (one+x)**2
       S2x   = sf_S2(x)
       res = CF*(half*TR)*(4 + 4*ln1mx + (10 - 4*(ln1mx - lnx) + 2*(-ln1mx&
            & + lnx)& 
           & **2 - (2*Pi**2)/3._dp)* pqg - lnx*(1 - 4*x) - lnx**2*(1 - 2&
           & *x) - 9*x) + CA*(half*TR)*(20.22222222222222_dp - 4*ln1mx + (&
           & -24.22222222222222_dp + 4*ln1mx - 2*ln1mx**2 + (44*lnx)/3._dp&
           & - lnx**2 + Pi**2/3._dp)*pqg + 2*pqgmx*S2x + 40/(9._dp*x) +&
           & (14*x)/9._dp - lnx**2*(2 + 8*x) + lnx*(-12.666666666666666_dp&
           & + (136*x)/3._dp))
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res + zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P1qg


  !======================================================================
  ! DIS version of the above
  function sf_P1qg_DIS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    
    write(0,*) 'sf_P1qg_DIS: DIS scheme splitting functions currently not supported' 
    stop
    res = sf_P1qg(y) - cf_CgF2MSbar(y)*(11*CA - 2*nf)/6.0_dp
  end function sf_P1qg_DIS


  !----------------------------------------------------------------------
  function sf_P1gq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lnx, ln1mx, pgq, pgqmx,S2x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnx = log(x); ln1mx = log(one - x)
       pgq   = (one + (one-x)**2)/x
       pgqmx = -(one + (one+x)**2)/x
       S2x   = sf_S2(x)
       res = CF*Tf*(-((2.2222222222222223_dp + (4*ln1mx)/3._dp)*pgq) - (4&
           & *x)/3._dp) + CF**2*(-2.5_dp - (3*ln1mx + ln1mx**2)*pgq - lnx&
           & **2*(1 - x/2._dp) - (7*x)/2._dp - 2*ln1mx*x + lnx*(2 + (7*x)&
           & /2._dp)) + CA*CF*(3.111111111111111_dp + pgq*(0.5_dp + (11&
           & *ln1mx)/3._dp + ln1mx**2 - 2*ln1mx*lnx + lnx**2/2._dp - Pi**2&
           & /6._dp) + pgqmx*S2x + (65*x)/18._dp + 2*ln1mx*x + (44*x**2)&
           & /9._dp + lnx**2*(4 + x) - lnx*(12 + 5*x + (8*x**2)/3._dp))
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res + zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_P1gq


  !======================================================================
  ! Alternative versions of the splitting functions, more similar
  ! to what is in the ESW book
  ! and some things for checking sum rules more easily
  function sf_P1qqBryan(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    res = sf_P1qqV(y) + sf_P1qqbarV(y) + two*nf * sf_P1qqS(y)
    !res = sf_P1qqBryanTyped(y)
  end function sf_P1qqBryan
  function sf_P1qgBryan(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    res = two*nf * sf_P1qg(y)
  end function sf_P1qgBryan
  function sf_P1minus(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    res = sf_P1qqV(y) - sf_P1qqbarV(y)
  end function sf_P1minus
  function sf_P1fromq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    res = sf_P1qqBryan(y) + sf_P1gq(y)
  end function sf_P1fromq
  function sf_P1fromg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    res = sf_P1gg(y) + sf_P1qgBryan(y)
  end function sf_P1fromg


  !----------------------------------------------------------------------
  ! Following was useful for debugging purposes...
  function sf_P1qqBryanTyped(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lnx, ln1mx, pqq, pqqmx, S2x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnx = log(x); ln1mx = log(one - x)
       pqq = two/(one-x) - one - x
       pqqmx = two/(one+x) - one + x
       S2x = sf_S2(x)
       res = CA*CF*((3.7222222222222223_dp + (11*lnx)/6._dp + lnx**2/2._dp&
            & - Pi**2/6._dp)*pqq - pqqmx*S2x + (14*(1 - x))/3._dp) +&
            &  CF**2*(-1 - ((3*lnx)/2._dp + 2*ln1mx*lnx)*pqq + 2*pqqmx*S2x + &
            &     lnx*(0.5_dp - (3*x)/2._dp) + x - (lnx**2*(1 + x))/2._dp) +&
            & CF*Tf*(-5.333333333333333_dp - &
            &     (1.1111111111111112_dp + (2*lnx)/3._dp)*pqq + 40/(9._dp*x) +&
            &     (40*x)/3._dp - (112*x**2)/9._dp - 2*lnx**2*(1 + x) +&
            &     lnx*(2 + 10*x + (16*x**2)/3._dp))

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       pqq = -two/(1-x)
       res = res + CA*CF*(3.7222222222222223_dp - Pi**2/6._dp)*pqq - (10*CF&
           & *pqq*Tf)/9._dp
    case(cc_DELTA)
       res = -(CF*(0.16666666666666666_dp + (2*Pi**2)/9._dp)*Tf) + CA&
           & *CF*(0.7083333333333334_dp + (11*Pi**2)/18._dp - 3*zeta3) +&
           & CF**2*(0.375_dp - Pi**2/2._dp + 6*zeta3)
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function sf_P1qqBryanTyped

  
  
  !======================================================================
  ! helper...
  function sf_S2(x) result(S2)
    use special_functions
    real(dp), intent(in) :: x
    real(dp)             :: S2, lnx
    lnx = log(x)
    S2 = -two*ddilog(-x) + half * lnx**2 - two*lnx*log(one+x) - pi**2&
         & /6.0_dp 
  end function sf_S2

  
  !======================================================================
  ! What follows is stuff for VFNS.
  ! It has been hacked in one way or another from the code provided
  ! by W.L. van Neerven (~/src/HEP/vannerven-vfns.f), and appropriate 
  ! citations are
  ! 
  !   M. Buza, Y. Matiounine, J. Smith, R. Migneron, W.L. van Neerven, 
  !   Nucl. Phys. B472 (1996) 611, hep-ph/9601302.
  !   M. Buza, Y. Matiounine, J. Smith, W.L. van Neerven,
  !   Eur. Phys. J. C1 (1998) 301. 
  !
  ! Pieces that are going to be needed for my purposes are those from that
  ! code which are free of logs of mu^2/m^2, since we will for the time being
  ! keep these two scales equal (and later on, if we do not then pieces
  ! with logs will be reconstructed appropriately?)
  !
  ! Thus the pieces needed, in the notation of the second of the above 
  ! papers, are:
  !
  ! For Delta (f_k + f_kbar):
  !   A^{2,NS}_{qq,H}    CODED
  !   A^{2,PS}_{qq,H}    <-- this is zero at O(as^2)
  !   A^{2,S}_{qg,H}     <-- this is zero at O(as^2)
  !
  ! For Delta (f_H + f_Hbar)
  !   A^{2,PS}_{Hq}      CODED
  !   A^{2,S}_{Hg}       CODED
  !
  ! For Delta (G)
  !   A^{2,S}_{gq,H}     CODED
  !   A^{2,S}_{gg,H}     CODED
  !
  ! NB: currently using cernlib version of WGPLG. This should be hacked
  !     out of CERNLIB and inserted explicitly into special functions.
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! This one hacked out of vanneerven-vfns.f90
  function sf_A2Sgq_H(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: z
    real(dp) :: ln1mz
    z = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       ln1mz = log(one - z)
       res = (four*(two/z-two+z)*ln1mz**2/3.0_dp+8.0_dp*(10.0_dp/z &
            &-10.0_dp+8.0_dp*z)*ln1mz/9.0_dp+(448.0_dp/z-448.0_dp+344.0_dp*z)&
            &/27.0_dp)*cf*tr
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res + zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * z
    !-- recall that above results are for (as/4pi)^2, whereas
    !   we actually use (as/2pi)^2
    res = res * 0.25_dp
  end function sf_A2Sgq_H

  !-----------------------------------------------------------------
  ! This one typed in by hand from article
  function sf_A2Sgg_H(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: z
    real(dp) :: lnz, lnz2, lnz3,ln1mz
    z = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnz = log(z); lnz2 = lnz*lnz; lnz3 = lnz2*lnz
       ln1mz = log(one - z)
       res = CF*TR*(four_thirds*(1+z)*lnz3+(6+10*z)*lnz2+(32+48*z)*lnz&
            &-8/z+80-48*z-24*z**2) + &
            &CA*TR*((four_thirds*(1+z)*lnz2+(52+88*z)*lnz/9.0_dp&
            &       -four_thirds*z*ln1mz)&
            &     + (224/(1-z)+556/z-628+548*z-700*z**2)/27.0_dp )
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - CA*TR*224.0_dp/(27.0_dp*(1-z))
    case(cc_DELTA)
       res = -15*CF*TR + CA*TR*10.0_dp/9.0_dp
    end select

    if (cc_piece /= cc_DELTA) res = res * z
    !-- recall that above results are for (as/4pi)^2, whereas
    !   we actually use (as/2pi)^2
    res = res * 0.25_dp
  end function sf_A2Sgg_H


  !-----------------------------------------------------------------
  ! This one partially hacked out of vanneerven-vfns.f90
  function sf_A2NSqq_H(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: z
    real(dp) :: lnz, lnz2, ln1mz
    z = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnz = log(z); lnz2 = lnz*lnz
       ln1mz = log(one - z)
       res = ((one+Z*Z)*(2.0D0*lnz2/3.0D0+20.0D0*lnz/9.0D0)&
            &/(one-Z)+8.0D0*(one-Z)*lnz/3.0D0+44.0D0/27.0D0 &
            &-268.0D0*Z/27.0D0 + 224.0_dp/(27.0_dp*(1-z)) )*CF*TR
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - (224.0_dp/(27.0_dp*(1-z)) )*CF*TR
    case(cc_DELTA)
       res = CF*TR*(-8*zeta3/3.0_dp + 40.0_dp*zeta2/9.0_dp+73.0_dp/18.0_dp)
    end select

    if (cc_piece /= cc_DELTA) res = res * z
    !-- recall that above results are for (as/4pi)^2, whereas
    !   we actually use (as/2pi)^2
    res = res * 0.25_dp
  end function sf_A2NSqq_H


  !-----------------------------------------------------------------
  ! This one largely hacked out of vanneerven-vfns.f90
  function sf_A2PShg(y) result(res)
    use special_functions
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: z
    real(dp) :: lnz, lnz2, lnz3, ln1mz, ln1mz2, ln1mz3
    real(dp) :: ln1pz, ln1pz2
    real(dp) :: S121MZ, S12MZ,S211MZ,S21MZ,S111MZ,S11MZ
    real(dp) :: A01, A02, B01, B02
    !complex(dp) :: WGPLG
    z = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       ! these will need to be sorted out properly
       S121MZ=WGPLG(1,2,1.0D0-Z) 
       S12MZ=WGPLG(1,2,-Z) 
       S211MZ=WGPLG(2,1,1.0D0-Z) 
       S21MZ=WGPLG(2,1,-Z) 
       S111MZ=WGPLG(1,1,1.0D0-Z) 
       S11MZ=WGPLG(1,1,-Z) 

       lnz = log(z); lnz2 = lnz*lnz; lnz3 = lnz2*lnz
       ln1mz = log(one - z)
       ln1mz2=ln1mz*ln1mz 
       ln1mz3=ln1mz2*ln1mz 
       ln1pz=log(1.0d0+z) 
       ln1pz2=ln1pz*ln1pz 

       ! C_F.T_r  PART                                                        
       A01=(1-2._dp*z+2._dp*z*z)*(8._dp*zeta3+4._dp*ln1mz3/3._dp &
            &-8._dp*ln1mz*s111mz+8._dp*zeta2*lnz-4._dp*lnz*ln1mz2+2._dp*lnz3 &
            &/3._dp-8._dp*lnz*s111mz+8._dp*s211mz-24._dp*s121mz)               
       A02=-(4._dp+96._dp*z-64._dp*z*z)*s111mz-(4._dp-48._dp*z &
            &+40._dp*z*z)*zeta2-(8._dp+48._dp*z-24._dp*z*z)*lnz*ln1mz &
            &+(4._dp+8._dp*z-12._dp*z*z)*ln1mz2-(1._dp+12._dp*z-20._dp*z*z) &
            &*lnz2-(52._dp*z-48._dp*z*z)*ln1mz-(16._dp+18._dp*z+48._dp*z*z) &
            &*lnz+26._dp-82._dp*z+80._dp*z*z+z*z*(-16._dp*zeta2*lnz &
            &+4._dp*lnz3/3._dp+ 16._dp*lnz*s111mz+ 32._dp*s121mz)              

       ! c_a.t_r  part                                                        
       B01=(1._dp-2._dp*z+2._dp*z*z)*(-4._dp*ln1mz3/3._dp+8._dp*ln1mz &
            &*s111mz-8._dp*s211mz)+(1._dp+2._dp*z+2._dp*z*z)*(-8._dp*zeta2 &
            &*ln1pz-16._dp*ln1pz*s11mz-8._dp*lnz*ln1pz2+&
            &4._dp*lnz2*ln1pz+8._dp*lnz &
            &*s11mz-8._dp*s21mz-16._dp*s12mz)+(16._dp+64._dp*z)*(2._dp*s121mz &
            &+lnz*s111mz)-(4._dp+8._dp*z)*lnz3/3._dp+(8._dp-32._dp*z &
            &+16._dp*z*z)*zeta3-(16._dp+64._dp*z)*zeta2*lnz                    
       B02=(16._dp*z+16._dp*z*z)*(s11mz+lnz*ln1pz)+(32._dp/z/3._dp+12._dp &
            &+64._dp*z-272._dp*z*z/3._dp)*s111mz-(12._dp+48._dp*z &
            &-260._dp*z*z/3._dp+32._dp/z/3._dp)*zeta2-4._dp*z*z*lnz*ln1mz &
            &-(2._dp+8._dp*z-10._dp*z*z)*ln1mz2+&
            &(2._dp+8._dp*z+46._dp*z*z/3._dp)&
            &*lnz2+(4._dp+16._dp*z-16._dp*z*z)*ln1mz-(56._dp/3._dp+172._dp*z &
            &/3._dp+1600._dp*z*z/9._dp)*lnz-448._dp/z/27._dp-4._dp/3._dp &
            &-628._dp*z/3._dp+6352._dp*z*z/27._dp                              

       res = TR*(CF*(A01+A02) + CA*(B01+B02))
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res 
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * z
    !-- recall that above results are for (as/4pi)^2, whereas
    !   we actually use (as/2pi)^2
    res = res * 0.25_dp
  end function sf_A2PShg

  !-----------------------------------------------------------------
  ! This one will use Vogts parameterisation
  function sf_A2PShg_vogt(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: z
    real(dp) :: A2HGA, A2HGC
    z = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = A2HGA(z)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res 
    case(cc_DELTA)
       res = A2HGC(zero)
    end select

    if (cc_piece /= cc_DELTA) res = res * z
    !-- recall that above results are for (as/4pi)^2, whereas
    !   we actually use (as/2pi)^2
    res = res * 0.25_dp
  end function sf_A2PShg_vogt


  !-----------------------------------------------------------------
  ! This one largely hacked out of vanneerven-vfns.f90
  function sf_A2PShq(y) result(res)
    use special_functions
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: z
    real(dp) :: lnz, lnz2, lnz3
    real(dp) :: S121MZ, S111MZ
    real(dp) :: A0
    !complex(dp) :: WGPLG
    z = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       ! these will need to be sorted out properly
       S121MZ=WGPLG(1,2,1.0D0-Z) 
       S111MZ=WGPLG(1,1,1.0D0-Z) 

       lnz = log(z); lnz2 = lnz*lnz; lnz3 = lnz2*lnz

       ! C_F.T_r  PART                                                        
       A0=(1._dp+Z)*(32._dp*S121MZ+16._dp*lnz*S111MZ-16._dp*ZETA2 &
            &*lnz-4._dp*lnz3/3._dp)+(32._dp/Z/3._dp+8._dp-8._dp*Z-32._dp &
            &*Z*Z/3._dp)*S111MZ+(-32._dp/Z/3._dp-8._dp+8._dp*Z &
            &+32._dp*Z*Z/3._dp)*ZETA2+(2._dp+10._dp*Z+16._dp*Z*Z/3._dp) &
            &*lnz2-(56._dp/3._dp+88._dp*Z/3._dp+448._dp*Z*Z/9._dp)*lnz &
            &-448._dp/Z/27._dp-4._dp/3._dp-124._dp*Z/3._dp+1600._dp*Z*Z &
            &/27._dp                                                           

       res = TR*CF*A0
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res 
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * z
    !-- recall that above results are for (as/4pi)^2, whereas
    !   we actually use (as/2pi)^2
    res = res * 0.25_dp
  end function sf_A2PShq
  

  !======================================================================
  ! polarized splitting functions...
  !======================================================================
  ! LO: from ESW.
  ! Original reference: Altarelli & Parisi NPB126 (1977) 298
  !----------------------------------------------------------------------
  function sf_DPgg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = ca*(two/(one-x)  - four*x + two)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - two*ca/(one-x)
    case(cc_DELTA)
       res = (11.0_dp*ca - four*nf*tr)/6.0_dp
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_DPgg

  !----------------------------------------------------------------------
  function sf_DPqq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = cf*(two/(1-x) - 1 - x)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - cf*two/(one-x)
    case(cc_DELTA)
       res = cf*three*half
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_DPqq

  !----------------------------------------------------------------------
  function sf_DPgq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = cf*(two-x)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_DPgq

  !----------------------------------------------------------------------
  function sf_DPqg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = tr*(two*x - one)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_DPqg

  !======================================================================
  ! NLO spin-dependent splitting functions

  ! THE CALCULATION OF THE TWO LOOP SPIN SPLITTING FUNCTIONS P(IJ)(1)(X).
  ! By R. Mertig (NIKHEF, Amsterdam), W.L. van Neerven (Leiden U.),. 
  ! INLO-PUB-6-95, NIKHEF-H-95-031, Jun 1995. 33pp.
  ! Published in Z.Phys.C70:637-654,1996
  ! e-Print Archive: hep-ph/9506451

  ! The spin dependent two loop splitting functions
  ! W. Vogelsang (Rutherford),. RAL-TR-96-020, Mar 1996. 25pp.
  ! Published in Nucl.Phys.B475:47-72,1996
  ! e-Print Archive: hep-ph/9603366

  ! Will use the Vogelsang paper for input. His convention coincides
  ! with the one used above with regards to alphas/two and derivative
  ! wrt ln Q^2. However it differs in use of nf factors. Will stay
  ! consistent with the unpolarized case (i.e. not include 2nf factors).
  !
  ! Tests carried out (NLL): comparison to omega=0 (N=1) momenta, eqs.54. 
  ! All work out (after fixing a line in wrong place), though since Pqg->0
  ! one cannot check normalisation of this one...
  !
  ! He also provides coefficient functions, should one 
  ! wish to implement them...

  !--------------------------------------------------
  ! identical to P1qqV
  function sf_DP1qqV(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    res = sf_P1qqV(y)
  end function sf_DP1qqV

  !--------------------------------------------------
  ! identical to -P1qqbarV
  function sf_DP1qqbarV(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    res = -sf_P1qqbarV(y)
  end function sf_DP1qqbarV
  
  !----------------------------------------------------------------------
  ! REMEMBER: 2nf factor NOT included!
  function sf_DP1qqS(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lnx
    x = exp(-y)
    res = zero
    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnx = -y
       res = two*CF*TR*((1-x) - (1-3*x)*lnx - (1+x)*lnx**2)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
    case(cc_DELTA)
       res = zero
    end select
    res = res * half ! Tf->TR accounts for nf; this accounts for 2

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_DP1qqS
  
  !----------------------------------------------------------------------
  ! REMEMBER: 2nf factor NOT included!
  function sf_DP1qg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lnx, ln1mx, dpqg, dpqgmx, S2x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnx = -y; ln1mx = log(one - x)
       dpqg   =  two*x - one
       dpqgmx = -two*x - one
       S2x   = sf_S2(x)
       res = CF*TR*(-22 + 27*x - 9*lnx + 8*(1-x)*ln1mx&
         &         + dpqg*(2*ln1mx**2 - 4*ln1mx*lnx &
         &                + lnx**2 - two/three*pisq))&
         & + CA*TR*((24-22*x) - 8*(1-x)*ln1mx + (2+16*x)*lnx  &
         &         - 2*(ln1mx**2-pisq/6.0_dp)*dpqg &
         &         - (2*S2x - 3*lnx**2) * dpqgmx)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res + zero
    case(cc_DELTA)
       res = zero
    end select
    res = res * half ! Tf->TR accounts for nf; this accounts for 2

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_DP1qg


  !----------------------------------------------------------------------
  function sf_DP1gq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lnx, ln1mx, dpgq, dpgqmx,S2x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnx = -y; ln1mx = log(one - x)
       dpgq   = 2 - x
       dpgqmx = 2 + x
       S2x   = sf_S2(x)
       res =   CF*Tf*(-four/9.0_dp*(x+4) - four/three*dpgq*ln1mx)&
            &+ CF**2*(-half - half*(4-x)*lnx - dpgqmx*ln1mx&
            &        + (-4 - ln1mx**2 + half*lnx**2)*dpgq)&
            &+ CF*CA*((4-13*x)*lnx + (10+x)*ln1mx/three + (41+35*x)/9.0_dp&
            &        + half*(-2*S2x + 3*lnx**2)*dpgqmx &
            &        + (ln1mx**2 - 2*ln1mx*lnx - pisq/6.0_dp)*dpgq)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res + zero
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x
  end function sf_DP1gq


  !----------------------------------------------------------------------
  ! Vogelsang says (p.10) that delta-function parts are just those
  ! of 
  function sf_DP1gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: lnx, ln1mx, S2x, dpgg, dpggmx
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       lnx = -y; ln1mx = log(one - x)
       dpgg   = one/(1-x) - 2*x + 1
       dpggmx = one/(1+x) + 2*x + 1
       S2x   = sf_S2(x)
       res = -CA*Tf*(4*(1-x) + four/three*(1+x)*lnx + 20.0_dp/9.0_dp*dpgg)&
            &-CF*Tf*(10*(1-x) + 2*(5-x)*lnx + 2*(1+x)*lnx**2)&
            &+CA**2*((29-67*x)*lnx/three - 9.5_dp*(1-x) + 4*(1+x)*lnx**2&
            &       - 2*S2x*dpggmx  + (67.0_dp/9.0_dp - 4*ln1mx*lnx &
            &                          + lnx**2 - Pi**2/3._dp)*dpgg)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       dpgg = -one/(1-x)
       res = res + CA**2*dpgg*(67.0_dp/9.0_dp - Pi**2/3._dp) - (20*CA&
           & *dpgg*Tf)/9._dp
    case(cc_DELTA)
       res = (-4*CA*Tf)/3._dp - CF*Tf + CA**2*(8.0_dp/three&
           & + 3*zeta3)
    end select

    if (cc_piece /= cc_DELTA) res = res * x
    
  end function sf_DP1gg

end module splitting_functions





