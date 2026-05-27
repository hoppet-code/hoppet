!HoppetSTARTHEADER
!
! Copyright (c) 2008-2026, Hoppet Authors
!
!----------------------------------------------------------------------
! This file is part of Hoppet <https://github.com/hoppet-code/hoppet>.
!
! Hoppet is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! Hoppet is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details
! <http://www.gnu.org/licenses/>.
! ----------------------------------------------------------------------
!
!HoppetENDHEADER
module xpns3s_2604_fit
  character(len=*), parameter :: name_xpns3s = "xpns3s_2604_fit"
contains
!
! ..The 4-loop MSbar splitting functions P_ns^(3)s, which contributes
!    to the N^3LO evolution of the total valence quark distribution. 
!    The expansion parameter is alpha_s/(4 pi), the scale mu_r = mu_f.
!
!    The expressions have been extracted from the Mathematica files
!    that form the ancilliary materials of 2604.09534. Permission to
!    extract and include the code in hoppet in the form of this
!    Fortran file has been explicitly granted by the authors of
!    2604.09534. The Fortran code and its correctness are the sole
!    responsibility of the hoppet authors.  
!
! =====================================================================
!
!
! ..The regular piece of P_ns^(3)s. 
!
  FUNCTION P3NSSA_2604_fit (Y, NF)
!
       IMPLICIT REAL*8 (A-Z)
       INTEGER nf
       DOUBLE PRECISION NC, zeta2, zeta3, zeta4, zeta5, zeta6
       PARAMETER (NC = 3.0D0)
       PARAMETER (zeta2 = 1.644934066848226D0)
       PARAMETER (zeta3 = 1.2020569031595942D0)
       PARAMETER (zeta4 = 1.0823232337111381D0)
       PARAMETER (zeta5 = 1.03692775514337D0)
       PARAMETER (zeta6 = 1.0173430619844492D0)
       DOUBLE PRECISION CA, CF, D33C, D4FF, D4FA
       PARAMETER (CA = NC, CF = (NC**2 - 1.0D0)/(2*NC))
       PARAMETER (D33C = CF * (NC**2 - 4.0D0)/(8.0D0*NC))
       PARAMETER (D4FF = CF * (NC**4 - 6.0D0*NC**2 + 18.0D0)/(48.0D0*NC**2))
       PARAMETER (D4FA = CF * (NC**2 + 6.0D0)/24.0D0*NC)
       
       X   = Y
       Y1  = 1.D0-Y
       DM  = 1.D0/Y1
       DL  = LOG (Y)
       DL1 = LOG (1.D0-Y)

       P3NSSA_2604_fit =-1.4222222222222223d0*D33C*NF**2&
            &*(3179.999999842691d0+2070.d0*DL+840.d0*DL**2&
            &+29.999999999999996d0*DL**3+14.999999999999998d0*DL**4&
            &+1.d0*DL**5-3437.5783515028093d0*x+1264.2848022004646d0&
            &*DL*x-310.27814976775085d0*DL**2*x-49.352612926730075d0&
            &*DL**3*x+8988.018403776245d0*x**2+3966.8267661463947d0&
            &*DL*x**2+913.6366713105419d0*DL**2*x**2&
            &+0.8680378274566417d0*DL**3*x**2-7559.778322888116d0*x&
            &**3+5287.284923979491d0*DL*x**3-599.7770423623881d0*DL&
            &**2*x**3+215.0377301522394d0*DL**3*x**3&
            &-1123.0291826790135d0*x**4+137.81857642873774d0*x**5&
            &-30.818711230975914d0*x**6+6.523460433077803d0*x**7&
            &-0.9464917377207254d0*x**8-49.91119182567704d0*DL1*Y1&
            &+8.692442624682197d0*DL1**2*Y1-0.00009798699021313083d0&
            &*DL1**3*Y1-88.79603482395439d0*DL1*Y1**2&
            &-12.707906039071162d0*DL1**2*Y1**2-0.7053534132693796d0&
            &*DL1**3*Y1**2-284.67730327004205d0*DL1*Y1**3&
            &+87.26475424064846d0*DL1**2*Y1**3-19.286702545333284d0&
            &*DL1**3*Y1**3-1635.d0*zeta2+389.99999999999994d0*DL&
            &*zeta2-540.d0*DL**2*zeta2-9.999999999999998d0*DL**3&
            &*zeta2+2250.d0*zeta3-959.9999999999999d0*DL*zeta3+180.d0&
            &*DL**2*zeta3-900.d0*zeta2*zeta3+389.99999999999994d0&
            &*zeta4+75.d0*DL*zeta4+1140.d0*zeta5)&
            &+2.1333333333333333d0*CF*D33C*NF*(8639.9999986435d0&
            &+6240.d0*DL+3720.d0*DL**2+340.d0*DL**3+132.5d0*DL**4&
            &-4.d0*DL**5+1.d0*DL**6+12467.856293687217d0*x&
            &+1678.7778453070252d0*DL*x+1304.1618636684277d0*DL**2*x&
            &+407.2392959093284d0*DL**3*x-73.98833128776009d0*DL**4*x&
            &-4074.719922261669d0*x**2+16683.6261152328d0*DL*x**2&
            &+8048.612799026701d0*DL**2*x**2+2547.100316344871d0*DL&
            &**3*x**2+281.0764383851888d0*DL**4*x**2&
            &-4093.867505435656d0*x**3-14934.4742079253d0*DL*x**3&
            &+6320.633206198991d0*DL**2*x**3-1497.809100697404d0*DL&
            &**3*x**3+229.1116674465202d0*DL**4*x**3&
            &+2683.6260320595343d0*x**4-611.3140033665157d0*x**5&
            &+165.74533523525392d0*x**6-38.53816747049772d0*x**7&
            &+6.4298898277845264d0*x**8+270.3849028263012d0*DL1*Y1&
            &+10.641927555198412d0*DL1**2*Y1+0.0001993854690794181d0&
            &*DL1**3*Y1+244.0207129826757d0*DL1*Y1**2&
            &+64.50576153009764d0*DL1**2*Y1**2+2.491086677195316d0&
            &*DL1**3*Y1**2+2247.8169296850415d0*DL1*Y1**3&
            &-728.7366310484178d0*DL1**2*Y1**3+76.43953599781163d0&
            &*DL1**3*Y1**3-17250.d0*zeta2-7080.d0*DL*zeta2-3975.d0*DL&
            &**2*zeta2+80.d0*DL**3*zeta2-95.d0*DL**4*zeta2-6540.d0&
            &*zeta3-5880.d0*DL*zeta3-240.d0*DL**2*zeta3-320.d0*DL**3&
            &*zeta3+2760.d0*zeta2*zeta3+3600.d0*DL*zeta2*zeta3&
            &+1800.d0*zeta3**2+12322.5d0*zeta4+900.d0*DL*zeta4&
            &+2355.d0*DL**2*zeta4-60.d0*zeta5+360.d0*DL*zeta5-237.5d0&
            &*zeta6)-3.2d0*CA*D33C*NF*(30573.333332461574d0&
            &+12033.333333333332d0*DL+2986.6666666666665d0*DL**2&
            &+460.d0*DL**3+73.33333333333333d0*DL**4&
            &-1.7777777777777777d0*DL**5+1.d0*DL**6&
            &+8858.35349762344d0*x+368.89406415869115d0*DL*x&
            &+529.3918203312218d0*DL**2*x+220.3674188539495d0*DL**3*x&
            &-24.65421913995769d0*DL**4*x-1285.121629737188d0*x**2&
            &+5186.062656756615d0*DL*x**2+1498.1714427487636d0*DL**2&
            &*x**2+682.6883864278088d0*DL**3*x**2+136.7178580600142d0&
            &*DL**4*x**2-1370.593556908863d0*x**3&
            &+1.2575674635955227d0*DL*x**3+707.4382501483999d0*DL**2&
            &*x**3-254.90622341238833d0*DL**3*x**3&
            &-35.14154037628444d0*DL**4*x**3-1538.1025557849434d0*x&
            &**4+605.3980492062203d0*x**5-178.85441124633255d0*x**6&
            &+39.48169708952239d0*x**7-4.672957967347179d0*x**8&
            &+183.18897754975504d0*DL1*Y1+3.1255594976234833d0*DL1**2&
            &*Y1+0.0005041553853719101d0*DL1**3*Y1&
            &+435.4188232704969d0*DL1*Y1**2+58.27606439133409d0*DL1&
            &**2*Y1**2+3.297423360635865d0*DL1**3*Y1**2&
            &+1410.573774627683d0*DL1*Y1**3-329.56589549458505d0*DL1&
            &**2*Y1**3+79.42559996133505d0*DL1**3*Y1**3-13760.d0&
            &*zeta2-6240.d0*DL*zeta2-1140.d0*DL**2*zeta2&
            &-22.22222222222222d0*DL**3*zeta2-45.d0*DL**4*zeta2&
            &-16270.d0*zeta3-5453.333333333333d0*DL*zeta3-1300.d0*DL&
            &**2*zeta3-280.d0*DL**3*zeta3+4760.d0*zeta2*zeta3+2480.d0&
            &*DL*zeta2*zeta3+1930.d0*zeta3**2-2753.333333333333d0&
            &*zeta4-550.d0*DL*zeta4-535.d0*DL**2*zeta4&
            &-2986.6666666666665d0*zeta5-3480.d0*DL*zeta5+362.5d0&
            &*zeta6)
       
       RETURN
     END FUNCTION P3NSSA_2604_FIT
!
! =================================================================av==
   end module xpns3s_2604_fit
