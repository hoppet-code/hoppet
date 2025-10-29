!******************************************************************************
!* This file is part of libome                                                *
!* Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
!* SPDX-License-Identifier: GPL-3.0-or-later                                  *
!******************************************************************************

program FortranExample
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    function ome_AqqQNSEven_reg(as, LM, NF, x) bind(C, name="ome_AqqQNSEven_reg") result(res)
      import c_double
      real(c_double), value, intent(in) :: as, LM, NF, x
      real(c_double) :: res
    end function
  end interface

  real(kind=kind(1d0)) :: pi, as, m2, mu2, NF, LM, x

  ! Pi
  pi = 3.1415926535897932d0
  ! strong coupling, as = \alpha_s/(4 \pi)
  as = 0.118d0/(4d0*pi)
  ! squared quark mass (on-shell scheme)
  m2 = 1.59d0**2
  ! squared renormalisation scale
  mu2 = 100.d0**2
  ! number of massless quark flavours
  NF = 3d0
  ! precompute the mass logarithm
  LM = log(m2/mu2)
  ! Bjorken x
  x = 0.2

  write (*,*) ome_AqqQNSEven_reg(as, LM, NF, x)
end program
