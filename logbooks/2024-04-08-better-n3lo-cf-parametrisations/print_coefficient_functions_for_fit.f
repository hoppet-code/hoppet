      program xc2ns3e_vs_xc2ns3p
      use xc2ns3e
      use xc2ns3p
      use xc3ns3e
      use xc3ns3p
      use xclns3e
      use xclns3p
      use n3locns
      use c2nsreg
      IMPLICIT REAL*8 (A - Z)
      COMPLEX*16 HC1, HC2, HC3, HC4, HC5 
      INTEGER NF2, N1, N2, NW
      PARAMETER ( N1 = -1, N2 = 1, NW = 5 ) 
      DIMENSION HC1(N1:N2),HC2(N1:N2,N1:N2),HC3(N1:N2,N1:N2,N1:N2), 
     ,     HC4(N1:N2,N1:N2,N1:N2,N1:N2), 
     ,     HC5(N1:N2,N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HR1(N1:N2),HR2(N1:N2,N1:N2),HR3(N1:N2,N1:N2,N1:N2), 
     ,     HR4(N1:N2,N1:N2,N1:N2,N1:N2), 
     ,     HR5(N1:N2,N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HRP1(N1:N2),HRP2(N1:N2,N1:N2),HRP3(N1:N2,N1:N2,N1:N2), 
     ,     HRP4(N1:N2,N1:N2,N1:N2,N1:N2), 
     ,     HRP5(N1:N2,N1:N2,N1:N2,N1:N2,N1:N2) 
      DIMENSION HI1(N1:N2),HI2(N1:N2,N1:N2),HI3(N1:N2,N1:N2,N1:N2), 
     ,     HI4(N1:N2,N1:N2,N1:N2,N1:N2), 
     ,     HI5(N1:N2,N1:N2,N1:N2,N1:N2,N1:N2)
      
!     implicit none
      integer :: nf, cc, ix
      character(100) :: x_char, nf_char, filename
      real * 8 :: logx, delx, x, dl, param, exact, large_x, param_new
      integer, parameter :: nbins=200
      real * 8, parameter :: logxmin = log(1d-12)!-30.0d0
      real * 8, parameter ::  logxmax=0.0d0
      real * 8, parameter :: logxomxmin = log(1d-7) ! smallest 1 - x value
      real * 8, parameter ::  logxomxmax= - log(1d-6) ! smallest x value
      integer i1,i2,i3,i4,i5

! Set reasonable default values. The exact and parametrised versions of
! the F2 and Fl coefficient functions differ at large x.

      cc = 1                   ! Do NC coefficient function
!      nf = 5

!      CALL SET_C2SOFT_N3LO(nf)  ! Needed to set some common pieces

!     C2 pieces ------------------------------------------------------
      logx = 0.0d0
      filename = 'coefficient_functions_for_fit.dat'
      delx = (logxomxmax - logxomxmin) / dble(nbins)
      open(unit = 99, file = trim(filename))
      write(99,*) '# x, param, large_x, exact'
      write(99,*) '# c2np (nf independent piece)'
      do ix = 1, nbins
         logx = logxomxmin + (ix - 0.5d0) * delx
         x = exp(-logx)/(1d0+ exp(-logx))
         dl = log(x)
         
         large_x = c2np3a_large_x(x,dl,0)
         exact = x2np3a(x,0,cc)
         param = c2np3a(x,dl,0,cc)
         
         write(99,*) x, param, large_x, exact
      enddo
      write(99,*) ''
      write(99,*) ''

      write(99,*) '# c2np (nf piece)'
      do ix = 1, nbins
         logx = logxomxmin + (ix - 0.5d0) * delx
         x = exp(-logx)/(1d0+ exp(-logx))
         dl = log(x)

         large_x = (c2np3a_large_x(x,dl,1) - c2np3a_large_x(x,dl,-1))
     $        * 0.5d0
         exact = (x2np3a(x,1,cc) - x2np3a(x,-1,cc)) * 0.5d0
         param = (c2np3a(x,dl,1,cc) - c2np3a(x,dl,-1,cc)) * 0.5d0
         
       write(99,*) x, param, large_x, exact
      enddo
      write(99,*) ''
      write(99,*) ''

      write(99,*) '# c2np (nf**2 piece)'
      do ix = 1, nbins
         logx = logxomxmin + (ix - 0.5d0) * delx
         x = exp(-logx)/(1d0+ exp(-logx))
         dl = log(x)

         large_x = (c2np3a_large_x(x,dl,1) + c2np3a_large_x(x,dl,-1) -
     $        2d0 * c2np3a_large_x(x,dl,0))* 0.5d0
         exact = (x2np3a(x,1,cc) + x2np3a(x,-1,cc) - 2d0*x2np3a(x,0,cc))
     $        * 0.5d0
         param = (c2np3a(x,dl,1,cc) + c2np3a(x,dl,-1,cc) - 2d0*c2np3a(x
     $        ,dl,0,cc))* 0.5d0
         
       write(99,*) x, param, large_x, exact
      enddo
      write(99,*) ''
      write(99,*) ''

!     END C2 pieces --------------------------------------------------

      
!     C3 pieces ------------------------------------------------------
      logx = 0.0d0
      cc = 0 ! Here cc = 0 means minus rather than valence which is what we want.
      filename = 'coefficient_functions_for_fit.dat'
      delx = (logxomxmax - logxomxmin) / dble(nbins)
      open(unit = 99, file = trim(filename))
!      write(99,*) '# x, param, large_x, exact'
      write(99,*) '# c3nm (nf independent piece)'
      do ix = 1, nbins
         logx = logxomxmin + (ix - 0.5d0) * delx
         x = exp(-logx)/(1d0+ exp(-logx))
         dl = log(x)
         
         large_x = c3nm3a_large_x(x,dl,0)
         exact = x3nm3a(x,dl,0,cc)
         param = c3nm3a(x,dl,0,cc)
         
         write(99,*) x, param, large_x, exact
      enddo
      write(99,*) ''
      write(99,*) ''

      write(99,*) '# c3nm (nf piece)'
      do ix = 1, nbins
         logx = logxomxmin + (ix - 0.5d0) * delx
         x = exp(-logx)/(1d0+ exp(-logx))
         dl = log(x)

         large_x = (c3nm3a_large_x(x,dl,1) - c3nm3a_large_x(x,dl,-1))
     $        * 0.5d0
         exact = (x3nm3a(x,dl,1,cc) - x3nm3a(x,dl,-1,cc)) * 0.5d0
         param = (c3nm3a(x,dl,1,cc) - c3nm3a(x,dl,-1,cc)) * 0.5d0
         
       write(99,*) x, param, large_x, exact
      enddo
      write(99,*) ''
      write(99,*) ''

      write(99,*) '# c3nm (nf**2 piece)'
      do ix = 1, nbins
         logx = logxomxmin + (ix - 0.5d0) * delx
         x = exp(-logx)/(1d0+ exp(-logx))
         dl = log(x)

         large_x = (c3nm3a_large_x(x,dl,1) + c3nm3a_large_x(x,dl,-1) -
     $        2d0 * c3nm3a_large_x(x,dl,0))* 0.5d0
         exact = (x3nm3a(x,dl,1,cc) + x3nm3a(x,dl,-1,cc) - 2d0*x3nm3a(x
     $        ,dl,0,cc))* 0.5d0
         param = (c3nm3a(x,dl,1,cc) + c3nm3a(x,dl,-1,cc) - 2d0*c3nm3a(x
     $        ,dl,0,cc))* 0.5d0
         
       write(99,*) x, param, large_x, exact
      enddo
      write(99,*) ''
      write(99,*) ''

!     END C3 pieces --------------------------------------------------

!     CL pieces ------------------------------------------------------
      logx = 0.0d0
      cc = 1 ! Charged current piece
      filename = 'coefficient_functions_for_fit.dat'
      delx = (logxomxmax - logxomxmin) / dble(nbins)
      open(unit = 99, file = trim(filename))
      write(99,*) '# x, param, large_x, exact'
      write(99,*) '# clnp (nf independent piece)'
      do ix = 1, nbins
         logx = logxomxmin + (ix - 0.5d0) * delx
         x = exp(-logx)/(1d0+ exp(-logx))
         dl = log(x)
         
         large_x = clnp3a_large_x(x,dl,0)
         exact = xlnp3a(x,0,cc)
         param = clnp3a(x,dl,0,cc)
         
         write(99,*) x, param, large_x, exact
      enddo
      write(99,*) ''
      write(99,*) ''

      write(99,*) '# clnp (nf piece)'
      do ix = 1, nbins
         logx = logxomxmin + (ix - 0.5d0) * delx
         x = exp(-logx)/(1d0+ exp(-logx))
         dl = log(x)

         large_x = (clnp3a_large_x(x,dl,1) - clnp3a_large_x(x,dl,-1))
     $        * 0.5d0
         exact = (xlnp3a(x,1,cc) - xlnp3a(x,-1,cc)) * 0.5d0
         param = (clnp3a(x,dl,1,cc) - clnp3a(x,dl,-1,cc)) * 0.5d0
         
       write(99,*) x, param, large_x, exact
      enddo
      write(99,*) ''
      write(99,*) ''

      write(99,*) '# clnp (nf**2 piece)'
      do ix = 1, nbins
         logx = logxomxmin + (ix - 0.5d0) * delx
         x = exp(-logx)/(1d0+ exp(-logx))
         dl = log(x)

         large_x = (clnp3a_large_x(x,dl,1) + clnp3a_large_x(x,dl,-1) -
     $        2d0 * clnp3a_large_x(x,dl,0))* 0.5d0
         exact = (xlnp3a(x,1,cc) + xlnp3a(x,-1,cc) - 2d0*xlnp3a(x,0,cc))
     $        * 0.5d0
         param = (clnp3a(x,dl,1,cc) + clnp3a(x,dl,-1,cc) - 2d0*clnp3a(x
     $        ,dl,0,cc))* 0.5d0
         
       write(99,*) x, param, large_x, exact
      enddo
      write(99,*) ''
      write(99,*) ''

!     END C2 pieces --------------------------------------------------


      end

      
