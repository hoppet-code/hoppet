      program xc2ns3e_vs_xc2ns3p
      use xc2ns3e
      use xc2ns3p
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
      real * 8 :: logx, delx, x, dl, param, exact, param_old
      integer, parameter :: nbins=200
      real * 8, parameter :: logxmin = log(1d-12)!-30.0d0
      real * 8, parameter ::  logxmax=0.0d0
      real * 8, parameter :: logxomxmin = log(1d-4) ! smallest 1 - x value
      real * 8, parameter ::  logxomxmax= - log(1d-6) ! smallest x value
      integer i1,i2,i3,i4,i5

      x = 0.9d0

      CALL HPLOG5 (X, NW, HC1,HC2,HC3,HC4,HC5, HR1,HR2,HR3,HR4,HR5,
     ,     HI1,HI2,HI3,HI4,HI5, N1, N2)

! Set reasonable default values. The exact and parametrised versions of
! the F2 and Fl coefficient functions differ at large x.

      cc = 1                   ! Do NC coefficient function
      nf = 5

      CALL SET_C2SOFT_N3LO(nf)  ! Needed to set some common pieces

      logx = 0.0d0
      filename = 'coefficient_functions_for_fit.dat'
      delx = (logxomxmax - logxomxmin) / dble(nbins)
      open(unit = 99, file = trim(filename))
      write(99,*)
     $     '# x,                          param                      exa
     $ct,                   exact - param'
      write(99,*) '# c2np (nf independent piece)'
      do ix = 1, nbins
         logx = logxomxmin + (ix - 0.5d0) * delx
!x = 1.0d0 - exp(logx)
         x = exp(-logx)/(1d0+ exp(-logx))
         dl = log(x)
         
         param = c2np3a_nopoly(x,dl,0,cc)
         exact = x2np3a(x,0,cc)
         
         write(99,*) x, param, exact, exact - param
      enddo
      write(99,*) ''
      write(99,*) ''

      write(99,*) '# c2np (nf piece)'
      do ix = 1, nbins
         logx = logxomxmin + (ix - 0.5d0) * delx
         !x = 1.0d0 - exp(logx)
         x = exp(-logx)/(1d0+ exp(-logx))
         dl = log(x)

         param = (c2np3a_nopoly(x,dl,1,cc) - c2np3a_nopoly(x,dl,-1,cc))
     $        * 0.5d0
         exact = (x2np3a(x,1,cc) - x2np3a(x,-1,cc)) * 0.5d0
         
       write(99,*) x, param, exact, exact - param
      enddo
      write(99,*) ''
      write(99,*) ''

      write(99,*) '# c2np (nf**2 piece)'
      do ix = 1, nbins
         logx = logxomxmin + (ix - 0.5d0) * delx
         !x = 1.0d0 - exp(logx)
         x = exp(-logx)/(1d0+ exp(-logx))
         dl = log(x)

         param = (c2np3a_nopoly(x,dl,1,cc) + c2np3a_nopoly(x,dl,-1,cc) -
     $        2d0 * c2np3a_nopoly(x,dl,0,cc))* 0.5d0
         exact = (x2np3a(x,1,cc) + x2np3a(x,-1,cc) - 2d0*x2np3a(x,0,cc))
     $        * 0.5d0
         
       write(99,*) x, param, exact, exact - param
      enddo
      write(99,*) ''
      write(99,*) ''

      ! New files with fit below
      
      logx = 0.0d0
      filename = 'coefficient_functions_new.dat'
      delx = (logxomxmax - logxomxmin) / dble(nbins)
      open(unit = 99, file = trim(filename))
      write(99,*)
     $     '# x,                          param                      ex
     $act,                   param/exact                    param_old/ex
     $act'
      do ix = 1, nbins
         logx = logxomxmin + (ix - 0.5d0) * delx
         !x = 1.0d0 - exp(logx)
         x = exp(-logx)/(1d0+ exp(-logx))
         dl = log(x)

         param = c2np3a_new(x,dl,nf,cc)
         param_old = c2np3a(x,dl,nf,cc)
         exact = x2np3a(x,nf,cc)
         
       write(99,*) x, param, exact, param/exact, param_old/exact 
      enddo
      write(99,*) ''
      write(99,*) ''
      
      end

      
