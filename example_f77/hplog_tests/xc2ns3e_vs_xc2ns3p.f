      program xc2ns3e_vs_xc2ns3p
      use xc2ns3e
      use xc2ns3p
      implicit none
      integer :: nf, cc, ix
      character(100) :: x_char, nf_char, filename
      real * 8 :: logx, delx, x, dl, param, exact
      integer, parameter :: nbins=200
      real * 8, parameter :: logxmin = -20.0d0
      real * 8, parameter ::  logxmax=0.0d0

! Set reasonable default values. The exact and parametrised versions of
! the F2 and Fl coefficient functions differ at large x.

      cc = 1                    ! Do NC coefficient function
      delx = (logxmax - logxmin) / dble(nbins)
      nf = 5

      CALL SET_C2SOFT_N3LO(nf)  ! Needed to set some common pieces

      logx = 0.0d0
      filename = 'xc2ns3e_vs_xc2ns3p.dat'
      open(unit = 99, file = trim(filename))
      write(99,*)
     $     '# 1-x,                          param                      e
     $xact,                   param/exact'
      do ix = 1, nbins
         logx = logxmin + (ix - 0.5d0) * delx
         x = 1.0d0 - exp(logx)
         dl = log(x)

         param = c2np3a(x,dl,nf,cc)
         exact = x2np3a(x,nf,cc)
         
       write(99,*) 1-x, param, exact, param/exact
      enddo
      
      end

      
