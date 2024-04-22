      program xc2ns3e_vs_xc2ns3p
      use xc2ns3e
      use xc2ns3p
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
      real * 8 :: logx, delx, x, dl, param, exact
      integer, parameter :: nbins=200
      real * 8, parameter :: logxmin = log(1d-12)!-30.0d0
      real * 8, parameter ::  logxmax=0.0d0
      integer i1,i2,i3,i4,i5

      x = 0.9d0

      CALL HPLOG5 (X, NW, HC1,HC2,HC3,HC4,HC5, HR1,HR2,HR3,HR4,HR5,
     ,     HI1,HI2,HI3,HI4,HI5, N1, N2)
      CALL HPOLY5 (X, NW, HRP1,HRP2,HRP3,HRP4,HRP5, N1, N2)

      do i1 = -1,1
         write(*,'(a3,i2,a1,E17.7,E17.7,F17.7)') 'H1(',i1,')', Hr1(i1),
     $        Hrp1(i1), Hrp1(i1)/Hr1(i1)
      enddo
      print*, ''
      do i1 = -1,1
         do i2 = -1,1
            write(*,'(a3,i2,i3,a1,E17.7,E17.7,F17.7)') 'H2(',i1,i2,')',
     $           Hr2(i1,i2),Hrp2(i1,i2), Hrp2(i1,i2)/Hr2(i1,i2)
         enddo
      enddo
      print*, ''
      do i1 = -1,1
         do i2 = -1,1
            do i3 = -1,1
               write(*,'(a3,i2,i3,i3,a1,E17.7,E17.7,F17.7)') 'H3(',i1,i2
     $              ,i3,')',Hr3(i1,i2,i3),Hrp3(i1,i2,i3), Hrp3(i1,i2,i3)
     $              /Hr3(i1 ,i2,i3)
            enddo
         enddo
      enddo
      print*, ''
      do i1 = -1,1
         do i2 = -1,1
            do i3 = -1,1
               do i4 = -1,1
                  write(*,'(a3,i2,i3,i3,i3,a1,E17.7,E17.7,F17.7)') 'H4('
     $                 ,i1,i2,i3,i4,')',Hr4(i1,i2,i3,i4),Hrp4(i1,i2,i3
     $                 ,i4), Hrp4(i1,i2,i3,i4)/Hr4(i1 ,i2,i3,i4)
               enddo
            enddo
         enddo
      enddo
      print*, ''
      do i1 = -1,1
         do i2 = -1,1
            do i3 = -1,1
               do i4 = -1,1
                  do i5 = -1,1
                     write(*,'(a3,i2,i3,i3,i3,i3,a1,E17.7,E17.7,F17.7)')
     $                    'H5(',i1,i2,i3,i4,i5,')',Hr5(i1,i2,i3,i4,i5)
     $                    ,Hrp5(i1,i2,i3,i4,i5), Hrp5(i1,i2,i3,i4,i5)
     $                    /Hr5(i1 ,i2,i3,i4,i5)
                  enddo
               enddo
            enddo
         enddo
      enddo
!      stop

! Set reasonable default values. The exact and parametrised versions of
! the F2 and Fl coefficient functions differ at large x.

      cc = 1                   ! Do NC coefficient function
      delx = (logxmax - logxmin) / dble(nbins)
      nf = 5

      CALL SET_C2SOFT_N3LO(nf)  ! Needed to set some common pieces

      logx = 0.0d0
      filename = 'xc2ns3e_vs_xc2ns3p_hplog.dat'
      open(unit = 99, file = trim(filename))
      write(99,*)
     $     '# 1-x,                          param                      e
     $xact,                   param/exact'
      do ix = 1, nbins
         logx = logxmin + (ix - 0.5d0) * delx
         x = 1.0d0 - exp(logx)
         dl = log(x)

         param = c2np3a(x,dl,nf,1) + c2np3a(x,dl,nf,0) 
         exact = x2np3a(x,nf,1) + x2np3a(x,nf,0)
         
       write(99,*) 1-x, param, exact, param/exact
      enddo
      
      end

      
