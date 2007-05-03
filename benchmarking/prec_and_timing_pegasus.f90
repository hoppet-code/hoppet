! Program for testing pegasus the same way as has been done with hoppet
! 
! ./prec_and_timing_pegasus -nrep 1 -nxQ 5000 -outputgrid [-ifast {0|1}]
!
! For this to compile you will need to manually add the set of pegasus
! files in the link part of the Makefile, (.../Pegasus/Objects2/*.o,
! once you've removed the "old" files from there; or even
! /Pegasus/Objects2/*[^old].o)
!
!======================================================================
module pegasus_communicator
  integer :: pegasus_ifast
end module pegasus_communicator

program prec_and_timing_pegasus
  use sub_defs_io
  use pegasus_communicator
  use types; use consts_dp
  implicit none
  !-------------------------------
  integer            :: order, order1, order2, nloop, i, nrep, nxQ,olnlnQ
  integer            :: hires
  real(dp)           :: dy, Qinit, Qmax, du, dlnlnQ
  real(dp)           :: ymax
  real(dp)           :: time_start, time_init_done, time_ev_done, time_end
  real(dp), pointer  :: vogt_init(:,:)
  logical            :: output, outputgrid, preev


  ymax  = dble_val_opt('-ymax',11.5_dp)
  Qinit = sqrt(two); Qmax = dble_val_opt('-Qmax',1e4_dp)
  pegasus_ifast = int_val_opt('-ifast',0)


  call cpu_time(time_start)
  CALL INITEVOL (2) ! start pegasus using the subroutine below
  call cpu_time(time_init_done)


  ! decide 
  nrep  = int_val_opt('-nrep',1)
  nxQ = int_val_opt('-nxQ',0)
  output = log_val_opt('-output') .or. log_val_opt('-outputgrid')
  outputgrid  = log_val_opt('-outputgrid')
  if (output) call output_info

  !-- security ----------------------
  if (.not. CheckAllArgsUsed(0)) stop
  !----------------------------------

  call cpu_time(time_ev_done)

  ! evolution & output
  do i = 1, nrep
     ! pegasus input initialisation from the subroutine below
     CALL INITINP (2)
     !CALL INITINP (1) ! init from file


     ! one form of output
     if (outputgrid) then
        call eval_output_grid()
     else
        call eval_output_lines()
     end if

  end do
  call cpu_time(time_end)
  write(0,'(a,4f10.5)') "Timings (init, preevln, evln) = ", &
       &   time_init_done-time_start, &
       &   time_ev_done-time_init_done, &
       &   (time_end-time_ev_done)/nrep
  write(6,'(a)',advance='no') "# "
  call time_stamp(6)
  call system("echo '# host:' `hostname`")
  ! record info about the cpu
  call system("grep -e name -e cache -e MHz /proc/cpuinfo | sed 's/^/# /'")
  if (output) write(6,'(a,4f10.5)') "# Timings (init, preevln, evln) = ", &
       &   time_init_done-time_start, &
       &   time_ev_done-time_init_done, &
       &   (time_end-time_ev_done)/nrep


contains

  !---------------------------------------------------------------
  !! some basic commentary
  subroutine output_info
    write(6,'(a)') '# '//trim(command_line())
    write(6,'(a)') '# Running Pegasus:'
    write(6,'(a,i5)')      '# ifast = ', pegasus_ifast
    write(6,'(a,4i10)')    '# dummy = '
    write(6,'(a,4f10.6)')  '# ymax  = '
    write(6,'(a,4i10)'  )  '# dummy = '
    write(6,'(a,4es10.2)') '# dummy = '

    write(6,'(a,f10.5,a,f10.5)') '# dummy = '
    write(6,'(a,f10.5,a,i5)')    '# dummy = '
  end subroutine output_info


  !-------------------------------------------------------------------
  !! output lines in the y,Q plane
  subroutine eval_output_lines()
    integer nn, j
    real(dp) :: y, Q, pdfval(-6:6),as
    nn = nxQ/4
    do j = 1, nn
       y = j*ymax/nn
       Q = Qmax - j*(Qmax-Qinit)/nn
       call PEGXPARTON(pdfval,as,exp(-y),Q**2)
       if (output .and. i==1) write(6,'(20es20.10)') y,Q,pdfval(0:4)
       !if (output .and. i==1) write(6,'(20es20.8)') y,Q,vogt_init(:,0:3).atx.(grid.with.exp(-y))
    end do

    if (output .and. i==1) write(6,*)
    if (output .and. i==1) write(6,*)
    do j = nn,1,-1
       y = j*ymax/nn
       Q = 4.0_dp + j*5.0_dp/nn
       call PEGXPARTON(pdfval,as,exp(-y),Q**2)
       if (output .and. i==1) write(6,'(20es20.10)') y,Q,pdfval(0:4)
    end do

    if (output .and. i==1) write(6,*)
    if (output .and. i==1) write(6,*)
    do j = nn,1,-1
       y = j*ymax/nn
       Q = Qmax*(1-j*0.2_dp/nn)
       call PEGXPARTON(pdfval,as,exp(-y),Q**2)
       if (output .and. i==1) write(6,'(20es20.10)') y,Q,pdfval(0:4)
    end do

    if (output .and. i==1) write(6,*)
    if (output .and. i==1) write(6,*)
    do j = nn,1,-1
       y = j*ymax/nn
       Q = sqrt(Qinit*Qmax)*(1+j*0.2_dp/nn)
       call PEGXPARTON(pdfval,as,exp(-y),Q**2)
       if (output .and. i==1) write(6,'(20es20.10)') y,Q,pdfval(0:4)
    end do
  end subroutine eval_output_lines


  !----------------------------------------------------------------
  !! output the results on a grid uniform in 
  !! zeta = ln 1/x + grid_a*(1-x)
  subroutine eval_output_grid()
    integer  :: nz, nQ, iz, iQ
    real(dp) :: zmax, zeta, y, Q, zQ, zQmax, pdfval(-6:6), as
    real(dp), parameter :: grid_a = 9, gridQ_a = 3

    nz = nint(sqrt(four*nxQ))
    nQ = nint(sqrt(0.25_dp*nxQ))-1
    !nz=100
    if (output .and. i == 1) write(0,*) 'nz, nQ = ', nz, nQ

    zmax = zeta_of_y(ymax, grid_a)
    zQmax = zeta_of_y(log(Qmax/Qinit), gridQ_a)
    do iQ = 0, nQ
       do iz = 1, nz
          zeta = iz * zmax/nz
          y    = y_of_zeta(zeta, grid_a)
          !zQ = (iQ) * zQmax / nQ
          zQ = (iQ+zeta/zmax) * zQmax / nQ
          Q = max(Qinit,min(Qmax,Qinit * exp(y_of_zeta(zQ, gridQ_a))))
          call PEGXPARTON(pdfval,as,exp(-y),Q**2)
          if (output .and. i == 1) write(6,'(20es20.10)') y,Q,pdfval(-5:5)
       end do
       if (output .and. i == 1) write(6,'(a)') 
    end do

  end subroutine eval_output_grid



  !-----------------------------------------------------------------
  !! return zeta = ln 1/x + a*(1-x)  (x = exp(-y))
  function zeta_of_y(y, a) result(zeta)
    real(dp), intent(in) :: y, a
    real(dp)             :: zeta
    zeta = y + a*(one - exp(-y))
  end function zeta_of_y


  !-----------------------------------------------------------------
  !! return inversion of zeta = ln 1/x + a*(1-x)  (x = exp(-y))
  function y_of_zeta(zeta, a) result(y)
    real(dp), intent(in) :: zeta, a
    real(dp)             :: y, x, diff_from_zero, deriv
    integer             :: iter
    real(dp), parameter :: eps = 1e-12_dp
    integer,  parameter :: maxiter = 100

    ! starting condition (and soln if a = 0)
    y = zeta 
    if (a /= zero) then
       do iter = 0, maxiter
          x = exp(-y);
          diff_from_zero = zeta - y - a*(one-x);
          ! we have found good solution
          if (abs(diff_from_zero) < eps) exit
          deriv = -one  - a*x;
          y = y - diff_from_zero / deriv;
       end do
    end if

    if (iter > maxiter) write(0,*) "y_of_zeta reached maxiter"

  end function y_of_zeta


  subroutine PEGXPARTON(pdfval,as,x,Q2)
    real(dp), intent(out) :: pdfval(-6:6),as
    real(dp), intent(in)  :: x, Q2
    call XPARTON(pdfval,as,x,Q2,-5,5,1)
    call swap(pdfval(1),pdfval(2))
    call swap(pdfval(-1),pdfval(-2))
  end subroutine PEGXPARTON
  
  subroutine swap(a,b)
    real(dp) :: a, b
    real(dp) :: tmp
    tmp = a
    a = b
    b = tmp
  end subroutine swap
  

end program prec_and_timing_pegasus


!
! ..File: usrinit.f
!
!
! ..Initialization parameters for speed, the heavy-flavour treatment,
!    the mode and order of the evolution, and the value of  FR2 = 
!    mu_f^2/mu_r^2  used by INITEVOL (IPAR) if called with  IPAR = 2.
!    See section 4.1 of the manual.
!
! =====================================================================
!
!
SUBROUTINE USRINIT (IFAST, IVFNS, NFF, IMODEV, NPORD, FR2)
  use pegasus_communicator
  ! 
  IMPLICIT DOUBLE PRECISION (A - Z)
  INTEGER IFAST, IVFNS, NFF, IMODEV, NPORD
  !
  IFAST  = pegasus_ifast
  !IFAST  = 1
  IVFNS  = 1
  NFF    = 4
  IMODEV = 1
  NPORD  = 2
  FR2    = 1.D0
  !
  RETURN
END SUBROUTINE USRINIT
!
! =================================================================av==

! ..File: usrinp.f  
!
!
! ..Input values for initial scale and alpha_s, the heavy-quark masses 
!    and the parton-distribution shapes used by  INITINP (IPAR)  if 
!    called with  IPAR = 2.  See section 4.2 of the manual.
!                   
! =====================================================================
!
!
SUBROUTINE USRINP (PUV, PDV, PLM, PLP, PSM, PSP, PGL, M20, ASI,&
     MC2, MB2, MT2, NFORM, IMOMIN, ISSIMP)
  ! 
  IMPLICIT DOUBLE PRECISION (A - Z)
  INTEGER IMOMIN, ISSIMP, NFORM 
  DIMENSION PUV(6), PDV(6), PLM(6), PLP(6), PSM(6), PSP(6), PGL(6)
  !
  M20 = 2D0
  ASI = 0.35000
  !
  ! ..The otherwise irrelevant shift of m_c ensures that the input is
  !   returned for m_c = mu_0 also in the VFNS. Useful for input checks.
  !
  MC2 = 2.0000000017732D0  ! use identical value to what hoppet has
  MB2 = 20.25D0
  MT2 = 3.0625D4
  !
  NFORM  = 1
  IMOMIN = 1
  ISSIMP = 0
  !
  PUV(1) =  2.0D0
  PUV(2) =  0.8D0
  PUV(3) =  3.0D0
  PUV(4) =  0.5D0
  PUV(5) =  0.0D0
  PUV(6) =  0.0D0
  !
  PDV(1) =  1.0D0
  PDV(2) =  0.8D0
  PDV(3) =  4.0D0
  PDV(4) =  0.5D0
  PDV(5) =  0.0D0
  PDV(6) =  0.0D0
  !
  PLM(1) =  0.193987D0
  PLM(2) =  0.9D0
  PLM(3) =  6.0D0
  PLM(4) =  0.5D0
  PLM(5) =  0.0D0
  PLM(6) =  0.0D0
  !           
  PLP(1) =  0.136565D0
  PLP(2) = -0.1D0
  PLP(3) =  6.0D0
  PLP(4) =  0.5D0
  PLP(5) =  0.0D0
  PLP(6) = -0.5D0
  !          
  PSM(1) =  0.0d0
  PSM(2) =  0.9D0
  PSM(3) =  6.0D0
  PSM(4) =  0.5D0
  PSM(5) =  0.0D0
  PSM(6) =  0.0D0
  !          
  PSP(1) =  0.027313d0
  PSP(2) = -0.1D0
  PSP(3) =  6.0D0
  PSP(4) =  0.5D0
  PSP(5) =  0.0D0
  PSP(6) = -0.5D0
  !
  PGL(1) =  1.0D0
  PGL(2) = -0.1D0
  PGL(3) =  5.0D0
  PGL(4) =  0.5D0
  PGL(5) =  0.0D0
  PGL(6) =  0.0D0
  !
  RETURN
END SUBROUTINE USRINP
!
! =================================================================av==

SUBROUTINE USRPINP
end SUBROUTINE USRPINP

SUBROUTINE USRPINIT
end SUBROUTINE USRPINIT
