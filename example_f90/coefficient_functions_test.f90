!! An example program using structure functions up to N3LO
!!
!!

program coefficient_functions_test
  use hoppet_v1
  implicit none

  call print()
contains 
  !----------------------------------------------------------------------
  ! fill the streamlined interface PDF table (possibly using hoppet's
  ! evolution)
  subroutine print()
    use xclns3e
    use xclns3p
    use XCDIFF3P
    use XCLPS3E
    use XCLSG3P
!    use XCLSG3PV
    implicit none
    integer, parameter :: nf_lc = 5, nbins=500
    real(dp), parameter :: logxmin = -12_dp, logxmax=zero
    integer :: ix
    real(dp) :: logx, delx, x, y
    character (len=20) :: filename

    delx = (logxmax - logxmin) / dble(nbins)
    
    filename = 'xclns3p.dat'
    open(unit = 99, file = trim(filename))
    logx = zero
    do ix = 1, nbins
       logx = logxmin + (ix - 0.5_dp) * delx
       x = exp(logx)
       y = -logx
       write(99,*) logx, CLNP3A(x, -y, nf_lc, 1), &
                         CLNP3A(x, -y, nf_lc, 0), &
                         CLNP3C(x, nf_lc) !, CLQ3DF(x, -y, nf_int, 0)
    enddo
    close(unit = 99)

    filename = 'xclns3e.dat'
    open(unit = 99, file = trim(filename))
    logx = zero
    do ix = 1, nbins
       logx = logxmin + (ix - 0.5_dp) * delx
       x = exp(logx)
       y = -logx
       write(99,*) logx, XLNP3A(x,nf_lc) 
    enddo
    close(unit = 99)

    filename = 'xclsg3p.dat'
    open(unit = 99, file = trim(filename))
    logx = zero
    do ix = 1, nbins
       logx = logxmin + (ix - 0.5_dp) * delx
       x = exp(logx)
       y = -logx
       write(99,*) logx, CLS3A(x, -y, nf_lc, 1), &
                         CLS3A(x, -y, nf_lc, 0), &
                         CLG3A(x, -y, nf_lc, 1), CLG3A(x, -y, nf_lc, 0)
    enddo
    close(unit = 99)

!    filename = 'xclsg3pv.dat'
!    open(unit = 99, file = trim(filename))
!    logx = zero
!    do ix = 1, nbins
!       logx = logxmin + (ix - 0.5_dp) * delx
!       x = exp(logx)
!       y = -logx
!       write(99,*) logx, CLS3AV(x, -y, nf_lc, 1), &
!                         CLS3AV(x, -y, nf_lc, 0), &
!                         CLG3AV(x, -y, nf_lc, 1), CLG3AV(x, -y, nf_lc, 0)
!    enddo
!    close(unit = 99)

    filename = 'xclps3e.dat'
    open(unit = 99, file = trim(filename))
    logx = zero
    do ix = 1, nbins
       logx = logxmin + (ix - 0.5_dp) * delx
       x = exp(logx)
       y = -logx
       write(99,*) logx, XLS3A(x,nf_lc) 
    enddo
    close(unit = 99)


  end subroutine print

end program coefficient_functions_test


