program compare2file
  use sub_defs_io
  use types; use consts_dp
  implicit none
  character(len=300) :: line1, line2
  integer :: idev1, idev2, iostat
  
  real(dp) :: max07_guds, max07_gudsc, err_guds, err_gudsc
  real(dp) :: max09_guds, max09_gudsc
  real(dp) :: res1(0:4), res2(0:4), err(0:4)
  real(dp) :: y,Q

  idev1 = idev_open_arg(1,status='old')
  idev2 = idev_open_arg(2,status='old')

  do
     read(idev1,'(a)',iostat=iostat) line1
     if (iostat /= 0) exit
     read(idev2,'(a)',iostat=iostat) line2
     read(line1,*,iostat=iostat) y,Q,res1
     if (iostat /= 0) cycle
     read(line2,*) y,Q,res2

     err(0:4) = abs(res1/res2-one)
     err_gudsc = maxval(err, mask=(res2 /= zero))
     err_guds  = maxval(err(0:3), mask=(res2(0:3)/=0))
     if (exp(-y) <= 0.9_dp) then
        max09_guds  = max(max09_guds,  err_guds)
        max09_gudsc = max(max09_gudsc, err_gudsc)
        !write(0,*) y,Q, max09_gudsc
     end if
     if (exp(-y) <= 0.7_dp) then
        max07_guds  = max(max07_guds,  err_guds)
        max07_gudsc = max(max07_gudsc, err_gudsc)
     end if
  end do

  write(6,"(4es20.10)") max07_guds, max07_gudsc, max09_guds, max09_gudsc
end program compare2file
