! program for comparing format produced by "small_fast_tab -outputgrid" 
!
! Usage:
!   progname file1 file2 {-summary | -channel ic [-minerr val] [-maxerr val]}
program compare2file_v2
  use sub_defs_io
  use types; use consts_dp
  implicit none
  character(len=300) :: line1, line2
  integer :: idev1, idev2, iostat
  integer, parameter :: maxn = 200

  real(dp) :: max07_guds, max07_all, err_guds, err_all
  real(dp) :: max09_guds, max09_all
  real(dp) :: res1(maxn,maxn,-5:5) = zero
  real(dp) :: res2(maxn,maxn,-5:5) = zero
  real(dp) :: err (maxn,maxn,-5:5)
  real(dp) :: yval(maxn), Qval(maxn,maxn), xval(maxn)
  logical  :: mask  (maxn,maxn,-5:5)=.true.
  logical  :: mask07(maxn,maxn,-5:5)=.false.
  logical  :: mask09(maxn,maxn,-5:5)=.false.
  real(dp) :: y,Q, minerr, maxerr, this_err, true_err
  integer  :: iy, iQ, ny, nQ, ic, idy, ml(3)

  idev1 = idev_open_arg(1,status='old')
  idev2 = idev_open_arg(2,status='old')

  iy = 1; iQ = 1
  do
     read(idev1,'(a)',iostat=iostat) line1
     if (iostat /= 0) exit
     read(idev2,'(a)',iostat=iostat) line2
     if (iostat /= 0) exit

     ! keep track of indicies
     if (trim(line1) == "") then
        ny = iy - 1
        iy = 1
        iQ = iQ + 1
     end if

     read(line1,*,iostat=iostat) yval(iy),Qval(iy,iQ),res1(iy,iQ,:)
     if (iostat /= 0) cycle
     read(line2,*) yval(iy),Qval(iy,iQ),res2(iy,iQ,:)
     iy = iy + 1
  end do
  nQ = iQ - 1

  ! establish a general protection mask -- if there is a sign change
  ! in the neighbourhood of the point then it should be eliminated.
  ! Work out dy neighbourhood
  idy = nint(0.4_dp / (yval(ny)-yval(ny-1)))
  if (log_val_opt('-protect')) then
     forall(iQ=1:nQ)
        forall(iy=1:ny)
           !mask(iy,iQ,:) = (abs(res2(iy,iQ,:)) &
           !     &               < 3.0_dp*min(abs(res2(max(iy-1,1), iQ,:)), &
           !     &                            abs(res2(min(iy+1,ny),iQ,:)))&
           !     & .or.      abs(res2(iy,iQ,:)) &
           !     &               < 3.0_dp*min(abs(res2(iy,max(iQ-1, 1),:)), &
           !     &                            abs(res2(iy,min(iQ+1,nQ),:))))
           ! check for neighbouring sign change
           !mask(iy,iQ,:) = res2(max(iy-1,1), iQ,:)*res2(min(iy+1,ny),iQ,:) >= 0
           ! check for sign change in x direction
           mask(iy,iQ,:) = mask(iy,iQ,:) .and. &
                &   minval(res2(max(iy-idy,1):min(iy+idy,ny), iQ,:),dim=1) *&
                &   maxval(res2(max(iy-idy,1):min(iy+idy,ny), iQ,:),dim=1) >= 0
           ! check for sign change in Q direction
           mask(iy,iQ,:) = mask(iy,iQ,:) .and. &
                &   minval(res2(iy,max(iQ-1,1):min(iQ+1,nQ),:),dim=1) *&
                &   maxval(res2(iy,max(iQ-1,1):min(iQ+1,nQ),:),dim=1) >= 0
        end forall
     end forall
  end if

  ! now get masks for the x range
  xval(:ny) = exp(-yval(:ny))
  mask07(:ny,:nQ,:) = mask(:ny,:nQ,:) .and. spread(spread(xval<0.7_dp,2,nQ),3,11)
  mask09(:ny,:nQ,:) = mask(:ny,:nQ,:) .and. spread(spread(xval<0.9_dp,2,nQ),3,11)

  ! now establish errors
  where (res2 /= zero)
     err = abs(res1/res2 - one)
  elsewhere
     err = zero
  end where
  
  
  max07_guds = maxval(err(:ny,:nQ,0:3),mask07(:ny,:nQ,0:3))
  max09_guds = maxval(err(:ny,:nQ,0:3),mask09(:ny,:nQ,0:3))
  max07_all  = maxval(err(:ny,:nQ,-5:5),mask07(:ny,:nQ,-5:5))
  max09_all  = maxval(err(:ny,:nQ,-5:5),mask09(:ny,:nQ,-5:5))

  ml = maxloc(err(:ny,:nQ,0:3),mask09(:ny,:nQ,0:3));
  !write(0,*) yval(ml(1)), Qval(ml(1),ml(2)), ml(3)-1

  if (log_val_opt("-summary")) then
     write(6,"(4es20.10)") max07_guds, max07_all, max09_guds, max09_all
  end if
  
  if (log_val_opt("-channel")) then
     ic = int_val_opt("-channel")
     minerr = dble_val_opt("-minerr",zero)
     maxerr = dble_val_opt("-maxerr",1e100_dp)
     do iQ = 1, nQ
        do iy = 1, ny
           if (ic == 11) then
              true_err = maxval(err(iy,iQ,-5:5))
              if (any(mask(iy,iQ,:))) then
                 this_err = maxval(err(iy,iQ,-5:5),mask(iy,iQ,-5:5))
              else
                 this_err = -one
              end if
           else 
              true_err = err(iy,iQ,ic)
              if (.not.mask(iy,iQ,ic)) then
                 this_err = -one
              else
                 this_err = true_err
              end if
           end if
           if (this_err >= minerr .and. this_err <= maxerr) then
              write(6,'(4es20.12)') yval(iy), Qval(iy,iQ), this_err, true_err
           end if
        end do
        write(6,*)
     end do
  end if

end program compare2file_v2
