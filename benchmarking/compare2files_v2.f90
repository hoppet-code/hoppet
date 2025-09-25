!! program for comparing format produced by "prec_and_timing [...] -outputgrid" 
!!
!! Usage:
!!   progname file1 file2 {-summary | -channel ic [-minerr val] [-maxerr val]} [-protect]
!!
!! Arguments:
!!
!! -summary: print a single-line summary of the relative differences between files
!!
!! -channel ic: the flavour channel to consider (ic = 11 for all flavours)
!!
!! -protect: apply a protection mask to eliminate points close to a sign change
!!           (hard-coded to be within 0.4 in y and within one grid point in Q)
!!
!! -minerr val: only print errors greater than val
!! -maxerr val: only print errors less than val
!!
program compare2file_v2
  use io_utils
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
  real(dp) :: y,Q, minerr, maxerr, this_err, true_err, Qmin
  integer  :: iy, iQ, ny, nQ, ic, icmax, idy, ml(3)
  real(dp) :: yval2, Qval2
  logical  :: do_help

  if (iargc() < 2) then
     do_help = .true.
  else
     do_help = trim(string_val_arg(1)) == "-h"
  end if
  if (do_help) then
       write(6,"(a)") "Usage: " // trim(command_line()) // " file1 file2 {-summary | -channel ic [-minerr val] [-maxerr val]} [-protect]"
       write(6,"(a)") ""
       write(6,"(a)") "Arguments:"
       write(6,"(a)") ""
       write(6,"(a)") "-summary: print a single-line summary of the relative differences between files"
       write(6,"(a)") ""
       write(6,"(a)") "-channel ic: print out an x,Q grid for the given flavour channel (ic = 11 for all flavours)"
       write(6,"(a)") ""
       write(6,"(a)") "-protect: apply a protection mask to eliminate points close to a sign change"
       write(6,"(a)") "          (hard-coded to be within 0.4 in y and within one grid point in Q)"
       write(6,"(a)") ""
       write(6,"(a)") "-minerr val: only print errors greater than val"
       write(6,"(a)") "-maxerr val: only print errors less than val"
       write(6,"(a)") "-Qmin val: only consider points with Q >= val (default 0.0)"
       stop
   end if

  idev1 = idev_open_arg(1,status='old')
  idev2 = idev_open_arg(2,status='old')
  Qmin = dble_val_opt("-Qmin",0.0_dp)

  write(6,"(a)") "# " // trim(command_line())
  write(6,"(a)") "# Comparing files " // trim(string_val_arg(1)) // " and " // trim(string_val_arg(2))

  iy = 1; iQ = 1
  outer: do
     read(idev1,'(a)',iostat=iostat) line1
     if (iostat /= 0) exit outer
     !do while(.true.)
     !   if (len_trim(line1) > 0 .and. line1(1:1) /= "") exit
     !end do
     read(idev2,'(a)',iostat=iostat) line2
     if (iostat /= 0) exit outer
     !do while(.true.)
     !   read(line2,*,iostat=iostat) yval2,Qval2,res2(iy,iQ,:)
     !   if (iostat == 0) exit
     !end do
     !read(idev2,'(a)',iostat=iostat) line2
     !if (iostat /= 0) exit

     ! keep track of indicies
     if (trim(line1) == "") then
        ny = iy - 1
        iy = 1
        iQ = iQ + 1
     end if

     read(line1,*,iostat=iostat) yval(iy),Qval(iy,iQ),res1(iy,iQ,:)
     if (iostat /= 0) cycle
     read(line2,*,iostat=iostat) yval2,Qval2,res2(iy,iQ,:)
     if (iostat /= 0) ERROR STOP "Second file has a different format than the first one"
     if (yval(iy) /= yval2 .or. Qval(iy,iQ) /= Qval2) then
        write(6,"(a)") "Mismatch in y or Q values in input files"
        write(6,"(a)") line1
        write(6,"(a)") line2
        stop
     end if
     iy = iy + 1
  end do outer
  nQ = iQ - 1

  write(6,"(a,f10.4,a,f20.4)") "# yrange ", minval(yval(1:ny)), "-", maxval(yval(1:ny))
  write(6,"(a,f10.4,a,f20.4)") "# Qrange ", minval(Qval(1:ny,1:nQ)), "-", maxval(Qval(1:ny,1:nQ))

  ! establish a general protection mask -- if there is a sign change
  ! in the neighbourhood of the point then it should be eliminated.
  ! Work out dy neighbourhood
  idy = nint(0.4_dp / (yval(ny)-yval(ny-1)))
  if (log_val_opt('-protect')) then
     forall(iQ=1:nQ)
        forall(iy=1:ny)
           ! check for sign change in x direction
           mask(iy,iQ,:) = mask(iy,iQ,:) .and. &
                &   minval(res2(max(iy-idy,1):min(iy+idy,ny), iQ,:),dim=1) *&
                &   maxval(res2(max(iy-idy,1):min(iy+idy,ny), iQ,:),dim=1) > 0
           ! check for sign change in Q direction
           mask(iy,iQ,:) = mask(iy,iQ,:) .and. &
                &   minval(res2(iy,max(iQ-1,1):min(iQ+1,nQ),:),dim=1) *&
                &   maxval(res2(iy,max(iQ-1,1):min(iQ+1,nQ),:),dim=1) > 0
        end forall
     end forall
  end if

  ! now get masks for the x range
  xval(:ny) = exp(-yval(:ny))
  mask  (:ny,:nQ,:) = mask(:ny,:nQ,:) .and. spread(Qval(:ny,:nQ)>=Qmin,3,11)
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
   write(6,"(a)") "# Summary of relative errors"
   
   write(6,"(a)") "# max rel error across: guds(x< 0.7), all(x< 0.7), guds(x< 0.9), all(x< 0.9)"
   write(6,"(4es20.10)") max07_guds, max07_all, max09_guds, max09_all
  end if
  
  if (log_val_opt("-channel")) then
     ic = int_val_opt("-channel")
     minerr = dble_val_opt("-minerr",zero)
     maxerr = dble_val_opt("-maxerr",1e100_dp)
     write(6,'(a,i3,a)') "# Flavour channel ", ic, &
        " (-5:5 means difference for that flavour; 11 means max difference across all flavours)"
     write(6,'(a)') "# Columns: y(=ln1/x) Q max_masked_difference max_difference maxflav_masked"
     do iQ = 1, nQ
        do iy = 1, ny
           if (ic == 11) then
              ! do all flavours
              true_err = maxval(err(iy,iQ,-5:5))
              if (any(mask(iy,iQ,:))) then
                 this_err = maxval(err(iy,iQ,-5:5),mask(iy,iQ,-5:5))
                 icmax = maxloc(err(iy,iQ,-5:5),mask=mask(iy,iQ,-5:5),dim=1) - 6 ! 1: -> -5:
              else
                 icmax = 11
                 this_err = -one
              end if
           else 
              icmax = ic
              true_err = err(iy,iQ,ic)
              if (.not.mask(iy,iQ,ic)) then
                 this_err = -one
              else
                 this_err = true_err
              end if
           end if
           if (this_err >= minerr .and. this_err <= maxerr) then
              write(6,'(4es20.12,i3,es20.12)') yval(iy), Qval(iy,iQ), this_err, true_err, icmax, merge(res1(iy,iQ,ic),zero,ic/=11)
           end if
        end do
        write(6,*)
     end do
  end if

  !-- security ----------------------
  if (.not. CheckAllArgsUsed(0)) stop
  !----------------------------------

end program compare2file_v2
