!
! First results: in general, need to watch out for placing 
!                first Q bin too low, because it "cuts" the evolution
!                and one doesn't really get the full result
!
! Worst channel for Q case is b+bbar
!
! Actually "best" way of testing Q evolution might be to 
! rerun the evolution from scratch each time, maybe even in a
! fixed-flavour number scheme? Answers might change by a factor of
! two roughly.
!
! Where first Qval is 5.0 GeV, in a VFNS, we're basically seeing
! NB: this is all wrong....

! du       err (always channel 11 = b+bbar, Q=5,y=0.1)
! 0.150    2.658068E-04
! 0.100    6.027373E-05
! 0.070    1.308940E-05
! 0.050    4.379481E-06
! 0.035    8.875060E-07
! 0.025    2.433153E-07
!
! The best channels seem to be uv,dv,ubar+dbar,s+sbar, about 6-10
! times better.
! 
! The above results correspond roughly to expectation of eps^4
module accuracy_helper
   use hoppet; use sort !(not included with hoppet)
  implicit none

  private

  integer, parameter, public :: accflv_n = 13
  integer, parameter, public :: accflv_g   = 1  ! gluon
  integer, parameter, public :: accflv_sng = 2  ! singlet
  integer, parameter, public :: accflv_uv  = 3  ! u - ubar
  integer, parameter, public :: accflv_dv  = 4  ! d - dbar
  integer, parameter, public :: accflv_umd = 5  ! u - d
  integer, parameter, public :: accflv_Lp  = 6  ! ubar + dbar
  integer, parameter, public :: accflv_Lm  = 7  ! ubar - dbar
  integer, parameter, public :: accflv_ss  = 8  ! s + sbar
  integer, parameter, public :: accflv_cs  = 9  ! c + cbar
  integer, parameter, public :: accflv_bs  = 10 ! b + bbar
  integer, parameter, public :: accflv_sv  = 11 ! s - sbar
  integer, parameter, public :: accflv_cv  = 12 ! c - cbar
  integer, parameter, public :: accflv_bv  = 13 ! b - bbar

  character(len=10), parameter, public :: accflv_names(accflv_n) = (/&
       "gluon    ", & ! accflv_g   = 1  !
       "singlet  ", & ! accflv_sng = 2  !
       "u-ubar   ", & ! accflv_uv  = 3  !
       "d-dbar   ", & ! accflv_dv  = 4  !
       "u-d      ", & ! accflv_umd = 5  !
       "ubar+dbar", & ! accflv_Lp  = 6  !
       "ubar-dbar", & ! accflv_Lm  = 7  !
       "s+sbar   ", & ! accflv_ss  = 8  !
       "c+cbar   ", & ! accflv_cs  = 9  !
       "b+bbar   ", & ! accflv_bs  = 10 !
       "s-sbar   ", & ! accflv_sv  = 11 !
       "c-cbar   ", & ! accflv_cv  = 12 !
       "b-bbar   "/)  ! accflv_bv  = 13 !

  public :: FillRefTable

  real(dp), parameter :: delta_y = 0.03_dp

  interface percentiles
     module procedure percentiles_1d, percentiles_2d, percentiles_3d
  end interface
  public :: percentiles

  public :: add_du_info

contains

  !======================================================================
  !! given a grid return a table of distributions, and optionally
  !! return a table of reference normalizations -- these are identical
  !! to the distributions except when there is a nearby sign change, in
  !! which case they are deduced from the average of the function
  !! at y+-delta_y.
  !!
  !! Currently best suited to cases where all(Qvals > Qinit)
  subroutine FillRefTable(grid, dh, yvals, Qvals, table, table_norm)
    type(grid_def),     intent(in)  :: grid
    type(dglap_holder), intent(in)  :: dh
    real(dp),           intent(in)  :: yvals(:), Qvals(:)
    real(dp),           intent(out) :: table(:,:,:)
    real(dp), optional, intent(out) :: table_norm(:,:,:)
    !--------------------------
    type(running_coupling) :: coupling
    integer :: ny, nQ, nflv, iQ, iy
    real(dp)          :: lastQ, slice(-6:6)
    real(dp)          :: accref_yplus(accflv_n), accref_yminus(accflv_n)
    real(dp), pointer :: pdf(:,:)
    real(dp)          :: Qinit

    ny   = assert_eq(size(yvals),size(table,1),'FillRefTable')
    nQ   = assert_eq(size(Qvals),size(table,3),'FillRefTable')
    nflv = assert_eq(accflv_n,   size(table,2),'FillRefTable')

    call AllocPDF(grid, pdf)
    pdf = unpolarized_dummy_pdf(xValues(grid))
    Qinit = sqrt(two)  ! Vogt starting scale

    call InitRunningCoupling(coupling,alfas=0.35_dp,Q=Qinit,nloop=3)

    lastQ = Qinit
    do iQ = 1, nQ
       call evolvePDF(dh, pdf, coupling, lastQ, Qvals(iQ))
       lastQ = Qvals(iQ)
       
       do iy = 1, ny
          slice = pdf(:,-6:6) .aty. (yvals(iy).with.grid)
          table(iy,:,iQ) = human2acc(slice)
          if (present(table_norm)) then
             accref_yplus  = human2acc(pdf(:,-6:6) .aty. (&
                  &           (yvals(iQ)+delta_y).with.grid))
             accref_yminus = human2acc(pdf(:,-6:6) .aty. (&
                  &           max(zero,yvals(iQ)-delta_y).with.grid))
             ! look out for sign changes
             where(accref_yplus*accref_yminus < zero)
                ! essentially take derivative as normalisation
                ! ( == we will look at error on position of zero).
                table_norm(iy,:,iQ) = &
                     &              (abs(accref_yplus - accref_yminus))&
                     &               / (two*delta_y)
             else where
                table_norm(iy,:,iQ) = abs(table(iy,:,iQ))
             end where
          end if
       end do
       
    end do

    ! eliminate problems associated with sign changes in Q...
    ! (recall we have already taken abs value).
    do iQ = 2, nQ-1
       where (table_norm(:,:,iQ) < table_norm(:,:,iQ+1) .and. &
            & table_norm(:,:,iQ) < table_norm(:,:,iQ-1))
          table_norm(:,:,iQ) = half * &
               &    (table_norm(:,:,iQ-1)+table_norm(:,:,iQ+1))
       end where
    end do
    

    call Delete(coupling)
    call Delete(pdf)
  end subroutine FillRefTable
  

  !======================================================================
  !! Given a pdf slice in "human" format, return one in the format
  !! needed for our accuracy tests (marginally different from those
  !! used with Vogt).
  function human2acc(slice) result(accslice)
    real(dp), intent(in) :: slice(ncompmin:)
    real(dp)             :: accslice(accflv_n)

    accslice(accflv_g  ) = slice(iflv_g)
    accslice(accflv_sng) = sum(slice(1:6)) + sum(slice(-6:-1))
    accslice(accflv_uv ) = slice(iflv_u   ) - slice(iflv_ubar)
    accslice(accflv_dv ) = slice(iflv_d   ) - slice(iflv_dbar)
    accslice(accflv_umd) = slice(iflv_u   ) - slice(iflv_d   )
    accslice(accflv_Lp ) = slice(iflv_ubar) + slice(iflv_dbar)
    accslice(accflv_Lm ) = slice(iflv_ubar) - slice(iflv_dbar)
    accslice(accflv_ss ) = slice(iflv_s   ) + slice(iflv_sbar)
    accslice(accflv_cs ) = slice(iflv_c   ) + slice(iflv_cbar)
    accslice(accflv_bs ) = slice(iflv_b   ) + slice(iflv_bbar)
    accslice(accflv_sv ) = slice(iflv_s   ) - slice(iflv_sbar)
    accslice(accflv_cv ) = slice(iflv_c   ) - slice(iflv_cbar)
    accslice(accflv_bv ) = slice(iflv_b   ) - slice(iflv_bbar)

  end function human2acc
  

  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolazrized evolution
  function unpolarized_dummy_pdf(xvals) result(dummy_pdf)
    real(dp), intent(in) :: xvals(:)
    real(dp)             :: dummy_pdf(size(xvals),ncompmin:ncompmax)
    real(dp) :: uv(size(xvals)), dv(size(xvals))
    real(dp) :: ubar(size(xvals)), dbar(size(xvals))
    !---------------------
    real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
    real(dp), parameter :: N_uv=5.107200_dp, N_dv = 3.064320_dp
    real(dp), parameter :: N_db = half*N_ls
  
    dummy_pdf = zero ! automatically in human rep

    !-- remember that these are all xvals*q(xvals)
    uv = N_uv * xvals**0.8_dp * (1-xvals)**3
    dv = N_dv * xvals**0.8_dp * (1-xvals)**4
    dbar = N_db * xvals**(-0.1_dp) * (1-xvals)**6
    ubar = dbar * (1-xvals)
    dummy_pdf(:,iflv_g) = N_g * xvals**(-0.1_dp) * (1-xvals)**5
        
    dummy_pdf(:,-iflv_s) = 0.2_dp*(dbar + ubar)
    dummy_pdf(:, iflv_s) = dummy_pdf(:,-iflv_s)
    dummy_pdf(:, iflv_u) = uv + ubar
    dummy_pdf(:,-iflv_u) = ubar
    dummy_pdf(:, iflv_d) = dv + dbar
    dummy_pdf(:,-iflv_d) = dbar
  end function unpolarized_dummy_pdf


  !======================================================================
  !! Return the percentiles of array as requested in percentiles_wanted,
  !! optionally masked.
  !!
  !! Solution used (sort compacted version of masked array) is OK if
  !! mask is one-off or there are many small masks, but suboptimal if
  !! one repeatedly investigates masks which are nearly always 1 [in
  !! that case a better solution would be to sort the array once and 
  !! for all and then work out the percentiles of subsets of the
  !! sorted array].
  !!
  !! Percentiles should be as a fraction of 1.
  function percentiles_1d(array,percentiles_wanted, mask) result(res)
    real(dp), intent(in), target   :: array(:), percentiles_wanted(:)
    logical,  intent(in), optional :: mask(:)
    real(dp)                       :: res(size(percentiles_wanted))
    !------------------------------------------------------
    real(dp), target  :: masked_array(size(array))
    real(dp), pointer :: ref_array(:)
    integer           :: n, i, ell, index_array(size(array))

    if (present(mask)) then
       ! make a compact copy of the masked array
       ell = 0
       do i = 1, size(array)
          if (mask(i)) then
             ell = ell + 1
             masked_array(ell) = array(i)
          end if
       end do
       ref_array => masked_array(1:ell)
    else
       ref_array => array(:)
    end if
    
    n = size(ref_array)
    forall(i = 1:n) index_array(i) = i
    call indexx(ref_array, index_array(1:n))

    forall(i = 1:size(percentiles_wanted)) 
       res(i) = ref_array(index_array(max(1,nint(percentiles_wanted(i)*n))))
    end forall
  end function percentiles_1d


  !======================================================================
  !! 2d version of percentiles
  function percentiles_2d(array, percentiles_wanted, mask) result(res)
    real(dp), intent(in), target   :: array(:,:), percentiles_wanted(:)
    logical,  intent(in), optional :: mask(:,:)
    real(dp)                       :: res(size(percentiles_wanted))
    !-----------------------------------------------------
    real(dp), target  :: masked_array(size(array))
    integer           :: i, j, ell
    if (present(mask)) then
       ! make a compact copy of the masked array
       ell = 0
       do j = 1, size(array,2)
          do i = 1, size(array,1)
             if (mask(i,j)) then
                ell = ell + 1
                masked_array(ell) = array(i,j)
             end if
          end do
       end do
       res = percentiles(masked_array(1:ell), percentiles_wanted)
    else
       masked_array = reshape(array,(/ size(array) /))
       res = percentiles(masked_array, percentiles_wanted)
    end if
  end function percentiles_2d
  

  !======================================================================
  !! 3d version of percentiles
  function percentiles_3d(array, percentiles_wanted, mask) result(res)
    real(dp), intent(in), target   :: array(:,:,:), percentiles_wanted(:)
    logical,  intent(in), optional :: mask(:,:,:)
    real(dp)                       :: res(size(percentiles_wanted))
    !-----------------------------------------------------
    real(dp), target  :: masked_array(size(array))
    integer           :: i, j, k, ell
    if (present(mask)) then
       ! make a compact copy of the masked array
       ell = 0
       do k = 1, size(array,3)
          do j = 1, size(array,2)
             do i = 1, size(array,1)
                if (mask(i,j,k)) then
                   ell = ell + 1
                   masked_array(ell) = array(i,j,k)
                end if
             end do
          end do
       end do
       res = percentiles(masked_array(1:ell), percentiles_wanted)
    else
       masked_array = reshape(array,(/ size(array) /))
       res = percentiles(masked_array, percentiles_wanted)
    end if
  end function percentiles_3d


  !======================================================================
  !! takes a grid info string and adds info on current du
  subroutine add_du_info(string,du)
    character(len=*), intent(inout) :: string
    real(dp),         optional      :: du
    !-------------------
    character(len=80) duinfo

    write(duinfo,'("|du=",f6.4)') default_or_opt(DefaultEvolutionDu(),du)
    string = trim(string)//trim(duinfo)
  end subroutine add_du_info
  
end module accuracy_helper


!======================================================================
!! Program options
!! -du duval :   establish quality of specified duval
!! -dytable  :   run with a table of ty
!! -big -div2 n -div3 n : loop over ymax vals, orders, dy
program determine_accuracies
   use hoppet; use io_utils
  use accuracy_helper
  implicit none
  real(dp)           :: ymax = 12.0_dp, dy_ref, dy, du_ref
  type(grid_def)     :: grid_ref, grid, gridarray(5)
  type(dglap_holder) :: dh_ref, dh
  integer            :: iy, iQ, i
  integer, parameter :: ny = 443, nQ = 6, ref_order=5
  ! leave quite a large spacing otherwise with coarse Q grid we
  ! learn very little...
  real(dp), parameter :: Qvals(nQ) = &
       &       (/5.0_dp,10.0_dp,30.0_dp,100.0_dp,1000.0_dp,1e4_dp/)
  real(dp), target :: yvals(ny), table(ny,accflv_n,nQ), du
  real(dp), target :: reftable(ny,accflv_n,nQ), reftable_norm(ny,accflv_n,nQ)
  real             :: test_run_time = 0, test_init_time = 0, start_time

  character(len=80) :: grid_string_ref, grid_string_test
  integer           :: outdev

  if (log_val_opt('-out')) then
     outdev = idev_open_opt("-out")
  else
     outdev = 6
  end if

  write(outdev,'(a)') '# '//trim(command_line())

  !dy_ref = 0.0125_dp
  dy_ref = dble_val_opt('-dy_ref',0.10_dp)
  du_ref = dble_val_opt('-du_ref',0.01_dp)
  write(0,*) 'dy_ref = ', dy_ref
  write(0,*) 'du_ref = ', du_ref
  call InitGridDef(gridarray(3),dy_ref/9.0_dp, 0.5_dp, ref_order)
  call InitGridDef(gridarray(2),dy_ref/3.0_dp, 2.0_dp, ref_order)
  call InitGridDef(gridarray(1),dy_ref,        ymax,   ref_order)
  call InitGridDef(grid_ref, gridarray(1:3), locked=.true.)
  call Delete(gridarray(1:3))

  call InitDglapHolder(grid_ref,dh_ref,nloop=3,nflo=3,nfhi=6)
  write(0,*) 'Reference initialisation done'

  ! uniform spacing (maybe non-uniform would be better?)
  forall(iy=1:ny) yvals(iy) = 0.1_dp + (iy-1)*(ymax-0.1_dp)/(ny-1)

  call SetDefaultEvolutionDu(du_ref)
  call GetGridInfoString(grid_ref,grid_string_ref)
  call add_du_info(grid_string_ref)

  call FillRefTable(grid_ref,dh_ref,yvals,Qvals,reftable,reftable_norm)

  if (log_val_opt('-du')) then
     du = dble_val_opt('-du')
     call SetDefaultEvolutionDu(du)
     write(outdev,*) "Evolution du = ", du
     grid = grid_ref
     dh   = dh_ref
     test_init_time = 0.0
     !call FillRefTable(grid,dh_ref,yvals,Qvals,table)
     call FillRefTable_with_timings()
     call printout_accuracies
  end if


  if (log_val_opt('-dytable')) then
     dy = 0.4
     do 
        if (dy < 2*dy_ref) exit
        call InitGridDef(gridarray(3),dy/9.0_dp, 0.5_dp, ref_order)
        call InitGridDef(gridarray(2),dy/3.0_dp, 2.0_dp, ref_order)
        call InitGridDef(gridarray(1),dy,        ymax,   ref_order)
        call InitGridDef(grid, gridarray(1:3), locked=.true.)
        call Delete(gridarray(1:3))

        call cpu_time(start_time)
        call InitDglapHolder(grid,dh,nloop=3,nflo=3,nfhi=6)
        call cpu_time(test_init_time); 
        test_init_time = test_init_time - start_time

        !call FillRefTable(grid,dh,yvals,Qvals,table)
        call FillRefTable_with_timings()
        call printout_accuracies

        call Delete(grid)
        call Delete(dh)
        dy = dy * 0.8_dp
     end do
  end if

  if (log_val_opt('-big')) call big_table

contains

  !======================================================================
  subroutine big_table
    real(dp) :: inner_ymax(1:3), lcl_ymax(1:3), dyarr(3)
    integer  :: order1,order2,order3, divisors(1:3), orders(3)
    integer, parameter :: inner_ymax_lims(2:3) = (/ 3.5_dp, 1.0_dp/)
    integer  :: ndone
    logical  :: dummy

    dummy = log_val_opt('-dummy')
    ndone = 0

    divisors(1) = 1
    divisors(2) = int_val_opt('-div2')
    divisors(3) = int_val_opt('-div3') * divisors(2)

    inner_ymax(1) = ymax

    ! ymax(3) loop
    inner_ymax(3) = 0.3_dp
    do 
       if (inner_ymax(3) > inner_ymax_lims(3)) exit

       ! ymax(2) loop
       inner_ymax(2) = 1.5_dp * inner_ymax(3)
       do
          if (inner_ymax(2) > inner_ymax_lims(2)) exit

          ! order loops
          do order1 = -6, 6
          do order2 = -6, 6
          do order3 = -6, 6

             orders = (/order1,order2,order3/)
             if (any(abs(orders) < 4)) cycle

             dy = 0.5
             do 
                if (dy < 2*dy_ref) exit

                ! adjust things to be exactly on-grid
                dyarr(1)  = ymax / nint(ymax/dy)
                dyarr(2:) = dyarr(1) / divisors(2:)
                
                lcl_ymax(1)  = ymax
                lcl_ymax(2:) = nint(inner_ymax(2:)/dyarr(2:))*dyarr(2:)

                ! check to make sure we have room for the orders
                ! we're asking for...
                !write(0,*) '-----------------------------'
                !write(0,*) dy
                !write(0,*) dyarr
                !write(0,*) lcl_ymax
                !write(0,*) orders
                !write(0,*) nint(lcl_ymax/dyarr)
                
                if (.not.any(nint(lcl_ymax/dyarr) < abs(orders)+1)) then
                   call InitGridDef(gridarray(3),dyarr(3), lcl_ymax(3), order3)
                   call InitGridDef(gridarray(2),dyarr(2), lcl_ymax(2), order2)
                   call InitGridDef(gridarray(1),dyarr(1), lcl_ymax(1), order1)
                   call InitGridDef(grid, gridarray(1:3), locked=.true.)
                   call Delete(gridarray(1:3))
                   
                   call cpu_time(start_time)
                   if (.not.dummy) &
                        &call InitDglapHolder(grid,dh,nloop=3,nflo=3,nfhi=6)
                   
                   call cpu_time(test_init_time); 
                   test_init_time = test_init_time - start_time
                   
                   !call FillRefTable(grid,dh,yvals,Qvals,table)
                   if (.not.dummy) then
                      call FillRefTable_with_timings()
                   else
                      table = zero
                   end if
                   !call printout_accuracies
                   if (.not.dummy .or. mod(ndone,1000)==0) call printout_accuracies
                   
                   call Delete(grid)
                   if (.not.dummy) call Delete(dh)
                   ndone = ndone + 1
                   write(0,*) ndone
                end if
                dy = dy * 0.8_dp
             end do ! dy loop
          end do
          end do
          end do    ! order loops

          inner_ymax(2) = inner_ymax(2) + 0.3_dp
       end do       ! inner_ymax(2) loop

       inner_ymax(3) = inner_ymax(3) + 0.1_dp
    end do          ! inner_ymax(3) loop
    
  end subroutine big_table



  !======================================================================
  subroutine FillRefTable_with_timings()
    real    :: ref_time
    integer :: nrun
    call cpu_time(ref_time); nrun = 0
    do
       call FillRefTable(grid,dh,yvals,Qvals,table)
       call cpu_time(test_run_time); nrun = nrun + 1
       test_run_time = test_run_time - ref_time
       if (test_run_time > 1.0) exit
    end do
    test_run_time = test_run_time / nrun
  end subroutine FillRefTable_with_timings


  !======================================================================
  !! A quick and dirty routine to print out accuracies, etc...
  subroutine printout_accuracies
    real(dp), pointer :: normtable(:,:,:)
    logical           :: ymask(ny,accflv_n,nQ)
    logical           :: flvmask(ny,accflv_n,nQ), flvmask2(ny,accflv_n,nQ)
    integer ::  iflv, posn(3,6), j
    !real(dp) :: percentls(6) = (/0.0_dp,0.5_dp, 0.9_dp, 0.99_dp, 0.999_dp, 1.0_dp /)
    real(dp) :: percentls(1) = (/0.99_dp/)
    real(dp) :: perc_res(size(percentls),6)

    forall (iy   = 1:ny) ymask(iy,:,:) = exp(-yvals(iy)) < 0.7_dp
    forall (iflv = 1:accflv_n) flvmask(:,iflv,:)  = iflv <=8
    forall (iflv = 1:accflv_n) flvmask2(:,iflv,:) = iflv <=5

    ! decide whether to use our special special table designed to
    ! avoid "dangerous zeros"
    if (log_val_opt('-norm')) then
       normtable => reftable_norm
    else
       normtable => reftable
    end if

    where (reftable /= zero)
       table = abs((table-reftable) / normtable)
    elsewhere
       table = zero
    end where

    posn(:,1) = maxloc(table)
    posn(:,2) = maxloc(table,ymask)            
    posn(:,3) = maxloc(table,flvmask)          
    posn(:,4) = maxloc(table,flvmask.and.ymask)
    posn(:,5) = maxloc(table,flvmask2)          
    posn(:,6) = maxloc(table,flvmask2.and.ymask)
    ! now get some percentiles
    perc_res(:,1) = percentiles(table,percentls)
    perc_res(:,2) = percentiles(table,percentls,ymask)            
    perc_res(:,3) = percentiles(table,percentls,flvmask)          
    perc_res(:,4) = percentiles(table,percentls,flvmask.and.ymask)
    perc_res(:,5) = percentiles(table,percentls,flvmask2)          
    perc_res(:,6) = percentiles(table,percentls,flvmask2.and.ymask)

    write(outdev,'(78("="))')
    write(outdev,'(a)') "Ref.: "//trim(grid_string_ref)
    call GetGridInfoString(grid,grid_string_test)
    call add_du_info(grid_string_test)
    write(outdev,'(a)') "Test: "//trim(grid_string_test)

    write(outdev,'(a,2es9.1,2f10.5)') & 
         & "Summary(max(4,5);99%(3);run/s;init/s):", &
         & max(table(posn(1,4),posn(2,4),posn(3,4)), &
         &     table(posn(1,5),posn(2,5),posn(3,5))),&
         & perc_res(:,3), &
         & test_run_time, &
         & test_init_time
    write(outdev,'(78("-"))')

    write(outdev,'("Accuracy/99%: ",es7.1," (chnl ",i2,")/",es7.1,&
         &" ; x<=0.7: ",es7.1," (chnl ",i2,")/",es7.1)') &
         & table(posn(1,1),posn(2,1),posn(3,1)), posn(2,1), perc_res(1,1),&
         & table(posn(1,2),posn(2,2),posn(3,2)), posn(2,2), perc_res(1,2)
    write(outdev,'("Acc(1-8)/99%: ",es7.1," (chnl ",i2,")/",es7.1,&
         &" ; x<=0.7: ",es7.1," (chnl ",i2,")/",es7.1)') &
         & table(posn(1,3),posn(2,3),posn(3,3)), posn(2,3), perc_res(1,3),&
         & table(posn(1,4),posn(2,4),posn(3,4)), posn(2,4), perc_res(1,4)
    write(outdev,'("Acc(1-5)/99%: ",es7.1," (chnl ",i2,")/",es7.1,&
         &" ; x<=0.7: ",es7.1," (chnl ",i2,")/",es7.1)') &
         & table(posn(1,5),posn(2,5),posn(3,5)), posn(2,5), perc_res(1,5),&
         & table(posn(1,6),posn(2,6),posn(3,6)), posn(2,6), perc_res(1,6)

    !write(outdev,'(a,es8.1,a,i2,a,es8.1,a,i2,a)') &
    !     &   'Worst acc(1-8):',table(posn(1,3),posn(2,3),posn(3,3)),&
    !     &   ' (channel ', posn(2,3),'); for x<=0.7:',&
    !     &   table(posn(1,4),posn(2,4),posn(3,4)),&
    !     &   ' (channel ', posn(2,4),')'
    !write(outdev,'(a,es8.1,a,i2,a,es8.1,a,i2,a)') &
    !     &   'Worst acc(1-5):',table(posn(1,5),posn(2,5),posn(3,5)),&
    !     &   ' (channel ', posn(2,5),'); for x<=0.7:',&
    !     &   table(posn(1,6),posn(2,6),posn(3,6)),&
    !     &   ' (channel ', posn(2,6),')'
    !do j = 1, size(percentls)
    !   write(outdev,'("percentile ",f6.3,":",6es8.1)') percentls(j), perc_res(j,:)
    !end do


    write(outdev,'(78("-"))')
    write(outdev,'(a2,a12,a8,a7,a8,a8," |",a8,a7,a8,a8)') &
         & "","channel","worst","y","Q","99%","x<0.7","y","Q","99%"
    write(outdev,'(78("-"))')
    do iflv = 1, accflv_n
       posn((/1,3/),1) = maxloc(table(:,iflv,:))
       posn((/1,3/),2) = maxloc(table(:,iflv,:),ymask(:,iflv,:))
       perc_res(:,1) = percentiles(table(:,iflv,:),percentls(:))
       perc_res(:,2) = percentiles(table(:,iflv,:),percentls(:),ymask(:,iflv,:))
       write(outdev,&
            &'(i2,a12,es8.1,f7.3,f8.1,es8.1," |",es8.1,f7.3,f8.1,es8.1)') &
            &  iflv, accflv_names(iflv), &
            &   table(posn(1,1),iflv,posn(3,1)),&  
            &   yvals(posn(1,1)), Qvals(posn(3,1)),&
            &   perc_res(:,1), &
            &   table(posn(1,2),iflv,posn(3,2)),&  
            &   yvals(posn(1,2)), Qvals(posn(3,2)),&
            &   perc_res(:,2)
       if (iflv == 5 .or. iflv == accflv_n) write(outdev,'(78("-"))')
    end do
  end subroutine printout_accuracies

end program determine_accuracies
