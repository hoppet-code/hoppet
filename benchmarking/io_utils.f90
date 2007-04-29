!! io_utils.f90
!!
!! version 1.1, GPS, 24 September 2003
!!
!======================================================================
!! All the various bits and pieces which could be useful regarding i/o,
!! such as opening the right files, getting numerical values for
!! command line arguments, etc...
!======================================================================


! unfortunate that we need an extra module...
module sub_defs_io_consts
  integer, parameter :: max_arg_len  = 180
  integer, parameter :: max_line_len = 600
end module sub_defs_io_consts



!----------------------------------------------------------------------
!! All the interfaces needed here and in the lcl_???.f90 routines
!----------------------------------------------------------------------
module sub_defs_io
  use sub_defs_io_consts

  interface
     subroutine open_arg(iarg, idev, ext, status, form, default_val)
       implicit none
       integer, intent(in) :: iarg, idev
       character(*), optional, intent(in) :: ext, status, form, default_val
     end subroutine open_arg
  end interface

  interface
     function idev_open_opt(opt, default_val, ext, status, form) result(idev)
       implicit none
       character(*), intent(in) :: opt
       character(*), optional, intent(in) :: default_val, ext, status, form
       integer  :: idev
     end function idev_open_opt
  end interface
  
  interface
     function idev_open_arg(iarg, ext, status, form) result(idev)
       implicit none
       integer             :: idev
       integer, intent(in) :: iarg
       character(*), optional, intent(in) :: ext, status, form
     end function idev_open_arg
  end interface

  interface
     function iargc_opt(opt)
       implicit none
       integer                      :: iargc_opt
       character(len=*), intent(in) :: opt
     end function iargc_opt
  end interface

  interface
     function command_line()
       use sub_defs_io_consts
       implicit none
       character(len=max_line_len) ::  command_line
     end function command_line
  end interface

  interface
     function string_val_opt(opt,default_val)
       use sub_defs_io_consts
       implicit none
       character(len=max_arg_len)             :: string_val_opt
       character(len=*), intent(in)           :: opt
       character(len=*), intent(in), optional :: default_val
     end function string_val_opt
  end interface
  
  interface
     function int_val_opt(opt,default_val)
       implicit none
       integer                      :: int_val_opt
       character(len=*), intent(in) :: opt
       integer, intent(in), optional :: default_val
     end function int_val_opt
  end interface

  interface
     function log_val_opt(opt,default_val)
       implicit none
       logical                       :: log_val_opt
       character(len=*), intent(in)  :: opt
       logical, intent(in), optional :: default_val
     end function log_val_opt
  end interface

  interface
     function dble_val_opt(opt,default_val)
       implicit none
       real(kind(1d0))                      :: dble_val_opt
       character(len=*), intent(in) :: opt
       real(kind(1d0)), intent(in), optional :: default_val
     end function dble_val_opt
  end interface

  interface
     integer function int_value(string)
        implicit none
        character(*),  intent(in)  ::   string
     end function int_value
  end interface

  interface
     integer function int_val_arg(argk, default_val)
       implicit none
       integer, intent(in)           :: argk
       integer, intent(in), optional :: default_val
     end function int_val_arg
  end interface

  interface
     real(kind(1d0)) function dble_val_arg(argk, default_val)
       implicit none
       integer,         intent(in)           :: argk
       real(kind(1d0)), intent(in), optional :: default_val
     end function dble_val_arg
  end interface

  interface
     function string_val_arg(argk, default_val)
       use sub_defs_io_consts
       implicit none
       character(len=max_arg_len)             :: string_val_arg
       integer,          intent(in)           :: argk
       character(len=*), intent(in), optional :: default_val
     end function string_val_arg
  end interface

  interface
     real(kind(1d0)) function dble_value(string)
       implicit none
       character(*),  intent(in)  ::   string
     end function dble_value
  end interface


  interface num2char
     function int2char(num,frmt)
       implicit none
       character(len=30) :: int2char
       integer, intent(in) :: num
       character(*), intent(in), optional :: frmt
     end function int2char
     function sp2char(num,frmt)
       implicit none
       character(len=30) :: sp2char
       real(kind(1.0)), intent(in) :: num
       character(*), intent(in), optional :: frmt
     end function sp2char
     function dp2char(num,frmt)
       implicit none
       character(len=30) :: dp2char
       real(kind(1d0)), intent(in) :: num
       character(*), intent(in), optional :: frmt
     end function dp2char
  end interface
  

  interface
     function get_new_device() result(dev)
       implicit none
       integer :: dev
     end function get_new_device
  end interface
  
  interface
     subroutine time_stamp(idev, string)
       implicit none
       integer,          intent(in),  optional :: idev
       character(len=*), intent(out), optional :: string
     end subroutine time_stamp
  end interface
  
  interface
     subroutine error_report(subname,line1,line2,line3,line4,action)
       implicit none
       character(*), intent(in) :: subname
       character(*), intent(in), optional :: line1, line2, line3, line4
       character(*), intent(in), optional :: action
     end subroutine error_report
  end interface
       
       
   !==================================================================
   ! ====
   ! Interfaces to the standard routines which provide the links to
   ! platform dependent calls
   !======================================================================

  interface
     subroutine lcl_write_2d(idev, array)
       implicit none
       integer,         intent(in) :: idev
       real(kind(1d0)), intent(in) :: array(1:, 1:)
     end subroutine lcl_write_2d
  end interface

  interface
     subroutine lcl_read_2d(idev, array, ifail)
       implicit none
       integer,           intent(in)  :: idev
       real(kind(1d0)),   intent(out) :: array(1:, 1:)
       integer, optional, intent(out) :: ifail
     end subroutine lcl_read_2d
  end interface

  interface
     integer function lcl_iargc()
     end function lcl_iargc
  end interface

  interface
     subroutine lcl_getarg(k, argk)
       implicit none
       integer,      intent(in)  :: k
       character(*), intent(out) :: argk
     end subroutine lcl_getarg
  end interface

  interface
     subroutine lcl_flush(idev)
       implicit none
       integer, intent(in) :: idev
     end subroutine lcl_flush
  end interface

  
  interface  
     function CheckAllArgsUsed(ErrDev) result(check)
       logical :: check
       integer, optional :: ErrDev
     end function CheckAllArgsUsed
  end interface
  interface  
     function CheckAllOptsUsed(ErrDev) result(check)
       logical :: check
       integer, optional :: ErrDev
     end function CheckAllOptsUsed
  end interface
  
end module sub_defs_io


!--------------------------------------------------------------
!! When an argument is used register the fact here. Then at the
!! end of program initialisation one can check to see whether all
!! arguments have been used
module track_args
  ! putting the "only" here ensures that this compiles under intel 7.1
  ! (not immediately clear why, though it could be related to fact
  ! that interface CheckAllArgsUsed is being used, but that the 
  ! corresponding subroutine uses track_args?)
  use sub_defs_io, only: lcl_iargc, lcl_getarg
  use sub_defs_io_consts
  implicit none
  private
  logical, allocatable, save ::argsused(:)
  logical :: started = .false.
  integer, save :: narg
  character(len=*), parameter :: sep = &
       &'=============================================================='

  public :: ta_RegisterUsed, ta_CheckAllArgsUsed, ta_CheckAllOptsUsed

contains

  subroutine ta_RegisterUsed(iarg)
    integer, intent(in)  :: iarg
    if (.not. started) then
       narg = lcl_iargc()
       allocate(argsused(narg))
       argsused = .false.
       started = .true.
    end if
    
    if (iarg > narg .or. iarg < 1) then
       write(0,*) 'In ta_RegisterUsed tried to set illegal arg', iarg
       stop
    end if
    argsused(iarg) = .true.
  end subroutine ta_RegisterUsed

  
  !--------------------------------------------------------
  !! ErrDev = device to which one should send the list of 
  !! unused args
  function ta_CheckAllOptsUsed(ErrDev) result(check)
    logical :: check
    integer, optional :: ErrDev
    integer :: i
    character(len=max_arg_len) :: argi
    logical :: writeit

    check = .true.
    writeit = present(ErrDev)
    do i = 1, narg
       if (argsused(i)) cycle
       call lcl_getarg(i, argi)
       if (argi(1:1) == "-") then
          if (check .and. writeit) then
             write(ErrDev,*) sep
             write(ErrDev,*) &
                  &'WARNING / ERROR: the following options were not recognized'
          end if
          check = .false.
          write(ErrDev,*) trim(argi)
       end if
    end do
    if ((.not. check) .and. writeit ) write(ErrDev,*) sep
  end function ta_CheckAllOptsUsed


  !--------------------------------------------------------
  !! ErrDev = device to which one should send the list of 
  !! unused args
  function ta_CheckAllArgsUsed(ErrDev) result(check)
    logical :: check
    integer, optional :: ErrDev
    integer :: i
    character(len=max_arg_len) :: argi

    if (narg > 0 ) then 
       check = all(argsused)
    else
       check = .true. 
       return
    end if
    if ((.not. check) .and. present(ErrDev)) then
       write(ErrDev,*) sep
       write(ErrDev,*) &
            &'WARNING / ERROR: the following args were not recognized'
       do i = 1, narg 
          if (.not. argsused(i)) then
             call lcl_getarg(i, argi)
             write(ErrDev,*) trim(argi)
          end if
       end do
       write(ErrDev,*) sep
    end if
  end function ta_CheckAllArgsUsed
  
  
end module track_args




!--------------------------------------------------------
!! ErrDev = device to which one should send the list of 
!! unused args
function CheckAllArgsUsed(ErrDev) result(check)
  use track_args
  logical :: check
  integer, optional :: ErrDev
  check = ta_CheckAllArgsUsed(ErrDev)
end function CheckAllArgsUsed
function CheckAllOptsUsed(ErrDev) result(check)
  use track_args
  logical :: check
  integer, optional :: ErrDev
  check = ta_CheckAllOptsUsed(ErrDev)
end function CheckAllOptsUsed


!----------------------------------------------------------------------
!! As open_arg -- except that it automatically allocates and returns a
!! device number
!----------------------------------------------------------------------
function idev_open_arg(iarg, ext, status, form) result(idev)
  use sub_defs_io, except => idev_open_arg
  implicit none
  integer             :: idev
  integer, intent(in) :: iarg
  character(*), optional, intent(in) :: ext, status, form
  !--------------------------------------------------------------------
  idev = get_new_device()
  call open_arg(iarg, idev, ext, status, form)
end function idev_open_arg


!---------------------------------------------------------------------
!! Open file indicated by command line option and return the 
!! device number.
!!
!! GPS 05/01/01
!---------------------------------------------------------------------
function idev_open_opt(opt, default_val, ext, status, form) result(idev)
  use sub_defs_io, except => idev_open_opt
  implicit none
  character(*), intent(in) :: opt
  character(*), optional, intent(in) :: default_val, ext, status, form
  integer  :: idev, iarg
  !--------------------------------------------------------------------
  idev = get_new_device()
  iarg = iargc_opt(opt) + 1 ! will return -100 if opt is not present
  if (iarg < 0) then
     write(0,*) 'ERROR in idev_open_opt:'
     write(0,*) 'command-line option "'//opt//'" missing'
     stop
  end if
  call open_arg(iarg, idev, ext, status, form, default_val)
end function idev_open_opt


!----------------------------------------------------------------------
!! Open file corresponding the name in command line argument 'iarg'
!! 
!! GPS 20/03/97
!----------------------------------------------------------------------
subroutine open_arg(iarg, idev, ext, status, form, default_val)
  use sub_defs_io, only : lcl_getarg, lcl_iargc, max_arg_len, time_stamp
  use track_args
  implicit none
  integer, intent(in) :: iarg, idev
  character(*), optional, intent(in) :: ext, status, form, default_val
  !----------------------------------------------------------------------
  character(len=max_arg_len) :: arg
  integer       :: l, le
  
  l = lcl_iargc()
  if (iarg > l .or. iarg < 0) then
     if (present(default_val)) then
        arg = trim(default_val)
     else
        write(0,*) 'Number of args is',l
        write(0,*) 'Could not open file corresponding to arg',iarg
        write(0,*) 'because that arg does not exist'
        call time_stamp(0)
        write(0,*) 'stop'
        stop
     end if
  else
     call lcl_getarg(iarg, arg)
     call ta_RegisterUsed(iarg)
  end if
  
  l = index(arg,' ') - 1
  if (present(ext)) then
     le = index(ext,' ') - 1
     if (le <= 0) le = len(ext)
     arg(l+1:l+le) = ext(1:le)
     l = l + le
  end if
  !-- the following stupid sequence should not be necessary, but the
  ! standard doesn't specify open as having `optional args'
  if (present(status)) then
     if (present(form)) then
        open(unit=idev,file=arg(1:l),status=status,form=form)
     else
        open(unit=idev,file=arg(1:l),status=status)
     end if
  else
     if (present(form)) then
        open(unit=idev,file=arg(1:l),form=form)
     else
        open(unit=idev,file=arg(1:l))
     end if
  end if

end subroutine open_arg


!----------------------------------------------------------------------
!! It can often be useful to have a record of the command line -- this 
!! returns a string which conatins the command name and all the arguments
!----------------------------------------------------------------------
function command_line()
  use sub_defs_io, except => command_line
  implicit none
  character(len=max_line_len) ::  command_line
  !----------------------------------------------------------------------
  character(len=max_arg_len) :: string
  integer :: i, n

  n = lcl_iargc()
  call lcl_getarg(0,command_line)
  do i = 1, n
     call lcl_getarg(i,string)
     command_line = trim(command_line)//' '//trim(string)
  end do
end function command_line

!----------------------------------------------------------------------
!! Return the index of the argument which matches the string opt -- 
!! Intended to be used to get hold of the presence of an option, and/or
!! to locate an option so that one can extract values that follow it
!!
!! Trailing spaces are ignored.
!!
!! Returns a negative number if the option is not found.
!!
!! GPS 24/01/98 
!----------------------------------------------------------------------
function iargc_opt(opt)
  use sub_defs_io, except => iargc_opt
  use track_args
  implicit none
  integer                      :: iargc_opt
  character(len=*), intent(in) :: opt
  !------------------------------------------------------------
  integer :: i, n
  character(len=max_arg_len) :: string

  n = lcl_iargc()
  do i = 1, n
     call lcl_getarg(i,string)
     if (trim(string) == trim(opt)) exit
  end do
  ! Check to see if argument found
  if (i > n) then
     !-- make sure that iargc_opt + reasonable number remains negative
     iargc_opt = -100
  else
     iargc_opt = i
     call ta_RegisterUsed(iargc_opt)
  end if
end function iargc_opt



!----------------------------------------------------------------------
!! Return the string value corresponding to the argument which follows
!! the command line option opt
!!
!! GPS 8/11/95 (CCN8 23)
!----------------------------------------------------------------------
function string_val_opt(opt,default_val)
  use sub_defs_io, except => string_val_opt
  use track_args
  implicit none
  character(len=max_arg_len)             :: string_val_opt
  character(len=*), intent(in)           :: opt
  character(len=*), intent(in), optional :: default_val
  integer :: i,n
  i = iargc_opt(opt)
  n = lcl_iargc()
  if (i >= n .or. i < 0) then
     if (present(default_val)) then
        string_val_opt = default_val
     else
        write(0,*) 'String value for option ',trim(opt),' has been&
             & requested'
        write(0,*) 'but that option is not present, and no default value&
             & was provided'
        stop
     end if
  else
     call lcl_getarg(i+1,string_val_opt)
     call ta_RegisterUsed(i+1)
  end if
end function string_val_opt

!----------------------------------------------------------------------
!! Return the integer value corresponding to the argument which follows
!! the command line option opt
!!
!! GPS 8/11/95 (CCN8 23)
!----------------------------------------------------------------------
function int_val_opt(opt,default_val)
  use sub_defs_io, except => int_val_opt
  implicit none
  integer                      :: int_val_opt
  character(len=*), intent(in) :: opt
  integer, intent(in), optional :: default_val
  !----------------------------------------------------------------------
  integer :: i, n

  i = iargc_opt(opt)
  n = lcl_iargc()
  if (i >= n .or. i < 0) then
     if (present(default_val)) then
        int_val_opt = default_val
     else
        write(0,*) 'Numerical value for option ',trim(opt),' has been&
             & requested'
        write(0,*) 'but that option is not present, and no default value&
             & was provided'
        stop
     end if
  else
     int_val_opt = int_val_arg(i+1,default_val)
  end if
end function int_val_opt
!----------------------------------------------------------------------
!! Similar to int_val_opt, but for dble value
!----------------------------------------------------------------------
function dble_val_opt(opt,default_val)
  use sub_defs_io, except => dble_val_opt
  implicit none
  real(kind(1d0))                      :: dble_val_opt
  character(len=*), intent(in) :: opt
  real(kind(1d0)), intent(in), optional :: default_val
  !----------------------------------------------------------------------
  integer :: i, n

  i = iargc_opt(opt)
  n = lcl_iargc()
  if (i >= n .or. i < 0) then
     if (present(default_val)) then
        dble_val_opt = default_val
     else
        write(0,*) 'Numerical value for option ',trim(opt),' has been&
             & requested'
        write(0,*) 'but that option is not present, and no default value&
             & was provided'
        stop
     end if
  else
     dble_val_opt = dble_val_arg(i+1,default_val)
  end if
end function dble_val_opt
!----------------------------------------------------------------------
!! Say opt is -xxx: if -xxx is present on the command line, then
!! return true; if -noxxx is present on the command line, then return false.
!! If neither is present on the command line, then return default_val,
!! which therefore is not optional. If both arguments are present, then
!! return default val.
!!
!! GPS 8/11/95 (CCN8 23)
!----------------------------------------------------------------------
function log_val_opt(opt,default_val)
  use sub_defs_io, except => log_val_opt
  implicit none
  logical                       :: log_val_opt
  character(len=*), intent(in)  :: opt
  logical, intent(in), optional :: default_val
  !----------------------------------------------------------------------
  character(len=30) :: noopt
  integer :: i, j

  i = len_trim(opt)
  noopt = opt(1:1)//'no'//opt(2:i)

  i = iargc_opt(trim(opt))
  j = iargc_opt(trim(noopt))
  
  if (i > 0 .neqv. j > 0) then
     log_val_opt = (i>0)
  else
     if (present(default_val)) then
        log_val_opt = default_val
     else
        log_val_opt = .false.
     end if
  end if
end function log_val_opt

  
!----------------------------------------------------------------------
!! Return the integer value corresponding to the argk^th command line
!! argument. If # args < argk then use default if it is present, otherwise
!! output an error message and stop
!!
!! GPS 8/11/95 (CCN8 23)
!----------------------------------------------------------------------
integer function int_val_arg(argk, default_val)
  use track_args
  use sub_defs_io, except => int_val_arg
  implicit none
  integer, intent(in)           :: argk
  integer, intent(in), optional :: default_val
  !------------------------------------------------------------
  integer             :: i
  character(len=max_arg_len) string

  i = lcl_iargc()
  if (i < argk .or. argk < 0) then
     if (present(default_val)) then
        int_val_arg = default_val
        return
     else
        write(0,*) 'Numerical value of arg',argk,' has been requested'
        write(0,*) 'but only',i,' args present, and no default value&
             & was provided'
        stop
     end if
  end if

  call lcl_getarg(argk, string)
  call ta_RegisterUsed(argk)
  if (trim(string) == '-' .and. present(default_val)) then
     int_val_arg = default_val
  else
     int_val_arg = int_value(string)
  end if

end function int_val_arg

      
!======================================================================
!! Routine to convert string into an integer
!! Added ability to cope with negative numbers
!! GPS 30/11/94+13/11/96
!======================================================================

integer function int_value(string) result(value)
      implicit none
      character(*),  intent(in)  ::   string
      !----------------------------------------------------------------
      character(11), parameter   ::  template = '0123456789-'
      character  ::  current
      integer    ::  i, j, sgn
      logical    ::  number_started

      value          = 0
      sgn = 1
      number_started = .false.

loop: do i = 1, len(string)
         current = string(i:i)
         j = index(template, current)
         if (j == 11) then
            if (number_started) then
               write(0,*) string(1:len(string)), ' is an invalid number'
               write(0,*) 'Stop.'
               stop
            end if
            sgn = -1
            number_started = .true.
            cycle
         end if
         if (j /= 0) number_started = .true.
         if (number_started) then
            if (j == 0) exit
            value = value * 10 + (j-1)
         endif
      end do loop
      value = value * sgn
end function int_value

!----------------------------------------------------------------------
!! Return the dble value corresponding to the argk^th command line
!! argument. If # args < argk then use default if it is present, otherwise
!! output an error message and stop
!!
!! GPS 8/11/95 (CCN8 23)
!----------------------------------------------------------------------
real(kind(1d0)) function dble_val_arg(argk, default_val)
  use track_args
  use sub_defs_io, except => dble_val_arg
  implicit none
  integer,         intent(in)           :: argk
  real(kind(1d0)), intent(in), optional :: default_val
  !------------------------------------------------------------
  integer             :: i
  character(80) string

  i = lcl_iargc()
  if (i < argk .or. argk < 0) then
     if (present(default_val)) then
        dble_val_arg = default_val
        return
     else
        write(0,*) 'Numerical value of arg',argk,' has been requested'
        write(0,*) 'but only',i,' args present, and no default value&
             & was provided'
        stop
     end if
  end if

  call lcl_getarg(argk, string)
  call ta_RegisterUsed(argk)
  if (trim(string) == '-' .and. present(default_val)) then
     dble_val_arg = default_val
  else
     dble_val_arg = dble_value(string)
  end if
end function dble_val_arg


!---------------------------------------------------------------
!! returns the string corresponding to argument argk.
!! If that argument is absent, return the default_val (if provided)
function string_val_arg(argk, default_val)
  use track_args
  use sub_defs_io, except => string_val_arg
  implicit none
  character(len=max_arg_len)             :: string_val_arg
  integer,          intent(in)           :: argk
  character(len=*), intent(in), optional :: default_val
  integer             :: i

  i = lcl_iargc()
  if (i < argk .or. argk < 0) then
     if (present(default_val)) then
        string_val_arg = trim(default_val)
        return
     else
        write(0,*) 'String value of arg',argk,' has been requested'
        write(0,*) 'but only',i,' args present, and no default value&
             & was provided'
        stop
     end if
  end if
  
  call lcl_getarg(argk, string_val_arg)
  call ta_RegisterUsed(argk)
end function string_val_arg
     

!======================================================================
!! Routine to convert string into an float
!! Added ability to cope with negative numbers. Cannot cope with
!! scientific notation or other such fancy things.
!!
!! Actually: current version just does a read on the string...
!!
!! GPS 30/11/94+13/11/96+18/11/96
!======================================================================
real(kind(1d0)) function dble_value(string)
      implicit none
      character(*),  intent(in)  ::   string
      !----------------------------------------------------------------
      character(12), parameter   ::  template = '0123456789-.'
      character  ::  current
      integer    ::  i, j, sgn
      logical    ::  number_started, dpoint
      real(kind(1d0)) :: jfact, pfact

      read(string,*) dble_value
      return
      !-- default sign is positive
      sgn = 1

      number_started = .false.
      dpoint = .false.
      dble_value          = 0
      jfact = 1d0
      pfact = 10d0

loop: do i = 1, len(string)
         current = string(i:i)
         j = index(template, current)
         if (j == 11) then
            if (number_started) call invalid_number
             sgn = -1
            number_started = .true.
            cycle
         end if
         if (j /= 0) number_started = .true.
         if (number_started) then
            if (j == 0) exit
            if (j == 12) then 
               if (dpoint) call invalid_number
               dpoint = .true.
               pfact = 1d0
               cycle
            end if
            if (dpoint) jfact = jfact * 0.1d0
            dble_value = dble_value * pfact + (j-1) * jfact
         endif
      end do loop
      dble_value = dble_value * sgn
contains

  subroutine invalid_number
    write(0,*) string(1:len(string)), ' is an invalid number'
    write(0,*) 'Stop.'
    stop
  end subroutine invalid_number
  
end function dble_value

!----------------------------------------------------------------------
!! A set of routines for converting a number to a string, with optional
!! format specifier. 
function int2char(num,frmt)
  implicit none
  character(len=30) :: int2char
  integer, intent(in) :: num
  character(*), intent(in), optional :: frmt
  if (present(frmt)) then
     write(int2char,frmt) num
  else
     write(int2char,*) num
  end if
  int2char = adjustl(int2char)
end function int2char
function sp2char(num,frmt)
  implicit none
  character(len=30) :: sp2char
  real(kind(1.0)), intent(in) :: num
  character(*), intent(in), optional :: frmt
  if (present(frmt)) then
     write(sp2char,frmt) num
  else
     write(sp2char,*) num
  end if
  sp2char = adjustl(sp2char)
end function sp2char
function dp2char(num,frmt)
  implicit none
  character(len=30) :: dp2char
  real(kind(1d0)), intent(in) :: num
  character(*), intent(in), optional :: frmt
  if (present(frmt)) then
     write(dp2char,frmt) num
  else
     write(dp2char,*) num
  end if
  dp2char = adjustl(dp2char)
end function dp2char


!======================================================================
!! Returns a new device if any are available, or else causes an error
!! dump
!!
!! GPS 13/04/98
!======================================================================
function get_new_device() result(dev)
  use sub_defs_io, except => get_new_device
  implicit none
  integer :: dev
  !--------------------------------------------------------------
  integer, parameter :: dev_low = 30
  logical :: exist, opened
  integer :: iostat
  character(len=*), parameter :: subname='get_new_device'

  dev = dev_low
  do
     inquire(unit=dev,iostat=iostat, exist=exist, opened=opened)
     if (iostat /= 0) call error_report(subname, 'iostat non zero')
     if (.not. exist) call error_report(subname, 'available devices exhausted')
     if (.not. opened) exit
     dev = dev + 1 
  end do

end function get_new_device


!==========================================================================
!! Routine to write `time stamp' to a given output device (idev). The stamp
!! also contains the name of the program that generates it (taken from arg0)
!!
!! GPS 4/12/96
!==========================================================================
subroutine time_stamp(idev,string)
  use sub_defs_io, excepts => time_stamp
  implicit none
  integer,          intent(in),  optional :: idev
  character(len=*), intent(out), optional :: string
  !----------------------------------------------------------------------
  character(79) :: prog_name
  integer       :: time(8)
  
  call lcl_getarg(0,prog_name)
  call date_and_time(values=time)
  if (present(idev)) then
     !write(idev,15) prog_name(1:index(prog_name,' ')-1),&
     !     & time(3),time(2),time(1),time(5),time(6),time(7)
     write(idev,15) trim(prog_name),&
          & time(3),time(2),time(1),time(5),time(6),time(7)
  end if
  if (present(string)) then
     !write(string,15) prog_name(1:index(prog_name,' ')-1),&
     !     & time(3),time(2),time(1),time(5),time(6),time(7)
     write(string,15) trim(prog_name),&
          & time(3),time(2),time(1),time(5),time(6),time(7)
  end if
  
15 format('Stamped by ',a,' on ',i2,'/',i2.2,'/',i4.4,' at ',i2.2,':'&
        & ,i2.2 ,':',i2.2) 
  call lcl_flush(idev)
end subroutine time_stamp

!======================================================================
!! Standard error reporting routine. Writes error to
!! stderr and then does one of the following (specified by contents of
!! "action")
!!
!!  core
!!  stop
!!  cont
!!
!! Default action is currently defined below, but maybe its definition
!! should be moved to a globally accessible location.
!!
!! GPS 19/08/97
!======================================================================
subroutine error_report(subname,line1,line2,line3,line4,action)
  use sub_defs_io, excepts => error_report; implicit none
  character(*), intent(in) :: subname
  character(*), intent(in), optional :: line1, line2, line3, line4
  character(*), intent(in), optional :: action
  !--------------------------------------------------------------------
  character(*), parameter :: default_action = 'core'
  character(10) :: laction
  real :: a, b
  
  if (present(action)) then
     laction = trim(action)
  else
     laction = default_action
  end if

  write(0,*)
  write(0,'(2a)') 'Error report from routine: ',subname
  write(0,*)
  if (present(line1)) then
     write(0,'(a)') 'Reason is:'
     write(0,'(a)') line1
  end if
  if (present(line2)) write(0,'(a)') line2
  if (present(line3)) write(0,'(a)') line3
  if (present(line4)) write(0,'(a)') line4
  write(0,*)
  call time_stamp(0)
  !------ allow for more information on exit.
  if (trim(laction) == 'core') then
     a = 1d0; b = 1d0
     write(0,'(a)') 'Error report routine is about to dump core'
     write(0,*)
     write(0,*) 1d0/sqrt(a-b)
     !-- just in case that did not work, do an explicit stop
     stop
  else if (trim(laction) /= 'cont') then
     write(0,'(a)') 'Program stopped'
     stop
  else
     write(0,'(a)') 'Execution continuing'
  end if
end subroutine error_report
