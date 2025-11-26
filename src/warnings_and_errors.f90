!======================================================================
!!
!! Routine for dealing with warnings that may crop up many times and 
!! which should be output only a limited number of times
!!
!! In f95 could have used a special type which was preinitialised...
!!
!! $Id: warnings_and_errors.f90,v 1.6 2004/02/26 19:02:19 salam Exp $
!======================================================================
module warnings_and_errors
  use types
  implicit none
  private

  integer, parameter :: base = 10000
  integer            :: n_warn_sources   = 0
  integer, parameter, public :: warn_id_INIT=0
  integer, parameter, public :: default_max_warn = 5

  interface wae_warn
     module procedure wae_warn_new, wae_warn_old
  end interface
  
  public :: wae_warn, wae_error, wae_setunit
  public :: wae_enable_cxx_exceptions

  integer, parameter :: stddev_in = 0
  integer            :: stddev = stddev_in

  logical, public, save :: throw_cxx_exceptions = .false.

  interface 
    subroutine hoppet_throw_runtime_error(cstr_ptr) bind(C,name="hoppet_throw_runtime_error")
      use iso_c_binding
      type(c_ptr), value :: cstr_ptr
    end subroutine hoppet_throw_runtime_error
  end interface
  public :: hoppet_throw_runtime_error



contains

  !! On some systems (linux, MacOS at the time of writing, 2025), 
  !! C++ exceptions can propagate through Fortran code correctly, and
  !! this flag, if set to true, will cause a C++ exception to be
  !! thrown instead of the Fortran error stop being invoked.
  !!
  !! Beware, this is not generally portable, nor guaranteed to be future proof.
  !! Use at your own risk.
  subroutine wae_enable_cxx_exceptions() bind(C, name='hoppetEnableCxxExceptionsFromFortran')
    use iso_c_binding
    use hoppet_term
    print '(a)', bg_white//red//bold//'*** Enabled C++ exceptions from hoppet Fortran; use with care: not portable ***'//reset
    throw_cxx_exceptions = .true.
  end subroutine wae_enable_cxx_exceptions
  !---------------------------------------------------------------------
  !! Routine to allow the output of a warning up to some maximum number of 
  !! times after which no such warning will be issues.
  !!
  !! On the first call for a given kind of warning, the warn_n should
  !! (must have save attribute in calling routine) should be the 
  !! maximum number of warnings that will be output
  !!
  subroutine wae_warn_new(warn_n, text, text2, text3, text4, intval, dbleval)
    integer, intent(inout)       :: warn_n
    character(len=*), intent(in) :: text
    character(len=*), intent(in), optional :: text2
    character(len=*), intent(in), optional :: text3
    character(len=*), intent(in), optional :: text4
    integer,          intent(in), optional :: intval
    real(kind(1d0)),  intent(in), optional :: dbleval
    !--------------------------------------

    if (warn_n > 0) then
       warn_n = warn_n - 1
       write(stddev,'(a)', advance='no') 'WARNING in '
       write(stddev,'(a)')  text
       if (present(text2)) write(stddev,'(a)') text2
       if (present(text3)) write(stddev,'(a)') text3
       if (present(text4)) write(stddev,'(a)') text4
       if (present(intval))  write(stddev,*) intval
       if (present(dbleval)) write(stddev,*) dbleval
       if (warn_n == 0) write(stddev,'(a)') &
            &'----- No more such warnings will be issued ------'
    end if
  end subroutine wae_warn_new


  !---------------------------------------------------------------------
  !! Routine to allow the output of a warning up to some maximum number of 
  !! times after which no such warning will be issues.
  !!
  !! On the first call for a given kind of warning, the warn_id should be 
  !! equal to warn_id_INIT; warn_id should have the SAVE attribute to
  !! allow wae_warn to keep track of the number of warnings of this type
  !!
  subroutine wae_warn_old(max_warn, warn_id, &
       &text, text2, text3, intval, dbleval)
    integer, intent(in)          :: max_warn
    integer, intent(inout)       :: warn_id
    character(len=*), intent(in) :: text
    character(len=*), intent(in), optional :: text2
    character(len=*), intent(in), optional :: text3
    integer,          intent(in), optional :: intval
    real(kind(1d0)),  intent(in), optional :: dbleval
    !--------------------------------------
    integer :: warn_index, nwarn

    !-- generate a new warn_id
    if (warn_id <= 0) then
       n_warn_sources = n_warn_sources + 1
       warn_id = n_warn_sources * base
    end if
    
    warn_index = warn_id / base
    nwarn      = warn_id - warn_index*base

    if (nwarn < max_warn) then
       if (max_warn > base-2) call wae_error('wae_warn',&
            & 'max_warn exceeded maximum allowed value; message was', text)
       !-- does this make any sense at all (GPS 8/1/03)?
       if (warn_id > huge(warn_id)) call wae_error('wae_warn',&
            & 'exceeded max capicity for distinct warnings; message was', text)
    
       warn_id = warn_id + 1
       write(stddev,'(a)', advance='no') 'WARNING in '
       write(stddev,'(a)')  text
       if (present(text2)) write(stddev,'(a)') text2
       if (present(text3)) write(stddev,'(a)') text3
       if (present(intval))  write(stddev,*) intval
       if (present(dbleval)) write(stddev,*) dbleval
       !-- if there is only an 1 warning to be written then
       !   avoid cluttering screen with too many messages
       if (nwarn == max_warn - 1 .and. max_warn > 1) write(stddev,'(a)') &
            &'----- No more such warnings will be issued ------'
    end if
  end subroutine wae_warn_old

  !======================================================================
  !! Report an error and then crash (by attempting floating point exception)
  !======================================================================
  subroutine wae_error(text1, text2, text3, text4, intval, dbleval)
    use hoppet_to_string
    use hoppet_term
    use hoppet_c_f_string_utils
    character(len=*), intent(in) :: text1
    character(len=*), intent(in), optional :: text2
    character(len=*), intent(in), optional :: text3
    character(len=*), intent(in), optional :: text4
    integer,          intent(in), optional :: intval
    real(kind(1d0)),  intent(in), optional :: dbleval
    !real :: a,b
    character(len=:), allocatable :: msg
    type(c_ptr) :: msg_cstr_ptr

    if (.not. throw_cxx_exceptions) then
      write(stddev,*)
      write(stddev,'(a)') '============================================================='
      write(stddev,'(a)', advance='no') red//bold//'FATAL ERROR in '//reset
    else
      msg = 'FATAL HOPPET ERROR in '
    end if

    msg = msg//trim(text1)
    if (present(text2))   msg = msg//' '//trim(text2)
    if (present(text3))   msg = msg//' '//trim(text3)
    if (present(text4))   msg = msg//' '//trim(text4)
    if (present(intval))  msg = msg//' (integer value: '//trim(to_string_int(intval))//')'
    if (present(dbleval)) msg = msg//' (double value: '//trim(to_string_dp(dbleval))//')'

    if (throw_cxx_exceptions) then
      msg_cstr_ptr = cstr_from_fortran_string(msg)
      !! deallocate the Fortran string after conversion because
      !! the C++ exceptions will not look after Fortran memory
      deallocate(msg)
      call hoppet_throw_runtime_error(msg_cstr_ptr)
    else
      write(stddev,'(a)') msg
      write(stddev,*)
      error stop 1
    end if

  end subroutine wae_error
  
  !!
  !! Set the unit for all output by the warning and error routines
  !!
  subroutine wae_setunit(unit)
    integer, intent(in) :: unit
    stddev = unit
  end subroutine wae_setunit
  
end module warnings_and_errors

