!======================================================================
! This contains a few useful bits and pieces, some of which are from NR.
! The parts from NR (assertions) are in the public domain (see section
! C1.2 of "Numerical Recipes in Fortran 90" by Press et al.)
!
! $Id: assertions.f90,v 1.1 2001/06/27 13:40:16 gsalam Exp $
!
!======================================================================
module assertions
  use types
  implicit none
  private
  
  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eq5,assert_eqn
  END INTERFACE
  public :: assert_eq

  interface default_or_opt
     MODULE PROCEDURE default_or_opt_dp, default_or_opt_sp,&
          & default_or_opt_int, default_or_opt_log
  end interface
  public :: default_or_opt

contains

  function default_or_opt_dp(xdef,xopt) result(x)
    real(dp), intent(in)           :: xdef
    real(dp), intent(in), optional :: xopt
    real(dp)                       :: x
    if (present(xopt)) then
       x = xopt
    else
       x = xdef
    end if
  end function default_or_opt_dp
  function default_or_opt_sp(xdef,xopt) result(x)
    real(sp), intent(in)           :: xdef
    real(sp), intent(in), optional :: xopt
    real(sp)                       :: x
    if (present(xopt)) then
       x = xopt
    else
       x = xdef
    end if
  end function default_or_opt_sp
  function default_or_opt_int(xdef,xopt) result(x)
    integer, intent(in)           :: xdef
    integer, intent(in), optional :: xopt
    integer                       :: x
    if (present(xopt)) then
       x = xopt
    else
       x = xdef
    end if
  end function default_or_opt_int
  function default_or_opt_log(xdef,xopt) result(x)
    logical, intent(in)           :: xdef
    logical, intent(in), optional :: xopt
    logical                       :: x
    if (present(xopt)) then
       x = xopt
    else
       x = xdef
    end if
  end function default_or_opt_log
  



  FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (0,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
  END FUNCTION assert_eq2
  !BL
  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (0,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3
  !BL
  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (0,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq4
  FUNCTION assert_eq5(n1,n2,n3,n4,n5,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4,n5
    INTEGER :: assert_eq5
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4 .and. n4 == n5) then
       assert_eq5=n1
    else
       write (0,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq5
  !BL
  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (0,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn
  !BL

end module assertions
