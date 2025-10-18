module hoppet_version
  implicit none
  private
  ! DO NOT EDIT the next line -- it is modified automatically by scripts/set-version.sh
  character(len=*), parameter :: private_version_string = "2.0.1"

public :: HoppetVersion

contains  

  !! return a string with the current version of HOPPET
  function HoppetVersion() result(vs)
    character(len=len(private_version_string)) :: vs
    vs = private_version_string
  end function HoppetVersion

  !! Takes a pointer to an array of C characters and if the version
  !! string length is < maxlen-1, then fills the C-array with a null
  !! terminated string with the version, and returns 0. If the version
  !! string length is >= maxlen-1, then does not fill the C-array,
  !! and returns the required length (including the null terminator).
  !!
  !! C-signature is 
  !!  int HoppetVersionC(char *cchar, int maxlen)
  function hoppetVersionC(cchar, maxlen) result (errcode_or_maxlen) bind(C, name="hoppetVersionC")
    use, intrinsic :: iso_c_binding
    implicit none
    ! Arguments
    type(c_ptr),    value        :: cchar          !< pointer to char buffer from C/C++
    integer(c_int), value        :: maxlen         !< max buffer length provided by caller
    integer(c_int)               :: errcode_or_maxlen
    ! Local variables
    integer :: verlen
    character(kind=c_char), pointer :: fbuf(:)

    ! Actual length of version string
    verlen = len_trim(private_version_string)

    if (verlen < maxlen) then
      ! Associate the C pointer with a Fortran char buffer of length maxlen
      call c_f_pointer(cchar, fbuf, [maxlen])

      fbuf(1:verlen) = transfer(private_version_string(1:verlen), fbuf(1:verlen))

      ! Null terminate
      fbuf(verlen+1) = c_null_char

      ! Return success (0)
      errcode_or_maxlen = 0
    else
       ! Not enough room: return required length including null terminator
       errcode_or_maxlen = verlen + 1
    end if
  end function HoppetVersionC

end module hoppet_version


subroutine HoppetWelcomeMessage
  use hoppet_version
  use hoppet_term
  implicit none
  character(len=:), allocatable :: welcome_str

  welcome_str = 'Welcome to HOPPET v. '//HoppetVersion()

  write(0,'(a)') '-----------------------------------------------------------'//bold
  !write(0,'(a)') '              Welcome to HOPPET v. '//trim(HoppetVersionString())
  write(0,'(a)') repeat(' ', (59-len(welcome_str))/2)//welcome_str
  write(0,'(a)') '    Higher Order Perturbative Parton Evolution Toolkit     '//reset
  write(0,'(a)') '-----------------------------------------------------------'
  !write(0,'(a)') ''
  write(0,'(a)') '                  Written (2001-2025) by                   '
  write(0,'(a)') '     Frederic Dreyer, Alexander Karlberg, Paolo Nason,     '
  write(0,'(a)') '      Juan Rojo, Gavin P. Salam and Giulia Zanderighi      '
  write(0,'(a)') ''
  write(0,'(a)') ' It is made available under the GNU public license.        '
  write(0,'(a)') ' If using it for scientific work, please cite              '
  write(0,'(a)') ' G.P. Salam & J. Rojo, CPC 180(2009)120 (arXiv:0804.3755). '
  write(0,'(a)') ' and A. Karlberg, P. Nason, G.P. Salam, G. Zanderighi      '
  write(0,'(a)') ' & F. Dreyer (arXiv:2510.09310).                           '
  write(0,'(a)') ' '
  write(0,'(a)') ' You are also encouraged to cite the original references,  '
  write(0,'(a)') ' for LO, NLO, NNLO, and N3LO splitting functions, the QCD  '
  write(0,'(a)') ' 1, 2, 3, and 4 loop beta functions and the PDF and        '
  write(0,'(a)') ' coupling mass threshold matching functions. Additionally  '
  write(0,'(a)') ' the DIS coefficient functions should be cited when used.  '
  write(0,'(a)') '-----------------------------------------------------------'
end subroutine HoppetWelcomeMessage
  
