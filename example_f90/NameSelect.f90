
!!
!! Module generated automatically by NameSelect.pl
!! Provides conversion from string codes to integer codes
!!
module NameSelect
  use dglap_choices
  
  implicit none
  private

  public :: CodeOfName
  public :: NameOfCode
  public :: code_val_opt

  integer, parameter :: maxlen_name = 40
  integer, parameter :: maxlen_longname = 120
  integer, parameter, public :: NameSelect_maxlen_name = maxlen_name
  integer, parameter, public :: NameSelect_maxlen_longname = maxlen_longname

contains

  !!
  !! Returns the code associated with 'name'; if status dummy arg
  !! is present a non-zero result indicates failure to find the name.
  !! Otherwise failure results in a hard error.
  !!
  function CodeOfName(name, status) result(code)
    character(len=*), intent(in)   :: name
    integer                        :: code
    integer, optional, intent(out) :: status

    if (present(status)) status = 0

    ! now follows code generated specifically for above used modules
    select case(trim(name))
    case('nnlo_splitting_exact')
      code = nnlo_splitting_exact
    case('nnlo_splitting_param')
      code = nnlo_splitting_param
    case('nnlo_splitting_Nfitav')
      code = nnlo_splitting_Nfitav
    case('nnlo_splitting_Nfiterr1')
      code = nnlo_splitting_Nfiterr1
    case('nnlo_splitting_Nfiterr2')
      code = nnlo_splitting_Nfiterr2
    case('nnlo_nfthreshold_exact')
      code = nnlo_nfthreshold_exact
    case('nnlo_nfthreshold_param')
      code = nnlo_nfthreshold_param
    case('factscheme_MSbar')
      code = factscheme_MSbar
    case('factscheme_DIS')
      code = factscheme_DIS
    case('factscheme_PolMSbar')
      code = factscheme_PolMSbar
    case default
      if (present(status)) then
        status = -1
      else
        write(0,*) 'ERROR. CodeOfName: unrecognized name "'//name//'"'
        stop
      endif
    end select
  end function CodeOfName


  ! standard code (independent of used modules)
  !!
  !! looks for the command-line argument 'option' and if it is
  !! present its value (optionally prefixed with 'prefix') is
  !! fed to CodeOfName
  !!
  function code_val_opt(option, default, prefix) result(code)
    use sub_defs_io
    character(len=*),           intent(in) :: option
    integer, optional,          intent(in) :: default
    character(len=*), optional, intent(in) :: prefix
    integer                                :: code
    !---------------------------------------
    character(len=maxlen_name) :: opt_val

    if (log_val_opt(option)) then
       opt_val = string_val_opt(option)
       if (present(prefix)) then
         code = CodeOfName(prefix//trim(opt_val))
       else
         code = CodeOfName(trim(opt_val))
       end if
    else
       if (present(default)) then
          code = default
       else
          write(0,*) 'Error in code_val_opt: command-line option '&
               &//option//' absent and no default provided'
       end if
    end if
  end function code_val_opt

  !!
  !! Returns the name of the given integer code which has
  !! the (optional) specified prefix
  !!
  function NameOfCode(code,prefix,longname) result(name)
    integer,          intent(in)            :: code
    character(len=*), intent(in),  optional :: prefix
    character(len=*), intent(out), optional :: longname
    character(len=maxlen_name)              :: name
    !----------------------------------------------
    integer :: nocc

    nocc = 0
    if (present(longname)) longname = ''

    ! code specific to used modules starts here
    if (PrefixMatches('nnlo_splitting_exact',prefix) .and. code == nnlo_splitting_exact) then
       nocc = nocc + 1; name = 'nnlo_splitting_exact'
    end if

    if (PrefixMatches('nnlo_splitting_param',prefix) .and. code == nnlo_splitting_param) then
       nocc = nocc + 1; name = 'nnlo_splitting_param'
    end if

    if (PrefixMatches('nnlo_splitting_Nfitav',prefix) .and. code == nnlo_splitting_Nfitav) then
       nocc = nocc + 1; name = 'nnlo_splitting_Nfitav'
    end if

    if (PrefixMatches('nnlo_splitting_Nfiterr1',prefix) .and. code == nnlo_splitting_Nfiterr1) then
       nocc = nocc + 1; name = 'nnlo_splitting_Nfiterr1'
    end if

    if (PrefixMatches('nnlo_splitting_Nfiterr2',prefix) .and. code == nnlo_splitting_Nfiterr2) then
       nocc = nocc + 1; name = 'nnlo_splitting_Nfiterr2'
    end if

    if (PrefixMatches('nnlo_nfthreshold_exact',prefix) .and. code == nnlo_nfthreshold_exact) then
       nocc = nocc + 1; name = 'nnlo_nfthreshold_exact'
    end if

    if (PrefixMatches('nnlo_nfthreshold_param',prefix) .and. code == nnlo_nfthreshold_param) then
       nocc = nocc + 1; name = 'nnlo_nfthreshold_param'
    end if

    if (PrefixMatches('factscheme_MSbar',prefix) .and. code == factscheme_MSbar) then
       nocc = nocc + 1; name = 'factscheme_MSbar'
    end if

    if (PrefixMatches('factscheme_DIS',prefix) .and. code == factscheme_DIS) then
       nocc = nocc + 1; name = 'factscheme_DIS'
    end if

    if (PrefixMatches('factscheme_PolMSbar',prefix) .and. code == factscheme_PolMSbar) then
       nocc = nocc + 1; name = 'factscheme_PolMSbar'
    end if

    ! common-code resumes
    if (nocc == 0) then
       if (present(prefix)) then
          write(0,*) 'Error in NameOfCode: could not find code',&
               &code,' with prefix "'//prefix//'"'
       else
          write(0,*) 'Error in NameOfCode: could not find code',&
               &code,' (without prefix)'
       end if
       stop
    end if

    if (nocc > 1) then
       if (present(prefix)) then
          write(0,*) 'Error in NameOfCode: several meanings for code',&
               &code,' with prefix "'//prefix//'"'
       else
          write(0,*) 'Error in NameOfCode: several meanings for code',&
               &code,' (without prefix)'
       end if
       stop
    end if

  end function NameOfCode



  !!
  !! for establishing whether a prefix matches a given string
  !!
  logical function PrefixMatches(string,prefix)
    character(len=*), intent(in)           :: string
    character(len=*), intent(in), optional :: prefix

    if (present(prefix)) then
       PrefixMatches = (index(string,prefix) == 1)
    else
       PrefixMatches = .true.
    end if
  end function PrefixMatches


end module NameSelect
