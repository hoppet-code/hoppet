module hoppet_version
  implicit none
  private
  !! when updating the version number, remember to do this also in 
  !! CMakeLists.txt and in the documentation.
  character(len=*), parameter :: private_version_string = "2.0.0-devel"

  public :: HoppetVersion

contains  

  !! return a string with the current version of HOPPET
  function HoppetVersion() result(vs)
    character(len=len(private_version_string)) :: vs
    vs = private_version_string
  end function HoppetVersion

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
  write(0,'(a)') '                   Written (2001-2025) by                  '
  write(0,'(a)') '     Frederic Dreyer, Alexander Karlberg, Paolo Nason,     '
  write(0,'(a)') '      Juan Rojo, Gavin P. Salam and Giulia Zanderighi      '
  write(0,'(a)') ''
  write(0,'(a)') ' It is made available under the GNU public license.        '
  write(0,'(a)') ' If using it for scientific work, please cite              '
  write(0,'(a)') ' G.P. Salam & J. Rojo, CPC 180(2009)120 (arXiv:0804.3755). '
  write(0,'(a)') ' and A. Karlberg, P. Nason, G.P. Salam, G. Zanderighi      '
  write(0,'(a)') ' & F. Dreyer (arXiv:2510.XXXXX).                           '
  write(0,'(a)') ' '
  write(0,'(a)') ' You are also encouraged to cite the original references,  '
  write(0,'(a)') ' for LO, NLO, NNLO, and N3LO splitting functions, the QCD  '
  write(0,'(a)') ' 1, 2, 3, and 4 loop beta functions and the PDF and        '
  write(0,'(a)') ' coupling mass threshold matching functions. Additionally  '
  write(0,'(a)') ' the DIS coefficient functions should be cited when used.  '
  write(0,'(a)') '-----------------------------------------------------------'
end subroutine HoppetWelcomeMessage
  
