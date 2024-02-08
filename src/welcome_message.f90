! keep this outside a module so that we can recompile (e.g. for
! changing version number) without modifying any dependences.
!
subroutine HoppetWelcomeMessage
  write(0,'(a)') '-----------------------------------------------------------'
  write(0,'(a)') '               Welcome to HOPPET v. 1.3.0-devel            '
  write(0,'(a)') '   Higher Order Perturbative Parton Evolution Toolkit      '
  write(0,'(a)') ''
  write(0,'(a)') '                   Written (2001-2023) by                  '
  write(0,'(a)') '     Frederic Dreyer, Alexander Karlberg, Paolo Nason,     '
  write(0,'(a)') '      Juan Rojo, Gavin P. Salam and Giulia Zanderighi      '
  write(0,'(a)') ''
  write(0,'(a)') ' It is made available under the GNU public license,        '
  write(0,'(a)') ' with the additional request that if you use it or any     '
  write(0,'(a)') ' derivative of it in scientific work then you should cite: '
  write(0,'(a)') ' G.P. Salam & J. Rojo, CPC 180(2009)120 (arXiv:0804.3755). '
  write(0,'(a)') ' and                                                       '
  write(0,'(a)') ' A. Karlberg, P. Nason, G.P. Salam, G. Zanderighi          '
  write(0,'(a)') ' & F. Dreyer (arXiv:2401.XXXXX).                           '
  write(0,'(a)') ' '
  write(0,'(a)') ' You are also encouraged to cite the original references,  '
  write(0,'(a)') ' for LO, NLO and NNLO splitting functions, the QCD         '
  write(0,'(a)') ' 1, 2, 3, and 4 loop beta functions and the PDF and        '
  write(0,'(a)') ' coupling mass threshold matching functions. Additionally  '
  write(0,'(a)') ' the DIS coefficient functions should be cited when used.  '
  write(0,'(a)') '-----------------------------------------------------------'
end subroutine HoppetWelcomeMessage
  
