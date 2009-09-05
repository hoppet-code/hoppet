! keep this outside a module so that we can recompile (e.g. for
! changing version number) without modifying any dependences.
!
subroutine HoppetWelcomeMessage
  write(0,'(a)') '-----------------------------------------------------------'
  write(0,'(a)') '             Welcome to HOPPET v. 1.1.2                    '
  write(0,'(a)') '   Higher Order Perturbative Parton Evolution Toolkit      '
  write(0,'(a)') ''
  write(0,'(a)') '        Written by Gavin P. Salam (2001-2008)'
  write(0,'(a)') '          with contributions from Juan Rojo'
  write(0,'(a)') ''
  write(0,'(a)') ' It is made available under the GNU public license,'
  write(0,'(a)') ' with the additional request that if you use it or any'
  write(0,'(a)') ' derivative of it in scientific work then you should cite:'
  write(0,'(a)') ' G.P. Salam and J. Rojo, arXiv:0804.3755.'
  write(0,'(a)') ' '
  write(0,'(a)') ' You are also encouraged to cite the original references,'
  write(0,'(a)') ' for LO, NLO and NNLO splitting functions, the QCD'
  write(0,'(a)') ' 1, 2 and 3 loop beta functions and the coupling and '
  write(0,'(a)') ' PDF and coupling mass threshold matching functions.'
  write(0,'(a)') '-----------------------------------------------------------'
end subroutine HoppetWelcomeMessage
  
