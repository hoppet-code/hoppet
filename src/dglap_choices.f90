!!
!! Parameters describing the range of possible "options" in the
!! DGLAP evolution. Options refer both to physical choices (such 
!! as schemes, polarization) as well as technical choices (choice of 
!! approximation used for NNLO splitting function).
!!
!! GPS: adapted from qcd.f90 18/09/2004
!!
module dglap_choices
  implicit none
  private

  integer, parameter, public :: nnlo_splitting_exact = -2
  integer, parameter, public :: nnlo_splitting_param = -1
  ! these three should keep their numerical values because
  ! of a correspondence with vogt imod values
  integer, parameter, public :: nnlo_splitting_Nfitav   =  0
  integer, parameter, public :: nnlo_splitting_Nfiterr1 =  1
  integer, parameter, public :: nnlo_splitting_Nfiterr2 =  2
  integer, public :: nnlo_splitting_variant = nnlo_splitting_param

  integer, parameter, public :: n3lo_splitting_exact = -2
  integer, parameter, public :: n3lo_splitting_param = -1
  ! these three should keep their numerical values because
  ! of a correspondence with vogt imod values
  integer, parameter, public :: n3lo_splitting_Nfitav   =  0
  integer, parameter, public :: n3lo_splitting_Nfiterr1 =  1
  integer, parameter, public :: n3lo_splitting_Nfiterr2 =  2
  integer, public :: n3lo_splitting_variant = n3lo_splitting_Nfitav
  ! As of 2024-04-16 there are several approximations available of the
  ! n3lo splitting functions. To maintain some backwards compatibility
  ! we have a switch below that allows the user to pick which set of
  ! approximations to use. The have progressively more moments.
  integer, parameter, public :: n3lo_splitting_approximation_up_to_2310_05744 = 100 !< Uses non-singlet of 1610.07477+1707.08315,
                                                                                    !< pure-singlet (qq) of 2302.07593,
                                                                                    !< qg of 2307.04158 and gq and gg of 2310.05744
  integer, parameter, public :: n3lo_splitting_approximation_up_to_2404_09701 = 101 !< Replaces gq with that of 2404.09701
  integer, parameter, public :: n3lo_splitting_approximation_up_to_2410_08089 = 102 !< Additionally replaces gg with that of 2410.08089
  integer, public :: n3lo_splitting_approximation =  n3lo_splitting_approximation_up_to_2410_08089 

  integer, parameter, public :: nnlo_nfthreshold_exact = -12
  integer, parameter, public :: nnlo_nfthreshold_param = -11
  integer, public :: nnlo_nfthreshold_variant = nnlo_nfthreshold_param

  integer, parameter, public :: n3lo_nfthreshold_on  = 1
  integer, parameter, public :: n3lo_nfthreshold_off = 0
  integer, public :: n3lo_nfthreshold = n3lo_nfthreshold_on

  integer, parameter, public :: factscheme_MSbar    = 1
  integer, parameter, public :: factscheme_DIS      = 2
  integer, parameter, public :: factscheme_PolMSbar = 3
  integer, parameter, public :: factscheme_FragMSbar = 4
  ! have these on split lines so that they are not caught
  ! by naming routines.
  integer, parameter,&
       & public :: factscheme_default = factscheme_MSbar
  integer, parameter, &
       &public :: factscheme_Poldefault = factscheme_PolMSbar

  public :: dglap_Set_nnlo_splitting
  public :: dglap_Set_n3lo_splitting
  public :: dglap_Set_n3lo_splitting_approximation
  public :: dglap_Set_nnlo_nfthreshold
  public :: dglap_Set_n3lo_nfthreshold
contains

  !-------- overkill ----------------------------------------
  subroutine dglap_Set_nnlo_splitting(variant)
    integer, intent(in) :: variant
    nnlo_splitting_variant = variant
  end subroutine dglap_Set_nnlo_splitting

  !-------- overkill ----------------------------------------
  subroutine dglap_Set_n3lo_splitting(variant)
    integer, intent(in) :: variant
    n3lo_splitting_variant = variant
  end subroutine dglap_Set_n3lo_splitting

  !-------- overkill ----------------------------------------
  subroutine dglap_Set_n3lo_splitting_approximation(variant)
    integer, intent(in) :: variant
    n3lo_splitting_approximation = variant
  end subroutine dglap_Set_n3lo_splitting_approximation

  !-------- overkill ----------------------------------------
  subroutine dglap_Set_nnlo_nfthreshold(variant)
    integer, intent(in) :: variant
    nnlo_nfthreshold_variant = variant
  end subroutine dglap_Set_nnlo_nfthreshold

  !-------- overkill ----------------------------------------
  subroutine dglap_Set_n3lo_nfthreshold(variant)
    integer, intent(in) :: variant
    n3lo_nfthreshold = variant
  end subroutine dglap_Set_n3lo_nfthreshold
end module dglap_choices
