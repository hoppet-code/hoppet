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

  integer, parameter, public :: nnlo_nfthreshold_exact = -12
  integer, parameter, public :: nnlo_nfthreshold_param = -11
  integer, public :: nnlo_nfthreshold_variant = nnlo_nfthreshold_param

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
  public :: dglap_Set_nnlo_nfthreshold
contains

  !-------- overkill ----------------------------------------
  subroutine dglap_Set_nnlo_splitting(variant)
    integer, intent(in) :: variant
    nnlo_splitting_variant = variant
  end subroutine dglap_Set_nnlo_splitting

  !-------- overkill ----------------------------------------
  subroutine dglap_Set_nnlo_nfthreshold(variant)
    integer, intent(in) :: variant
    nnlo_nfthreshold_variant = variant
  end subroutine dglap_Set_nnlo_nfthreshold
end module dglap_choices
