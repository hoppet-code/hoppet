!----------------------------------------------------------------------
!
! This routine stores info about numbers of components used, 
! labelling of those components, and more importantly to a general
! used it contains routines for converting between the evolution storage
! format and the "human" storage format
!
!
! $Id: pdf_representation.f90,v 1.5 2003/11/22 00:45:39 salam Exp $
!----------------------------------------------------------------------
module pdf_representation
  use types; use consts_dp; use assertions; use warnings_and_errors
  implicit none
  !----------------------------------------------------------------------
  ! As of 6/01/2003, ncompmax and iflv_max must be the same (actually
  ! not anymore since we added extra representation iflv_info)
  integer, parameter, public :: iflv_max=6
  integer, parameter, public :: iflv_min=-iflv_max
  integer, parameter, public :: ncomponents=iflv_max
  !------------------------------------------------------
  integer, parameter, public :: iflv_g=0, iflv_sigma=1, iflv_V=-1
  integer, parameter, public :: iflv_d = 1, iflv_u = 2, iflv_s = 3, iflv_c = 4
  integer, parameter, public :: iflv_b = 5, iflv_t = 6
  integer, parameter, public :: iflv_dbar = -iflv_d
  integer, parameter, public :: iflv_ubar = -iflv_u
  integer, parameter, public :: iflv_sbar = -iflv_s
  integer, parameter, public :: iflv_cbar = -iflv_c
  integer, parameter, public :: iflv_bbar = -iflv_b
  integer, parameter, public :: iflv_tbar = -iflv_t
  integer, parameter, public :: iflv_info =  iflv_max + 1
  integer, parameter, public :: ncompmin = iflv_min
  integer, parameter, public :: ncompmax = iflv_info

  integer, parameter, public :: pdfr_Human = -1000001
  !integer, parameter, public :: pdfr_Evln  = 1

  integer, parameter :: default_ibase = 1


  type pdf_rep
     integer :: nf, ibase
  end type pdf_rep
  
  public :: pdf_rep

  interface CopyHumanPdfToEvln
     module procedure pdfr_HumanToEvln_sc, pdfr_HumanToEvln_1d
     module procedure pdfr_HumanToEvln_nf_sc, pdfr_HumanToEvln_nf_1d
  end interface
  
  interface CopyEvlnPdfToHuman
     module procedure pdfr_EvlnToHuman_sc, pdfr_EvlnToHuman_1d
     module procedure pdfr_EvlnToHuman_nf_sc, pdfr_EvlnToHuman_nf_1d
  end interface

  public :: CopyHumanPdfToEvln, CopyEvlnPdfToHuman
  public :: GetPdfRep, LabelPdfAsRep
  public :: LabelPdfAsHuman
  public :: DefaultEvlnRep

contains

  !======================================================================
  !! return the pdf_rep object corresponding to the default Evln
  !! representation for the specified value of nf.
  function DefaultEvlnRep(nf_lcl) 
    integer, intent(in) :: nf_lcl
    type(pdf_rep)       :: DefaultEvlnRep
    DefaultEvlnRep%nf    = nf_lcl
    DefaultEvlnRep%ibase = default_ibase
  end function DefaultEvlnRep


  
  !-------------------------------------------------------------------
  !! Take a "human format" pdf set and convert it to a format
  !! suitable for evolution, as described near the beginning of
  !! dglap_objects.
  !!
  !! What is called "k" there, here is known as prep%ibase
  !!
  !! Human format is that produced by the pdf_general module.
  !! i.e. (tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t)
  pure subroutine pdfr_HumanToEvln_sc(prep, qh, qe)
    type(pdf_rep), intent(in)  :: prep
    real(dp),      intent(in)  :: qh(ncompmin:)
    real(dp),      intent(out) :: qe(ncompmin:)
    !----------------------------------------------
    real(dp) :: tmp
    real(dp) :: qplus_base, qminus_base
    integer  :: i, j

    qe(iflv_g) = qh(iflv_g)
    qe(iflv_sigma) = sum(qh(1:prep%nf))
    tmp          = sum(qh(-prep%nf:-1))
    qe(iflv_V)     = qe(iflv_sigma) - tmp
    qe(iflv_sigma) = qe(iflv_sigma) + tmp
    
    qplus_base  = qh(prep%ibase) + qh(-prep%ibase)
    qminus_base = qh(prep%ibase) - qh(-prep%ibase)
    do i = 2, prep%nf
       if (i > prep%ibase) then
          j = i
       else
          j = i-1
       end if
       qe( i) = qh(j) + qh(-j) - qplus_base
       qe(-i) = qh(j) - qh(-j) - qminus_base
    end do
    !-- keep the rest clean...
    do i = prep%nf+1, ncomponents
       !qe(i) = zero
       !qe(-i) = zero
       qe(i) = qh(i)
       qe(-i) = qh(-i)
    end do
    
  end subroutine pdfr_HumanToEvln_sc

  !------------------------------------------------
  !! a vector version of the above
  subroutine pdfr_HumanToEvln_1d(prep, qh, qe)
    type(pdf_rep), intent(in)  :: prep
    real(dp),      intent(in)  :: qh(:,ncompmin:)
    real(dp),      intent(out) :: qe(:,ncompmin:)
    integer :: n, i

    n = assert_eq(size(qh,dim=1),size(qe,dim=1),'pdfr_HumanToEvln_1d')
    if (GetPdfRep(qh) /= pdfr_Human) call wae_error('pdfr_HumanToEvln_1d',&
         &'qh is not in "Human" format')

    do i = 1, n
       call pdfr_HumanToEvln_sc(prep, qh(i,:), qe(i,:))
    end do
    !call LabelPdfAsRep(qe,pdfr_Evln)
    call LabelPdfAsRep(qe,prep%nf)
  end subroutine pdfr_HumanToEvln_1d
  

  !---------------------------------------------------------
  !! Take a pdf set in "evolution format" and convert it back
  !! to "human format".
  subroutine pdfr_EvlnToHuman_sc(prep, qe, qh)
    type(pdf_rep), intent(in)  :: prep
    real(dp),      intent(in)  :: qe(ncompmin:)
    real(dp),      intent(out) :: qh(ncompmin:)
    !---------------------------------------------
    integer  :: i, j
    real(dp) :: tmp
    
    qh(iflv_g) = qe(iflv_g)
    !-- first construct individual + and - distributions
    qh( prep%ibase) = (qe(iflv_sigma) - sum(qe(2:prep%nf))   )/prep%nf
    qh(-prep%ibase) = (qe(iflv_V)     - sum(qe(-prep%nf:-2)) )/prep%nf

    do i = 2, prep%nf
       if (i > prep%ibase) then
          j = i
       else
          j = i-1
       end if
       qh( j) = qe( i) + qh( prep%ibase)
       qh(-j) = qe(-i) + qh(-prep%ibase)
    end do

    !-- then go from + and - to q and qbar
    do i = 1, prep%nf
       tmp = qh(-i)
       qh(-i) = half*(qh(i) - tmp)
       qh( i) = half*(qh(i) + tmp)
    end do
    
    !-- make sure the rest is zero
    do i = prep%nf+1, iflv_max
       !qh( i) = zero
       !qh(-i) = zero
       qh( i) = qe(i)
       qh(-i) = qe(-i)
    end do
    
  end subroutine pdfr_EvlnToHuman_sc
  

  !------------------------------------------------
  !! a vector version of the above
  subroutine pdfr_EvlnToHuman_1d(prep, qe, qh)
    type(pdf_rep), intent(in)  :: prep
    real(dp),      intent(in)  :: qe(:,ncompmin:)
    real(dp),      intent(out) :: qh(:,ncompmin:)
    integer :: n, i
    n = assert_eq(size(qh,dim=1),size(qe,dim=1),'pdfr_EvlnToHuman_1d')
    !if (GetPdfRep(qe) /= pdfr_Evln) &
    if (GetPdfRep(qe) /= prep%nf) &
         &call wae_error('pdf_EvlnToHuman_1d', 'qe is not in correct "Evln" format')
    do i = 1, n
       call pdfr_EvlnToHuman_sc(prep, qe(i,:), qh(i,:))
    end do
    call LabelPdfAsRep(qh,pdfr_Human)
  end subroutine pdfr_EvlnToHuman_1d
  

  !=== next follow overloaded versions with integer (nf) spec of rep
  subroutine pdfr_HumanToEvln_nf_sc(nf_lcl, qh, qe)
    integer ,      intent(in)  :: nf_lcl
    real(dp),      intent(in)  :: qh(:)
    real(dp),      intent(out) :: qe(:)
    call CopyHumanPdfToEvln(DefaultEvlnRep(nf_lcl), qh, qe)
  end subroutine pdfr_HumanToEvln_nf_sc
  !
  subroutine pdfr_HumanToEvln_nf_1d(nf_lcl, qh, qe)
    integer ,      intent(in)  :: nf_lcl
    real(dp),      intent(in)  :: qh(:,:)
    real(dp),      intent(out) :: qe(:,:)
    call CopyHumanPdfToEvln(DefaultEvlnRep(nf_lcl), qh, qe)
  end subroutine pdfr_HumanToEvln_nf_1d
  !
  subroutine pdfr_EvlnToHuman_nf_sc(nf_lcl, qe, qh)
    integer ,      intent(in)  :: nf_lcl
    real(dp),      intent(in)  :: qe(:)
    real(dp),      intent(out) :: qh(:)
    call CopyEvlnPdfToHuman(DefaultEvlnRep(nf_lcl), qe, qh)
  end subroutine pdfr_EvlnToHuman_nf_sc
  !
  subroutine pdfr_EvlnToHuman_nf_1d(nf_lcl, qe, qh)
    integer ,      intent(in)  :: nf_lcl
    real(dp),      intent(in)  :: qe(:,:)
    real(dp),      intent(out) :: qh(:,:)
    call CopyEvlnPdfToHuman(DefaultEvlnRep(nf_lcl), qe, qh)
  end subroutine pdfr_EvlnToHuman_nf_1d
  


  !--------------------------------------------------------------
  !! Label the pdf with a "key" corresponding to the representation
  !! Make part of Evln the key random so that if we subtract two pdfs
  !! in the Evln representation they will not accidentally give
  !! something in the human representation.
  subroutine LabelPdfAsRep(q,irep)
    use random
    real(dp), intent(inout) :: q(0:,ncompmin:)
    integer,  intent(in)    :: irep
    
    if (ubound(q,dim=2) /= ncompmax) call wae_error('LabelPdfAsRep',&
         &'upper bound of q does not correspond to ncompmax; it is:',&
         &intval=ubound(q,dim=2))

    if (ubound(q,dim=1) < 4) call wae_error('LabelPdfAsRep',&
         &'grid is too small to hold pdf flavour representation info; size is:',&
         &intval=ubound(q,dim=1))

    ! very wasteful labelling, but in f90 it is hard to see 
    ! what else can be done...
    select case(irep)
    case(pdfr_Human)  
       q(:,iflv_info) = zero
    case(1:iflv_max)
       !q(0,iflv_info)  = pi*1e-1_dp
       !q(1,iflv_info)  = one+ran()
       !q(2:,iflv_info) = zero
       q(0,iflv_info)   = pi*1e-1_dp
       q(1,iflv_info)   = one+ran()
       q(2:3,iflv_info) = irep*q(0:1,iflv_info)
       q(4:,iflv_info) = zero
    case default
       call wae_error('LabelPdfAsRep','Unrecognized irep:',intval=irep)
    end select
  end subroutine LabelPdfAsRep

  
  !======================================================================
  !! Label the PDF as being in the human representation
  subroutine LabelPdfAsHuman(q)
    real(dp), intent(inout) :: q(0:,ncompmin:)
    call LabelPdfAsRep(q,pdfr_human)
  end subroutine LabelPdfAsHuman
  
  !-------------------------------------------------------------
  !! This tells us what representation it is in.
  !! 
  !! If it is in hte human representation, then it will return
  !! pdfr_Human, otherwise it will return the number of flavours
  !! associated with the current representation
  function GetPdfRep(q) result(irep)
    real(dp), intent(in) :: q(0:,ncompmin:)
    integer              :: irep
    real(dp)             :: drep
    real(dp), parameter  :: rep_tolerance = 1e-7_dp

    if (ubound(q,dim=2) /= ncompmax) call wae_error('GetPdfRep',&
         &'upper bound of q does not correspond to ncompmax; it is:',&
         &intval=ubound(q,dim=2))

    !if (q(0,iflv_info) == zero .neqv. q(1,iflv_info) == zero) then
    !   call wae_error('GetPdfRep', 'Inconsistent behaviour in iflv_info')
    !else 
    if (q(0,iflv_info) == zero .and. q(1,iflv_info) == zero) then
       irep = pdfr_Human
    else
       !irep = pdfr_Evln
       ! we now get the number of flavours as a ratio of different
       ! entries in index ncompmax, to help distinguish the different
       ! representation for each separate nf value
       drep = (abs(q(2,iflv_info))+abs(q(3,iflv_info)))/&
            & (abs(q(0,iflv_info))+abs(q(1,iflv_info)))
       irep = nint(drep)
       if (abs(drep - irep) > rep_tolerance) then
          call wae_error('GetPdfRep',&
               &'representation seems to be inconsistent (non-integer):',&
               &dbleval=drep)
       end if
    end if
  end function GetPdfRep
  
  
end module pdf_representation
