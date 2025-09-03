!-----------------------------------------------------------
!! Routines related to allocation and initialisation of parton 
!! distributions.
module pdf_general
  use types; use consts_dp
  !-- do some recycling (actually only the 2d versions will
  !   be of any use, but hopefully that will not cause problems)
  use convolution
  use pdf_representation
  implicit none
  private

  interface AllocPDF
     module procedure pdfgen_AllocPDF_0d, pdfgen_AllocPDF_1d
  end interface
  public :: AllocPDF


  interface InitPDF
     module procedure pdfgen_InitPDF_, pdfgen_InitPDF_a,&
          & pdfgen_InitPDF_ai
  end interface
  public :: InitPDF

  interface InitPDFSub
     module procedure pdfgen_InitPDFSub_, pdfgen_InitPDFSub_a,&
          & pdfgen_InitPDFSub_ai
  end interface
  public :: InitPDFSub

  public :: InitPDF_LHAPDF

  interface AllocInitPDF
     module procedure pdfgen_AllocInitPDF_, pdfgen_AllocInitPDF_a,&
          & pdfgen_AllocInitPDF_ai
  end interface
  public :: AllocInitPDF

  interface AllocInitPDFSub
     module procedure pdfgen_AllocInitPDFSub_, &
          & pdfgen_AllocInitPDFSub_a, pdfgen_AllocInitPDFSub_ai
  end interface
  public :: AllocInitPDFSub

  interface operator(.anti.)
     module procedure pdfgen_anti_0d, pdfgen_anti_1d
  end interface
  public :: operator(.anti.)
contains
  !-------------------------------------------------------
  ! allocate a parton distribution; dimensionality refers
  ! to extra directions over and above (x,nf).
  subroutine pdfgen_AllocPDF_0d(grid,q)
    type(grid_def), intent(in) :: grid
    real(dp),       pointer    :: q(:,:)

    call AllocGridQuant(grid,q,ncompmin, ncompmax)
    call LabelPdfAsRep(q,pdfr_Human)
  end subroutine pdfgen_AllocPDF_0d

  !-------------------------------------------------------
  subroutine pdfgen_AllocPDF_1d(grid,q, nl, nh)
    type(grid_def), intent(in) :: grid
    real(dp),       pointer    :: q(:,:,:)
    integer ,       intent(in) :: nl, nh
    integer :: i

    call AllocGridQuant(grid,q,ncompmin, ncompmax, nl, nh)
    do i = nl, nh
       call LabelPdfAsRep(q(:,:,i),pdfr_Human)
    end do
  end subroutine pdfgen_AllocPDF_1d



  !----------------------------------------------------------------------
  ! All of these routine just redirect to the official routine, with
  ! the added features that if it's a standard QCD PDF, then they
  ! redirect only the components iflv_min:iflv_max.
  !
  ! They also label the result as being in the "Human" representation
  recursive subroutine pdfgen_InitPDF_(grid, gq, func)
    real(dp),         intent(inout) :: gq(0:,ncompmin:)
    type(grid_def),   intent(in)    :: grid
    interface
       function func(x,n)
         use types; implicit none
         real(dp), intent(in) :: x
         integer , intent(in) :: n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------
    if (ubound(gq,dim=2) == ncompmax) then
       ! if it's a standard QCD-only PDF, then pass just -6:6 to the subroutine
       call InitGridQuant(grid, gq(:,iflv_min:iflv_max), func)
    else
       ! otherwise pass the whole thing; the subroutine must be aware that
       ! index 7 is a special ("representation") index
       call InitGridQuant(grid, gq(:,iflv_min:iflv_max), func)
    end if
    call LabelPdfAsRep(gq(:,ncompmin:ncompmax),pdfr_Human)
  end subroutine pdfgen_InitPDF_

  !----------------------------------------------------------------------
  ! version intended for use when there is an extra argument whose
  ! value is fixed and needs to be passed to func
  !
  ! updated for multi
  recursive subroutine pdfgen_InitPDF_a(grid, gq, func, axtra)
    real(dp),         intent(inout) :: gq(0:,ncompmin:)
    type(grid_def),   intent(in)    :: grid
    real(dp),         intent(in)    :: axtra
    interface
       function func(x,axtra,n)
         use types; implicit none
         real(dp), intent(in) :: x,axtra
         integer , intent(in) :: n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------
    if (ubound(gq,dim=2) == ncompmax) then
       ! if it's a standard QCD-only PDF, then pass just -6:6 to the subroutine
       call InitGridQuant(grid, gq(:,iflv_min:iflv_max), func, axtra)
    else
       ! otherwise pass the whole thing; the subroutine must be aware that
       ! index 7 is a special ("representation") index
       call InitGridQuant(grid, gq(:,iflv_min:iflv_max), func, axtra)
    end if
    call LabelPdfAsRep(gq(:,ncompmin:ncompmax),pdfr_Human)
  end subroutine pdfgen_InitPDF_a

  !----------------------------------------------------------------------
  ! version intended for use when there is an extra argument whose
  ! value is fixed and needs to be passed to func
  !
  ! updated for multi
  recursive subroutine pdfgen_InitPDF_ai(grid, gq, func, axtra, ixtra)
    real(dp),         intent(inout) :: gq(0:,ncompmin:)
    type(grid_def),   intent(in)    :: grid
    real(dp),         intent(in)    :: axtra
    integer,          intent(in)    :: ixtra
    interface
       function func(x,axtra,ixtra,n)
         use types; implicit none
         real(dp), intent(in) :: x,axtra
         integer , intent(in) :: ixtra,n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------
    if (ubound(gq,dim=2) == ncompmax) then
       ! if it's a standard QCD-only PDF, then pass just -6:6 to the subroutine
       call InitGridQuant(grid, gq(:,iflv_min:iflv_max), func, axtra, ixtra)
    else
       ! otherwise pass the whole thing; the subroutine must be aware that
       ! index 7 is a special ("representation") index
       call InitGridQuant(grid, gq(:,iflv_min:iflv_max), func, axtra, ixtra)
    end if
    call LabelPdfAsRep(gq(:,ncompmin:ncompmax),pdfr_Human)
  end subroutine pdfgen_InitPDF_ai


  !----------------------------------------------------------------------
  ! updated for multi
  recursive subroutine pdfgen_InitPDFSub_(grid, gq, sub)
    real(dp),         intent(inout) :: gq(0:,ncompmin:)
    type(grid_def),   intent(in)    :: grid
    interface
       subroutine sub(y,res)
         use types; implicit none
         real(dp), intent(in)  :: y
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    call InitGridQuantSub(grid, gq(:,iflv_min:iflv_max), sub)
    call LabelPdfAsRep(gq,pdfr_Human)
  end subroutine pdfgen_InitPDFSub_

  recursive subroutine pdfgen_InitPDFSub_a(grid, gq, sub, axtra)
    real(dp),         intent(inout) :: gq(0:,ncompmin:)
    type(grid_def),   intent(in)    :: grid
    real(dp),         intent(in)    :: axtra
    interface
       subroutine sub(y, axtra, res)
         use types; implicit none
         real(dp), intent(in)  :: y, axtra
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    call InitGridQuantSub(grid, gq(:,iflv_min:iflv_max), sub, axtra)
    call LabelPdfAsRep(gq,pdfr_Human)
  end subroutine pdfgen_InitPDFSub_a

  recursive subroutine pdfgen_InitPDFSub_ai(grid, gq, sub, axtra, ixtra)
    real(dp),         intent(inout) :: gq(0:,ncompmin:)
    type(grid_def),   intent(in)    :: grid
    real(dp),         intent(in)    :: axtra
    integer,          intent(in)    :: ixtra
    interface
       subroutine sub(y, axtra, ixtra, res)
         use types; implicit none
         real(dp), intent(in)  :: y, axtra
         integer,  intent(in)  :: ixtra
         real(dp), intent(out) :: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------
    call InitGridQuantSub(grid, gq(:,iflv_min:iflv_max), sub, axtra, ixtra)
    call LabelPdfAsRep(gq,pdfr_Human)
  end subroutine pdfgen_InitPDFSub_ai



  !======================================================================
  !! Initialise the subroutine using an LHAPDF style subroutine
  !!
  !! It expects LHAsub to set components -6:6 with plain QCD
  !! evolution, but the full size of the flavour dimension if the
  !! upper limit is anything other than ncompax=7 (e.g. for QED
  !! evolution it should go from -6 to ncompmaxLeptons=11).
  subroutine InitPDF_LHAPDF(grid, gq, LHAsub, Q)
    real(dp),         intent(inout) :: gq(0:,ncompmin:)
    type(grid_def),   intent(in)    :: grid
    real(dp),         intent(in)    :: Q
    interface
       subroutine LHAsub(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine LHAsub
    end interface
    !-------------------------------------------
    integer :: ncompmax_lcl

    ! if the upper bound of gq(:,:) is the usual ncompmax
    ! then assume the user's routine will fill -6:6.
    ! Otherwise the user should fill -6:ubound(gq,dim=2)
    ncompmax_lcl = ubound(gq,dim=2)
    if (ncompmax_lcl == ncompmax) ncompmax_lcl = iflv_max

    !call InitGridQuantLHAPDF(grid, gq(:,iflv_min:iflv_max), LHAsub, Q)
    ! PN & GZ: the following modification is needed for the streamlined interface with QED evolution
    call InitGridQuantLHAPDF(grid, gq(:,iflv_min:ncompmax_lcl), LHAsub, Q)        
    call LabelPdfAsRep(gq,pdfr_Human)
  end subroutine InitPDF_LHAPDF
  



  !------------------------------------------------------
  ! allocate and initialise in various forms!
  subroutine pdfgen_AllocInitPDF_(grid,q, func)
    type(grid_def), intent(in) :: grid
    real(dp),       pointer    :: q(:,:)
    interface
       function func(y,n)
         use types; implicit none
         real(dp), intent(in) :: y
         integer , intent(in) :: n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------

    call AllocPDF(grid,q)
    call InitPDF(grid, q, func)
  end subroutine pdfgen_AllocInitPDF_

  !--------------------------------------------------
  subroutine pdfgen_AllocInitPDF_a(grid,q, func, axtra)
    type(grid_def), intent(in) :: grid
    real(dp),       pointer    :: q(:,:)
    real(dp),       intent(in) :: axtra
    interface
       function func(y,axtra,n)
         use types; implicit none
         real(dp), intent(in) :: y,axtra
         integer , intent(in) :: n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------

    call AllocPDF(grid,q)
    call InitPDF(grid, q, func, axtra)
  end subroutine pdfgen_AllocInitPDF_a

  !-------------------------------------------------------------
  subroutine pdfgen_AllocInitPDF_ai(grid,q, func, axtra, ixtra)
    type(grid_def), intent(in) :: grid
    real(dp),       pointer    :: q(:,:)
    real(dp),       intent(in) :: axtra
    integer,        intent(in) :: ixtra
    interface
       function func(y,axtra,ixtra,n)
         use types; implicit none
         real(dp), intent(in) :: y,axtra
         integer , intent(in) :: ixtra,n
         real(dp)             :: func(n)
       end function func
    end interface
    !-----------------------------------------

    call AllocPDF(grid,q)
    call InitPDF(grid, q, func, axtra, ixtra)
  end subroutine pdfgen_AllocInitPDF_ai

  !------------------------------------------------------
  ! allocate and initialise in various forms!
  subroutine pdfgen_AllocInitPDFSub_(grid,q, sub)
    type(grid_def), intent(in) :: grid
    real(dp),       pointer    :: q(:,:)
    interface
       subroutine sub(y,res)
         use types; implicit none
         real(dp), intent(in) :: y
         real(dp), intent(out):: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------

    call AllocPDF(grid,q)
    call InitPDFSub(grid, q, sub)
  end subroutine pdfgen_AllocInitPDFSub_

  !--------------------------------------------------
  subroutine pdfgen_AllocInitPDFSub_a(grid,q, sub, axtra)
    type(grid_def), intent(in) :: grid
    real(dp),       pointer    :: q(:,:)
    real(dp),       intent(in) :: axtra
    interface
       subroutine sub(y,axtra,res)
         use types; implicit none
         real(dp), intent(in) :: y,axtra
         real(dp), intent(out):: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------

    call AllocPDF(grid,q)
    call InitPDFSub(grid, q, sub, axtra)
  end subroutine pdfgen_AllocInitPDFSub_a

  !-------------------------------------------------------------
  subroutine pdfgen_AllocInitPDFSub_ai(grid,q, sub, axtra, ixtra)
    type(grid_def), intent(in) :: grid
    real(dp),       pointer    :: q(:,:)
    real(dp),       intent(in) :: axtra
    integer,        intent(in) :: ixtra
    interface
       subroutine sub(y,axtra,ixtra,res)
         use types; implicit none
         real(dp), intent(in) :: y,axtra
         integer , intent(in) :: ixtra
         real(dp), intent(out):: res(:)
       end subroutine sub
    end interface
    !-----------------------------------------

    call AllocPDF(grid,q)
    call InitPDFSub(grid, q, sub, axtra, ixtra)
  end subroutine pdfgen_AllocInitPDFSub_ai

  !---------------------------------------------------
  ! write it out in full rather than using f90 array
  ! notation shortcuts to reduce need for assumptions about
  ! relation between iflv_min and ncompmin
  function pdfgen_anti_0d(q) result(antiq)
    real(dp), intent(in) :: q(:,ncompmin:)
    real(dp)             :: antiq(size(q,1),lbound(q,2):ubound(q,2))
    integer :: i

    !-- try to catch stupid and illegal uses...
    if (ubound(q,2) /= ncompmax) call wae_error('pdfgen_anti_0d',&
         &'ubound of second dimension of q should be ncompmax, instead it is',&
         &intval=ubound(q,2))

    do i = ncompmin, ncompmax
       if (i >= iflv_min .and. i <= iflv_max) then
          antiq(:,i) = q(:,-i)
       else
          antiq(:,i) = q(:,i)
       end if
    end do
  end function pdfgen_anti_0d
  
  !---------------------------------------------------
  function pdfgen_anti_1d(q) result(antiq)
    real(dp), intent(in) :: q(:,ncompmin:,:)
    real(dp)             :: antiq(size(q,1),lbound(q,2):ubound(q,2),size(q,3))
    integer :: i
    do i = 1, size(q,dim=3)
       antiq(:,:,i) = pdfgen_anti_0d(q(:,:,i))
    end do
  end function pdfgen_anti_1d
  
end module pdf_general
