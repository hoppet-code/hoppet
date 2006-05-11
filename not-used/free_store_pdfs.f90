module free_store_pdfs
  use types; use consts_dp
  use convolution; use pdf_general
  use warnings_and_errors
  implicit none
  private

  type pdf_holder
     real(dp), pointer :: pdf(:,:) => null()
  end type pdf_holder
  ! pointers to the actual pdfs
  type(pdf_holder), pointer, save :: pdfs(:)
  ! where to look in pdfs to find an available slot
  integer         , pointer, save :: free_queue(:)
  ! where to look in free_queue for next pdf
  integer                  , save :: position_on_queue

  public :: free_store_pdfs_init
  public :: new_handle
  public :: set_pointer_to_pdf
  public :: release_pdf

contains

  !----------------------------------------------
  !! initialise the free store of pdfs
  subroutine free_store_pdfs_init
    integer, parameter :: initial_size = 100
    integer            :: i
    allocate(pdfs(initial_size))
    allocate(free_queue(initial_size))
    forall(i = 1:size(free_queue)) free_queue(i) = i
    position_on_queue = 1
  end subroutine free_store_pdfs_init

  !-----------------------------------------
  !! return a new pdf handle
  integer function new_handle(gd)
    type(grid_def), intent(in) :: gd
    
    if (position_on_queue > size(free_queue)) call resize_free_store

    new_handle = free_queue(position_on_queue)
    position_on_queue = position_on_queue + 1

    call pdfgen_AllocPDF(gd,pdfs(new_handle)%pdf)
  end function new_handle
  
  !--------------------------------------------
  !! Returns a pointer to the pdf...
  subroutine  set_pointer_to_pdf(handle, pdf) 
    integer, intent(in) :: handle
    real(dp), pointer   :: pdf(:,:)

    call validate(handle, 'get_pdf')
    pdf => pdfs(handle)%pdf
  end subroutine set_pointer_to_pdf
  
  

  !-------------------------------------------
  !! declares the pdf referred to by the handle to be no longer of any
  !! use (so memory can be deallocated).
  subroutine release_pdf(handle)
    integer, intent(in) :: handle
    real(dp), pointer   :: pdf(:,:)

    call validate(handle, 'release_pdf')
    
    ! release the memory
    deallocate(pdfs(handle)%pdf)

    ! update the free_queue array to indicate that this handle
    ! is once again available
    position_on_queue = position_on_queue - 1
    free_queue(position_on_queue) = handle
  end subroutine release_pdf
  

  !-------------------------------------------------------------
  !! Increase the size of the free store to allow for extra pdfs
  subroutine resize_free_store
    type(pdf_holder), pointer :: new_pdfs(:)
    integer,          pointer :: new_free_queue(:)
    integer :: old_size, new_size, i

    old_size = size(pdfs)
    new_size = 2*old_size

    ! allocate new arrays
    allocate(new_pdfs(new_size))
    allocate(new_free_queue(new_size))
    
    ! copy old information on free_queue + get new information
    new_free_queue(1:old_size) = free_queue  ! maybe not necessary?
    forall(i=old_size+1:new_size) new_free_queue(i) = i

    ! copy pointers to pdfs
    do i = 1, old_size
       new_pdfs(i) = pdfs(i)
    end do

    ! deallocate old structures & replace with new
    deallocate(pdfs,free_queue)
    pdfs => new_pdfs
    free_queue => new_free_queue
    
  end subroutine resize_free_store
  

  !--------------------------------------------
  !! check that a handle is legitimate
  subroutine validate(handle,source)
    integer,          intent(in) :: handle
    character(len=*), intent(in) :: source
    if (handle > 0 .and. handle <= size(pdfs)) then
       if (associated(pdfs(handle)%pdf)) return
    end if
    call wae_error('validate ('//source//')',&
         &         'handle value does not correspond to a known pdf',&
         &         intval = handle)
  end subroutine validate
  
end module free_store_pdfs
