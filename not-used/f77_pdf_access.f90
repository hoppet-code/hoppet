!======================================================================
!! Module containing the locally used information for the f77 pdf acess
!! routines 
module f77_pdf_access
  use types; use consts_dp
  use convolution; use pdf_general; use conv_objects
  use holders; use pdf_general
  use free_store_pdfs
  implicit none

  
  type(grid_def),     save :: gd, gdarray(3)
  type(sigma_holder), save :: sh

  interface 
     function pdfconv_new_pdf(pdf_subroutine,Q) result(pdf_handle)
       use types
       implicit none
       real(dp), intent(in) :: Q
       integer              :: pdf_handle
       interface
          subroutine pdf_subroutine(x,Q,res)
            use types; implicit none
            real(dp), intent(in)  :: x,Q
            real(dp), intent(out) :: res(*)
          end subroutine pdf_subroutine
       end interface
     end function pdfconv_new_pdf
  end interface
  
end module f77_pdf_access



!======================================================================
!! do all necessary initialisation
subroutine pdfconv_start(dy,nf)
  use f77_pdf_access
  implicit none
  real(dp), intent(in) :: dy
  integer,  intent(in) :: nf
  !-------------------------------------
  real(dp)          :: ymax
  integer           :: order

  ! initialise our grids
  ymax = 12.0_dp
  order = 5
  call conv_InitGridDef(gdarray(3),dy/9.0_dp,0.5_dp, order=order)
  call conv_InitGridDef(gdarray(2),dy/3.0_dp,2.0_dp, order=order)
  call conv_InitGridDef(gdarray(1),dy,       ymax  ,order=order)
  call conv_InitGridDef(gd,gdarray(1:3),locked=.true.)

  ! initialise the pdf free store
  call free_store_pdfs_init

  ! initialise splitting-function holder
  call holder_InitSigma(gd,sh,nloop=2,nflo=nf,nfhi=nf)
  call holder_SetNf(sh,nf)
end subroutine pdfconv_start


!======================================================================
!! returns an integer handle to the pdf corresponding to calling
!! pdf_subroutine(x,Q,f) at scale Q
function pdfconv_new_pdf(pdf_subroutine,Q) result(pdf_handle)
  use f77_pdf_access
  implicit none
  real(dp), intent(in) :: Q
  integer              :: pdf_handle
  interface
     subroutine pdf_subroutine(x,Q,res)
       use types; implicit none
       real(dp), intent(in)  :: x,Q
       real(dp), intent(out) :: res(*)
     end subroutine pdf_subroutine
  end interface
  !---------------------------------------------
  real(dp), pointer :: qdist(:,:)

  pdf_handle = new_handle(gd)
  call set_pointer_to_pdf(pdf_handle, qdist)
  call pdfgen_InitPDF_LHAPDF(gd,qdist,pdf_subroutine,Q)
end function pdfconv_new_pdf


!======================================================================
!! returns a handle (conv_handle) to the result of the convolution
!!
!!       P_LO \otimes pdf_handle
!!
!! where P_LO is the leading order splitting function matrix with the
!! normalisation such that the LO DGLAP equation is
!!
!!       d pdf / dln Q^2 = alphas/2pi * P_LO \otimes pdf
!!
function pdfconv_P_LO(pdf_handle) result(conv_handle)
  use f77_pdf_access
  implicit none
  integer, intent(in) :: pdf_handle
  integer             :: conv_handle
  !----------------------------------
  real(dp), pointer :: qdist(:,:), conv_qdist(:,:)

  call set_pointer_to_pdf(pdf_handle, qdist)
  conv_handle = new_handle(gd)
  call set_pointer_to_pdf(conv_handle, conv_qdist)
  conv_qdist = sh%P .conv. qdist
end function pdfconv_P_LO


!======================================================================
!! returns a handle (conv_handle) to the result of the convolution
!!
!!       P_NLO \otimes pdf_handle
!!
!! where P_NLO is the leading order splitting function matrix with the
!! normalisation such that the NLO DGLAP equation is
!!
!!       d pdf / dln Q^2 = alphas/2pi * P_LO \otimes pdf
!!                            + (alphas/2pi)^2 * P_NLO \otimes pdf
function pdfconv_P_NLO(pdf_handle) result(conv_handle)
  use f77_pdf_access
  implicit none
  integer, intent(in) :: pdf_handle
  integer             :: conv_handle
  !----------------------------------
  real(dp), pointer :: qdist(:,:), conv_qdist(:,:)

  call set_pointer_to_pdf(pdf_handle, qdist)
  conv_handle = new_handle(gd)
  call set_pointer_to_pdf(conv_handle, conv_qdist)
  conv_qdist = sh%P1 .conv. qdist
end function pdfconv_P_NLO



!======================================================================
!! Places in f(-6:6) the quantity x*f(x) for each of the flavours
!! ranging from -6:6
subroutine pdfconv_eval_pdf(handle, x, f)
  use f77_pdf_access
  implicit none
  integer,  intent(in)  :: handle
  real(dp), intent(in)  :: x
  real(dp), intent(out) :: f(-6:6)
  !---------------------------------------------
  real(dp), pointer :: qdist(:,:)

  call set_pointer_to_pdf(handle, qdist)
  f = qdist .atx. (x.with.gd)
end subroutine pdfconv_eval_pdf


!======================================================================
!! Declare the PDF referred to via the supplied handle to be no longer
!! needed, allowing the memory to be released back to the system. (It
!! is important to make use of this if you are examining many
!! different PDFs, one after another, otherwise memory consumption
!! will grow to be very large).
subroutine pdfconv_release_pdf(handle)
  use f77_pdf_access
  implicit none
  integer, intent(in) :: handle
  call release_pdf(handle)
end subroutine pdfconv_release_pdf
