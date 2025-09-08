!!
!!

! include the Fortran Template Library Macros, which define the CAT
! macro in a way that works across compilers
#include "ftlMacros.inc"

#define EvalPdfTable_yQ_orderNNNN CAT(EvalPdfTable_yQ_order,__HOPPET_InterpOrder__)
#define fill_interp_weightsNNNN   CAT(fill_interp_weights,__HOPPET_InterpOrder__)

  ! 2025-09-08 attempt to explore potential for speedup in PDF
  ! evaluation when we hard-code and inline as much as possible
  subroutine EvalPdfTable_yQ_orderNNNN(tab,y,Q,res) 
    use interpolation_coeffs; use convolution
    type(pdf_table), intent(in), target :: tab
    real(dp),        intent(in)         :: y, Q
    real(dp),        intent(out)        :: res(iflv_min:)
    !----------------------------------------
    integer  :: iflv
    integer, parameter :: NN = __HOPPET_InterpOrder__, halfNN=(NN-1)/2
    real(dp) :: y_wgts(0:NN), lnlnQ_wgts(0:NN), lnlnQ
    integer  :: i_nf, iylo, ilnlnQ, igd, iy_offset, iQ
    real(dp) :: ynorm, lnlnQ_norm
    type(grid_def),   pointer :: grid, subgd
    type(pdfseginfo), pointer :: seginfo


    !----- get the info for the y interpolation
    call tab_get_grid_ptr(tab, y, subgd, iy_offset)

    ynorm = y / subgd%dy
    iylo = int(ynorm) - halfNN
    if (iylo < 0) iylo = 0
    if (iylo + NN > subgd%ny) iylo = subgd%ny - NN
    call fill_interp_weightsNNNN(ynorm - iylo, y_wgts)
    iylo = iylo + iy_offset

    !----- next deal with the Q interpolation
    lnlnQ = lnln(tab,Q)
    if (lnlnQ < tab%lnlnQ_min) then
      lnlnQ = tab%lnlnQ_min
    else if (lnlnQ > tab%lnlnQ_max) then
      call wae_error("EvalPdfTable_yQ_orderNNNN","Q was too large",dbleval=Q)
    endif
    if (.not. tab%nf_info_associated) call wae_error("EvalPdfTable_yQ_orderNNNN",&
         &"tab%nf_info_associated was not set")
    do i_nf = tab%nflo, tab%nfhi-1
      if (lnlnQ < tab%seginfo(i_nf)%lnlnQ_lo) exit
    end do

    seginfo => tab%seginfo(i_nf)
    lnlnQ_norm = (lnlnQ - seginfo%lnlnQ_lo) * seginfo%inv_dlnlnQ
    if (seginfo%ilnlnQ_hi - seginfo%ilnlnQ_lo < NN) then
      call wae_error("EvalPdfTable_yQ_orderNNNN","not enough points in Q segment")
    end if
    ilnlnQ = int(lnlnQ_norm) + seginfo%ilnlnQ_lo - halfNN
    if (ilnlnQ    < seginfo%ilnlnQ_lo) ilnlnQ = seginfo%ilnlnQ_lo
    if (ilnlnQ+NN > seginfo%ilnlnQ_hi) ilnlnQ = seginfo%ilnlnQ_hi - NN
    call fill_interp_weightsNNNN(lnlnQ_norm - (ilnlnQ-seginfo%ilnlnQ_lo), lnlnQ_wgts)


    ! now do the interpolation
    ! first do the flavours we know we will need (this loop is easier to unroll)
    do iflv = iflv_min, iflv_max
      res(iflv) = sum(tab%tab(iylo:iylo+NN, iflv,ilnlnQ  ) * y_wgts) * lnlnQ_wgts(0)
      do iQ = 1, NN
        res(iflv) = res(iflv) + sum(tab%tab(iylo:iylo+NN, iflv,ilnlnQ+iQ) * y_wgts) * lnlnQ_wgts(iQ)
      end do
      !res(iflv) = sum(tab%tab(iylo:iylo+NN, iflv,ilnlnQ  ) * y_wgts) * lnlnQ_wgts(0) &
      !    + sum(tab%tab(iylo:iylo+NN, iflv,ilnlnQ+1) * y_wgts) * lnlnQ_wgts(1) &
      !    + sum(tab%tab(iylo:iylo+NN, iflv,ilnlnQ+2) * y_wgts) * lnlnQ_wgts(2)
    end do
    
    ! and then do any remaining flavours (separating things gains us a couple of ns)
    do iflv = iflv_max+1, tab%tab_iflv_max
      res(iflv) = sum(tab%tab(iylo:iylo+NN, iflv,ilnlnQ  ) * y_wgts) * lnlnQ_wgts(0)
      do iQ = 1, NN
        res(iflv) = res(iflv) + sum(tab%tab(iylo:iylo+NN, iflv,ilnlnQ+iQ) * y_wgts) * lnlnQ_wgts(iQ)
      end do
    end do
  end subroutine EvalPdfTable_yQ_orderNNNN


#undef EvalPdfTable_yQ_orderNNNN
#undef fill_interp_weightsNNNN