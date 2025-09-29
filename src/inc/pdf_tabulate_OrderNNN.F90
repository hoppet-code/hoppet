!! code that is included in pdf_tabulate.f90, with __HOPPET_InterpOrder__
!! acting as a kind of template parameter, so as to help generate
!! efficient code (e.g. where hard-coded loop bounds facilitate
!! loop unrolling)
!!
!! This code gets included within pdf_tabulate module and so automatically
!! has access to all the module's variables and procedures.
!!

!! To test macro expansion try
!!   gfortran -E -D__HOPPET_InterpOrder__=4 ../src/inc/pdf_tabulate_OrderNNN.F90 -o /tmp/a.i
!! and inspect /tmp/a.i

! include the Fortran Template Library Macros, which define the CAT
! macro in a way that works across compilers
#include "ftlMacros.inc"
#include "pdf_tabulate.inc"

#ifndef __HOPPET_InterpOrderY__
#define __HOPPET_InterpOrderY__ __HOPPET_InterpOrder__
#endif
#ifndef __HOPPET_InterpOrderQ__
#define __HOPPET_InterpOrderQ__ __HOPPET_InterpOrder__
#endif


#if __HOPPET_InterpOrderY__ > __HOPPET_tab_min_subgrid_ny__
#error "Interpolation order is too high for the current grid configuration: __HOPPET_InterpOrder__ > __HOPPET_tab_min_subgrid_ny__"
#endif
#if __HOPPET_InterpOrderQ__ > __HOPPET_tab_min_nQ__
#error "Interpolation order is too high for the current grid configuration: __HOPPET_InterpOrder__ > __HOPPET_tab_min_nQ__"
#endif



! this definition refers to functions in the interpolation_coeffs module (interpolation.f90)
! It will expand, e.g., as fill_interp_weights3 if __HOPPET_InterpOrder__ is 3
#define fill_interp_weightsNNY   CAT(fill_interp_weights,__HOPPET_InterpOrderY__)
#define fill_interp_weightsNNQ   CAT(fill_interp_weights,__HOPPET_InterpOrderQ__)
! the next two definitions are shorthands for functions defined below.
! Again, NNNN -> __HOPPET_InterpOrder__
#define EvalPdfTable_get_weights_orderNNNN CAT3(EvalPdfTable_get_weights_order,__HOPPET_InterpOrderY__,__HOPPET_InterpOrderQ__)
#define EvalPdfTable_yQ_orderNNNN CAT3(EvalPdfTable_yQ_order,__HOPPET_InterpOrderY__,__HOPPET_InterpOrderQ__)
#define EvalPdfTable_yQf_orderNNNN CAT3(EvalPdfTable_yQf_order,__HOPPET_InterpOrderY__,__HOPPET_InterpOrderQ__)
#define EvalPdfTable1D_yQ_orderNNNN CAT3(EvalPdfTable1D_yQ_order,__HOPPET_InterpOrderY__,__HOPPET_InterpOrderQ__)

subroutine EvalPdfTable_get_weights_orderNNNN(tab,y,Q,y_wgts, lnlnQ_wgts, iylo, ilnlnQ)
  use interpolation_coeffs; use convolution
  type(pdf_table), intent(in), target :: tab
  real(dp),        intent(in)         :: y, Q
  real(dp),        intent(out)        :: y_wgts(0:__HOPPET_InterpOrderY__), lnlnQ_wgts(0:__HOPPET_InterpOrderQ__)
  integer,         intent(out)        :: iylo, ilnlnQ
  !----------------------------------------
  integer, parameter :: NNY = __HOPPET_InterpOrderY__, halfNNY=(NNY)/2 !halfNNY=(NNY-1)/2
  integer, parameter :: NNQ = __HOPPET_InterpOrderQ__, halfNNQ=(NNQ)/2 !halfNNQ=(NNQ-1)/2
  real(dp) :: lnlnQ
  type(grid_def),   pointer :: grid, subgd
  type(pdfseginfo), pointer :: seginfo
  integer  :: i_nf, igd, iy_offset, iQ
  real(dp) :: ynorm, lnlnQ_norm

  !----- get the info for the y interpolation
  ! (NB 2025-09-09, tried manually inlining the routine on M2Pro-gfortran15-O3,
  ! but made no difference to timing)
  call tab_get_grid_ptr(tab, y, subgd, iy_offset)

  ynorm = y / subgd%dy
  iylo = int(ynorm) - halfNNY
  if (iylo < 0) iylo = 0
  if (iylo + NNY > subgd%ny) iylo = subgd%ny - NNY
  call fill_interp_weightsNNY(ynorm - iylo, y_wgts)
  iylo = iylo + iy_offset

  !----- next deal with the Q interpolation
  lnlnQ = lnln(tab,Q)

  if (lnlnQ < tab%lnlnQ_min) then
    lnlnQ = tab%lnlnQ_min
    ilnlnQ = lbound(tab%lnlnQ_vals,1)
    lnlnQ_wgts = zero
    if (tab%freeze_at_Qmin) then
      lnlnQ_wgts(0) = one
    endif
    return
  else if (lnlnQ > tab%lnlnQ_max) then
    call wae_error("EvalPdfTable_yQ_orderNNNN","Q was too large",dbleval=Q)
  endif

  if (.not. tab%nf_info_associated) then
    seginfo => tab%seginfo_no_nf
  else
    do i_nf = tab%nflo, tab%nfhi-1
      if (lnlnQ < tab%seginfo(i_nf)%lnlnQ_hi) exit
    end do
    seginfo => tab%seginfo(i_nf)
  end if
  if (seginfo%ilnlnQ_hi == seginfo%ilnlnQ_lo) then
    lnlnQ_wgts = zero
    lnlnQ_wgts(0) = 1.0_dp
    return
  endif

  lnlnQ_norm = (lnlnQ - seginfo%lnlnQ_lo) * seginfo%inv_dlnlnQ
  if (seginfo%ilnlnQ_hi - seginfo%ilnlnQ_lo < NNQ) then
    call wae_error("EvalPdfTable_yQ_orderNNNN","not enough points in Q segment")
  end if
  ilnlnQ = int(lnlnQ_norm) + seginfo%ilnlnQ_lo - halfNNQ
  if (ilnlnQ    < seginfo%ilnlnQ_lo) ilnlnQ = seginfo%ilnlnQ_lo
  if (ilnlnQ+NNQ > seginfo%ilnlnQ_hi) ilnlnQ = seginfo%ilnlnQ_hi - NNQ
  call fill_interp_weightsNNQ(lnlnQ_norm - (ilnlnQ-seginfo%ilnlnQ_lo), lnlnQ_wgts)
end subroutine EvalPdfTable_get_weights_orderNNNN


! 2025-09-08 attempt to explore potential for speedup in PDF
! evaluation when we hard-code and inline as much as possible
subroutine EvalPdfTable_yQ_orderNNNN(tab,y,Q,res) 
  use interpolation_coeffs; use convolution
  type(pdf_table), intent(in), target :: tab
  real(dp),        intent(in)         :: y, Q
  real(dp),        intent(out)        :: res(iflv_min:)
  !----------------------------------------
  integer, parameter :: NNY = __HOPPET_InterpOrderY__
  integer, parameter :: NNQ = __HOPPET_InterpOrderQ__
  real(dp) :: y_wgts(0:NNY), lnlnQ_wgts(0:NNQ)
  integer  :: iylo, ilnlnQ, iQ, iflv

  call EvalPdfTable_get_weights_orderNNNN(tab, y, Q, y_wgts, lnlnQ_wgts, iylo, ilnlnQ)

  ! now do the interpolation.
  ! first do the flavours we know we will need (this loop is easier to unroll)
  do iflv = iflv_min, iflv_max
    res(iflv) = sum(tab%tab(iylo:iylo+NNY, iflv,ilnlnQ  ) * y_wgts) * lnlnQ_wgts(0)
    do iQ = 1, NNQ
      res(iflv) = res(iflv) + sum(tab%tab(iylo:iylo+NNY, iflv,ilnlnQ+iQ) * y_wgts) * lnlnQ_wgts(iQ)
    end do
    !res(iflv) = sum(tab%tab(iylo:iylo+NNY, iflv,ilnlnQ  ) * y_wgts) * lnlnQ_wgts(0) &
    !    + sum(tab%tab(iylo:iylo+NNY, iflv,ilnlnQ+1) * y_wgts) * lnlnQ_wgts(1) &
    !    + sum(tab%tab(iylo:iylo+NNY, iflv,ilnlnQ+2) * y_wgts) * lnlnQ_wgts(2)
  end do

  ! and then do any remaining flavours (separating things gains us a couple of ns)
  do iflv = iflv_max+1, tab%tab_iflv_max
    res(iflv) = sum(tab%tab(iylo:iylo+NNY, iflv,ilnlnQ  ) * y_wgts) * lnlnQ_wgts(0)
    do iQ = 1, NNQ
      res(iflv) = res(iflv) + sum(tab%tab(iylo:iylo+NNY, iflv,ilnlnQ+iQ) * y_wgts) * lnlnQ_wgts(iQ)
    end do
  end do
end subroutine EvalPdfTable_yQ_orderNNNN

subroutine EvalPdfTable1D_yQ_orderNNNN(tab,y,Q,res) 
  use interpolation_coeffs; use convolution
  type(pdf_table), intent(in), target :: tab(:)
  real(dp),        intent(in)         :: y, Q
  real(dp),        intent(out)        :: res(iflv_min:,:)
  !----------------------------------------
  integer, parameter :: NNY = __HOPPET_InterpOrderY__
  integer, parameter :: NNQ = __HOPPET_InterpOrderQ__
  real(dp) :: y_wgts(0:NNY), lnlnQ_wgts(0:NNQ)
  integer  :: iylo, ilnlnQ, iQ, iflv, itab
  integer :: iflv_max

  iflv_max = tab(1)%tab_iflv_max

  call EvalPdfTable_get_weights_orderNNNN(tab(1), y, Q, y_wgts, lnlnQ_wgts, iylo, ilnlnQ)

  do itab = 1, size(tab)
    ! now do the interpolation.
    ! first do the flavours we know we will need (this loop is easier to unroll)
    do iflv = iflv_min, iflv_max
      res(iflv, itab) = sum(tab(itab)%tab(iylo:iylo+NNY, iflv,ilnlnQ  ) * y_wgts) * lnlnQ_wgts(0)
      do iQ = 1, NNQ
        res(iflv, itab) = res(iflv, itab) + sum(tab(itab)%tab(iylo:iylo+NNY, iflv,ilnlnQ+iQ) * y_wgts) * lnlnQ_wgts(iQ)
      end do
      !res(iflv) = sum(tab%tab(iylo:iylo+NNY, iflv,ilnlnQ  ) * y_wgts) * lnlnQ_wgts(0) &
      !    + sum(tab%tab(iylo:iylo+NNY, iflv,ilnlnQ+1) * y_wgts) * lnlnQ_wgts(1) &
      !    + sum(tab%tab(iylo:iylo+NNY, iflv,ilnlnQ+2) * y_wgts) * lnlnQ_wgts(2)
    end do

    ! and then do any remaining flavours (separating things gains us a couple of ns)
    do iflv = iflv_max+1, iflv_max
      res(iflv, itab) = sum(tab(itab)%tab(iylo:iylo+NNY, iflv,ilnlnQ  ) * y_wgts) * lnlnQ_wgts(0)
      do iQ = 1, NNQ
        res(iflv, itab) = res(iflv, itab) + sum(tab(itab)%tab(iylo:iylo+NNY, iflv,ilnlnQ+iQ) * y_wgts) * lnlnQ_wgts(iQ)
      end do
    end do

  end do ! itab
end subroutine EvalPdfTable1D_yQ_orderNNNN


function EvalPdfTable_yQf_orderNNNN(tab,y,Q,iflv) result(res)
  use interpolation_coeffs; use convolution
  type(pdf_table), intent(in), target :: tab
  real(dp),        intent(in)         :: y, Q
  integer,         intent(in)         :: iflv
  real(dp)                            :: res
  !----------------------------------------
  integer, parameter :: NNY = __HOPPET_InterpOrderY__
  integer, parameter :: NNQ = __HOPPET_InterpOrderQ__
  real(dp) :: y_wgts(0:NNY), lnlnQ_wgts(0:NNQ)
  integer  :: iylo, ilnlnQ, iQ

  call EvalPdfTable_get_weights_orderNNNN(tab, y, Q, y_wgts, lnlnQ_wgts, iylo, ilnlnQ)

  ! diagnostics
  !write(6,*) "orderNNNN: iylo = ", iylo, "y_wgts = ", y_wgts
  !write(6,*) "orderNNNN: ilnlnQ = ", ilnlnQ, "lnlnQ_wgts = ", lnlnQ_wgts
  !write(6,*) "orderNNNN: Qvals = ", tab%Q_vals(ilnlnQ:ilnlnQ+NN)

  ! now do the interpolation.
  ! first do the flavours we know we will need (this loop is easier to unroll)
  res = sum(tab%tab(iylo:iylo+NNY, iflv,ilnlnQ  ) * y_wgts) * lnlnQ_wgts(0)
  do iQ = 1, NNQ
    res = res + sum(tab%tab(iylo:iylo+NNY, iflv,ilnlnQ+iQ) * y_wgts) * lnlnQ_wgts(iQ)
  end do
end function EvalPdfTable_yQf_orderNNNN

! clear up all definitions, both local and those that were supplied
#undef EvalPdfTable_yQ_orderNNNN
#undef EvalPdfTable_yQf_orderNNNN
#undef fill_interp_weightsNNNN
#undef __HOPPET_InterpOrderY__
#undef __HOPPET_InterpOrderQ__
#undef __HOPPET_InterpOrder__