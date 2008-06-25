!!
!!

program resummation_example
  use hoppet_v1
  use resum_kfact
  use resum_grids_out ! Module to output grids for LO, NLO, NLOres SFs
  use interpolation ! (Why we need to call this module explicitely? )
  use convolution

  implicit none

  logical :: iread
  integer :: ix,ialpha,i,nf_in,k,jx,i_interp
  integer, parameter ::  n_interp=5
  integer, parameter ::  n_interp2=6
  integer ::  ix_interp(0:n_interp),ix_interp2(0:n_interp2),nk
  integer  :: order,nx

  real(dp) ::  xPgg_resum(nalphas,nxgrid,4)
  real(dp) ::  xPqq_resum(nalphas,nxgrid,4)
  real(dp) ::  xPqg_resum(nalphas,nxgrid,4)
  real(dp) ::  xPgq_resum(nalphas,nxgrid,4)
  real(dp) ::  x, split_interp, split_interp2
  real(dp) ::  x_interp(0:n_interp),weights(0:n_interp)
  real(dp) ::  x_interp2(0:n_interp2),weights2(0:n_interp2)
  real(dp) ::  resum_kfact1, resum_kfact2 ,ymax,dy,y
  type(grid_def) :: grid, gdarray(4), grid2
  real(dp), pointer :: xpdf(:),xpdf_conv_sr1(:),xpdf_conv_sr2(:)
  real(dp), pointer :: xpdf_conv_qq(:),xpdf_conv_qg(:)
  real(dp), pointer :: xpdf_conv_gg(:),xpdf_conv_gq(:)
  real(dp),pointer :: xvals(:)

  type(grid_conv) :: xPgg, xPqg,xPqq, xPgq
  type(split_mat) :: P_leading
  type(dglap_holder) :: dglap_h
  integer :: nloop, nf_lcl
  real(dp) :: nf_lcl_dp
  integer, parameter :: factscheme = factscheme_MSbar

  double complex :: N
  double complex :: P0NS, P0SG(2,2)

  !! Read resummed tables
  call resum_kfact_read

  !! Set number of active quarks
  nf_in=4
  call  qcd_SetNf(nf_in)

  !! Output resummed structure functions grids
  call resum_splitfun_grid(xPgg_resum,xPqq_resum,xPgq_resum,xPqg_resum) 
  
  ialpha=25
  
  open(unit=10,status="unknown",&
    & file="splresggnlo.res")
  open(unit=11,status="unknown",&
    & file="splresqqnlo.res")
   open(unit=12,status="unknown",&
    & file="splresgqnlo.res")
  open(unit=13,status="unknown",&
    & file="splresqgnlo.res")

  !! Save resummed structure function grids
  do ix=1,nxgrid
     
     write(10,"(6(f15.7,1x))") ,1.0_dp/xgrid(ix),&
      & alphas(ialpha), (xPgg_resum(ialpha,ix,i),i=1,4)
     write(11,"(6(f15.7,1x))") ,1.0_dp/xgrid(ix),&
      & alphas(ialpha), (xPqq_resum(ialpha,ix,i),i=1,4)
     write(12,"(6(f15.7,1x))") ,1.0_dp/xgrid(ix),&
      & alphas(ialpha), (xPgq_resum(ialpha,ix,i),i=1,4)
     write(13,"(6(f15.7,1x))") ,1.0_dp/xgrid(ix),&
      & alphas(ialpha), (xPqg_resum(ialpha,ix,i),i=1,4)

  enddo

  close(10)
  close(11)
  close(12)
  close(13)

  !! Obtain the weights of the general interpolation
  !! for each of the various grids

  open(unit=10,status="unknown",file="splresggnlo-interp.res")
  open(unit=11,status="unknown",file="kfatct-ggnlo-interp.res")
  nk=10000
  do k=1,nk

     ! Select the point in x
     x = 1e-8 * ( 1 / 1e-8 )**(dfloat(k)/nk)
          
     call interpolation_grid(n_interp, x, xgrid, ix_interp, x_interp)
     call general_interpolation_weights(x, x_interp, weights)

     ! Compute interpolated resummed splitting function for gg
     split_interp = 0_dp
     resum_kfact1 = 0_dp
     do i_interp = 0, n_interp
        split_interp =  split_interp + weights(i_interp) * &
        & xPgg_resum(ialpha,ix_interp(i_interp),1)
        resum_kfact1 =  resum_kfact1 + weights(i_interp) * &
        & xPgg_resum(ialpha,ix_interp(i_interp),4)
     enddo
     
     call interpolation_grid(n_interp2, x, xgrid, ix_interp2, x_interp2)
     call general_interpolation_weights(x, x_interp2, weights2)

     ! Compute interpolated resummed splitting function for gg
     split_interp2=0
     resum_kfact2=0
     do i_interp = 0, n_interp2
        split_interp2 =  split_interp2 + weights2(i_interp) * &
        & xPgg_resum(ialpha,ix_interp2(i_interp),1)
        resum_kfact2 =  resum_kfact2 + weights2(i_interp) * &
        & xPgg_resum(ialpha,ix_interp2(i_interp),4)
     enddo
     
     
     write(10,*) 1.0_dp/x,split_interp, &
      &  abs( (split_interp - split_interp2 )/ &
      &  split_interp2 )

     write(11,*) 1.0_dp/x,resum_kfact2, &
      &  abs( (resum_kfact2 - resum_kfact1 )/ &
      &  resum_kfact2 )

   !  write(6,*) 1.0_dp/x,jx,split_interp, &
   !   &  xPgg_resum(ialpha,jx,1) 

  enddo
  close(10)
  close(11)

 write(6,*) " "
 write(6,*) " Checks interpolated resummed SFs "
 write(6,*) " "

open(unit=11,status="unknown",file="kfactres-ggnlo-interp.res")
nk=10000
cc_piece = CC_REAL ! Real radiation piece
do k=1,nk

   ! Select the point in x
   x = 1e-8 * ( 1 / 1e-8 )**(dfloat(k)/nk)
   y=log(1.0_dp/x)
   
   write(11,"(E12.7,1x,4(f10.5,1x))") 1.0_dp/x,sf_Pres_gg(y),sf_Pres_qg(y), &
        & sf_Pres_qq(y),sf_Pres_gq(y)
   
enddo
close(11)


!
!----------------------------
!
! Define splitting function object with the interpolated
! splitting function and check its moments

write(6,*) " "
write(6,*) " Checks MSR with interpolated splitting functions "
write(6,*) " "
  
! set up the grid itself -- we use 4 nested subgrids
! order = 6
! ymax = 16.0_dp
! dy   = 0.05_dp
! call InitGridDef(gdarray(4),dy/27.0_dp,0.2_dp, order=order)
! call InitGridDef(gdarray(3),dy/9.0_dp,0.5_dp, order=order)
! call InitGridDef(gdarray(2),dy/3.0_dp,2.0_dp, order=order)
! call InitGridDef(gdarray(1),dy,       ymax  ,order=order)
! call InitGridDef(grid,gdarray(1:4),locked=.true.)

! Single fine grid (problem with xvals with nested grids)
 order = 6
 ymax = 20_dp
 dy   = 0.1_dp
call InitGridDef(grid,dy,       ymax  ,order=order)

! Initialize DGLAP holder
nloop=1 ! L0
nf_lcl=4 ! nf active
call qcd_SetNf(nf_lcl)
call InitDglapHolder(grid, dglap_h, factscheme, nloop)

! Allocation
call AllocGridConv(grid,xPgg)
call AllocGridConv(grid,xPqg)
call AllocGridConv(grid,xPgq)
call AllocGridConv(grid,xPqq)
call AllocSplitMat(grid,P_leading,nf_lcl)

call AllocGridQuant(grid,xpdf)
call AllocGridQuant(grid,xpdf_conv_sr1)
call AllocGridQuant(grid,xpdf_conv_sr2)
call AllocGridQuant(grid,xpdf_conv_qq)
call AllocGridQuant(grid,xpdf_conv_gg)
call AllocGridQuant(grid,xpdf_conv_qg)
call AllocGridQuant(grid,xpdf_conv_gq)


! Get x values
call AllocGridQuant(grid,xvals)
xvals = xValues(grid)
! Suitable function for Mellin moments
xpdf(:) = 1.0_dp/xvals(:)

! Strange values from xvals if grid is used !!!
! Problems with xvals with nested grids?
! call AllocGridQuant(grid,xvals)
! xvals = xValues(grid)
nx=size(xvals)



write(6,*) " "
write(6,*) " Compute LO an dim and Momentum sum rules "
write(6,*) " "

!----------------------------------
!
! Initialize LO splittingh functions
call InitSplitMat( P_leading, dglap_h%P_LO )
! The two calls should lead to identical results ???
call InitGridConv(grid,xPgg, sf_Pgg )
!call InitGridConv(grid,xPgg, P_leading%gg )
call InitGridConv(grid,xPqg, sf_Pqg )
! Add missing 2*nf term
nf_lcl_dp = nf_lcl
call Multiply( xPqg, 2.0_dp * nf_lcl_dp ) 
call InitGridConv(grid,xPgq, sf_Pgq )
call InitGridConv(grid,xPqq, sf_Pqq )
!
!----------------------------------

! Check LO momentum sum rule and N=2 anomalous dimensions
xpdf_conv_sr2 = (xPgg.conv.xpdf) + (xPqg.conv.xpdf)
xpdf_conv_sr1 = (xPqq.conv.xpdf) + (xPgq.conv.xpdf)
xpdf_conv_qq = (xPqq.conv.xpdf) 
xpdf_conv_gg = (xPgg.conv.xpdf) 
xpdf_conv_qg = (xPqg.conv.xpdf) 
xpdf_conv_gq = (xPgq.conv.xpdf) 

N=(2,0)
CALL evolinit  
call ANDIM_LO(N,nf_lcl,P0NS,P0SG)
!write(6,*) real(P0NS)/2_dp

! Why not exact result for gg?

do ix=nx-1,nx-1

   y = log(1.0_dp/xvals(ix))
   write(6,'(i4,1x,E10.2,4(1x,f8.5,1x,f8.5))') &
        & ix,xvals(ix), &
        & xvals(ix)*EvalGridQuant(grid,xpdf_conv_qq,y), real(P0SG(1,1))/2_dp,  &
        & xvals(ix)*EvalGridQuant(grid,xpdf_conv_gg,y), real(P0SG(2,2))/2_dp,  &
        & xvals(ix)*EvalGridQuant(grid,xpdf_conv_gq,y), &
        & real(P0SG(2,1))/2_dp,  &
        & xvals(ix)*EvalGridQuant(grid,xpdf_conv_qg,y),&
        & real(P0SG(1,2))/2_dp 

   write(6,'(i4,1x,E10.2,2(1x,f8.5,1x,f8.5))') &
        & ix,xvals(ix), &
        & xvals(ix)*EvalGridQuant(grid,xpdf_conv_sr1,y), &
        & real(P0SG(1,1))/2_dp + real(P0SG(2,1))/2_dp,  &
        & xvals(ix)*EvalGridQuant(grid,xpdf_conv_sr2,y), &
        & real(P0SG(2,2))/2_dp + real(P0SG(1,2))/2_dp 

enddo


write(6,*) " "
write(6,*) " Compute NLO an dim and Momentum sum rules "
write(6,*) " "

!----------------------------------
!
! Initialize NLO splitting functions
call InitGridConv(grid,xPgg, sf_P1gg )
call InitGridConv(grid,xPqg, sf_P1qg )
call Multiply( xPqg, 2.0_dp * nf_lcl_dp ) 
!
!----------------------------------

! Check NLO momentum sum rule and N=2 anomalous dimensions
xpdf_conv_sr2 = (xPgg.conv.xpdf) + (xPqg.conv.xpdf)
!xpdf_conv_sr1 = (xPqq.conv.xpdf) + (xPgq.conv.xpdf)
!xpdf_conv_qq = (xPqq.conv.xpdf) 
xpdf_conv_gg = (xPgg.conv.xpdf) 
xpdf_conv_qg = (xPqg.conv.xpdf) 
!xpdf_conv_gq = (xPgq.conv.xpdf) 

N=(2,0)  
call ANDIM_NLO(N,nf_lcl,P0NS,P0SG)

open(unit=10,status="unknown",file="msr-nlo.res")
do ix=1,nx-1

   y = log(1.0_dp/xvals(ix))
  
   write(10,'(i4,1x,E10.2,2(1x,f10.5,1x,f8.5))') &
        & ix,1.0_dp/xvals(ix), &
        & xvals(ix)*EvalGridQuant(grid,xpdf_conv_sr1,y), &
        & real(P0SG(1,1))/2_dp + real(P0SG(2,1))/2_dp,  &
        & xvals(ix)*EvalGridQuant(grid,xpdf_conv_sr2,y), &
        & real(P0SG(2,2))/2_dp + real(P0SG(1,2))/2_dp 

enddo
close(10)


write(6,'(i4,1x,E10.2,2(1x,f8.5,1x,f10.5))') &
     & nx-1,xvals(nx-1), &
     & xvals(nx-1)*EvalGridQuant(grid,xpdf_conv_gg,y), real(P0SG(2,2))/4_dp, &
     & xvals(nx-1)*EvalGridQuant(grid,xpdf_conv_qg,y),&
     & real(P0SG(1,2))/4_dp 

write(6,'(i4,1x,E10.2,2(1x,f8.5,1x,f8.5))') &
     & nx-1,xvals(nx-1), &
     & xvals(nx-1)*EvalGridQuant(grid,xpdf_conv_sr1,y), &
     & real(P0SG(1,1))/2_dp + real(P0SG(2,1))/2_dp,  &
     & xvals(nx-1)*EvalGridQuant(grid,xpdf_conv_sr2,y), &
     & real(P0SG(2,2))/2_dp + real(P0SG(1,2))/2_dp 



write(6,*) " "
write(6,*) " Compute NLOres  Momentum sum rule "
write(6,*) " "

!----------------------------------
!
! Initialize NLOres splitting functions
call InitGridConv(grid,xPgg, sf_Pres_gg )
call InitGridConv(grid,xPqg, sf_Pres_qg )
call InitGridConv(grid,xPqq, sf_Pres_qq )
call InitGridConv(grid,xPgq, sf_Pres_gq )
!----------------------------------

! Check NLOres momentum sum rule 
xpdf_conv_sr2 = (xPgg.conv.xpdf) + (xPqg.conv.xpdf)
xpdf_conv_sr1 = (xPqq.conv.xpdf) + (xPgq.conv.xpdf)
xpdf_conv_qq = (xPqq.conv.xpdf) 
xpdf_conv_gg = (xPgg.conv.xpdf) 
xpdf_conv_qg = (xPqg.conv.xpdf) 
xpdf_conv_gq = (xPgq.conv.xpdf) 

open(unit=10,status="unknown",file="msr-nlores.res")
do ix=1,nx-1

   y = log(1.0_dp/xvals(ix))
   
   write(10,'(i4,1x,E10.2,2(1x,f8.5,1x,f8.5))') &
        & ix,1.0_dp/xvals(ix), &
        & xvals(ix)*EvalGridQuant(grid,xpdf_conv_sr1,y), &
        & real(P0SG(1,1))/2_dp + real(P0SG(2,1))/2_dp, & 
        & xvals(ix)*EvalGridQuant(grid,xpdf_conv_sr2,y), &
        & real(P0SG(2,2))/2_dp + real(P0SG(1,2))/2_dp 

enddo
close(10)

ix=nx-1
write(6,'(i4,1x,E10.2,2(1x,f8.5,1x,f8.5))') &
     & ix,xvals(ix), &
      & xvals(ix)*EvalGridQuant(grid,xpdf_conv_sr1,y), &
      & real(P0SG(1,1))/2_dp + real(P0SG(2,1))/2_dp ,&
     & xvals(ix)*EvalGridQuant(grid,xpdf_conv_sr2,y), &
     & real(P0SG(2,2))/2_dp + real(P0SG(1,2))/2_dp 


end program resummation_example


