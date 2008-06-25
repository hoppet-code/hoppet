
module resum_grids_out
  use types; use consts_dp; use  splitting_functions;
  use convolution_communicator;
  use resum_kfact;
  implicit none

  public :: resum_splitfun_grid

contains

 subroutine resum_splitfun_grid(xPgg_resum,xPqq_resum,xPgq_resum,xPqg_resum) 

   real(dp), intent(out)  :: xPgg_resum(nalphas,nxgrid,4)
   real(dp), intent(out)  :: xPqq_resum(nalphas,nxgrid,4)
   real(dp), intent(out)  :: xPqg_resum(nalphas,nxgrid,4)
   real(dp), intent(out)  :: xPgq_resum(nalphas,nxgrid,4)
   
   integer :: ix,ialpha,nf_local
   real(dp) :: y
  

   ! Construct the resummed NLO splitting functions
   ! at the points of the grid where the resummation
   ! additive corrections have been tabulated

   ! Real part of splitting functions
   CC_PIECE = cc_REAL
   nf_local=4

   do ialpha=1,nalphas

      do ix=1,nxgrid

         !Change to variable y
         y=log( 1_dp / xgrid(ix) )

         ! Construct resummed splitting functions
         ! HOPPET return xP
         ! gg
         xPgg_resum(ialpha,ix,1) = alphas(ialpha)*sf_Pgg(y)/(2_dp*pi) &
          & + ( alphas(ialpha)/(2_dp*pi) )**2_dp * sf_P1gg(y) &
          & + alphas(ialpha)* kfact_splresggnlo(ialpha,ix)
         ! qq
         ! Again missing factor of nf
         xPqq_resum(ialpha,ix,1) = alphas(ialpha)*sf_Pqq(y)/(2_dp*pi) &
          & + 2_dp * nf_local * &
          & ( alphas(ialpha)/(2_dp*pi) )**2_dp * sf_P1qqS(y) &
          & + alphas(ialpha)* kfact_splresqqnlo(ialpha,ix)
         ! qg
         ! Missing factor of 2Nf in qg splitting functions
         xPqg_resum(ialpha,ix,1) = 2_dp*nf_local &
           & *alphas(ialpha)*sf_Pqg(y)/(2_dp*pi) &
          & + 2_dp*nf_local &
          & * ( alphas(ialpha)/(2_dp*pi) )**2_dp * sf_P1qg(y) &
          & + alphas(ialpha)* kfact_splresqgnlo(ialpha,ix)
         ! gq
         xPgq_resum(ialpha,ix,1) = alphas(ialpha)*sf_Pgq(y)/(2_dp*pi) &
          & + ( alphas(ialpha)/(2_dp*pi) )**2_dp * sf_P1gq(y) &
          & + alphas(ialpha)* kfact_splresgqnlo(ialpha,ix)

!!$         if(ialpha.eq.25) &
!!$          & write(6,*)alphas(ialpha), xgrid(ix),kfact_splresqqnlo(ialpha,ix)

         ! As a cross-check add also output of
         ! unresummed splitting functions
         ! LO
          xPgg_resum(ialpha,ix,2) = alphas(ialpha)*sf_Pgg(y)/(2_dp*pi)
          xPqq_resum(ialpha,ix,2) = alphas(ialpha)*sf_Pqq(y)/(2_dp*pi)
          xPqg_resum(ialpha,ix,2) = 2_dp*nf_local & 
               & * alphas(ialpha)*sf_Pqg(y)/(2_dp*pi)
          xPgq_resum(ialpha,ix,2) = alphas(ialpha)*sf_Pgq(y)/(2_dp*pi)
         ! NLO
          xPgg_resum(ialpha,ix,3) = alphas(ialpha)*sf_Pgg(y)/(2_dp*pi) &
          & + ( alphas(ialpha)/(2_dp*pi) )**2_dp * sf_P1gg(y)
          ! Note various qq splitting functions at NLO
          ! depending on quark flavour combination
          ! Resummation correction applies to Singlet splitting function
          xPqq_resum(ialpha,ix,3) = alphas(ialpha)*sf_Pqq(y)/(2_dp*pi) &
          & + 2_dp * nf_local * &
          & ( alphas(ialpha)/(2_dp*pi) )**2_dp * sf_P1qqS(y)
          ! Factors of 2NF
          xPqg_resum(ialpha,ix,3) = 2_dp*nf_local &
          & * alphas(ialpha)*sf_Pqg(y)/(2_dp*pi) &
          & + 2_dp*nf_local &
          & * ( alphas(ialpha)/(2_dp*pi) )**2_dp * sf_P1qg(y)
          xPgq_resum(ialpha,ix,3) = alphas(ialpha)*sf_Pgq(y)/(2_dp*pi) &
          & + ( alphas(ialpha)/(2_dp*pi) )**2_dp * sf_P1gq(y)

         !write(6,*) xgrid(ix),alphas(ialpha), Pgg_resum(ialpha,ix)

          ! Differences between resummed and NLO unresummed P's
          ! gg
         xPgg_resum(ialpha,ix,4) = kfact_splresggnlo(ialpha,ix)
         ! qq
         xPqq_resum(ialpha,ix,4) = kfact_splresqqnlo(ialpha,ix)
         ! qg
         xPqg_resum(ialpha,ix,4) = kfact_splresqgnlo(ialpha,ix)
         ! gq
         xPgq_resum(ialpha,ix,4) = kfact_splresgqnlo(ialpha,ix)
       
      enddo
      
   enddo

   end subroutine resum_splitfun_grid

  !-----------------------------------

end module resum_grids_out
