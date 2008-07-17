!======================================================================
!
! Read and save ABF resummed K factors
!
!-----------------------------------------------

module resum_kfact
  use types; use consts_dp
  use convolution_communicator;
  implicit none

  public :: resum_kfact_read

  ! Set parameters of the resummation correction grids
  integer, parameter :: nalphas = 57
  integer, parameter :: nxgrid = 59

  real(dp) :: alphas(nalphas),xgrid(nxgrid)

  real(dp) :: kfact_splresggnlo(nalphas,nxgrid)
  real(dp) :: kfact_splresgqnlo(nalphas,nxgrid)
  real(dp) :: kfact_splresqgnlo(nalphas,nxgrid)
  real(dp) :: kfact_splresqqnlo(nalphas,nxgrid)

  real(dp) :: kfact_coresggnlo(nalphas,nxgrid)
  real(dp) :: kfact_coresgqnlo(nalphas,nxgrid)
  real(dp) :: kfact_coresqgnlo(nalphas,nxgrid)
  real(dp) :: kfact_coresqqnlo(nalphas,nxgrid)


contains
  !----------------------------------------------------------------------
  subroutine resum_kfact_read
    ! Read resummed K factors

    integer :: ix,ialpha
    character*23 ::  filepath

    filepath ="../resummed-tables/nf4/"
  
    ! Corrections for resummed splitting functions
   open(unit=10,status="unknown",&
    & file="../resummed-tables/nf4/splresggnlo.dat")
   open(unit=11,status="unknown",&
    & file="../resummed-tables/nf4/splresqgnlo.dat")
   open(unit=12,status="unknown",&
    & file="../resummed-tables/nf4/splresgqnlo.dat")
   open(unit=13,status="unknown",&
    & file="../resummed-tables/nf4/splresqqnlo.dat")

   ! Corrections for resummed coefficient functions
    open(unit=14,status="unknown",&
         & file="../resummed-tables/nf4/coresggnlo.dat")
    open(unit=15,status="unknown",&
         & file="../resummed-tables/nf4/coresqgnlo.dat")
    open(unit=16,status="unknown",&
         & file="../resummed-tables/nf4/coresgqnlo.dat")
    open(unit=17,status="unknown",&
         & file="../resummed-tables/nf4/coresqqnlo.dat")
    
   do ialpha=1,nalphas

      do ix=1,nxgrid

         read(10,*) alphas(ialpha),xgrid(ix),kfact_splresggnlo(ialpha,ix)
         read(11,*) alphas(ialpha),xgrid(ix),kfact_splresqgnlo(ialpha,ix)
         read(12,*) alphas(ialpha),xgrid(ix),kfact_splresgqnlo(ialpha,ix)
         read(13,*) alphas(ialpha),xgrid(ix),kfact_splresqqnlo(ialpha,ix)

       !  if(ialpha.eq.25) &
       !   & write(*,"(i2,1x,f5.2,1x,i2,1x,f10.8,1x,f6.3)") &
       !    &ialpha,alphas(ialpha),ix,xgrid(ix),kfact_splresggnlo(ialpha,ix)

         read(14,*) alphas(ialpha),xgrid(ix),kfact_coresggnlo(ialpha,ix)
         read(15,*) alphas(ialpha),xgrid(ix),kfact_coresqgnlo(ialpha,ix)
         read(16,*) alphas(ialpha),xgrid(ix),kfact_coresgqnlo(ialpha,ix)
         read(17,*) alphas(ialpha),xgrid(ix),kfact_coresqqnlo(ialpha,ix)


      enddo
      
   enddo
   write(6,*) "Reading grids done!"

   close(10)
   close(11)
   close(12)
   close(13)

    
  end subroutine resum_kfact_read

!-------------------------------------


end module resum_kfact


