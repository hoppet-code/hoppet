
! Common location for an user-defined table
module cteq_table_module

   use hoppet_v1
   implicit none
   type(pdf_table) :: tab
   type(dglap_holder) :: dh
   type(running_coupling) :: coupling
   type(grid_def) :: grid

 end module cteq_table_module

program test_hoppet
! Use module in which the user-defined table is defined
use cteq_table_module
! Use the main module
use hoppet_v1

implicit none 

real(dp) :: pdf0(-6,6)
real(dp) :: Q0

integer nloop, factscheme,nf_lcl
   
call InitGridDef(grid,dy=0.2_dp,ymax=8.0_dp,order=2)
   
! Initialize DGLAP holder
nloop=3
factscheme = factscheme_MSbar
call qcd_SetNf(4)
call InitDglapHolder(grid,dh,factscheme,nloop)

Q0=1.41
! Preevolution
call PreEvolvePdfTable(tab,Q0,dh,coupling)
! Evolution
call EvolvePdfTable(tab,pdf0)





end program test_hoppet

program test_hoppet
! Use module in which the user-defined table is defined
use cteq_table_module
! Use the main module
use hoppet_v1

implicit none 


call Subrout1
! Evolution
call SubroutB





end program test_hoppet

!-----------------------------------------------------------


