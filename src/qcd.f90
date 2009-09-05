! $Id: qcd.f90,v 1.13 2004/09/18 14:39:29 salam Exp $
module qcd
  use types; use consts_dp

  integer,  parameter :: nf_def = 5
  real(dp), parameter :: ca_def = 3, cf_def = four/three, tr_def = half
  real(dp), parameter :: tf_def = tr_def * nf_def

  !-- the following are all modifiable, but have default values
  real(dp), public :: ca = ca_def
  real(dp), public :: cf = cf_def
  real(dp), public :: tr = tr_def

  integer,  public :: nf_int = nf_def
  integer,  public :: nf_u = nf_def/2
  integer,  public :: nf_d = (nf_def+1)/2
  real(dp), public :: nf = nf_def, tf = tf_def
  ! HACK TO GET GLUON + SINGLET EVOLUTION in 0 & 1
  ! For it to work, quark component in 2 must be zero.
  !integer, parameter, public :: nf_u = nf_def, nf_d = 0

  ! beta2 is from Tarasov JINR P2-82-900 (Dubna 1982)
  ! Larin & Vermaseren NIKHEF-H/92-17 hep-ph/9302208
  real(dp), public :: beta0 = (11*ca_def - four*tf_def)/(12*pi)
  real(dp), public :: twopi_beta0 = (11*ca_def - four*tf_def)/6.0_dp
  real(dp), public :: beta1 = (17*ca_def**2 - tf_def*(10*ca_def + 6*cf_def))&
       &                        / (24*pisq)
  real(dp), public :: beta2 = (2857*ca_def**3  &
       & + (54*cf_def**2-615*cf_def*ca_def-1415*ca_def**2)*two*tf_def&
       & + (66*cf_def + 79*ca_def)*(two*tf_def)**2)/(3456*pi**3)

  !--- nomenclature is: 
  !  alpha_s(nf+1,mu) = alpha_s(nf,mu) + 
  !               sum_{mn} alphastep_mn * (alphas)^(m+1) ln^n(mu^2/m(mu^2)^2)
  !
  ! taken from hep-ph/9411260; should also refer to 
  ! Bernreuther, Wetzel, NPB197 (1982) 228
  ! Bernreuther, Annals of Phys. 151 (1983) 127
  !
  ! Only use this if (mass_steps_on)?
  logical,  public :: mass_steps_on = .true.
  !logical,  public, parameter :: mass_steps_on = .false.
  !---- the following is for the MSbar mass
  ! taken from hep-ph/9411260; should also refer to 
  ! Bernreuther, Wetzel, NPB197 (1982) 228
  ! Bernreuther, Annals of Phys. 151 (1983) 127
  real(dp), public :: alphastep11 = four*tr_def/(12*pi)
  real(dp), public :: alphastep22 = (four*tr_def/(12*pi))**2
  real(dp), public :: alphastep21 = tr_def*(10*ca_def + 6*cf_def)/(24*pisq)
  real(dp), public :: alphastep20_msbar=&
       &(13.0_dp/48.0_dp*cf_def-two/9.0_dp*ca_def)&
       &                            *tr_def/pisq
  !-- for the pole mass, take this from hep-ph/9706430 
  !   (Chetyrkin, Kniehl and Steinhauser), PRL 79 (1997) 2184
  !   though it is not the original referece, just a "container"
  !-- expression not known for variable colour factors?
  !   NB: they express nf-1 in terms of nf...
  !real(dp), public :: alphastep20_pole  = -11.0_dp/72.0_dp/pisq
  !real(dp), public :: alphastep20_pole  =   7.0_dp/24.0_dp/pisq
  real(dp), public :: alphastep20_pole  = &
       &(15.0_dp/16.0_dp*cf_def-two/9.0_dp*ca_def) * tr_def/pisq

  real(dp), public :: cmw_K = ca_def*(67.0_dp/18.0_dp - pisq/6.0_dp) &
       &                      - tf_def * 10.0_dp/9.0_dp
  
  !!!! Taken from Moch, Vermaseren & Vogt: NB TR dependence not in...
  real(dp), parameter :: cmw_K2_def = &
       &  ca_def**2 * ( 245._dp/24._dp - 67._dp/9._dp*zeta2 &
       &             + 11.0_dp/6._dp * zeta3 + 11.0_dp/5._dp * zeta2**2)&
       &+ two*cf_def*tf_def * (-55._dp/24._dp + 2*zeta3)&
       &+ two*ca_def*tf_def * (-209._dp/108._dp + 10._dp/9._dp*zeta2 &
       &                   - 7._dp/3._dp * zeta3)&
       &- four*tf_def**2 /27._dp
  real(dp), public :: cmw_K2  = cmw_K2_def
  real(dp), public :: mvv_A3  = 16*cf_def*cmw_K2_def
  real(dp), public :: mvv_A3G = 16*ca_def*cmw_K2_def


  !--- is it useful to have a fake upper entry?
  ! put charm quark mass just above sqrt(2)
  real(dp), public, parameter :: &
       &        quark_masses_def(6) = (/1e-10_dp,1e-10_dp,0.15_dp,&
       &                           1.414213563_dp, 4.5_dp, 175.0_dp/)


  public :: qcd_SetNf, qcd_SetGroup!, qcd_SetVogtImod
  
contains

  !----------------------------------------------------------------------
  subroutine qcd_SetNf(nf_in)
    integer :: nf_in
    nf_int = nf_in
    nf     = nf_int
    nf_u   = nf_int/2
    nf_d   = (nf_int+1)/2
    tf     = tr * nf
    call qcd_set_beta0
  end subroutine qcd_SetNf
  
  !----------------------------------------------------------------------
  subroutine qcd_SetGroup(ca_in,cf_in,tr_in)
    real(dp), intent(in) :: ca_in,cf_in,tr_in
    ca = ca_in
    cf = cf_in
    tr = tr_in
    tf = nf * tr
    call qcd_set_beta0
  end subroutine qcd_SetGroup

  !----------------------------------------------------------------------
  subroutine qcd_set_beta0
    beta0       = (11*ca - four*tf)/(12*pi)
    twopi_beta0 = twopi * beta0
    beta1       = (17*ca**2 - tf*(10*ca + 6*cf)) / (24*pisq)
    beta2       = (2857*ca**3 + (54*cf**2-615*CF*CA-1415*ca**2)*two*tf&
         &           + (66*cf + 79*ca)*(two*tf)**2)/(3456*pi**3)
    !-- matching thresholds -----------------
    alphastep11 = four*tr/(12*pi)
    alphastep22 = (four*tr/(12*pi))**2
    alphastep21 = tr*(10*ca + 6*cf)/(24*pisq)
    alphastep20_msbar = (13.0_dp/48.0_dp*cf-two/9.0_dp*ca)*tr/pisq    
    !  -11.0_dp/72.0_dp/pisq    
    !-- expression not know for variable colour factors
    !alphastep20_pole  = 7.0_dp/24.0_dp/pisq
    alphastep20_pole  = (15.0_dp/16.0_dp*cf-two/9.0_dp*ca)*tr/pisq
    cmw_K       = ca*(67.0_dp/18.0_dp - pisq/6.0_dp) - tf*10.0_dp/9.0_dp
    cmw_K2 = &
       &  ca**2 * ( 245._dp/24._dp - 67._dp/9._dp*zeta2 &
       &             + 11.0_dp/6._dp * zeta3 + 11.0_dp/5._dp * zeta2**2)&
       &+ two*cf*tf * (-55._dp/24._dp + 2*zeta3)&
       &+ two*ca*tf * (-209._dp/108._dp + 10._dp/9._dp*zeta2 &
       &                   - 7._dp/3._dp * zeta3)&
       &- four*tf**2 /27._dp
    mvv_A3  = 16*cf*cmw_K2
    mvv_A3G = 16*ca*cmw_K2
  end subroutine qcd_set_beta0
  
!!$  !-------- overkill ----------------------------------------
!!$  subroutine qcd_SetVogtImod(imod)
!!$    integer, intent(in) :: imod
!!$    vogt_imod = imod
!!$  end subroutine qcd_SetVogtImod
  
end module qcd
