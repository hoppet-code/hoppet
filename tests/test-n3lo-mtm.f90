! program to test our incorporation of the N3LO mass-threshold matching
! code from Johannes Blumlein
program test_blumlein
  use types; use consts_dp
  use blumlein_interfaces
  use mass_thresholds_n3lo
  use vfnsmod
  use daind_fake
  use convolution_communicator
  use splitting_functions
  !use A3PSqq_He
  !use A3Sqg_He
  !use A3PSqq_Hp
  !use A3Sqg_Hp
  use io_utils
  use streamlined_interface
  use dglap_choices
  use vfnsmod
  implicit none
  real(dp) :: GGPLU_a2, GGPLU_a3
  real(dp) :: a4pi(2)= (/0.1_dp, 1.0_dp/), ggplu_arr(2), coeffs(2:3), a2pi
  real(dp) :: x = 0.1_dp, lnmQmu2 = zero, res_aps1, res_aps2, res1, res2, res3
  integer  :: i, momN
  integer  :: nflcl = 4, nf_light, nf_heavy
  real(dp) :: nf_light_dp, nf_heavy_dp
  type(comparator) :: comparators(11)
  type(JB) :: A2Sgg_H
  type(JB) :: A3Sgg_H
  type(JB) :: A3NSm_H
  type(grid_conv) :: gc_A2Sgg_H_JB, gc_A2Sgg_H_HP, gc_A3Sgg_H, gc_A3NSm_H
  real(dp), pointer :: powx(:), conv_res(:)
  !--------------------------
  real(dp) :: dy
  integer  :: nloop
  real(dp) :: ymax, Qmin, Qmax, dlnlnQ
  integer  :: order
  real(dp) :: PSL_res(2:3)
  integer  :: idev
  dy = 0.05_dp
  nloop = 4
  ymax = dble_val_opt("-ymax",22.0d0) ! go to quite high y, to avoid finite-y effects in moments
  Qmin = 1.0d0
  Qmax = 28000d0 ! twice LHC c.o.m.
  dlnlnQ = dy/4.0_dp  ! min(0.5_dp*dy,0.07_dp)
  order = -6

  idev = idev_open_opt("-o")
  x = dble_val_opt("-x", x)
  a4pi(1) = dble_val_opt("-a4pi", a4pi(1))
  nflcl = int_val_opt("-nf", nflcl)
  if (.not.CheckAllArgsUsed(0)) stop

  call dglap_Set_nnlo_nfthreshold(nnlo_nfthreshold_exact)
  call hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order,factscheme_MSbar)

  !write(6,*) AGGPLU(0.1_dp, 4, 1.0_dp, lnmQmu2)
  !write(6,*) AGGPLU(0.1_dp, 4, 2.0_dp, lnmQmu2) / 8.0_dp
  

  !call hoppetStart(0.05_dp, 4)
  ! make sure SetNf comes after hoppetStart (which resets it)
  call qcd_SetNf(nflcl)
  call SetNfDglapHolder(dh, nflcl)
  nf_heavy = nf_int
  nf_light = nf_heavy - 1
  nf_light_dp = real(nf_light, dp)
  nf_heavy_dp = real(nf_heavy, dp)

  write(6,*) "-------------- Running tests that came with VFNS.f ---------------"
  ! The moment tests are designed to work with the following parameters
  ! nf,as,LL,   3.0000000000000000       0.10000000000000001        10.000000000000000
  call VFNS3sub(z_in=0.1_dp, nflightdp_in=3.0_dp, a4pi_in=0.1_dp, LL_in=10.0_dp)

  write(6,*) "-------------- comparing Bluemlein α^2 with Hoppet's α^2 ---------------"
  write(6,*) "nf_light = ", nf_light, " nf_heavy = ", nf_heavy
  comparators( 1) = comparator("GGPLU", GGPLU, sf_A2Sgg_H , cc_VIRT    )
  comparators( 2) = comparator("GGREG", GGREG, sf_A2Sgg_H , cc_REALVIRT)
  comparators( 3) = comparator("GGDEL", GGDEL, sf_A2Sgg_H , cc_DELTA   )
  comparators( 4) = comparator("NSPLU", NSPLU, sf_A2NSqq_H, cc_VIRT    )
  comparators( 5) = comparator("NSREG", NSREG, sf_A2NSqq_H, cc_REALVIRT)
  comparators( 6) = comparator("NSDEL", NSDEL, sf_A2NSqq_H, cc_DELTA   )

  comparators( 7) = comparator("QGL  ", QGL  , null()     , CC_REAL    )
  comparators( 8) = comparator("PSL  ", PSL  , null()     , CC_REAL    )
  comparators( 9) = comparator("GQ   ", GQ   , sf_A2Sgq_H , CC_REAL    )
  comparators(10) = comparator("QG   ", QG   , sf_A2PShg  , CC_REAL    )
  comparators(11) = comparator("PS   ", PS   , sf_A2PShq  , CC_REAL    )

  write(6,'(a)') "   piece     x        Blumlein's code           hoppet code                 difference"
  do i = 1, size(comparators)
    call comparators(i)%compare(x)
  end do

  !write(6,*) separate_orders(x, nf_light, QGL), AQG(x, nf_light, one, zero)
  write(6,*) "AQG = ", AQG(x, nf_light_dp, a4pi(1), zero)
  !write(6,*) A3Sqg_H_exact(x, nf_light)
  write(6,*) PSL(x, nf_light_dp, a4pi(1), zero)

  write(6,*) "Results including the a^3 factor:"
  res_aps1 = APS1(x, nf_light_dp, a4pi(1), 0.0_dp)
  res_aps2 = APS2(x, nf_light_dp, a4pi(1), zero)
  write(6,*) "PS separate orders ", separate_orders(x, nf_light, PS) * a4pi(1)**(/2,3/)
  write(6,*) "APS1, APS2: ", res_aps1, res_aps2
  !write(6,*) "A3PSqq_H_exact: ", A3PSqq_H_exact(x, nf_light) * a4pi(1)**3

  ! AK additions:
  ! I think PSL corresponds to PSqq_H and QGL to Sqg_H
  ! So this should be the comparison
  write(6,*) '----------------Comparison between AKs functions and Blumleins------------'
  write(6,*) "PSL separate orders ", separate_orders(x, nf_light, PSL) * a4pi(1)**(/2,3/)
  !write(6,*) "A3PSqq_H_exact: ", A3PSqq_H_exact(x, nf_light) * a4pi(1)**3
  !write(6,*) "A3PSqq_H_exact: ", A3PSqq_H_param(x, nf_light) * a4pi(1)**3
  write(6,*) "QGL separate orders ", separate_orders(x, nf_light, QGL) * a4pi(1)**(/2,3/)
  !write(6,*) "A3Sqg_H_exact: ", A3Sqg_H_exact(x, nf_light) * a4pi(1)**3
  !write(6,*) "A3Sqg_H_param: ", A3Sqg_H_param(x, nf_light) * a4pi(1)**3
  ! Here we print out the functions, without the a2pi pieces
  do i = 1, 100
     x = 0.01_dp ** (i/100.0_dp)
!     PSL_res = separate_orders(x, nf_light, QG)
!     write(idev,*) x, PSL_res(2), PSL_res(3) + AQGtotal(x, nf_light_dp, a4pi(2), zero)
     PSL_res = separate_orders(x, nf_light, GGDEL)
     write(idev,*) x, PSL_res(2), PSL_res(3) + AGGDEL(x, nf_light_dp, a4pi(2), zero)
  end do
  write(6,*) '--------------------------------------------------------------------------'
  
  write(6,*) "AGGREG0, AGGREG1: ", AGGREG0(x, nf_light_dp, a4pi(1), zero), AGGREG1(x, nf_light_dp, a4pi(1), zero)

  A3Sgg_H = JB(GGREG, GGPLU, GGDEL, 3) !!<< this is not the full answer, but the part without the A
  A2Sgg_H = JB(GGREG, GGPLU, GGDEL, order=2)
  cc_piece = cc_DELTA
  res1 = A2Sgg_H%integrand(log(1/x))
  res2 = sf_A2Sgg_H(log(1/x))
  write(6,*) "A2Sgg_H: ", res1, res2, res1-res2

  
  !call InitGridConv(grid, gc_A2Sgg_H_JB, A2Sgg_H)
  !call InitGridConv(grid, gc_A2Sgg_H_HP, sf_A2Sgg_H)
  !call InitGridConv(grid, gc_A2Sgg_H_JB, JB(GGREG,  GGPLU,  GGDEL,  order=2))
  !call InitGridConv(grid, gc_A3Sgg_H   , JB(GGREG,  GGPLU,  GGDEL,  order=3))
  !call AddWithCoeff(      gc_A3Sgg_H   , JB(AGGREG, AGGPLU, AGGDEL, order=3))
  ! call InitGridConv(grid, gc_A2Sgg_H_JB, JB(null() , null(),  GGDEL,  order=2))
  ! call InitGridConv(grid, gc_A3Sgg_H   , JB(null() , null(),  GGDEL,  order=3))
  ! call AddWithCoeff(      gc_A3Sgg_H   , JB(null() , null(), AGGDEL,  order=3))
  ! temporarily misuse some variable names
  !call InitGridConv(grid, gc_A2Sgg_H_JB, JB(reg=QG, a=AQG, order=2))
  !call InitGridConv(grid, gc_A3Sgg_H   , JB(reg=QG, a=AQG, order=3))

  write(6,*) "------------------------- preparing MTM_N3LO -------------------------"
  !call InitMTMN3LO(grid, dh%MTM_N3LO)
  write(6,*) "------------------------- done preparing MTM_N3LO -------------------------"

  call AllocGridQuant(grid, powx)
  call AllocGridQuant(grid, conv_res)
  a2pi = a4pi(1)*2

  write(6,*) "----------------------- AQG(z) -----------------------"
  do i = 1, 100
    x = 0.01_dp ** (i/100.0_dp)
    write(6,*) x, AQGtotal(x, nf_light_dp, a4pi(1), zero), "AQGtotal"
  end do

  ! do non-singlet moment with N=1 (other moments are infinite)
  write(6,*) "------------------------------------------------"
  momN = 1
  write(6,*) "momN = ", momN
  call compare_mom("NSqq_H", momN, dh%MTM_NNLO%NSqq_H, dh%MTM_N3LO%NSqq_H, MOMNS , .true.)

  do momN = 2, 6, 2

    write(6,*) "------------------------------------------------"
    write(6,*) "momN = ", momN

    call compare_mom("PShq  ", momN, dh%MTM_NNLO%PShq  , dh%MTM_N3LO%PShq  , MOMPS , .false.,MOMAPS)
    call compare_mom("PShg  ", momN, dh%MTM_NNLO%PShg  , dh%MTM_N3LO%PShg  , MOMQG , .false.)
    call compare_mom("PSqq_H", momN, dh%MTM_NNLO%PSqq_H, dh%MTM_N3LO%PSqq_H, MOMPSL, .false.)
    call compare_mom("NSqq_H", momN, dh%MTM_NNLO%NSqq_H, dh%MTM_N3LO%NSqq_H, MOMNS , .true.)
    call compare_mom("Sgg_H ", momN, dh%MTM_NNLO%Sgg_H , dh%MTM_N3LO%Sgg_H , MOMGG , .true.)
    call compare_mom("Sgq_H ", momN, dh%MTM_NNLO%Sgq_H , dh%MTM_N3LO%Sgq_H , MOMGQ , .false.,MOMAGQ)
    call compare_mom("Sqg_H ", momN, dh%MTM_NNLO%Sqg_H , dh%MTM_N3LO%Sqg_H , MOMQGL, .false.)
  end do

  ! set up a non-singlet object
  A3NSm_H = JB(OREG_znfasLL, OPLUS_znfasLL, ODEL_znfasLL, 203)
  call InitGridConv(grid, gc_A3NSm_H, A3NSm_H)
  do momN = 1, 7, 2
    write(6,*) "Odd moment momN = ", momN
    !                                this is dummy=0
    if (momN/=1) call compare_mom("NSm_H(hard-coded)", momN, dh%MTM_NNLO%NSqq_h, gc_A3NSm_H, ODDMOMNS_NxnfasLL , .true.)
    call compare_mom("NSm_H(integrated)", momN, dh%MTM_NNLO%NSqq_h, gc_A3NSm_H, ODDMOMNSint_NxnfasLL , .true.)
  end do


contains

  function ODDMOMNS_NxnfasLL(momN, x, nf_light_dp, a4pi, LL) result(res)
    integer   :: momN
    real(dp)  :: x, nf_light_dp, a4pi, LL
    !------------------------------
    real(dp) :: res
    res = ODDMOMNS(momN,a4pi,LL,nf_light_dp)
  end function ODDMOMNS_NxnfasLL

  function ODDMOMNSint_NxnfasLL(momN, x, nf_light_dp, a4pi, LL) result(res)
    integer   :: momN
    real(dp)  :: x, nf_light_dp, a4pi, LL
    !------------------------------
    real(dp) :: res
    res = MOMREG(a4pi,LL,nf_light_dp,momN) + ODEL(a4pi,LL,nf_light_dp)
    if (momN /= 1) res = res + MOMPLUS(a4pi,LL,nf_light_dp,momN)
  end function ODDMOMNSint_NxnfasLL


  subroutine compare_mom(name, momN, MTM2, MTM3, momfn, delta_one, momfn2)
    character(len=*), intent(in) :: name
    integer,          intent(in) :: momN
    type(grid_conv),  intent(in) :: MTM2, MTM3
    procedure(JBmom)             :: momfn
    logical,          intent(in) :: delta_one
    procedure(JBmom), optional   :: momfn2
    !------------------------------
    real(dp) :: res2,res3,HPsum, JBsums(1:2),JBcoeffs_a2pi(2:3), LL
    integer  :: i
    character(len=5) :: bad
    real(dp), parameter :: tolerance = 1e-4

    LL = zero
    do i = 1,2
      JBsums(i) = momfn(momN, x, nf_light_dp, a4pi(i), LL)
      if (present(momfn2)) JBsums(i) = JBsums(i) + momfn2(momN, x, nf_light_dp, a4pi(i), LL)
    end do
    JBcoeffs_a2pi = separate23(a4pi, JBsums, delta_one) / two**[2,3]

    res2 = gc_moment(MTM2, one * momN)
    res3 = gc_moment(MTM3, one * momN)
    HPsum = merge(1,0,delta_one) + a2pi**2*res2 + a2pi**3 * res3
    bad = ""
    if (any(abs(JBcoeffs_a2pi - (/res2,res3/)) > tolerance)) bad = "*BAD*"
    if (abs(HPsum - JBsums(1))                 > tolerance)  bad = "*BAD*"

    write(6,*) name, " JB : res2,3,sum:", JBcoeffs_a2pi, JBsums(1)
    write(6,*) name, " MTM: res2,3,sum:", res2, res3, HPsum, bad
  end subroutine compare_mom
end program test_blumlein
