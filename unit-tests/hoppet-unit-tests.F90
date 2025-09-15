!! a block of code to help with timing
!! It will run only if an existing variable timing_name coincides
!! with the NAME that is set
#define TIME_START(NAME,NREP) block;real(dp)::t_start,t_end;integer::irep,nrep=(NREP);\
character(len=*),parameter::name=NAME
#define TIME_LOOP  if (trim(timing_name)==name) then; call cpu_time(t_start); do irep = 1, nrep
! OUT is there to make sure we accumulate some kind of result and print
! is, so the result is not optimized away
#define TIME_END(OUT) end do;call cpu_time(t_end);print *,"TIME:",name,": ",(t_end-t_start)/nrep*1e9_dp,"ns",OUT;end if;end block


module hoppet_unit_tests_setup
  implicit none
  
contains
  subroutine hoppet_setup()
    use streamlined_interface
    use pdfs_for_benchmarks

    real(dp) :: dy, ymax, dlnlnQ, Qmin, Qmax, muR_Q
    real(dp) :: asQ, Q0alphas, Q0pdf
    real(dp) :: mc,mb,mt
    integer  :: order, nloop

    order = -6
    ymax  = 12.0_dp
    dy    = 0.1_dp
    ! and parameters for Q tabulation
    Qmin=1.0_dp
    Qmax=28000.0_dp
    dlnlnQ = dy/4.0_dp
    ! and number of loops to initialise!
    nloop = 3
    call hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,&
       &         order,factscheme_MSbar)

    ! Set heavy flavour scheme
    mc = 1.414213563_dp   ! sqrt(2.0_dp) + epsilon
    mb = 4.5_dp
    mt = 175.0_dp
    call hoppetSetVFN(mc, mb, mt)

    ! Streamlined evolution
    
    ! Set parameters of running coupling
    asQ = 0.35_dp
    Q0alphas = sqrt(2.0_dp)
    muR_Q = 1.0_dp
    Q0pdf = sqrt(2.0_dp) ! The initial evolution scale
    call hoppetEvolve(asQ, Q0alphas, nloop,muR_Q, benchmark_pdf_unpolarized_lha, Q0pdf)
  end subroutine hoppet_setup

end module hoppet_unit_tests_setup

!!!!!!!!!!!!!!!!!!!!!!!
program hoppet_unit_tests
  use types
  use unit_tests
  use hoppet_unit_tests_setup
  use interpolation
  use interpolation_coeffs
  implicit none

  ! decide which timing variant to run
  character(len=100):: timing_name = "fill_interp_weights6:x4"


  print '(a)', "Running HOPPET unit tests"
  call hoppet_setup()

  call test_interpolation_coeffs()
  call test_tab_eval()

  if (unit_test_failures > 0) then
    ! print a message in red
    print '(a)', red//trim(to_string(unit_test_failures))//' unit tests failed.'//reset
    stop 1
  else
    ! print a message in green
    print '(a)', green//'All '//trim(to_string(unit_test_successes))//' unit tests passed.'//reset
  end if

contains

  !! checks equivalence between interpolation weights from
  !! uniform_interpolation_weights and fill_interp_weightsN
  subroutine test_interpolation_coeffs()
    real(dp) :: x
    real(dp) :: weights(0:4)
    integer   :: i, ix
    integer, parameter :: nlo = 1, nhi = 6
    real(dp) :: weights1(0:nhi), weights2(0:nhi)
    real(dp), parameter :: xvals(*) = [0.2_dp,1.5_dp,3.4_dp, 1.0_dp]

    do ix = 1, size(xvals)
      x = xvals(ix)
      do i = nlo, nhi
        weights1 = 0.0_dp
        call uniform_interpolation_weights(x, weights1(0:i))
        select case (i)
          case (1)
            call fill_interp_weights1(x, weights2)
          case (2)
            call fill_interp_weights2(x, weights2)
          case (3)
            call fill_interp_weights3(x, weights2)
          case (4)
            call fill_interp_weights4(x, weights2)
          case (5)
            call fill_interp_weights5(x, weights2)
          case (6)
            call fill_interp_weights6(x, weights2)
          case Default
            call fail("Unsupported order "//trim(to_string(i))//" in test_interpolation_coeffs")
            continue
        end select

        call check_approx_eq("fill_interp_weightsN, order = "//trim(to_string(i))//", x="//trim(to_string(x)), &
                              weights2(0:i), weights1(0:i), tol_abs = 1e-10_dp)
      end do
    end do

  ! some basic timing code
  ! NB: it only runs if timing_name=="fill_interp_weights6:x4"
  TIME_START("fill_interp_weights6:x4",10**5)
  real(dp) :: weights_sum(0:6) = 0.0_dp
  TIME_LOOP
    ! cover a range of x values to include special cases & not
    do ix = 1, size(xvals)
      call fill_interp_weights6(xvals(ix), weights2)
      weights_sum = weights_sum + weights2
    end do
  TIME_END(weights_sum)

  end subroutine test_interpolation_coeffs


  !! a number of tests of the evaluation of the PDF tables
  subroutine test_tab_eval()
    use pdf_tabulate
    use pdf_representation
    use streamlined_interface
    real(dp), parameter :: xvals(*) = [1e-5_dp, 0.12_dp, 0.9_dp]
    real(dp), parameter :: Qvals(*) = [1.0_dp, sqrt(2.0_dp), 3.0_dp, 7.0_dp, 100.0_dp, 700.0_dp]
    !real(dp), parameter :: xvals(*) = [1e-5_dp]
    !real(dp), parameter :: Qvals(*) = [3.0_dp]
    real(dp) :: x, Q, xpdf(iflv_min:tables(0)%tab_iflv_max), xpdff
    integer  :: iflv, ix, iQ

    ! first check that single-flavour evaluation gives the same answer
    ! as all-flavour evaluation
    do ix = 1, size(xvals)
      x = xvals(ix)
      do iQ = 1, size(Qvals)
        Q = Qvals(iQ)
        call EvalPdfTable_xQ(tables(0), x, Q, xpdf)
        do iflv = iflv_min, iflv_max
          xpdff = EvalPdfTable_xQf(tables(0), x, Q, iflv)
          call check_approx_eq_0d("EvalPdfTable_xQf, x="//trim(to_string(x))//&
                               ", Q="//trim(to_string(Q))//", iflv="//trim(to_string(iflv)), &
                               xpdff, xpdf(iflv), tol_abs = 1e-10_dp)
        end do
      end do
    end do 

    ! Add tests that the general table-evaluation routine and the 
    ! order-specific ones agree
    block
      real(dp) :: xpdf_ord(iflv_min:tables(0)%tab_iflv_max)
      procedure(EvalPdfTable_yQ_interface ), pointer :: EvalPdfTable_yQ_order
      procedure(EvalPdfTable_yQf_interface), pointer :: EvalPdfTable_yQf_order
      integer, parameter :: orders(*) = [2,3,4]
      integer            :: iord

      do iord = 1, size(orders)
        select case (orders(iord))
          !case (1)
          !  EvalPdfTable_yQ_order => EvalPdfTable_yQ_order1
          !  EvalPdfTable_yQf_order => EvalPdfTable_yQf_order1
          case (2)
            EvalPdfTable_yQ_order => EvalPdfTable_yQ_order22
            EvalPdfTable_yQf_order => EvalPdfTable_yQf_order22
          case (3)
            EvalPdfTable_yQ_order => EvalPdfTable_yQ_order33
            EvalPdfTable_yQf_order => EvalPdfTable_yQf_order33
          case (4)
            EvalPdfTable_yQ_order => EvalPdfTable_yQ_order44
            EvalPdfTable_yQf_order => EvalPdfTable_yQf_order44
          case Default
            call fail("Unsupported order "//trim(to_string(orders(iord)))//" in test_tab_eval")
            cycle
        end select
        call PdfTableOverrideInterpOrders(orders(iord), orders(iord))
        do ix = 1, size(xvals)
          x = xvals(ix)
          do iQ = 1, size(Qvals)
            Q = Qvals(iQ)
            call EvalPdfTable_yQ(tables(0), -log(x), Q, xpdf)
            call EvalPdfTable_yQ_order(tables(0), -log(x), Q, xpdf_ord)
            call check_approx_eq("EvalPdfTable_yQ, order="//trim(to_string(orders(iord)))//&
                                 ", x="//trim(to_string(x))//", Q="//trim(to_string(Q)), &
                                 xpdf_ord, xpdf, tol_abs = 1e-10_dp)
            do iflv = iflv_min, iflv_max
              xpdff = EvalPdfTable_yQf_order(tables(0), -log(x), Q, iflv)
              call check_approx_eq_0d("EvalPdfTable_yQf, order="//trim(to_string(orders(iord)))//&
                                   ", x="//trim(to_string(x))//", Q="//trim(to_string(Q))//", iflv="//trim(to_string(iflv)), &
                                   xpdf_ord(iflv), xpdff, tol_abs = 1e-10_dp)
            end do
          end do
        end do
      end do
    end block
  end subroutine test_tab_eval

end program hoppet_unit_tests
