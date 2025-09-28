#include "timing.inc"

!!!!!!!!!!!!!!!!!!!!!!!
program hoppet_unit_tests
  use types
  use unit_tests
  use test_interpolation
  use test_alphas
  use io_utils
  implicit none

  if (log_val_opt("-h", .false.) .or. log_val_opt("--help", .false.)) then
    print *, "Usage: hoppet-unit-tests [-h|--help] [-time TIMING_NAME] [-list-timing]"
    print *, "  -h, --help          Print this help message"
    print *, "  -time TIMING_NAME   Run timing for the specified test (default: none)"
    print *, "                      (anything whose initial characters match TIMING_NAME will be timed)"
    print *, "  -list-timing        List all available timing names"
    stop 0
  end if
  print '(a)', "Running HOPPET unit tests"
  call hoppet_setup()  
  timing_name = string_val_opt("-time","none") ! get timing name from command line
  list_timing = log_val_opt("-list-timing",.false.)
  if (.not. CheckAllArgsUsed(0)) stop 1

  call test_interpolation_coeffs()
  call test_tab_eval()
  call unit_tests_alphas()

  if (unit_test_failures > 0) then
    ! print a message in red
    print '(a)', red//trim(to_string(unit_test_failures))//' unit tests failed.'//reset
    stop 1
  else
    ! print a message in green
    print '(a)', green//'All '//trim(to_string(unit_test_successes))//' unit tests passed.'//reset
  end if

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


  !! a number of tests of the evaluation of the PDF tables
  !! NB: when trying to include this in test_interpolation.F90, it gave 
  !! an internal compiler error with gfortran 15.1.0
  subroutine test_tab_eval()
    use interpolation
    use interpolation_coeffs
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
      real(dp) :: xpdf_ord(iflv_min:tables(0)%tab_iflv_max), xpdf1D_ord(iflv_min:tables(0)%tab_iflv_max,1)
      procedure(EvalPdfTable_yQ_interface ), pointer :: EvalPdfTable_yQ_order
      procedure(EvalPdfTable_yQf_interface), pointer :: EvalPdfTable_yQf_order
      procedure(EvalPdfTable1D_yQ_interface), pointer :: EvalPdfTable1D_yQ_order
      integer, parameter :: orders(*) = [22,33,44,54,64]
      integer            :: iord, order_y, order_Q

      do iord = 1, size(orders)
        order_y = orders(iord)/10
        order_Q = mod(orders(iord),10)
        call PdfTableSetInterpPointers(order_y, order_Q, &
                  EvalPdfTable_yQ_order, EvalPdfTable_yQf_order, EvalPdfTable1D_yQ_order)
        if (     .not. associated(EvalPdfTable_yQ_order) &
            .or. .not.associated(EvalPdfTable_yQf_order)&
            .or. .not.associated(EvalPdfTable1D_yQ_order)) then
          call fail("Unsupported order "//trim(to_string(orders(iord)))//" in test_tab_eval")
          cycle
        end if

        ! this makes sure that the call below to the any_order routine takes into account
        ! the override
        call PdfTableOverrideInterpOrders(order_y, order_Q)

        do ix = 1, size(xvals)
          x = xvals(ix)
          do iQ = 1, size(Qvals)
            Q = Qvals(iQ)
            ! call the "any_order" routine to make sure we don't have the automatic
            ! bypass to one of the hard-coded routines
            call EvalPdfTable_yQ_any_order(tables(0), -log(x), Q, xpdf)
            ! then call our own pointer to a hard-coded routine
            call EvalPdfTable_yQ_order(tables(0), -log(x), Q, xpdf_ord)
            call check_approx_eq("EvalPdfTable_yQ, order="//trim(to_string(orders(iord)))//&
                                 ", x="//trim(to_string(x))//", Q="//trim(to_string(Q)), &
                                 xpdf_ord, xpdf, tol_abs = 1e-10_dp)

            call EvalPdfTable1D_yQ_order(tables(0:0), -log(x), Q, xpdf1D_ord)
            call check_approx_eq("EvalPdfTable1D_yQ, order="//trim(to_string(orders(iord)))//&
                                 ", x="//trim(to_string(x))//", Q="//trim(to_string(Q)), &
                                 xpdf1D_ord(:,1), xpdf, tol_abs = 1e-10_dp)

            ! then flavour-by-flavour
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
