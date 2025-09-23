module split_helper
  use types;
  use consts_dp;
  use convolution
  use convolution_communicator
  use qcd
  implicit none
contains

  real(dp) function piecewise_fn(x) result(res)
    real(dp), intent(in) :: x
    if (x < 0) then
      res = -x
    else
      res = x**2
    end if
  end function piecewise_fn

  !----------------------------------------------------------------------
  function broken_Pqq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = cf*(one+x**2)/(one-x)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - cf*two/(one-x)
    case(cc_DELTA)
       res = cf*three*half
    end select

    if (x < half) res = two * res 

    if (cc_piece /= cc_DELTA) res = res * x
  end function broken_Pqq

  !----------------------------------------------------------------------
  function zerohalf_Pqq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: plus
    x = exp(-y)
    res = zero

    if (x > half .or. cc_piece == cc_DELTA) return

    plus = cf*(one+x**2)/(one-x)
    !plus = -two * log(x/(one-x)) * plus

    !! REALVIRT will be zero, because it adds REAL and VIRT
    !! which are equal and opposite
    select case(cc_piece)
    case(cc_REAL)
      res = plus
    case(cc_VIRT)
      res = -plus
    end select

    res = res * x
  end function zerohalf_Pqq

  !----------------------------------------------------------------------
  function halfone_Pqq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    real(dp) :: plus
    x = exp(-y)
    res = zero

    if (x < half .or. cc_piece == cc_DELTA) return

    plus = cf*(one+x**2)/(one-x)
    !plus = -two * log(x/(one-x)) * plus

    !! REALVIRT will be zero, because it adds REAL and VIRT
    !! which are equal and opposite
    select case(cc_piece)
    case(cc_REAL)
      res = plus
    case(cc_VIRT)
      res = -plus
    end select

    res = res * x
  end function halfone_Pqq


end module split_helper

program main
  use types
  use consts_dp
  use split_helper
  use integrator
  use io_utils
  use convolution
  use splitting_functions

  implicit none
  real(dp) :: eps = 1e-6_dp
  real(dp) :: yn, ys, yn1, yn2, ys1, ys2
  integer :: i, n
  logical :: all_good = .true.
  !print *, y

  yn = ig_LinWeight(piecewise_fn, -1.0_dp, 2.0_dp, one, one, eps)
  ys = ig_LinWeight(piecewise_fn, -1.0_dp, 2.0_dp, one, one, eps, (/0.0_dp/))
  call check("LinWeight(1.0,1.0)", yn, ys)

  yn = ig_LinWeight(piecewise_fn, -1.0_dp, 2.0_dp, half, one, eps)
  ys = ig_LinWeight(piecewise_fn, -1.0_dp, 2.0_dp, half, one, eps, (/0.0_dp/))
  call check("LinWeight(0.5,2.0)", yn, ys)

  yn1 = ig_PolyWeight(piecewise_fn, -one, two, (/-one,two/), 1, eps)
  yn2 = ig_PolyWeight(piecewise_fn, -one, two, (/-one,two/), 2, eps)
  ys1 = ig_PolyWeight(piecewise_fn, -one, two, (/-one,two/), 1, eps, split=(/0.0_dp/))
  ys2 = ig_PolyWeight(piecewise_fn, -one, two, (/-one,two/), 2, eps, split=(/0.0_dp/))
  call check("PolyWeight node 1", yn1, ys1)
  call check("PolyWeight node 2", yn2, ys2)
  call check("PolyWeight sum   ", half*yn1+yn2, half*ys1+ys2)

  ! yn1 = ig_PolyWeight_expand(piecewise_fn, -one, two, (/-one,two/), 1, eps)
  ! yn2 = ig_PolyWeight_expand(piecewise_fn, -one, two, (/-one,two/), 2, eps)
  ! ys1 = ig_PolyWeight_expand(piecewise_fn, -one, two, (/-one,two/), 1, eps, split=(/0.0_dp/))
  ! ys2 = ig_PolyWeight_expand(piecewise_fn, -one, two, (/-one,two/), 2, eps, split=(/0.0_dp/))
  ! call check("PolyWeight_expand node 1", yn1, ys1)
  ! call check("PolyWeight_expand node 2", yn2, ys2)
  ! call check("PolyWeight_expand sum   ", half*yn1+yn2, half*ys1+ys2)

  call convolution_checks()

  n = int(dble_val_opt("-n", 0.0_dp))
  print *, "Now trying to integrate ", n, " times."
  do i = 1, n
    ys = ys + ig_LinWeight(piecewise_fn, -1.0_dp, 2.0_dp, half, one, eps, (/0.0_dp/))
    !y = y + ig_LinWeight(piecewise_fn, -1.0_dp, 0.0_dp, half, one, eps)
  end do


  if (.not. all_good) stop 1
contains

  subroutine check(label, a,b, limit)
    use assertions
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: a, b
    real(dp), intent(in), optional :: limit
    real(dp) :: limit_lcl
    print *, label, " NO   split:", a
    print *, label, " WITH split:", b

    limit_lcl = default_or_opt(eps * (one + abs(a) + abs(b)), limit)
    if (abs(a-b) > limit_lcl) then
      all_good = .false.
      print *, "Error: |a-b| =  |", a - b , "| > ", limit_lcl
    end if
  end subroutine check


  subroutine convolution_checks()
    real(dp) :: dy, ymax = 10.0_dp
    type(grid_def) :: grid
    type(grid_conv) :: gc_broken_pqq, gc_broken_pqq_split, gc_pqq
    type(grid_conv) :: gc_zerohalf_pqq, gc_halfone_pqq
    type(grid_conv) :: gc_zerohalf_pqq_split, gc_halfone_pqq_split
    real(dp), pointer :: pdf(:), conv_res1(:), conv_res2(:), xvals
    real(dp) :: dx=0.001_dp, x
    integer  :: ix
    real(dp) :: split(1) = (/log(two)/)


    dy = dble_val_opt("-dy", 0.1_dp)
    dx = dble_val_opt("-dx", 0.1_dp)
    call InitGridDefDefault(grid, dy, ymax)

    call InitGridConv(grid, gc_broken_pqq, broken_Pqq)
    call InitGridConv(grid, gc_broken_pqq_split, broken_Pqq, split=split)

    call InitGridConv(grid, gc_zerohalf_pqq, zerohalf_Pqq)
    call InitGridConv(grid, gc_halfone_pqq , halfone_Pqq )
    call InitGridConv(grid, gc_zerohalf_pqq_split, zerohalf_Pqq, split=split)
    call InitGridConv(grid, gc_halfone_pqq_split , halfone_Pqq , split=split)

    !call AddWithCoeff(gc_broken_pqq_split, broken_Pqq, split=split)
    call InitGridConv(grid, gc_pqq, sf_Pqq)
    !call InitGridConv(grid, gc_broken_pqq_split)
    !call AddWithCoeff(gc_broken_pqq_split, sf_Pqq, split=(/log(two)/))

    call AllocGridQuant(grid, pdf)
    call AllocGridQuant(grid, conv_res1)
    call AllocGridQuant(grid, conv_res2)

    pdf = (one - xValues(grid))**3 * xValues(grid)**(-0.2_dp)
    !conv_res1 = gc_halfone_pqq       * pdf
    !conv_res2 = gc_halfone_pqq_split * pdf
    !conv_res1 = gc_zerohalf_pqq       * pdf
    !conv_res2 = gc_zerohalf_pqq_split * pdf
    !conv_res1 = gc_zerohalf_pqq       * pdf + gc_halfone_pqq       * pdf - gc_pqq * pdf
    !conv_res2 = gc_zerohalf_pqq_split * pdf + gc_halfone_pqq_split * pdf - gc_pqq * pdf

    ! for the unit test: compare the evaluation with complete splitting function
    ! to the evaluation with the two halves of the splitting function
    conv_res1 = gc_pqq * pdf
    conv_res2 = gc_zerohalf_pqq_split * pdf + gc_halfone_pqq_split * pdf 
    !conv_res2 = gc_zerohalf_pqq * pdf + gc_halfone_pqq * pdf 

    if  (maxval(abs(conv_res1-conv_res2)) > eps) then
      print *, "Error: maxval(abs(conv_res1-conv_res2)) = ", maxval(abs(conv_res1-conv_res2))
      all_good = .false.
    end if
    

    do ix = max(1,int(0.02/dx)), int(0.9/dx)
      x = ix * dx
      ! XX makes it easy to grep out the result
      print *, "XX", x, conv_res1.atx.(x.with.grid), conv_res2.atx.(x.with.grid)
    end do

  end subroutine convolution_checks
end program main
