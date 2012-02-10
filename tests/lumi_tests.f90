!! program to test parton luminosities
!!
!!
module lumi_elements
  use hoppet_v1
  use integrator
  use interpolation
  implicit none

  real(dp), pointer :: gq_copy(:), gq2_copy(:)
  type(grid_def) :: grid_copy
  real(dp) :: y_copy

contains

  function lumi_dgauss(grid, gq1, gq2, y) result(res)
    type(grid_def), intent(in) :: grid
    real(dp),        intent(in), target  :: gq1(0:), gq2(0:)
    real(dp), intent(in) :: y
    real(dp) :: res
    gq_copy => gq1
    gq2_copy => gq2
    grid_copy = grid
    y_copy = y

    res = ig_LinWeight(lumi_dgauss_integrand, zero, y, one, one, 1e-10_dp)
  end function lumi_dgauss
  

  function lumi_dgauss_integrand(y) result(res)
    real(dp), intent(in) :: y
    real(dp) :: res

    res = EvalGridQuant(grid_copy, gq_copy, y) &
         &       * EvalGridQuant(grid_copy, gq2_copy, y_copy-y)
  end function lumi_dgauss_integrand
  


  recursive subroutine conv_AddGridConv_alt(gc, gq)
    type(grid_conv), intent(inout) :: gc
    real(dp),        intent(in), target    :: gq(0:)

    gq_copy => gq
    grid_copy = gc%grid
    call AddWithCoeff(gc, conv_AddGridConv_integrand)
  end subroutine conv_AddGridConv_alt
  
  function conv_AddGridConv_integrand(y) result(res)
    use convolution_communicator
    real(dp), intent(in) :: y
    real(dp)             :: res
    
    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = EvalGridQuant(grid_copy, gq_copy, y)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       res = res - zero
    case(cc_DELTA)
       res = zero
    end select
  end function conv_AddGridConv_integrand
  

  subroutine polypoly_test()
    real(dp) :: nodesA(5), nodesB(5), result
    integer :: i, j
    nodesA = (/ 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp/)
    nodesB = (/ 1.0_dp, 0.0_dp,-1.0_dp, -2.0_dp,-3.0_dp/)
    result = 0.0_dp
    do i = 1, 5
      do j = 1, 5
        result = result + ig_PolyPoly(zero, one, nodesA, i, nodesB, j)
      end do
    end do
    
    write(0,*) result

  end subroutine polypoly_test
  
  
  !! function that integrates the product of two interpolation
  !! polynomials, specified by nodesA and nodesB and the "1" nodes
  !! inodeA_one and inodeB_one
  !!
  Recursive FUNCTION ig_PolyPoly(A,B,nodesA,inodeA_one,&
       &                             nodesB,inodeB_one) result(cgauss64)
    real(dp), intent(in) :: A,B,nodesA(:),nodesB(:)
    integer,  intent(in) :: inodeA_one, inodeB_one
    real(dp) :: zero_nodesA(size(nodesA)-1)
    real(dp) :: zero_nodesB(size(nodesB)-1)
    real(dp) :: norm_nodes
    integer  :: i, j
    REAL(dp) :: AA,BB,U,C1,C2,S8,S16,H, CGAUSS64, pmult,mmult,Const
    real(dp), parameter :: z1 = 1, hf = half*z1, cst = 5*Z1/1000
    real(dp) :: X(12), W(12)
    CHARACTER(len=*), parameter ::  NAME = 'cgauss64'
    
    DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/ 
    DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/ 
    DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/ 
    DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/ 
    DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/ 
    DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/ 
    DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/ 
    DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/ 
    DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/ 
    DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/ 
    DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/ 
    DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/ 

    !-- first set up the structure for the polynomial interpolation--
    j = 0
    do i = 1, size(nodesA)
       if (i /= inodeA_one) then
          j = j + 1
          zero_nodesA(j) = nodesA(i)
       end if
    end do
    j = 0
    do i = 1, size(nodesB)
       if (i /= inodeB_one) then
          j = j + 1
          zero_nodesB(j) = nodesB(i)
       end if
    end do
    norm_nodes = 1.0_dp / product(nodesA(inodeA_one) - zero_nodesA)
    norm_nodes = norm_nodes / product(nodesB(inodeB_one) - zero_nodesB)

    H=0 
    IF(B .EQ. A) GO TO 99 
    CONST=CST/ABS(B-A) 
    BB=A 
1   AA=BB 
    BB=B 
2   C1=HF*(BB+AA) 
    C2=HF*(BB-AA) 
    ! Gaussian integration with 8 points should work up to
    ! order 16, so we only need the 8-point version if the
    ! polynomial we're integrating isn't of too high order
    if (size(nodesA)+size(nodesB) < 15) then
      S8=0 
      DO I = 1,4 
        U=C2*X(I) 
        pmult = product(c1+u - zero_nodesA)*product(c1+u - zero_nodesB)
        mmult = product(c1-u - zero_nodesA)*product(c1-u - zero_nodesB)
        S8=S8+W(I)*(pmult+mmult) * norm_nodes
      end do
      H = S8*C2
    else
      S16=0 
      DO I = 5,12 
        U=C2*X(I) 
        pmult = product(c1+u - zero_nodesA)*product(c1+u - zero_nodesB)
        mmult = product(c1-u - zero_nodesA)*product(c1-u - zero_nodesB)
        S16=S16+W(I)*(pmult+mmult) * norm_nodes
      end do
      H = C2*S16 
    end if
    !write(0,*) C2*S8, C2*S16
99  cgauss64=H 
  end function ig_PolyPoly
  

  recursive subroutine conv_AddGridConv_plain(gc, gq)
    type(grid_conv), intent(inout) :: gc
    real(dp),        intent(in)    :: gq(0:)
    integer :: ny, isub
    ny = assert_eq(gc%grid%ny, ubound(gq,1),"conv_AddGridConv_gq")
    if (gc%grid%nsub == 0) then
      gc%conv(:,1) = gq * gc%grid%dy
    else
      do isub = 1, gc%grid%nsub
        call conv_AddGridConv_gq(gc%subgc(isub), &
                  &              gq(gc%grid%subiy(isub):gc%grid%subiy(isub+1)-1))
      end do
    end if
  end subroutine conv_AddGridConv_plain
  
  recursive subroutine conv_AddGridConv_gq(gc, gq)
    type(grid_conv), intent(inout) :: gc
    real(dp),        intent(in)    :: gq(0:)
    real(dp), allocatable :: matrix(:,:), nodesA(:), nodesB(:)
    integer :: ny, isub, i, j, k, order, n

    ny = assert_eq(gc%grid%ny, ubound(gq,1),"conv_AddGridConv_gq")
    if (gc%grid%nsub == 0) then
      ! now prepare an interpolating setup
      order = abs(gc%grid%order)
      allocate(matrix(order+1, order+1))
      allocate(nodesA(order+1))
      allocate(nodesB(order+1))
      nodesA = (/ (j,  j=0, order) /) * gc%grid%dy
      nodesB = (/ (1-j,j=0, order) /) * gc%grid%dy
      do i = 1, order+1
        do j = 1, order+1
          matrix(i,j) = ig_PolyPoly(zero, gc%grid%dy, nodesA, i, nodesB, j)
        end do
      end do
      gc%conv(:,1) = zero
      do i = 0, ny-1
        !     6-5-4-3-2-1-0        PDF with which we will later convolute
        !           ==========      <- interpolation region (point where add to conv)
        !           ---             <- integration region             |
        !     =========             <- interpolation region           ^
        !         0-1-2-3-4-5-6    PDF converted to splitting fn      k
        do j = 1, order+1
          k = i + j-1          ! entry in conv
          if (k > ny) exit
          n = min(order+1,i+2) ! number of gq entries to use
          gc%conv(k,1) = gc%conv(k,1) + sum(matrix(j,1:n) * gq(i+1:i-n+2:-1))
          !write(0,*) ny, i, j, k, n, gq(i), sum(matrix(:,:)), sum(matrix(j,1:n) * gq(i:i-n+1:-1))
        end do
      end do
      !gc%conv(:,1) = gq * gc%grid%dy
      !write(0,*) gc%conv(20,1), gq(20) * gc%grid%dy
      !write(0,*) 'matrix: ', matrix
    else
      do isub = 1, gc%grid%nsub
        call conv_AddGridConv_gq(gc%subgc(isub), &
                  &              gq(gc%grid%subiy(isub):gc%grid%subiy(isub+1)-1))
      end do
    end if
    
  end subroutine conv_AddGridConv_gq
  
  !======================================================================
  !> attempt at using interpolation to double range of acceptable 
  !> results at large x from the straightforward sum over bins
  function lumi_multi(grid, gq1, gq2) result(lumi)
    type(grid_def), intent(in) :: grid
    real(dp), intent(in) :: gq1(0:), gq2(0:)
    real(dp)             :: lumi(0:ubound(gq1,1))
    !-----------------------
    type(grid_def), pointer :: fine, coarse
    integer  :: isub

    ! get our first guess for the luminosity function
    lumi = lumi_plain(grid, gq1, gq2)

    ! if the grid doesn't have the right structure, we just give up
    if (grid%nsub == 0 .or. .not. grid%locked) return

    ! patch things up order by order
    do isub = 1, grid%nsub-1
      fine => grid%subgd(isub)
      coarse => grid%subgd(isub+1)
      call lumi_multi_do_sub(fine,&
           &                 gq1(grid%subiy(isub):grid%subiy(isub)+fine%ny),&
           &                 gq2(grid%subiy(isub):grid%subiy(isub)+fine%ny),&
           &                 coarse,&
           &                 gq1(grid%subiy(isub+1):grid%subiy(isub+1)+coarse%ny),&
           &                 gq2(grid%subiy(isub+1):grid%subiy(isub+1)+coarse%ny),&
           &                 lumi(grid%subiy(isub+1):grid%subiy(isub+1)+coarse%ny))
    end do
    
  end function lumi_multi
  
  !----------------------------------------------------------------------
  subroutine lumi_multi_do_sub(grid_fine, gq1_fine, gq2_fine,&
       &                       grid_coarse, gq1_coarse, gq2_coarse, lumi_coarse)
    type(grid_def), intent(in) :: grid_fine, grid_coarse
    real(dp), intent(in) :: gq1_fine(0:), gq2_fine(0:)
    real(dp), intent(in) :: gq1_coarse(0:), gq2_coarse(0:)
    real(dp), intent(inout) :: lumi_coarse(:)
    !--------------------------------------
    integer  :: ny_fine, ny_coarse, max_ny_coarse
    real(dp) :: gq1_tmp(0:2*ubound(gq1_fine,1))
    real(dp) :: gq2_tmp(0:2*ubound(gq2_fine,1))
    integer  :: order, range_up, range_down, range_down_last
    integer, parameter :: nmax = 9 ! max order, as uniform_interpolation_weights
    !real(dp) :: weights(0:nmax,nint(grid_coarse%dy / grid_fine%dy)-1)
    real(dp) :: weights(0:nmax)!,nint(grid_coarse%dy / grid_fine%dy)-1)
    integer  :: i, step, iy_coarse

    step = nint(grid_coarse%dy / grid_fine%dy)

    ! make sure the ny we use is a multiple of step
    ny_coarse = (grid_fine%ny/step)
    ny_fine   = ny_coarse * step

    ! first copy the high-resolution piece across to a temporary array
    gq1_tmp(0:ny_fine) = gq1_fine(0:ny_fine)
    gq2_tmp(0:ny_fine) = gq2_fine(0:ny_fine)
    
    ! next get take the low resolution part and interpolate it
    ! first establish the appropriate range
    max_ny_coarse = min(2*ny_coarse, grid_coarse%ny)
    
    ! decide where we will interpolate
    order = grid_coarse%order
    if (grid_coarse%ny < order) call wae_error("lumi_multi_do_sub",&
         &"grid_coarse%ny < order, grid_fine%dy = ", dbleval=grid_fine%dy)

    ! first do the trivial points
    gq1_tmp(ny_fine+step:max_ny_coarse*step:step) = &
         &                      gq1_coarse(ny_coarse+1:max_ny_coarse)
    gq2_tmp(ny_fine+step:max_ny_coarse*step:step) = &
         &                      gq2_coarse(ny_coarse+1:max_ny_coarse)


    do iy_coarse = ny_coarse, max_ny_coarse-1
      range_down = max(0, iy_coarse-(order)/2)
      range_up   = min(grid_coarse%ny, range_down + order)
      range_down = range_up - order  ! our order check earlier on ensures saftey here
      
      ! now get the weights
      do i = 1, step-1
        ! NB: with a tiny bit more work, could avoid repetitive determination
        !     of interpolation weights; but for now ignore this 
        call uniform_interpolation_weights(i*(one/step)+iy_coarse-range_down,&
             &                             weights(0:order))
        gq1_tmp(iy_coarse*step+i) = sum(weights(0:order)*gq1_coarse(range_down:range_up))
        gq2_tmp(iy_coarse*step+i) = sum(weights(0:order)*gq2_coarse(range_down:range_up))
      end do
      
    end do

  end subroutine lumi_multi_do_sub
  
  

  !----------------------------------------------------------------------
  function lumi_plain(grid, gq1, gq2) result(lumi)
    type(grid_def), intent(in) :: grid
    real(dp), intent(in) :: gq1(0:), gq2(0:)
    real(dp)             :: lumi(0:ubound(gq1,1))
    type(grid_conv) :: gc
    call InitGridConv(grid,gc) ! sets it to zero
    call conv_AddGridConv_plain(gc, gq1)
    lumi = gc .conv. gq2
  end function lumi_plain
  

  function lumi_simple(grid, gq1, gq2) result(lumi)
    type(grid_def), intent(in) :: grid
    real(dp), intent(in) :: gq1(0:), gq2(0:)
    real(dp)             :: lumi(0:ubound(gq1,1))
    integer :: i, ny

    type(grid_conv) :: gc
    call InitGridConv(grid,gc) ! sets it to zero
    call conv_AddGridConv_gq(gc, gq1)
    !call conv_AddGridConv_alt(gc, gq1)
    !call conv_AddGridConv_plain(gc, gq1)
    !gc%conv(:,1) = gq1*grid%dy
    lumi = gc .conv. gq2

    ! write(0,*) grid%ny, ubound(gq1,1),ubound(gq2,1)
    ! ny = assert_eq(grid%ny, ubound(gq1,1),ubound(gq2,1),"lumi_simple")
    ! if (grid%nsub == 0) then
    !   do i = 0, ny
    !     lumi(i) = sum(gq1(0:i)*gq2(i:0:-1))*grid%dy
    !   end do
    ! end if
  end function lumi_simple
    
  
end module lumi_elements


program lumi_tests
  use hoppet_v1; use sub_defs_io
  use lumi_elements
  implicit none
  !-- for PDFs
  type(grid_def) :: grid, gdarray(4)
  real(dp), pointer  :: pdf(:,:), lumi(:), lumi_p(:)
  real(dp) :: dy, ymax, sum, element, lumi1, lumi2, lumi_p_aty
  integer  :: i, n, order, j

  call polypoly_test()

  dy = dble_val_opt('-dy',0.1_dp)
  ymax = 10.0_dp
  order = int_val_opt("-order",-5)
  if (log_val_opt("-single-grid")) then
    write(0,*) "grid with single spacing"
    call InitGridDef(grid,ymax=10.0_dp, dy=dy, order=order)
  else
    write(0,*) "grid with multiple spacings"
    call InitGridDef(gdarray(4),dy/27.0_dp, 0.2_dp, order=order)
    call InitGridDef(gdarray(3),dy/9.0_dp,  0.5_dp, order=order)
    call InitGridDef(gdarray(2),dy/3.0_dp,  2.0_dp, order=order)
    call InitGridDef(gdarray(1),dy,         ymax  ,order=order)
    call InitGridDef(grid,gdarray(1:4),locked=.true.)
  end if

  call AllocPDF(grid, pdf)
  call AllocGridQuant(grid, lumi)
  call AllocGridQuant(grid, lumi_p)

  pdf = unpolarized_dummy_pdf(xValues(grid))

  n = ubound(pdf, 1)
  sum = 0.0_dp
  do i = 0, n
    element = pdf(i,iflv_g)*pdf(n-i,iflv_u)
    !write(6,*) i*dy, element
    sum = sum + element
  end do
  sum = sum * dy

  ! do i = 0, n, nint(0.2_dp/grid%dy)
  !   write(6,*) i*grid%dy, lumi(i)
  ! end do
  do j = 1, 1
    lumi = lumi_simple(grid, pdf(:,iflv_g), pdf(:,iflv_g))
  end do
  lumi_p = lumi_plain(grid, pdf(:,iflv_g), pdf(:,iflv_g))

  do i = 0, 50
    do j = 1, 1
      lumi2 = lumi_dgauss(grid, pdf(:,iflv_g), pdf(:,iflv_g), i*0.2_dp) 
    end do
    lumi1 = EvalGridQuant(grid, lumi, i*0.2_dp)
    lumi_p_aty = EvalGridQuant(grid, lumi_p, i*0.2_dp)
    write(6,*) i*0.2_dp, lumi_p_aty, lumi1, lumi2
  end do

contains
  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolazrized evolution
  function unpolarized_dummy_pdf(xvals) result(pdf)
    real(dp), intent(in) :: xvals(:)
    real(dp)             :: pdf(size(xvals),ncompmin:ncompmax)
    real(dp) :: uv(size(xvals)), dv(size(xvals))
    real(dp) :: ubar(size(xvals)), dbar(size(xvals))
    !---------------------
    real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
    real(dp), parameter :: N_uv=5.107200_dp, N_dv = 3.064320_dp
    real(dp), parameter :: N_db = half*N_ls
  
    pdf = zero
    ! clean method for labelling as PDF as being in the human representation
    call LabelPdfAsHuman(pdf)

    !-- remember that these are all xvals*q(xvals)
    uv = N_uv * xvals**0.8_dp * (1-xvals)**3
    dv = N_dv * xvals**0.8_dp * (1-xvals)**4
    dbar = N_db * xvals**(-0.1_dp) * (1-xvals)**6
    ubar = dbar * (1-xvals)
    pdf(:,iflv_g) = N_g * xvals**(-0.1_dp) * (1-xvals)**5
        
    pdf(:,-iflv_s) = 0.2_dp*(dbar + ubar)
    pdf(:, iflv_s) = pdf(:,-iflv_s)
    pdf(:, iflv_u) = uv + ubar
    pdf(:,-iflv_u) = ubar
    pdf(:, iflv_d) = dv + dbar
    pdf(:,-iflv_d) = dbar
  end function unpolarized_dummy_pdf

end program lumi_tests


