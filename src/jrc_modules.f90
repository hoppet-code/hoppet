!------------------------------------------------------

module jrc_modules
use types ! So that we can use the predefined type double precision
use qcd
implicit none
public

! Variable to pass in convolutions
real(dp) :: xx

! Count number of calls of a given routine 
integer :: number_of_calls

contains
    
  ! Subroutine to define our pdf
  subroutine example_gluon(y,g)
     
    real(dp),intent(in) :: y
    real(dp),intent(out) :: g(:)
    real(dp) ::x

    x=exp(-y)
    
    g = 1.0_dp * x**(0.3_dp) * ( 1_dp -x )**( 3.0_dp )
    
  end subroutine example_gluon

  ! Function to define our pdf
  real(dp) function toy_gluon(y)

    real(dp), intent(in) :: y
    real(dp) :: z

    z=exp(-y)

    toy_gluon = 10.0_dp * z**(-0.8_dp) * (1_dp -z)**(5_dp)

  end function toy_gluon


! Our evolution kernel including virtual and delta corrections
! in the NS case
  real(dp) function GammaNS_func(y)

    use convolution_communicator ! To communicate real, virtual and delta terms
    implicit none
  
    integer :: i
    real(dp), intent(in) :: y
    real(dp) :: x,q2i,q2f,gammax_nsvu
    real(dp) :: zmax
    complex(dp) :: efnns(4),efnns24(2),efnsg(2,2),zn

    x=exp(-y)

    q2i = 2_dp
    q2f = 10000_dp

    GammaNS_func=zero
    
    !Compute delta function contribution
    zn = (2d0,0d0)    
    call ZFUNC(zn,x,q2i,q2f,EFNNS,EFNNS24,EFNSG)

    ! We want to evolve the u_v distribution
    !Delta function contribution
    if(cc_piece == cc_DELTA) then
       GammaNS_func=real(efnns(3))
       !print *,"x, gm = ",x,real(efnns(3))
    else 
       ! Real contribution
       if(cc_piece == cc_REAL .or. cc_piece == cc_REALVIRT) &
            & GammaNS_func = gammax_nsvu(x,q2i,q2f)/x
       ! Virtual contribution
       if(cc_piece == cc_VIRT .or. cc_piece == cc_REALVIRT) &
            & GammaNS_func = GammaNS_func - x * gammax_nsvu(x,q2i,q2f)
    endif
    
  end function GammaNS_func

  ! Initial NonSinglet pdf
  real(dp) function input_pdf_ns(y)

    real(dp), intent(in) :: y
    real(dp) :: z
    real(dp) :: pdf(13)

    z=exp(-y)

    call PDFINLH(z,pdf)

    ! pdf(3) -> uv pdf
    input_pdf_ns = pdf(3)

  end function input_pdf_ns

! gg splitting function including contribution from
! delta function and plus distributions
  real(dp) function Pgg_func(y)
    use types
    use convolution_communicator ! Provides cc_piece, cc_REAL ...
    use qcd 
    implicit none
    real(dp),intent(in) :: y
    real(dp)            :: x

    x=exp(-y); Pgg_func = zero
    
    ! Delta function contribution to the splitting function
    if(cc_piece == cc_DELTA) then
       Pgg_func = (11*CA)/6.0_dp
    else
       if(cc_piece == cc_REAL .or. cc_piece == cc_REALVIRT) &
           & Pgg_func = 2_dp*CA*( x/(one-x) + (one-x)/x + x*(one-x) )
       if(cc_piece == cc_VIRT .or. cc_piece == cc_REALVIRT) &
            & Pgg_func = Pgg_func - 2_dp*CA*( one/(one-x) )
       ! Multiply by a factor of x
       Pgg_func = Pgg_func * x
    end if
    
  end function Pgg_func

! Exact integration of Pgg_func using dgauss
  real(dp) function Pgg_func_int(y)
    
    implicit none
    real(dp),intent(in) :: y
    real(dp) :: eps,DGAUSS

    !external toy_gluon
    
    eps=0.0001_dp
    !Pgg_func_int = zero
    Pgg_func_int=DGAUSS(toy_gluon,0.1_dp,0.8_dp,eps)
    
  end function Pgg_func_int

! Exact integration of Pgg_func using dgauss
  real(dp) function unit(y)
    use types

    implicit none
    real(dp),intent(in) :: y

    unit= one

  end function unit

! The evolved u_v pdf with the exact NNPDF evolution
! kernel and dgauss integration
  real(dp) function evolved_uv(y)

    implicit none
    real(dp),intent(in) :: y
    real(dp) :: eps,DGAUSS,x
    real(dp) :: pdf0(1:13)
    real(dp) :: q2i,q2f,gamma
    real(dp) :: pdfevol
    complex(dp) :: efnns(4),efnns24(2),efnsg(2,2),zn
    ! Do not declare the functions which go as argument in
    ! dgauss unless they velong to another module!

    x=exp(-y)
    
    evolved_uv=0

    ! Set accuracy
    eps=0.0000001_dp

    ! call initial pdf
    call PDFINLH(x,pdf0)

    ! Set evolution scale
    Q2i=2_dp
    Q2f=10000_dp

    ! Compute gamma
    zn = (2d0,0d0)    
    call ZFUNC(zn,x,Q2i,Q2f,EFNNS,EFNNS24,EFNSG)
    
    !Set internal variable
    xx=x
    pdfevol=dgauss(pdfevolxint,x,1.0_dp,eps) 

    gamma=real(efnns(3))-DGAUSS(gammax_wrp,0.0_dp,x,eps)
    evolved_uv=pdfevol+gamma*pdf0(3)
    
  end function evolved_uv

  real(dp) function gammax_wrp(y)

    implicit none
    real(dp),intent(in) :: y
    real(dp) :: q2i,q2f,gammax_nsvu

    q2i=2_dp
    q2f=10000_dp

    gammax_wrp = gammax_nsvu(y,Q2i,Q2f)

  end function gammax_wrp

  real(dp) function pdfevolxint(y)

    implicit none
    real(dp),intent(in) :: y
    real(dp) :: q2i,q2f,gammax_nsvu
    real(dp) :: pdf0(1:13),pdf1(1:13)

    q2i=2_dp
    q2f=10000_dp

    CALL PDFINLH(xx,pdf0)
    CALL PDFINLH(xx/y,pdf1)

!    print *,xx

    if(xx>y) then
       print *,"Problem =",xx," > ",y
       call exit(-10)
    endif
    
    pdfevolxint = (gammax_nsvu(y,Q2i,Q2f)/y)*( pdf1(3)/y - y*pdf0(3) )
    ! write(6,'(5f15.8)')y,pdfevolxint

  end function pdfevolxint

! Subroutine with the same interface al LHAPDF
  subroutine LHAsub(x,Q,res)
    implicit none

    real(dp),intent(in) :: x,Q
    real(dp),intent(out) :: res(*)
    real(dp) :: pdfin(1:13)

    call PDFINLH(x,pdfin)
    call PDFCOMP(pdfin,res)

  end subroutine LHAsub

subroutine LHAsub2(x,res2)
    implicit none

    real(dp),intent(in) :: x
    real(dp),intent(out) :: res2
    real(dp) :: pdfin(1:13),res(-6:6)

    call PDFINLH(x,pdfin)
    call PDFCOMP(pdfin,res)

    res2=res(0)

  end subroutine LHAsub2
  
  
end module jrc_modules

!--------------------------------------------------
