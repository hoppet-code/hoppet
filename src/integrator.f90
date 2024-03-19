!======================================================================
! Need a (set of?) special integrator(s) for going from 
! splitting/coefficient functions to splitting/coefficient arrays
! $Id: integrator.f90,v 1.2 2001/07/17 13:51:01 gsalam Exp $
module integrator
  use types; use consts_dp
  use warnings_and_errors
  implicit none
  private

  public :: ig_LinWeight, ig_LinWeightSing, ig_PolyWeight, ig_PolyWeight_expand

contains

  !======================================================================
  ! Function which integrates F weighted with the linear function which has
  ! values AMult and BMult at A & B respectively.
  !
  ! Result evaluated roughly with precision EPS (rel or abs 
  ! whichever is largest).
  !
  ! Perhaps in general one is better off writing a version of this routine
  ! which has two function, which get multiplied. Then this routine
  ! would just set parameters for one of those functions. This would
  ! allow the easy generalisation to the case with more complex weight
  ! functions.
  ! 
  Recursive FUNCTION ig_LinWeight(F,A,B,AMult,BMult,EPS, split) result(cgauss64)
    real(dp), intent(in) :: A,B,AMult,BMult,EPS
    real(dp), intent(in), optional :: split(:)
    real(dp), allocatable :: edges(:)
    integer               :: i, n
    real(dp) :: lmult, rmult
    REAL(dp) :: AA,BB,U,C1,C2,S8,S16,H, CGAUSS64, pmult,mmult,Const
    real(dp), parameter :: z1 = 1, hf = half*z1, cst = 5*Z1/1000
    real(dp) :: X(12), W(12)
    interface
       function f(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
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

    if (present(split)) then
      ! allocation cost is about 10ns, to be compared to integration
      ! time for a single simple function of c. 31ns. But without 
      ! splitting, integration time would be much higher (900ns versus
      ! a total of 73 ns with splitting), as well as less reliable.
      allocate(edges(size(split)+2))
      call split_limits(A,B,split,edges, n)
      cgauss64 = zero
      do i = 1, n
         lmult = ((edges(i  ) - A)/(B-A) * (BMult -AMult) + AMult)
         rmult = ((edges(i+1) - A)/(B-A) * (BMult -AMult) + AMult)
         cgauss64 = cgauss64 + ig_LinWeight(F,edges(i),edges(i+1),lmult,rmult,EPS)
      end do
      return
    end if

    H=0 
    IF(B .EQ. A) GO TO 99 
    CONST=CST/ABS(B-A) 
    BB=A 
1   AA=BB 
    BB=B 
2   C1=HF*(BB+AA) 
    C2=HF*(BB-AA) 
    S8=0 
    DO I = 1,4 
       U=C2*X(I) 
       pmult = ((c1+u) - A)/(B-A) * (BMult -AMult) + AMult
       mmult = ((c1-u) - A)/(B-A) * (BMult -AMult) + AMult
       S8=S8+W(I)*(F(C1+U)*pmult+F(C1-U)*mmult) 
    END DO
    S16=0 
    DO I = 5,12 
       U=C2*X(I) 
       pmult = ((c1+u) - A)/(B-A) * (BMult -AMult) + AMult
       mmult = ((c1-u) - A)/(B-A) * (BMult -AMult) + AMult
       S16=S16+W(I)*(F(C1+U)*pmult+F(C1-U)*mmult) 
    END DO
    S16=C2*S16 
    IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN 
       H=H+S16 
       IF(BB .NE. B) GO TO 1 
    ELSE 
       BB=C1 
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2 
       H=0 
       !CALL MTLPRT(NAME,'D113.1','TOO HIGH ACCURACY REQUIRED') 
       write(0,*) NAME,'D113.1','TOO HIGH ACCURACY REQUIRED'
       GO TO 99 
    END IF
99  cgauss64=H 
  end function ig_LinWeight


  !-------------------------------------------------------
  ! Try to improve convergence on nasty integrals
  Recursive FUNCTION ig_LinWeightSing(F,A_in,B_in,AMult,BMult,EPS) &
       &result(cgauss64)
    real(dp), intent(in) :: A_in,B_in,AMult,BMult,EPS
    REAL(dp) :: AA,BB,U,C1,C2,S8,S16,H, CGAUSS64, pmult,mmult,Const
    real(dp), parameter :: z1 = 1, hf = half*z1, cst = 5*Z1/1000
    real(dp) :: X(12), W(12), A, B, xp, xm, wp, wm
    interface
       function f(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
    integer :: i
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


    ! change variables specifically for case of singularity close to zero
    A = sqrt(A_in); B = sqrt(B_in)
    
    
    H=0 
    IF(B .EQ. A) GO TO 99 
    CONST=CST/ABS(B-A) 
    BB=A 
1   AA=BB 
    BB=B 
2   C1=HF*(BB+AA) 
    C2=HF*(BB-AA) 
    S8=0 
    DO I = 1,4 
       U=C2*X(I)
       !-- get original variables...
       xp = (c1+u)**2; xm = (c1-u)**2
       pmult = (xp - A_in)/(B_in-A_in) * (BMult -AMult) + AMult
       mmult = (xm - A)/(B-A) * (BMult -AMult) + AMult
       !-- account for change of weight due two change of var
       pmult = pmult * two*(c1+u)
       mmult = mmult * two*(c1-u)
       S8=S8+W(I)*(F(xp)*pmult+F(xm)*mmult) 
    END DO
    S16=0 
    DO I = 5,12 
       U=C2*X(I) 
       !-- get original variables...
       xp = (c1+u)**2; xm = (c1-u)**2
       pmult = (xp - A_in)/(B_in-A_in) * (BMult -AMult) + AMult
       mmult = (xm - A)/(B-A) * (BMult -AMult) + AMult
       !-- account for change of weight due two change of var
       pmult = pmult * two*(c1+u)
       mmult = mmult * two*(c1-u)
       S16=S16+W(I)*(F(xp)*pmult+F(xm)*mmult) 
    END DO
    S16=C2*S16 
    IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN 
       H=H+S16 
       IF(BB .NE. B) GO TO 1 
    ELSE 
       BB=C1 
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2 
       H=0 
       !CALL MTLPRT(NAME,'D113.1','TOO HIGH ACCURACY REQUIRED') 
       write(0,*) NAME,'D113.1','TOO HIGH ACCURACY REQUIRED'
       GO TO 99 
    END IF
99  cgauss64=H 
  end function ig_LinWeightSing
  

  !----------------------------------------------------------------------
  ! 
  ! Integrate F with a polynomial weight function which is zero at all
  ! nodes except that indicated by the index inode_one
  ! 
  ! If const is present then it is added to the weight function
  Recursive FUNCTION ig_PolyWeight(F,A,B,nodes,inode_one,EPS,wgtadd, split) result(cgauss64)
    real(dp), intent(in) :: A,B,nodes(:),EPS
    integer,  intent(in) :: inode_one
    real(dp), intent(in), optional :: wgtadd
    real(dp), intent(in), optional :: split(:)
    real(dp), allocatable :: edges(:)
    real(dp) :: zero_nodes(size(nodes)-1), norm_nodes, lcl_wgtadd
    integer  :: i, j, n
    REAL(dp) :: AA,BB,U,C1,C2,S8,S16,H, CGAUSS64, pmult,mmult,Const
    real(dp), parameter :: z1 = 1, hf = half*z1, cst = 5*Z1/1000
    real(dp) :: X(12), W(12)
    interface
       function f(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
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
    do i = 1, size(nodes)
       if (i /= inode_one) then
          j = j + 1
          zero_nodes(j) = nodes(i)
       end if
    end do
    norm_nodes = 1.0_dp / product(nodes(inode_one) - zero_nodes)
    if (present(wgtadd)) then
       lcl_wgtadd = wgtadd
    else
       lcl_wgtadd = zero
    end if

    H=0 
    IF(B .EQ. A) GO TO 99 
    CONST=CST/ABS(B-A) 
    BB=A 
1   AA=BB 
    BB=B 
2   C1=HF*(BB+AA) 
    C2=HF*(BB-AA) 
    S8=0 
    DO I = 1,4 
       U=C2*X(I) 
       pmult = product(c1+u - zero_nodes) * norm_nodes + lcl_wgtadd
       mmult = product(c1-u - zero_nodes) * norm_nodes + lcl_wgtadd
       S8=S8+W(I)*(F(C1+U)*pmult+F(C1-U)*mmult) 
    END DO
    S16=0 
    DO I = 5,12 
       U=C2*X(I) 
       pmult = product(c1+u - zero_nodes) * norm_nodes + lcl_wgtadd
       mmult = product(c1-u - zero_nodes) * norm_nodes + lcl_wgtadd
       S16=S16+W(I)*(F(C1+U)*pmult+F(C1-U)*mmult) 
    END DO
    S16=C2*S16 
    IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN 
       H=H+S16 
       IF(BB .NE. B) GO TO 1 
    ELSE 
       BB=C1 
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2 
       H=0 
       !CALL MTLPRT(NAME,'D113.1','TOO HIGH ACCURACY REQUIRED') 
       write(0,*) NAME,'D113.1','TOO HIGH ACCURACY REQUIRED'
       GO TO 99 
    END IF
99  cgauss64=H 
  end function ig_PolyWeight

  
  !----------------------------------------------------------------------
  ! 
  ! Integrate F with a polynomial weight function which is zero at all
  ! nodes except that indicated by the index inode_one
  ! 
  ! If const is present then it is added to the weight function
  !
  ! This version expands around the series in powers of
  ! v=(x-nodes(inode_one))~0 when v is small, to avoid numerical issues
  ! with convolution of real pieces
  Recursive FUNCTION ig_PolyWeight_expand(F,A,B,nodes,inode_one,EPS,wgtadd) result(cgauss64)
    use warnings_and_errors
    real(dp), intent(in) :: A,B,nodes(:),EPS
    integer,  intent(in) :: inode_one
    real(dp), intent(in), optional :: wgtadd
    real(dp) :: zero_nodes(size(nodes)-1), norm_nodes, lcl_wgtadd, expand1, expand2
    integer  :: i, j, k, l
    REAL(dp) :: AA,BB,U,C1,C2,S8,S16,H, CGAUSS64, pmult,mmult,Const
    real(dp), parameter :: z1 = 1, hf = half*z1, cst = 5*Z1/1000
    real(dp) :: X(12), W(12)
    ! the expansion threshold determines how close to inode we get
    ! before we switch to the expansion around inode. The
    ! switch location is about epsilon**(1/3)*dy, corresponding
    ! to the fact that we include terms up to (x-inode)**2
    ! (we assume nodes have relatively uniform spacing). 
    real(dp), parameter :: expansion_threshold_fraction = 3e-6_dp
    real(dp) :: expansion_threshold
    interface
       function f(x)
         use types; implicit none
         real(dp), intent(in) :: x
         real(dp)             :: f
       end function f
    end interface
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
    do i = 1, size(nodes)
       if (i /= inode_one) then
          j = j + 1
          zero_nodes(j) = nodes(i)
       end if
    end do
    norm_nodes = 1.0_dp / product(nodes(inode_one) - zero_nodes)
    
    ! use the spacing between the first two nodes to decide the actual
    ! expansion threshold
    expansion_threshold = expansion_threshold_fraction * abs(nodes(2) - nodes(1))
    expand1 = zero
    expand2 = zero
    ! GPS 2024-02-06: I think that the code is only correct when nodes(inode_one) is zero
    ! -- otherwise one should replace -1/zero_nodes(i) with 1/(nodes(inode_one) - zero_nodes(i))
    if (nodes(inode_one) /= zero) call wae_error("ig_PolyWeight_expand: nodes(inode_one) /= zero but is instead", &
                               &   dbleval = nodes(inode_one))
    do i = 1, size(zero_nodes)
       expand1 = expand1 - 1.0_dp / zero_nodes(i)
       do k=i+1, size(zero_nodes)
         expand2 = expand2 + one/(zero_nodes(i)*zero_nodes(k))
       end do
    end do

    if (present(wgtadd)) then
       lcl_wgtadd = wgtadd
    else
       lcl_wgtadd = zero
    end if

    H=0 
    IF(B .EQ. A) GO TO 99 
    CONST=CST/ABS(B-A) 
    BB=A 
1   AA=BB 
    BB=B 
2   C1=HF*(BB+AA) 
    C2=HF*(BB-AA) 
    S8=0 
    DO I = 1,4 
       U=C2*X(I)
       if (abs(c1+u) < expansion_threshold) then
          ! NB: in situations where lcl_wgtadd is exactly -one, (one +
          ! lcl_wgtadd) should be evaluated before anything else, giving
          ! exactly zero, and ensuring that we have high accuracy for
          ! the remaining terms.
          pmult = (one + lcl_wgtadd) + (c1+u)*expand1 + ((c1+u)**2)*expand2
       else
          pmult = product(c1+u - zero_nodes) * norm_nodes + lcl_wgtadd
       endif
       if (abs(c1-u) < expansion_threshold) then
          mmult = (one + lcl_wgtadd) + (c1-u)*expand1 + ((c1-u)**2)*expand2
       else
          mmult = product(c1-u - zero_nodes) * norm_nodes + lcl_wgtadd
       endif
       ! pmult = product(c1+u - zero_nodes) * norm_nodes + lcl_wgtadd
       ! mmult = product(c1-u - zero_nodes) * norm_nodes + lcl_wgtadd
       S8=S8+W(I)*(F(C1+U)*pmult+F(C1-U)*mmult) 
    END DO
    S16=0 
    DO I = 5,12 
       U=C2*X(I) 
       if (abs(c1+u) < expansion_threshold) then
          pmult = (one + lcl_wgtadd) + (c1+u)*expand1 + ((c1+u)**2)*expand2
       else
          pmult = product(c1+u - zero_nodes) * norm_nodes + lcl_wgtadd
       endif
       if (abs(c1-u) < expansion_threshold) then
          mmult = (one + lcl_wgtadd) + (c1-u)*expand1 + ((c1-u)**2)*expand2
       else
          mmult = product(c1-u - zero_nodes) * norm_nodes + lcl_wgtadd
       endif
       ! pmult = product(c1+u - zero_nodes) * norm_nodes + lcl_wgtadd
       ! mmult = product(c1-u - zero_nodes) * norm_nodes + lcl_wgtadd
       S16=S16+W(I)*(F(C1+U)*pmult+F(C1-U)*mmult)
    END DO
    S16=C2*S16 
    IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN 
       H=H+S16 
       IF(BB .NE. B) GO TO 1 
    ELSE 
       BB=C1 
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2 
       H=0 
       !CALL MTLPRT(NAME,'D113.1','TOO HIGH ACCURACY REQUIRED') 
       write(0,*) NAME,'D113.1','TOO HIGH ACCURACY REQUIRED'
       GO TO 99 
    END IF
99  cgauss64=H 
  end function ig_PolyWeight_expand
  

  !! given limits a,b and points at which to split the integration
  !! (split), set edges to the limits of the subintervals over which to integrate 
  !! and n to the number of non-zero sub-intervals. 
  !!
  !! Requirements are that a <= b and split(:) points should be ordered.
  !!
  !! The edges array should be at least of size size(split)+2
  subroutine split_limits(a, b, split, edges, n)
     real(dp), intent(in) :: a, b, split(:)
     real(dp), intent(out) :: edges(:)
     integer,  intent(out) :: n
     integer :: i
     if (a > b) call wae_error("integrator::split_limits, a > b, which is not allowed")
     if (size(edges) < size(split)+2) call wae_error("integrator::split_limits, edges array too small")
     edges(1) = a
     n = 0
     do i = 1, size(split)
        if (split(i) > a .and. split(i) < b) then
           n = n + 1
           edges(n+1) = split(i)
        end if
     end do
     n = n + 1
     edges(n+1) = b
  end subroutine split_limits 
end module integrator
