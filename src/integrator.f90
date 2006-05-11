!======================================================================
! Need a (set of?) special integrator(s) for going from 
! splitting/coefficient functions to splitting/coefficient arrays
! $Id: integrator.f90,v 1.2 2001/07/17 13:51:01 gsalam Exp $
module integrator
  use types; use consts_dp
  implicit none
  private

  public :: ig_LinWeight, ig_LinWeightSing, ig_PolyWeight
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
  Recursive FUNCTION ig_LinWeight(F,A,B,AMult,BMult,EPS) result(cgauss64)
    real(dp), intent(in) :: A,B,AMult,BMult,EPS
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
                                                                        
    H=0 
    IF(B .EQ. A) GO TO 99 
    CONST=CST/ABS(B-A) 
    BB=A 
1   AA=BB 
    BB=B 
2   C1=HF*(BB+AA) 
    C2=HF*(BB-AA) 
    S8=0 
    DO 3 I = 1,4 
       U=C2*X(I) 
       pmult = ((c1+u) - A)/(B-A) * (BMult -AMult) + AMult
       mmult = ((c1-u) - A)/(B-A) * (BMult -AMult) + AMult
3      S8=S8+W(I)*(F(C1+U)*pmult+F(C1-U)*mmult) 
    S16=0 
    DO 4 I = 5,12 
       U=C2*X(I) 
       pmult = ((c1+u) - A)/(B-A) * (BMult -AMult) + AMult
       mmult = ((c1-u) - A)/(B-A) * (BMult -AMult) + AMult
4      S16=S16+W(I)*(F(C1+U)*pmult+F(C1-U)*mmult) 
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
    DO 3 I = 1,4 
       U=C2*X(I)
       !-- get original variables...
       xp = (c1+u)**2; xm = (c1-u)**2
       pmult = (xp - A_in)/(B_in-A_in) * (BMult -AMult) + AMult
       mmult = (xm - A)/(B-A) * (BMult -AMult) + AMult
       !-- account for change of weight due two change of var
       pmult = pmult * two*(c1+u)
       mmult = mmult * two*(c1-u)
3      S8=S8+W(I)*(F(xp)*pmult+F(xm)*mmult) 
    S16=0 
    DO 4 I = 5,12 
       U=C2*X(I) 
       !-- get original variables...
       xp = (c1+u)**2; xm = (c1-u)**2
       pmult = (xp - A_in)/(B_in-A_in) * (BMult -AMult) + AMult
       mmult = (xm - A)/(B-A) * (BMult -AMult) + AMult
       !-- account for change of weight due two change of var
       pmult = pmult * two*(c1+u)
       mmult = mmult * two*(c1-u)
4      S16=S16+W(I)*(F(xp)*pmult+F(xm)*mmult) 
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
  Recursive FUNCTION ig_PolyWeight(F,A,B,nodes,inode_one,EPS,wgtadd) result(cgauss64)
    real(dp), intent(in) :: A,B,nodes(:),EPS
    integer,  intent(in) :: inode_one
    real(dp), intent(in), optional :: wgtadd
    real(dp) :: zero_nodes(size(nodes)-1), norm_nodes, lcl_wgtadd
    integer  :: i, j
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
    DO 3 I = 1,4 
       U=C2*X(I) 
       pmult = product(c1+u - zero_nodes) * norm_nodes + lcl_wgtadd
       mmult = product(c1-u - zero_nodes) * norm_nodes + lcl_wgtadd
3      S8=S8+W(I)*(F(C1+U)*pmult+F(C1-U)*mmult) 
    S16=0 
    DO 4 I = 5,12 
       U=C2*X(I) 
       pmult = product(c1+u - zero_nodes) * norm_nodes + lcl_wgtadd
       mmult = product(c1-u - zero_nodes) * norm_nodes + lcl_wgtadd
4      S16=S16+W(I)*(F(C1+U)*pmult+F(C1-U)*mmult) 
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

end module integrator
