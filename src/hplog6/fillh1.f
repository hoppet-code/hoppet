      subroutine fillh1(y,H1,HY1,Hi1,n1,n2) 
** fillh1 evaluates the 1dhpl's of weight 1 
      implicit double precision (a-h,o-z) 
      complex*16 H1 
      dimension H1(n1:n2) 
      dimension HY1(n1:n2) 
      dimension Hi1(n1:n2) 
      parameter (pi   = 3.14159265358979324d0) 
      if ( n1.eq.-1) then 
        if ( y.ge.-1.d0 ) then 
          HY1(-1) = log(1.d0+y) 
          Hi1(-1) = 0.d0 
        elseif ( y.lt.-1.d0 ) then 
          HY1(-1) = log(-1.d0-y) 
          Hi1(-1) = 1.d0 
        endif 
        H1(-1) = dcmplx(HY1(-1),pi*Hi1(-1)) 
      endif 
      if ( y.ge.0.d0 ) then 
        HY1(0) = log(y) 
*        Hi1(0) = 0.d0 
      elseif ( y.lt.0.d0 ) then 
        HY1(0) = log(-y) 
        Hi1(0) = 1.d0 
      endif 
      H1(0) = dcmplx(HY1(0),pi*Hi1(0)) 
      if ( n2.eq.1 ) then 
        if ( y.ge.1.d0 ) then 
          HY1(1) = - log(-1.d0+y) 
          Hi1(1) = 1.d0 
        elseif ( y.lt.1.d0 ) then 
          HY1(1) = - log(1.d0-y) 
          Hi1(1) = 0.d0 
        endif 
        H1(1) = dcmplx(HY1(1),pi*Hi1(1)) 
      endif 
      return 
      end 
************************************************************************ 
