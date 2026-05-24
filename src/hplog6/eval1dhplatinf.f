      subroutine eval1dhplatinf(y,nw,H1,H2,H3,H4,H5,H6, 
     $                          HY1,HY2,HY3,HY4,HY5,HY6,
     $                          Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2)
** evaluates 1dhpl's in the inf-range  (r2+1) < abs(y) 
** evaluating first the H(..,x=1/y) by calling eval1dhplat0(x)  
** and then expressing H(..,y=1/x) in terms of H(..,x) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4,H5,H6 
      dimension H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          H5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2), 
     $          H6(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HY5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2),
     $          HY6(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi6(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      dimension HX1(n1:n2),HX2(n1:n2,n1:n2),HX3(n1:n2,n1:n2,n1:n2), 
     $          HX4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HX5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2), 
     $          HX6(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      parameter (pi   = 3.14159265358979324d0) 
** the nw = 1 case 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
** the nw > 1 case 
      x = 1.d0/y 
*      print*,' eval1dhplatinf: y = ',y,', x = ',x 
      call fillirr1dhplat0(x,nw,HX1,HX2,HX3,HX4,HX5,HX6,n1,n2) 
** fillirr1dhplatinf takes care automatically of all the immaginary 
** parts as well as of the jump across y=1 
      call fillirr1dhplatinf(x,nw,HX1,HX2,HX3,HX4,HX5,HX6, 
     $                            HY1,HY2,HY3,HY4,HY5,HY6, 
     $                            Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4,H5,H6, 
     $                     HY1,HY2,HY3,HY4,HY5,HY6,
     $                     Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 
      return 
      end 
************************************************************************ 
