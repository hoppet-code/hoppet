      subroutine eval1dhplat0(y,nw,H1,H2,H3,H4,H5,H6, 
     $                          HY1,HY2,HY3,HY4,HY5,HY6,
     $                          Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 
** evaluates 1dhpl's in the 0-range  -(r2-1) < y <= (r2-1) 
** by direct series expansion (Bernoulli-accelerated) 
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
** evaluate the irreducible 1dHPL's first 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
      call fillirr1dhplat0(y,nw,HY1,HY2,HY3,HY4,HY5,HY6,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4,H5,H6, 
     $                     HY1,HY2,HY3,HY4,HY5,HY6,
     $                     Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 
      return 
      end 
************************************************************************ 
