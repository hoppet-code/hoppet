      subroutine eval1dhplin1(y,nw,H1,H2,H3,H4,H5,H6, 
     $                          HY1,HY2,HY3,HY4,HY5,HY6,
     $                          Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2)
** evaluates 1dhpl's for y=1 (explicit values are tabulated)
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
      parameter (pi   = 3.14159265358979324d0) 
** evaluate the irreducible 1dHPL's first 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
      call fillirr1dhplin1(y,nw,HY1,HY2,HY3,HY4,HY5,HY6,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4,H5,H6, 
     $                     HY1,HY2,HY3,HY4,HY5,HY6,
     $                     Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 
      if (n2.eq.0) return
** correct the ill-defined entries
      HY2(1,0) = - HY2(0,1) 
      Hi2(1,0) = 0d0 
      H2(1,0) = dcmplx(HY2(1,0),Hi2(1,0)*pi)
      if ( nw.eq.2 ) return
      HY3(1,0,0) = HY3(0,0,1) 
      Hi3(1,0,0) = 0d0 
      H3(1,0,0) = dcmplx(HY3(1,0,0),Hi3(1,0,0)*pi)
      if ( nw.eq.3 ) return
      HY4(1,0,0,0) = -HY4(0,0,0,1) 
      Hi4(1,0,0,0) = 0d0 
      H4(1,0,0,0) = dcmplx(HY4(1,0,0,0),Hi4(1,0,0,0)*pi)
      if ( nw.eq.4 ) return
      HY5(1,0,0,0,0) = HY5(0,0,0,0,1) 
      Hi5(1,0,0,0,0) = 0d0 
      H5(1,0,0,0,0) = dcmplx(HY5(1,0,0,0,0),Hi5(1,0,0,0,0)*pi)
      if ( nw.eq.5 ) return
      HY6(1,0,0,0,0,0) = -HY6(0,0,0,0,0,1) 
      Hi6(1,0,0,0,0,0) = 0d0 
      H6(1,0,0,0,0,0) = dcmplx(HY6(1,0,0,0,0,0),Hi6(1,0,0,0,0,0)*pi)
      return 
      end 
************************************************************************ 
