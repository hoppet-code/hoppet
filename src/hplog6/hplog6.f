      subroutine hplog6(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5,Hc6, 
     $                       Hr1,Hr2,Hr3,Hr4,Hr5,Hr6,
     $                       Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 
****** 
** x is the argument of the 1dHPL's (1 dimensional Harmonic PolyLogarithms) 
**   to be evaluated; 
** nw is the maximum weight of the required 1dHPL's; 
**    the maximum allowed value of nw of this implementation is 6; 
** Hc1,Hc2,Hc3,Hc4,Hc5,Hc6 are the complex*16 values of the 1dHPL; 
**    they must all be supplied in the arguments even if some of them 
**    are not to be evaluated; 
** Hr1,Hr2,Hr3,Hr4,Hr5,Hr6 are the double precision real parts of 
**    Hc1,Hc2,Hc3,Hc4,Hc5,Hc6; 
** Hi1,Hi2,Hi3,Hi4,Hi5,Hi6 are the double precision immaginary parts of 
**    Hc1,Hc2,Hc3,Hc4,Hc5,Hc6 divided by pi=3.114159.... 
** n1,n2 is the required range of indices, the allowed ranges are 
**    (0,1), (-1,0), (-1,1) ; 
** current version: must be (-1,1)
****** 
      implicit double precision (a-h,o-z) 
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5,Hc6 
      dimension Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), 
     $          Hc4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hc5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hc6(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), 
     $          Hr4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hr5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hr6(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2),
     $          Hi6(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      common /fillred/infilldim,infill(3) 
      parameter (r2   = 1.4142135623730950488d0) 
** check on the weight nw 
      if ( (nw.lt.1).or.(nw.gt.6) ) then 
        print*, ' illegal call of eval1dhpl with second argument', 
     $          ' (the weight) = ',nw 
        print*, ' the allowed values of the weight are 1,..,6 ' 
        stop
      endif
** check on the range n1:n2 
      if ( (n1.eq.-1).and.(n2.eq.0) ) then 
        infilldim =  2 
        infill(1) =  0 
        infill(2) = -1  
        print*, ' illegal call of hplog6 with index range'
        print*, ' only (-1,1) is supported ' 
        stop
      elseif ( (n1.eq.0).and.(n2.eq.1) ) then 
        infilldim =  2 
        infill(1) =  0 
        infill(2) =  1  
        print*, ' illegal call of hplog6 with index range'
        print*, ' only (-1,1) is supported ' 
        stop
      elseif ( (n1.eq.-1).and.(n2.eq.1) ) then 
        infilldim =  3 
        infill(1) =  0 
        infill(2) = -1  
        infill(3) =  1  
      else 
        print*, ' illegal call of eval1dhpl with the two last ', 
     $          'arguments = (',n1,',',n2,')' 
        print*, ' the allowed values are (-1,0), (0,1), (-1,1) ' 
        stop 
      endif 
** setting the immaginary parts equal to zero 
      call setzero(nw,Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 
** looking at the range of the argument 
*      r2 = sqrt(2.d0) 
      r2m1 = r2 - 1 
      r2p1 = r2 + 1 
      if ( ( x.gt.-r2m1 ).and.( x.le.r2m1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplat0 ' 
        call eval1dhplat0(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5,Hc6, 
     $                         Hr1,Hr2,Hr3,Hr4,Hr5,Hr6,
     $                         Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 
        return 
      elseif ( x.eq.1d0 ) then
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplin1 ' 
        call eval1dhplin1(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5,Hc6, 
     $                         Hr1,Hr2,Hr3,Hr4,Hr5,Hr6,
     $                         Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2)
        return 
      elseif ( ( x.gt.r2m1 ).and.( x.le.r2p1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplat1 ' 
        call eval1dhplat1(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5,Hc6, 
     $                         Hr1,Hr2,Hr3,Hr4,Hr5,Hr6,
     $                         Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 
        return
      elseif ( ( x.gt.r2p1 ) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatinf ' 
        call eval1dhplatinf(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5,Hc6, 
     $                           Hr1,Hr2,Hr3,Hr4,Hr5,Hr6,
     $                           Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 
        return   
      elseif ( ( x.le.-r2p1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatminf ' 
        call eval1dhplatminf(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5,Hc6, 
     $                            Hr1,Hr2,Hr3,Hr4,Hr5,Hr6,
     $                            Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 
        return 
      elseif ( x.eq.-1d0 ) then
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplinm1 ' 
        call eval1dhplinm1(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5,Hc6, 
     $                          Hr1,Hr2,Hr3,Hr4,Hr5,Hr6,
     $                          Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 
        return 
      elseif ( ( x.gt.-r2p1 ).and.( x.le.-r2m1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatm1 ' 
        call eval1dhplatm1(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5,Hc6, 
     $                          Hr1,Hr2,Hr3,Hr4,Hr5,Hr6,
     $                          Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 
        return 
      endif 
** 
      end 
************************************************************************ 
! AK: This is in fact just a little wrapper since hplog and hplog6
! clash. It just calls hplog6 after checking that nw is less-than-equal
! 4.
      subroutine hplog(x,nw,Hc1,Hc2,Hc3,Hc4, Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3
     $     ,Hi4,n1,n2) 
****** 
** x is the argument of the 1dHPL's (1 dimensional Harmonic PolyLogarithms) 
**   to be evaluated; 
** nw is the maximum weight of the required 1dHPL's; 
**    the maximum allowed value of nw of this implementation is 4; 
** Hc1,Hc2,Hc3,Hc4 are the complex*16 values of the 1dHPL; 
**    they must all be supplied in the arguments even if some of them 
**    are not to be evaluated; 
** Hr1,Hr2,Hr3,Hr4 are the double precision real parts of 
**    Hc1,Hc2,Hc3,Hc4; 
** Hi1,Hi2,Hi3,Hi4 are the double precision immaginary parts of 
**    Hc1,Hc2,Hc3,Hc4 divided by pi=3.114159.... 
** n1,n2 is the required range of indices, the allowed ranges are 
**    (0,1), (-1,0), (-1,1) ; 
****** 
      implicit double precision (a-h,o-z) 
      complex*16 Hc1,Hc2,Hc3,Hc4,Hc5,Hc6 
      dimension Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), 
     $          Hc4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hc5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hc6(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), 
     $          Hr4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hr5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hr6(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2),
     $          Hi6(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
** check on the weight nw 
      if ( (nw.lt.1).or.(nw.gt.4) ) then 
        print*, ' illegal call of hplog with second argument', 
     $          ' (the weight) = ',nw 
        print*, ' the allowed values of the weight are 1,..,4 ' 
        stop
      endif
      
      call hplog6(x,nw,Hc1,Hc2,Hc3,Hc4,Hc5,Hc6, Hr1,Hr2,Hr3,Hr4,Hr5,Hr6,
     $     Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 

      
      end 
