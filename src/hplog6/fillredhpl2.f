      subroutine FILLREDHPL2(iflag,H1,H2,i1,i2,na,nb) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2 
      dimension H1(i1:i2),H2(i1:i2,i1:i2) 
*23456789012345678901234567890123456789012345678901234567890123456789012 
* must be called with ordered indices na <= nb 
*      print*,' FILLREDHPL2, iflag =',iflag 
      if ( na.eq.nb ) then 
        H2(na,na) = 1.d0/2*( H1(na) )**2 
      else 
        H2(nb,na) = + H1(na)*H1(nb) - H2(na,nb) 
        if ( iflag.eq.1 ) then 
          call printer2(na,nb) 
        endif 
      endif 
      return 
      end 
************************************************************************ 
