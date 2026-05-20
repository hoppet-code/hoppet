      subroutine FILLREDHPL3(iflag,H1,H2,H3,i1,i2,ia,ib,ic) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3 
      dimension H1(i1:i2),H2(i1:i2,i1:i2),H3(i1:i2,i1:i2,i1:i2) 
* must be called with "properly ordered" indices 
* note in particular the remapping of, say, (ia,ib,ic) into 
* (na,na,nb) of ReducerTest.out 
      na = ia 
      if ( (ia.eq.ib).and.(ib.eq.ic) ) then 
* case (na,na,na) 
        H3(na,na,na) = 1.d0/6*( H1(na) )**3 
* ic cannot be anymore equal to ia 
      else if ( ic.eq.ia ) then 
        print*,' FILLREDHPL3, error 1, called with arguments ' 
        print*,'               ',ia,ib,ic 
        stop 
      else if ( ia.eq.ib ) then 
* case (na,na,nb) 
        nb = ic 
        if ( iflag.eq.1 ) then 
          call printer3(na,na,nb) 
        endif 
        H3(na,nb,na) = + H1(na)*H2(na,nb) - 2*H3(na,na,nb) 
        H3(nb,na,na) = + 1.d0/2*H1(na)*H1(na)*H1(nb) 
     $                 - H1(na)*H2(na,nb) + H3(na,na,nb) 
* ib cannot be anymore equal to ia 
      else if ( ib.eq.ia ) then 
        print*,' FILLREDHPL3, error 2, called with arguments ' 
        print*,'               ',ia,ib,ic 
        stop 
      else if ( ib.eq.ic ) then 
* case (na,nb,nb) 
        nb = ib 
        if ( iflag.eq.1 ) then 
          call printer3(na,nb,nb) 
        endif 
        H3(nb,na,nb) = + H1(nb)*H2(na,nb) - 2*H3(na,nb,nb) 
        H3(nb,nb,na) = + 1.d0/2*H1(na)*H1(nb)*H1(nb) 
     $                 - H1(nb)*H2(na,nb) + H3(na,nb,nb) 
* no need to protect against ic.eq.ib 
* when arriving here all indices are different 
      else 
* case (na,nb,nc)    all indices are different 
        nb = ib 
        nc = ic 
        if ( iflag.eq.1 ) then 
          call printer3(na,nb,nc) 
          call printer3(na,nc,nb) 
        endif 
        H3(nb,na,nc) = + H1(nb)*H2(na,nc) 
     $                 - H3(na,nb,nc) - H3(na,nc,nb) 
        H3(nb,nc,na) = + H1(na)*H2(nb,nc) 
     $                 - H1(nb)*H2(na,nc) + H3(na,nc,nb) 
        H3(nc,na,nb) = + H1(nc)*H2(na,nb) 
     $                 - H3(na,nb,nc) - H3(na,nc,nb) 
        H3(nc,nb,na) = + H1(na)*H1(nb)*H1(nc) - H1(na)*H2(nb,nc) 
     $                 - H1(nc)*H2(na,nb) + H3(na,nb,nc) 
      endif 
*23456789012345678901234567890123456789012345678901234567890123456789012 
      return 
      end 
************************************************************************ 
