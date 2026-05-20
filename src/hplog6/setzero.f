      subroutine setzero(nw,Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2) 
** initializes with 0 the elements of the arrays 
      implicit double precision (a-h,o-z) 
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi5(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2), 
     $          Hi6(n1:n2,n1:n2,n1:n2,n1:n2,n1:n2,n1:n2) 
      do k1=n1,n2 
        Hi1(k1) = 0.d0 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            Hi2(k1,k2) = 0.d0 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                Hi3(k1,k2,k3) = 0.d0 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    Hi4(k1,k2,k3,k4) = 0.d0 
                    if ( nw.gt.4 ) then 
                      do k5=n1,n2 
                        Hi5(k1,k2,k3,k4,k5) = 0.d0 
                        if ( nw.gt.5 ) then 
                           do k6=n1,n2 
                              Hi6(k1,k2,k3,k4,k5,k6) = 0.d0 
                           enddo
                        endif
                      enddo 
                    endif 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      return 
      end 
************************************************************************ 
