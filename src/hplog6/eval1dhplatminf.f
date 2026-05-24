      subroutine eval1dhplatminf(y,nw,H1,H2,H3,H4,H5,H6, 
     $                          HY1,HY2,HY3,HY4,HY5,HY6,
     $                          Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2)
** evaluates 1dhpl's in the (-1)-range y  <= -(r2+1) 
** evaluating first the H(..,-y) by calling eval1dhplatinf(-y), 
** and then expressing H(..,y) in terms of H(..,-y) 
      implicit double precision (a-h,o-z) 
      complex*16 H1,H2,H3,H4,H5,H6 
      complex*16 G1,G2,G3,G4,G5,G6  
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
      dimension G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G5(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G6(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY5(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY6(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      dimension Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi5(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi6(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
** 
      common /fillred/infilldim,infill(3) 
      dimension istorfill(3)
      dimension nphase(-1:1) 
      data nphase/-1,1,-1/ 
      parameter (pi   = 3.14159265358979324d0) 
*      print*,' eval1dhplatm1: y = ',y 
      if (infilldim.eq.2) then
         do i=1,2
            istorfill(i) = infill(i)
            infill(i) = -istorfill(i)
         enddo
      endif
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,Gi5,Gi6,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplatinf(-y,nw,G1,G2,G3,G4,G5,G6, 
     $                  GY1,GY2,GY3,GY4,GY5,GY6,
     $                  Gi1,Gi2,Gi3,Gi4,Gi5,Gi6,-n2,-n1) 
      if (infilldim.eq.2) then
         do i=1,2
            infill(i) = istorfill(i)
         enddo
      endif
** fill the arrays H's 
      do k1=n1,n2 
        nph1 = nphase(k1) 
        HY1(k1) =   nph1*GY1(-k1) 
        Hi1(k1) = - nph1*Gi1(-k1) 
        H1(k1)  = dcmplx(HY1(k1),Hi1(k1)*pi) 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            nph2 = nph1*nphase(k2) 
            HY2(k1,k2) =   nph2*GY2(-k1,-k2) 
            Hi2(k1,k2) = - nph2*Gi2(-k1,-k2) 
            H2(k1,k2)  = dcmplx(HY2(k1,k2),Hi2(k1,k2)*pi) 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                nph3 = nph2*nphase(k3) 
                HY3(k1,k2,k3) =   nph3*GY3(-k1,-k2,-k3) 
                Hi3(k1,k2,k3) = - nph3*Gi3(-k1,-k2,-k3) 
                H3(k1,k2,k3)  = dcmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi) 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    nph4 = nph3*nphase(k4) 
                    HY4(k1,k2,k3,k4) =   nph4*GY4(-k1,-k2,-k3,-k4) 
                    Hi4(k1,k2,k3,k4) = - nph4*Gi4(-k1,-k2,-k3,-k4) 
                    H4(k1,k2,k3,k4)  = 
     $                    dcmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi) 
                if ( nw.gt.4 ) then 
                  do k5=n1,n2 
                    nph5 = nph4*nphase(k5) 
                    HY5(k1,k2,k3,k4,k5) = nph5*GY5(-k1,-k2,-k3,-k4,-k5) 
                    Hi5(k1,k2,k3,k4,k5) =-nph5*Gi5(-k1,-k2,-k3,-k4,-k5) 
                    H5(k1,k2,k3,k4,k5)  = 
     $                dcmplx(HY5(k1,k2,k3,k4,k5),Hi5(k1,k2,k3,k4,k5)*pi) 
                if ( nw.gt.5 ) then 
                  do k6=n1,n2 
                    nph6 = nph5*nphase(k6) 
              HY6(k1,k2,k3,k4,k5,k6) = nph6*GY6(-k1,-k2,-k3,-k4,-k5,-k6) 
              Hi6(k1,k2,k3,k4,k5,k6) =-nph6*Gi6(-k1,-k2,-k3,-k4,-k5,-k6) 
              H6(k1,k2,k3,k4,k5,k6)  = 
     $         dcmplx(HY6(k1,k2,k3,k4,k5,k6),Hi6(k1,k2,k3,k4,k5,k6)*pi) 
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
