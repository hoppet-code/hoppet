      subroutine fillred1dhpl(nw,H1,H2,H3,H4,H5,H6, 
     $                          HY1,HY2,HY3,HY4,HY5,HY6,
     $                          Hi1,Hi2,Hi3,Hi4,Hi5,Hi6,n1,n2)
* fills the reducible 1dhpl from the irreducible set
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
      common /fillred/infilldim,infill(3) 
      parameter (pinv = 0.318309886183790672d0) 
      parameter (pi   = 3.14159265358979324d0) 
** combining real and immaginary into the complex value 
      do k1=n1,n2 
      do k2=n1,n2 
        H2(k1,k2) = dcmplx(HY2(k1,k2),Hi2(k1,k2)*pi) 
        if ( nw.gt.2 ) then 
          do k3=n1,n2 
            H3(k1,k2,k3) = dcmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi) 
            if ( nw.gt.3 ) then 
              do k4=n1,n2 
                H4(k1,k2,k3,k4) = 
     $               dcmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi) 
            if ( nw.gt.4 ) then 
              do k5=n1,n2 
                H5(k1,k2,k3,k4,k5) = 
     $             dcmplx(HY5(k1,k2,k3,k4,k5),Hi5(k1,k2,k3,k4,k5)*pi) 
            if ( nw.gt.5 ) then 
              do k6=n1,n2 
                H6(k1,k2,k3,k4,k5,k6) = 
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
      enddo 
** evaluating the reduced HPL's 
** iflag = 0 to suppress auxiliary printings of FILLREDHPLx 
      iflag = 0 
      do ia =  1,infilldim 
      do ib = ia,infilldim 
        call FILLREDHPL2(iflag,H1,H2,n1,n2,infill(ia),infill(ib)) 
        if ( nw.gt.2 ) then 
          do ic = ib,infilldim 
            call FILLREDHPL3(iflag,H1,H2,H3,n1,n2, 
     $                          infill(ia),infill(ib),infill(ic)) 
            if ( nw.gt.3 ) then 
              do id = ic,infilldim 
                call FILLREDHPL4(iflag,H1,H2,H3,H4,n1,n2, 
     $               infill(ia),infill(ib),infill(ic),infill(id)) 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
      if (nw.gt.4) call FILLREDHPL5(iflag,H1,H2,H3,H4,H5,n1,n2)
      if (nw.gt.5) call FILLREDHPL6(iflag,H1,H2,H3,H4,H5,H6,n1,n2)
** extractin real and immaginary parts from the complex value 
      do k1=n1,n2 
      do k2=n1,n2 
        HY2(k1,k2) =  dble(H2(k1,k2)) 
        Hi2(k1,k2) = dimag(H2(k1,k2))*pinv 
        if ( nw.gt.2 ) then 
          do k3=n1,n2 
            HY3(k1,k2,k3) =  dble(H3(k1,k2,k3)) 
            Hi3(k1,k2,k3) = dimag(H3(k1,k2,k3))*pinv 
            if ( nw.gt.3 ) then 
              do k4=n1,n2 
                HY4(k1,k2,k3,k4) =  dble(H4(k1,k2,k3,k4)) 
                Hi4(k1,k2,k3,k4) = dimag(H4(k1,k2,k3,k4))*pinv 
            if ( nw.gt.4 ) then 
              do k5=n1,n2 
                HY5(k1,k2,k3,k4,k5) =  dble(H5(k1,k2,k3,k4,k5)) 
                Hi5(k1,k2,k3,k4,k5) = dimag(H5(k1,k2,k3,k4,k5))*pinv 
            if ( nw.gt.5 ) then 
              do k6=n1,n2 
              HY6(k1,k2,k3,k4,k5,k6) =  dble(H6(k1,k2,k3,k4,k5,k6)) 
              Hi6(k1,k2,k3,k4,k5,k6) = dimag(H6(k1,k2,k3,k4,k5,k6))*pinv 
              enddo 
            endif 
              enddo 
            endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
      return 
      end 
************************************************************************ 
