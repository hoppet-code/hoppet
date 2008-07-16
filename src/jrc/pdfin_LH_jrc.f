**********************************************************************
*
*     File pdfin_LH.f
*
*     Input PDFs in x-space, for the Les Houches comparison
*     hep-ph/0204316
*
*     AG: 20.03.2007
*
***********************************************************************

      SUBROUTINE PDFINLH(x,pdf)
*
      IMPLICIT none
*
      DOUBLE PRECISION x,pdf(13)
*
      DOUBLE PRECISION uv,dv,ubar,dbar,ss,sbar
*
      DOUBLE PRECISION a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12

*     PDF parameters
      a1=1.700000d0  
      a2=-1.1d0  
      a3=5d0	
      a4=5.1072d0    
      a5=-0.2d0  
      a6=3d0 	
      a7=3.06432d0   
      a8=-0.2d0  
      a9=4d0	
      a10=0.1939875d0 
      a11=-1.1d0  
      a12=6d0 
*
      uv = a4*(x**a5)*(1d0-x)**a6
      dv = a7*(x**a8)*(1d0-x)**a9

      dbar = a10*(x**a11)*(1d0-x)**a12
      ubar = (1.d0-x)*dbar

      ss   = 0.2d0*(ubar+dbar)
      sbar = ss
*
*     Sigma=uv+dv+2(ubar+dbar)+ s+sbar
*     with LH assumptions it reduces to uv+dv+2.4*(2-x)*dbar

      pdf(1)  = uv + dv + 2.d0*(ubar+dbar) + (ss + sbar)

*     Gluon

      pdf(2)  = a1*(x**a2)*(1d0-x)**a3
*
*     Non Singlet distributions
*
      pdf(3)  = uv
      pdf(4)  = dv
      pdf(5)  = 0.d0
      pdf(6)  = 0.d0
      pdf(7)  = 0.d0
      pdf(8)  = 0.d0
*     T3= uv - dv -2*x*dbar with LH assumptions
      pdf(9)  = uv + 2.d0*ubar - (dv + 2.d0*dbar)
      pdf(10) = uv + 2.d0*ubar + dv + 2.d0*dbar - 2.d0 * (ss + sbar)
      pdf(11) = pdf(1)
      pdf(12) = pdf(1)
      pdf(13) = pdf(1)

*
      RETURN
      END
*
*=======================================================================
