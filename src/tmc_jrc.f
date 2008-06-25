****************************************************************
*
*     Target Mass Correction computation.
*     The TMC factor is multiplied by the coefficient functions 
*     in the wcoeff.f subroutine
*
*     TMC(1) -> F2 TMC
*     TMC(2) -> FL TMC
*     TMC(3) -> F3 TMC
*
*     MU: 1/8/2007
*
****************************************************************

      SUBROUTINE COMPUTE_TMC(N,X,Q2,TMC)
*
      IMPLICIT none
*
      include '../../commons/tmc.h'
      include '../../commons/expptkin.h'
      include '../../commons/observ.h'
*
      DOUBLE COMPLEX N
      DOUBLE PRECISION X,Q2,MTGT,TAU
      DOUBLE COMPLEX TMC(3)
*
      DOUBLE PRECISION SQTAU,SQTAUI
*
      MTGT=MP
*
      TAU = 1.d0 + (4.d0*MTGT**2*xpt**2)/Q2
*      
      SQTAU  = DSQRT(TAU)
      SQTAUI = 1.d0/SQTAU
*

      TMC(1) = (1.d0+SQTAU)**2/(4*TAU**1.5d0)
     #     *(1.d0+(3.d0*(1.d0-SQTAUI))/(N+1.d0))
      TMC(2) = (1D0-TAU)*(1.d0+SQTAU)**2/(4*TAU**1.5d0)
     #     *(1.d0-((3.d0-TAU)*(1.d0+SQTAUI))/
     #     (4D0*(TAU**2)*(N+1.d0)))
      TMC(3) = (1.d0+SQTAU)/(2*TAU)
     #     *(1.d0+(2.d0*(1.d0-SQTAUI))/(N+1.d0))
*
      RETURN
      END

*---------------------------------------------------
