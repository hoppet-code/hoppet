
****************************************************************
*
*     File: alphas.f
*     
*     Routines related to the computation of alpha_s.
*
*     FFN: returns the value of alpha_s computed for Fixed
*          Flavour Number.
*     VFN: returns the value of alpha_s computed for Variable
*          Flavour Number.
*     
*     Input:
*     Arguments: q2th,mz2,asz,qq2,ipt
*     Functions: FUNC that evaluates alphas
*
*     Output: alphas
*
****************************************************************
*
      SUBROUTINE FFN(q2,alphas,ipt,FUNC,nf)
*
      IMPLICIT none

*
      include "../../commons/consts.h"
      include "../../commons/alphas.h"

      INTEGER nf
      INTEGER ipt
      DOUBLE PRECISION q2,qq2
      DOUBLE PRECISION alphas,func,alphasref,qq2ref

      alphasref=asref
      qq2ref=q2ref
      qq2=q2

      alphas=FUNC(nf,qq2ref,alphasref,qq2,ipt)

      RETURN
      END

*
      SUBROUTINE VFN(qq2,alphas,ipt,FUNC)
*
      IMPLICIT none
*
      include "../../commons/consts.h"
      include "../../commons/alphas.h"
*
      INTEGER nfi,nff,dnf,snf,ipt
      double precision q2,qq2,FUNC,alphasref,qq2ref
      double precision alphas,c2,asi
*     c2 is the same coefficient used in eq. 2.42
*     of hep-ph/0408244 and obtained in eq. 10 
*     of hep-ph/9706430. In the following it is divided
*     by (4*pi)^2 to match the notations
      parameter(c2=14d0/3d0)
      external FUNC

      q2=qq2

      if(q2.gt.q2th(6))then
         nff=6
      elseif(q2.gt.q2th(5))then
         nff=5
      elseif(q2.gt.q2th(4))then
         nff=4
      else
         nff=3
      endif
      if(q2ref.gt.q2th(6))then
         nfi=6
      elseif(q2ref.gt.q2th(5))then
         nfi=5
      elseif(q2ref.gt.q2th(4))then
         nfi=4
      else
         nfi=3
      endif

      alphasref = asref
      qq2ref = q2ref

 10   if(nff.eq.nfi) then
         alphas = FUNC(nfi,qq2ref,alphasref,q2,ipt)
         return
      else
         if(nff.gt.nfi)then
            dnf=1
            snf=1
         else
            dnf=-1
            snf=0
         endif
         asi = FUNC(nfi,qq2ref,alphasref,q2th(nfi+snf),ipt)
         if(ipt.ge.2)then
            if(nff.gt.nfi) asi=asi+(c2/(4d0*pi)**2d0)*asi**3d0
            if(nff.lt.nfi) asi=asi-(c2/(4d0*pi)**2d0)*asi**3d0
         endif
         alphasref = asi
         qq2ref = q2th(nfi+snf)
         nfi = nfi+dnf
         goto 10
      endif   
      end
*
*     Routines for the computation of alpha_s with fixed number of 
*     flavours
* 
*     - AlphaMZ: computes alpha_s as function of alpha_s at a given
*                refernce scale.
*
*     - AlphaExactBeta: Exact solution of the QCD beta function 
*                       equation using fourth order Runge-Kutta 
*                       algorithm.
*                       (Used in PEGASUS and Les Houches tables)
*
*     - AlphaLambda: computes alpha_s as a function of Lambda_QCD
*      
      FUNCTION AlphaMZ(nfi,mz2,asz,q2,ipt)

      IMPLICIT none

      include "../../commons/consts.h"
      include "../../commons/beta.h"

      integer nfi,ipt
      double precision q2ref,q2,asi,asref
      double precision alo,t,as,den,mz2,asz,AlphaMZ

      q2ref=mz2
      asref=asz/4d0/pi

      asi=asref
      t=log(q2/q2ref)
      den=1+beta0(nfi)*asi*t
      alo=asi/den
*
*     LO
*
      as=alo
*
*     NLO
*
      if(ipt.ge.1)as=alo*(1-b1(nfi)*alo*log(den))
*
*     NNLO
*
      if(ipt.ge.2)then
         as=alo*(1d0+(alo*(alo-asi)*(b2(nfi)-b1(nfi)**2d0)
     #        +as*b1(nfi)*dlog(as/asi)))
      endif
*
      AlphaMZ=4d0*pi*as
*
      RETURN
      END
*
*
*
      FUNCTION AlphaLambda(nf,q2,ipt)
*
      IMPLICIT none
*
      include "../../commons/beta.h"
      include "../../commons/consts.h"
*
      integer nf,ipt
      double precision AlphaLambda
      double precision qcdlam2,q2,as,as0
      double precision logt,t,xlam
      parameter(qcdlam2=2.6908239d-2)
*
      xlam = 0.d0
      if(nf.eq.4)xlam=.294d0
      if(nf.eq.5)xlam=.210d0
*
      as=0d0
      t=log(q2/xlam**2d0)
      logt=dlog(t)
*
*     LO
*
      as0=1d0/(beta0(nf)*t)
      as=as0
*
*     NLO
*
      if(ipt.eq.1) then
         as=(1d0-b1(nf)*logt/(beta0(nf)*t))/beta0(nf)/t
      endif
*      
*     NNLO
*
      if(ipt.eq.2) then
         as=(1d0-b1(nf)*logt/(beta0(nf)*t)+((b1(nf)**2d0*(dlog(t))**2d0
     #        -b1(nf)**2d0*dlog(t)+b2(nf)-b1(nf)**2d0)/
     #        (beta0(nf)**2d0*t**2d0)))
     #        /beta0(nf)/t
      endif
*
      AlphaLambda=4d0*pi*as

      RETURN
      END
*
*
*

*
      FUNCTION AlphaExactBeta(nf,r20,as0,r2,ipt)
*
      include "../../commons/beta.h" 
      include "../../commons/consts.h"
      
      integer NFMIN, NFMAX, NF, NSTEP, K1,ipt,nnf
      double precision fbeta
      double precision as,as0,r20,r2
      double precision dlr,lrrat,sxth
      double precision xk0,xk1,xk2,xk3

      PARAMETER (NFMIN = 3, NFMAX = 6)
      parameter(NSTEP=50)
      PARAMETER (SXTH = 0.16666 66666 66666 D0 )
*     
*     ..The beta functions FBETAn at N^nLO for n = 1, 2, and 3
*     
      
      NNF = NF

      AS0 = AS0/4d0/pi
      
      AS = AS0
      LRRAT = dLOG (R2/R20)
      DLR = LRRAT / NSTEP
*     
*     ..Solution of the evolution equation depending on  NAORD
*   (fourth-order Runge-Kutta beyond the leading order)
*     
      IF (IPT.EQ.0) THEN
*     
         AS = AS0 / (1d0+ BETA0(NF) * AS0 * LRRAT)
*     
      ELSE IF (IPT.EQ.1) THEN
*     
         DO 2 K1 = 1, NSTEP
            XK0 = DLR * FBETA (AS,NNF,IPT)
            XK1 = DLR * FBETA (AS + 0.5d0 * XK0,NNF,IPT)
            XK2 = DLR * FBETA (AS + 0.5d0 * XK1,NNF,IPT)
            XK3 = DLR * FBETA (AS + XK2,NNF,IPT)
            AS = AS + SXTH * (XK0 + 2d0* XK1 + 2d0* XK2 + XK3)
  2       CONTINUE
*     
      ELSE IF (IPT.EQ.2) THEN
*     
         DO 3 K1 = 1, NSTEP
            XK0 = DLR * FBETA (AS,nnf,ipt)
            XK1 = DLR * FBETA (AS + 0.5d0 * XK0,nnf,ipt)
            XK2 = DLR * FBETA (AS + 0.5d0 * XK1,nnf,ipt)
            XK3 = DLR * FBETA (AS + XK2,nnf,ipt)
            AS = AS + SXTH * (XK0 + 2d0* XK1 + 2d0* XK2 + XK3)
 3       CONTINUE
*     
      END IF

      alphaexactbeta = as*(4d0*pi)

      RETURN
      END
*     
*     =================================================================av==


      
      function fbeta(a,nf,ipt)
      
      implicit none
      include "../../commons/beta.h"
      double precision fbeta,a
      integer nf,ipt
      
      if(ipt.eq.1)then
         FBETA = - A**2d0 * ( BETA0(NF) + A * BETA1(NF) )
      elseif(ipt.eq.2)then
         FBETA = - A**2d0 * ( BETA0(NF) + A * ( BETA1(NF)
     ,     + A * BETA2(NF) ) )
      endif

      return
      end
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCcccccccc
