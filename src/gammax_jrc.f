****************************************************************
*
*     File: gammax.f
*
*     Contains the functions to compute x*Gamma(x,Q2i,Q2f) as
*     Mellin inversions of the evolution operators computed 
*     in N-space.
*
*     AG: 26.01.07
*
****************************************************************
*
      FUNCTION gammax_nsvu(x,Q2i,Q2f)
      IMPLICIT none
      DOUBLE PRECISION gammax_nsvu
      DOUBLE PRECISION x,Q2i,Q2f
      DOUBLE PRECISION evpdf,A,alphas,AlphaExactBeta

      external AlphaExactBeta

      include "../../commons/alphas.h"
      include "../../commons/evflags.h"
      include "../../commons/evscale.h"
      include "../../commons/consts.h"
      include "../../commons/observ.h"

!     Initialize common blocks
      open(unit=10,status="old",file="../../commons/evol.par")
      open(unit=11,status="old",file="../../commons/alphas.par")
      read(10,*) Q20
      read(10,*) IPT
      read(10,*) IMODEV      
      read(10,*) IVFN
      read(10,*) ITMC
      read(10,*) A,A,A,A
      read(10,*) A
      read(10,*) A
      read(10,*) A
      read(11,*) ASREF,Q2REF
      read(11,*) Q2TH(4),Q2TH(5),Q2TH(6)
      !write(6,*) Q2TH(4),Q2TH(5),Q2TH(6)
      close(10)
      close(11)

!     Initialize constants
      CALL evolinit          ! Initialization of the evolution constants

!     Initialze coupling constant with same
!     values as the LH parameters
      ASREF=0.35d0
      Q2REF=2.0d0
      CALL FFN(q2f,alphas,ipt,AlphaExactBeta,5)
      ASQ = alphas/(4.d0*PI)
      !write(6,*) q2f,alphas

      !write(6,*) x,q2i,q2f

!     Set observable
      OBSERV=0

!     Compute the x-space nonsinglet evolution kernel
      call fixedtalbot_nsvu(evpdf,x,Q2i,Q2f)
!      gammax_nsvu = evpdf
!     Subtrcation at N=2
      gammax_nsvu = x*evpdf

      RETURN
      END
*
*-----------------------------------------------------------------
*
      FUNCTION gammax_nsm(x,Q2i,Q2f)
      IMPLICIT none
      DOUBLE PRECISION gammax_nsm
      DOUBLE PRECISION x,Q2i,Q2f
      DOUBLE PRECISION evpdf

      call fixedtalbot_nsm(evpdf,x,Q2i,Q2f)
      gammax_nsm = x*evpdf

      RETURN
      END
*
*-----------------------------------------------------------------
*
      FUNCTION gammax_nsp(x,Q2i,Q2f)
      IMPLICIT none
      DOUBLE PRECISION gammax_nsp
      DOUBLE PRECISION x,Q2i,Q2f
      DOUBLE PRECISION evpdf

      call fixedtalbot_nsp(evpdf,x,Q2i,Q2f)
      gammax_nsp = x*evpdf

      RETURN
      END
*
*-----------------------------------------------------------------
*
      FUNCTION gammax_nsvd(x,Q2i,Q2f)
      IMPLICIT none
      DOUBLE PRECISION gammax_nsvd
      DOUBLE PRECISION x,Q2i,Q2f
      DOUBLE PRECISION evpdf

      call fixedtalbot_nsvd(evpdf,x,Q2i,Q2f)

      gammax_nsvd = x*evpdf

      RETURN
      END
*
*-----------------------------------------------------------------
*

      FUNCTION gammax_ns24q(x,Q2i,Q2f)
      IMPLICIT none
      DOUBLE PRECISION gammax_ns24q
      DOUBLE PRECISION x,Q2i,Q2f
      DOUBLE PRECISION evpdf

      call fixedtalbot_ns24q(evpdf,x,Q2i,Q2f)

      gammax_ns24q = x*evpdf

      RETURN
      END
*
*-----------------------------------------------------------------
*
      FUNCTION gammax_ns24g(x,Q2i,Q2f)
      IMPLICIT none
      DOUBLE PRECISION gammax_ns24g
      DOUBLE PRECISION x,Q2i,Q2f
      DOUBLE PRECISION evpdf

      call fixedtalbot_ns24g(evpdf,x,Q2i,Q2f)

      gammax_ns24g = x*evpdf

      RETURN
      END
*
*-----------------------------------------------------------------
*
      FUNCTION gammax_qq(x,Q2i,Q2f)
      IMPLICIT none
      DOUBLE PRECISION gammax_qq
      DOUBLE PRECISION x,Q2i,Q2f
      DOUBLE PRECISION evpdf

      call fixedtalbot_qq(evpdf,x,Q2i,Q2f)

      gammax_qq = x*evpdf

      RETURN
      END
*-----------------------------------------------------------------
*
      FUNCTION gammax_qg(x,Q2i,Q2f)
      IMPLICIT none
      DOUBLE PRECISION gammax_qg
      DOUBLE PRECISION x,Q2i,Q2f
      DOUBLE PRECISION evpdf

      call fixedtalbot_qg(evpdf,x,Q2i,Q2f)

      gammax_qg = x*evpdf

      RETURN
      END
*
*-----------------------------------------------------------------
*
      FUNCTION gammax_gq(x,Q2i,Q2f)
      IMPLICIT none
      DOUBLE PRECISION gammax_gq
      DOUBLE PRECISION x,Q2i,Q2f
      DOUBLE PRECISION evpdf

      call fixedtalbot_gq(evpdf,x,Q2i,Q2f)

      gammax_gq = x*evpdf

      RETURN
      END
*
*-----------------------------------------------------------------
*
      FUNCTION gammax_gg(x,Q2i,Q2f)
      IMPLICIT none
      DOUBLE PRECISION gammax_gg
      DOUBLE PRECISION x,Q2i,Q2f
      DOUBLE PRECISION evpdf

      call fixedtalbot_gg(evpdf,x,Q2i,Q2f)

      gammax_gg = x*evpdf

      RETURN
      END

*
* ------- 1) Fixed Talbot 
*
      SUBROUTINE fixedtalbot_nsp(xfunc,x,Q2i,Q2f)
      IMPLICIT none

      double precision x,xfunc
      double precision Q2i,Q2f
      integer m
      double precision theta,pi,sigma,t,tmp,r
      double complex s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2)

      pi=acos(-1d0)
      m = 16
      t=-dlog(x)

      tmp  = 0d0
      tmp2 = 0d0

      r=2d0*m/5d0/t

      do theta=pi/m,(m-1d0)*pi/m,pi/m
        sigma=theta+(theta/tan(theta)-1d0)/tan(theta)
        s=r*theta*DCMPLX(1d0/tan(theta),1d0)
        s1 = s + (1d0,0d0)
        call ZFUNC(s1,x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG)
        tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfns(1)
        tmp=tmp+dreal(tmp2)
      enddo

      call zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG)

      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfns(1)))
      
      RETURN
      END
*
*
*
      SUBROUTINE fixedtalbot_nsm(xfunc,x,Q2i,Q2f)
      IMPLICIT none

      double precision x,xfunc
      double precision Q2i,Q2f
      integer m
      double precision theta,pi,sigma,t,tmp,r
      double complex s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2)

      pi=acos(-1d0)
      m = 16
         
      t=-dlog(x)

      tmp  = 0d0
      tmp2 = 0d0

      r=2d0*m/5d0/t

      do theta=pi/m,(m-1d0)*pi/m,pi/m
        sigma=theta+(theta/tan(theta)-1d0)/tan(theta)
        s=r*theta*DCMPLX(1d0/tan(theta),1d0)
        s1 = s + (1d0,0d0)
        call zfunc(s1,x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG)
        tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfns(2)
        tmp=tmp+dreal(tmp2)
      enddo
      call zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG)
      
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfns(2)))
      
      RETURN
      END
*
*
*
      SUBROUTINE fixedtalbot_nsvu(xfunc,x,Q2i,Q2f)
      IMPLICIT none

      double precision x,xfunc
      double precision Q2i,Q2f
      integer m
      double precision theta,pi,sigma,t,tmp,r
      double complex s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2)

      pi=acos(-1d0)
      m = 16
         
      t=-dlog(x)

      tmp  = 0d0
      tmp2 = 0d0

      r=2d0*m/5d0/t

      do theta=pi/m,(m-1d0)*pi/m,pi/m
        sigma=theta+(theta/tan(theta)-1d0)/tan(theta)
        s=r*theta*DCMPLX(1d0/tan(theta),1d0)
        s1 = s + (1d0,0d0)
        call zfunc(s1,x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG)
        tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfns(3)
        tmp=tmp+dreal(tmp2)
      enddo

      call zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG)
      
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfns(3)))
      
      RETURN
      END
*
*
*
      SUBROUTINE fixedtalbot_nsvd(xfunc,x,Q2i,Q2f)

      IMPLICIT none

      double precision x,xfunc
      double precision Q2i,Q2f
      integer m
      double precision theta,pi,sigma,t,tmp,r
      double complex s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2)

      pi=acos(-1d0)
      m = 16
         
      t=-dlog(x)

      tmp  = 0d0
      tmp2 = 0d0

      r=2d0*m/5d0/t

      do theta=pi/m,(m-1d0)*pi/m,pi/m
        sigma=theta+(theta/tan(theta)-1d0)/tan(theta)
        s=r*theta*DCMPLX(1d0/tan(theta),1d0)
        s1 = s + (1d0,0d0)
        call zfunc(s1,x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG)
        tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfns(4)
        tmp=tmp+dreal(tmp2)
      enddo
      call zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG)
      
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfns(4)))
      
      RETURN
      END
*
*
*
      SUBROUTINE fixedtalbot_ns24q(xfunc,x,Q2i,Q2f)
      IMPLICIT none

      double precision x,xfunc
      double precision Q2i,Q2f
      integer m
      double precision theta,pi,sigma,t,tmp,r
      double complex s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2)

      pi=acos(-1d0)
      m=16
         
      t=-dlog(x)

      tmp  = 0d0
      tmp2 = 0d0

      r=2d0*m/5d0/t

      do theta=pi/m,(m-1d0)*pi/m,pi/m
        sigma=theta+(theta/tan(theta)-1d0)/tan(theta)
        s=r*theta*DCMPLX(1d0/tan(theta),1d0)
        s1 = s + (1d0,0d0)
        call ZFUNC(s1,x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG)
        tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfns24(1)
        tmp=tmp+dreal(tmp2)
      enddo
      call zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG)
      
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfns24(1)))

      RETURN
      END
*
*
*
      SUBROUTINE fixedtalbot_ns24g(xfunc,x,Q2i,Q2f)

      IMPLICIT none

      double precision x,xfunc
      double precision Q2i,Q2f
      integer m
      double precision theta,pi,sigma,t,tmp,r
      double complex s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2)

      pi=acos(-1d0)
      m=16
         
      t=-dlog(x)

      tmp  = 0d0
      tmp2 = 0d0

      r=2d0*m/5d0/t

      do theta=pi/m,(m-1d0)*pi/m,pi/m
        sigma=theta+(theta/tan(theta)-1d0)/tan(theta)
        s=r*theta*DCMPLX(1d0/tan(theta),1d0)
        s1 = s + (1d0,0d0)
        call ZFUNC(s1,x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG)
        tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfns24(2)
        tmp=tmp+dreal(tmp2)
      enddo
      call zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,ZPDFNS,ZPDFNS24,ZPDFSG)
      
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfns24(2)))
      
      RETURN
      END
*
*
*
      SUBROUTINE fixedtalbot_qq(xfunc,x,Q2i,Q2f)
      implicit none

      double precision x,xfunc
      double precision Q2i,Q2f
      integer m
      double precision theta,pi,sigma,t,tmp,r
      double complex s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2)

      pi=acos(-1d0)
      m=16
         
      t=-dlog(x)

      tmp  = 0d0
      tmp2 = 0d0

      r=2d0*m/5d0/t

      do theta=pi/m,(m-1d0)*pi/m,pi/m
        sigma=theta+(theta/tan(theta)-1d0)/tan(theta)
        s=r*theta*DCMPLX(1d0/tan(theta),1d0)
        s1 = s + (1d0,0d0)
        call zfunc(s1,x,q2i,q2f,zpdfns,ZPDFNS24,zpdfsg)
        tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfsg(1,1)
        tmp=tmp+dreal(tmp2)
      enddo
      call zfunc(DCMPLX(r+1d0,0d0),x,q2i,q2f,zpdfns,ZPDFNS24,zpdfsg)
      
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfsg(1,1)))
      
      return
      end
*
*
*
      SUBROUTINE fixedtalbot_qg(xfunc,x,Q2i,Q2f)
      implicit none

      double precision x,xfunc
      double precision Q2i,Q2f
      integer m
      double precision theta,pi,sigma,t,tmp,r
      double complex s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2)

      pi=acos(-1d0)
      m=16
         
      t=-dlog(x)

      tmp  = 0d0
      tmp2 = 0d0

      r=2d0*m/5d0/t

      do theta=pi/m,(m-1d0)*pi/m,pi/m
        sigma=theta+(theta/tan(theta)-1d0)/tan(theta)
        s=r*theta*DCMPLX(1d0/tan(theta),1d0)
        s1 = s + (1d0,0d0)
        call zfunc(s1,x,q2i,q2f,zpdfns,ZPDFNS24,zpdfsg)
        tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfsg(1,2)
        tmp=tmp+dreal(tmp2)
      enddo
      call zfunc(DCMPLX(r+1d0,0d0),x,q2i,q2f,zpdfns,ZPDFNS24,zpdfsg)
      
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfsg(1,2)))
      
      return
      end
*
*
      SUBROUTINE fixedtalbot_gq(xfunc,x,Q2i,Q2f)
      implicit none

      double precision x,xfunc
      double precision Q2i,Q2f
      integer m
      double precision theta,pi,sigma,t,tmp,r
      double complex s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2)

      pi=acos(-1d0)
      m=16
         
      t=-dlog(x)

      tmp  = 0d0
      tmp2 = 0d0

      r=2d0*m/5d0/t

      do theta=pi/m,(m-1d0)*pi/m,pi/m
        sigma=theta+(theta/tan(theta)-1d0)/tan(theta)
        s=r*theta*DCMPLX(1d0/tan(theta),1d0)
        s1 = s + (1d0,0d0)
        call zfunc(s1,x,q2i,q2f,zpdfns,ZPDFNS24,zpdfsg)
        tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfsg(2,1)
        tmp=tmp+dreal(tmp2)
      enddo
      call zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,zpdfns,ZPDFNS24,zpdfsg)
      
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfsg(2,1)))
      
      return
      end
*
*
*
      SUBROUTINE fixedtalbot_gg(xfunc,x,Q2i,Q2f)
      implicit none

      double precision x,xfunc
      double precision Q2i,Q2f

      integer m
      double precision theta,pi,sigma,t,tmp,r
      double complex s,s1,tmp2,zpdfns(4),zpdfsg(2,2),ZPDFNS24(2)

      pi=acos(-1d0)
      m=16
         
      t=-dlog(x)

      tmp  = 0d0
      tmp2 = 0d0

      r=2d0*m/5d0/t

      do theta=pi/m,(m-1d0)*pi/m,pi/m
        sigma=theta+(theta/tan(theta)-1d0)/tan(theta)
        s=r*theta*DCMPLX(1d0/tan(theta),1d0)
        s1 = s + (1d0,0d0)
        call zfunc(s1,x,Q2i,Q2f,zpdfns,ZPDFNS24,zpdfsg)
        tmp2=exp(t*s1)*DCMPLX(1d0,sigma)*zpdfsg(2,2)
        tmp=tmp+dreal(tmp2)
      enddo
      call zfunc(DCMPLX(r+1d0,0d0),x,Q2i,Q2f,zpdfns,ZPDFNS24,zpdfsg)
      
      xfunc=r/m*(tmp+.5d0*exp((r+1d0)*t)*dreal(zpdfsg(2,2)))
      
      return
      end


*-------------------------------------------------------------
