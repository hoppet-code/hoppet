C----------------------------------------------------------------------
C
C  An f77 example program using a tabulation. It outputs a subset of
C  table 15 of hep-ph/0511119 and this output should be identical
C  to the contents of the file tabulation_example.default_output
C 
C  NB: for the full functionality used in generating the HeraLHC and
C      Les Houches comparison tables, see ../benchmarks/benchmarks.f90
C      and carefully read the comments at the beginning. Subtleties
C      exist in particular wrt the treatment of scales muF/=muR.
C
C  This program can be compiled with a plain f77 compiler (e.g. g77),
C  but then the appropriate f95 libraries must be included at
C  link-time. The Makefile contains illustrations of which libraries are
C  needed, for a range of compilers
C
C----------------------------------------------------------------------
      program tabulation_example
      implicit none
      !--- variables defining the grid and evolution parameters
      double precision dy, ymax, Qmin, Qmax, dlnlnQ
      integer          nloop, order
      !---------------
      external         heralhc_init
      double precision Q, ourpdf(-6:6)
      double precision asQ0, Q0, xmu
      integer          ix, scheme
      double precision heralhc_xvals(9)
      data heralhc_xvals/1d-5,1d-4,1d-3,1d-2,0.1d0,
     $                   0.3d0,0.5d0,0.7d0,0.9d0/

      ! start the dglap evolution/convolution package 
      ymax  = 12d0      ! max value of ln 1/x
      dy    = 0.1d0     ! the internal grid spacing (smaller->higher accuarcy)
                        ! 0.1 should provide at least 10^{-3} accuracy 
      Qmin  = 1.0d0     ! smallest Q value in tabulation
      Qmax  = 1d4       ! largest Q value in tabulation
      dlnlnQ = dy/4     ! tabulation spacing in dlnlnQ (dy/4 recommended)
      nloop  = 3        ! the number of loops to initialise (max=3!)
      order  = -6       ! numerical interpolation order (-6 is a good choice)
      scheme = 1        ! 1=unpol-MSbar, 2=unpol-DIS, 3=Pol-MSbar

      ! call this once at the beginning of your program
      call hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order,
     $                         scheme)
      write(6,'(a)') "Splitting functions initialised!"

      ! tell hoppet to use a variable flavour number scheme with
      ! the following c,b,t quark masses
      call hoppetSetVFN(1.414213563d0,4.5d0,175d0)

      ! set the initial scale and the coupling there (in general the PDF
      ! and coupling may be specified at different scales -- this is not
      ! done here) and 
      Q0   = sqrt(2d0)
      asQ0 = 0.35d0
      ! the ratio xmu = mu_F/mu_R to be used in the evolution.
      xmu  = 1.0d0
      
      ! carry out the evolution to create a tabulation, corresponding
      ! to the initial condition heralhc_init(...) given below
      call hoppetEvolve(asQ0,Q0,nloop,xmu,heralhc_init,Q0)
      ! alternatively if you need to repeat the same evolution on very
      ! many pdf sets, used a cached evolution (set up once, use many times
      ! and gain a factor 3-4 in speed.
      !call hoppetPreEvolve(asQ0,Q0,nloop,xmu,Q0)
      !call hoppetCachedEvolve(heralhc_init) 
      write(6,'(a)') "Evolution done!"

      ! print out some results
      Q = 100.0d0
      write(6,'(a)')
      write(6,'(a,f8.3,a)')"           Evaluating PDFs at Q = ",Q," GeV"
      write(6,'(a5,2a12,a14,a10,a12)') "x",
     $     "u-ubar","d-dbar","2(ubr+dbr)","c+cbar","gluon"
      do ix = 1, 9
         call hoppetEval(heralhc_xvals(ix),Q,ourpdf)
         write(6,'(1p,e7.1,5e12.4)') heralhc_xvals(ix),
     $        ourpdf(2)-ourpdf(-2), 
     $        ourpdf(1)-ourpdf(-1), 
     $        2*(ourpdf(-1)+ourpdf(-2)),
     $        (ourpdf(-4)+ourpdf(4)),
     $        ourpdf(0)
      end do

      ! free all memory that was allocated by hoppet
      call hoppetDeleteAll()
      
      end


      !--------------------------------------------------------------
      ! The initial condition used in hep-ph/0511119.
      !
      ! The subroutine returns x*q(x) in the array pdf(-6:6), containing
      ! flavours (tbar,bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b,t)
      !
      ! The variable Q is set by the calling routine to be equal to the
      ! scale at which the PDF is being requested. Below it is unused
      ! because the initial condition is here provided only for scale
      ! Q0=sqrt(2)GeV.
      !
      subroutine heralhc_init(x,Q,pdf)
      double precision x,Q,pdf(-6:6)
      double precision uv, dv
      double precision ubar, dbar
      !---------------------
      double precision N_g, N_ls, N_uv, N_dv, N_db
      data N_g,N_ls,N_uv,N_dv/1.7d0, 0.387975d0, 5.107200d0,3.064320d0/

      N_db = 0.5d0 * N_ls
      uv = N_uv * x**0.8d0 * (1d0-x)**3
      dv = N_dv * x**0.8d0 * (1d0-x)**4
      dbar = N_db * x**(-0.1d0) * (1d0-x)**6
      ubar = dbar * (1d0-x)

      pdf(0) = N_g * x**(-0.1d0) * (1d0-x)**5
      pdf(-3) = 0.2d0*(dbar + ubar)
      pdf( 3) = pdf(-3)
      pdf( 2) = uv + ubar
      pdf(-2) = ubar
      pdf( 1) = dv + dbar
      pdf(-1) = dbar

      pdf( 4) = 0
      pdf( 5) = 0
      pdf( 6) = 0
      pdf(-4) = 0
      pdf(-5) = 0
      pdf(-6) = 0
      end 
