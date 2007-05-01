C----------------------------------------------------------------------
C
C Program that takes an initial condition from LHAPDF at some scale,
C creates a hoppet tabulation by evolving it from that scale, and
C then compares the result at some other scale.
C
C This provdes a way of cross checking the evolution implicitly
C contained in public PDF sets
C 
C----------------------------------------------------------------------
      program compare_lhapdf_hoppet
      implicit none
      !--- variables defining the grid and evolution parameters
      double precision dy, ymax, Qmin, Qmax, dlnlnQ
      integer          nloop, order
      !---------------
      external         evolvePDF
      double precision x, Q, lhapdf(-6:6), ourpdf(-6:6)
      double precision asmz, Q0, xmu
      integer          iloop,  nf, iflv, i, scheme

      ! initialise an LHAPDF set
      call InitPDFsetByName("cteq61.LHgrid")
      call InitPDF(0)

      ! start the dglap evolution/convolution package 
      ymax  = 14d0      ! max value of ln 1/x
      dy    = 0.1d0     ! the internal grid spacing (smaller->higher accuarcy)
                        ! 0.1 should provide at least 10^{-3} accuracy 
      Qmin  = 1.0d0     ! smallest Q value in tabulation
      Qmax  = 1d4       ! largest Q value in tabulation
      dlnlnQ = 0.025d0  ! tabulation spacing in dlnlnQ (dy/4 recommended)
      nloop = 2         ! the number of loops to initialise (max=3!)
      order = -6        ! numerical interpolation order (-6 is a good choice)
      scheme = 1        ! 1=unpol-MSbar, 2=unpol-DIS, 3=Pol-MSbar

      ! call this once at the beginning of your program
      call hoppetStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order,
     $                         scheme)

      ! initialise our PDF using the LHAPDF subroutine for PDF-access
      ! at the starting scale Q0, assuming as(mz)=0.118 and evolve to 
      ! the other scales. You can replace the LHAPDF routine with your
      ! own one (with the same interface)
      !
      ! Call this routine each time you want to change 
      ! input PDFs
      asmz = 0.118d0
      Q0   = 20.0d0
      xmu  = 1.0d0
      !call hoppetSetVFN(1.5d0,4.7d0,175d0)
      call hoppetEvolve(asmz,91.2d0,nloop,xmu,evolvePDF,Q0)
      ! alternatively if you need to repeat the same evolution on very
      ! many pdf sets, used a cached evolution (set up once, use many times
      ! and gain a factor 3-4 in speed.
      !call hoppetPreEvolve(asmz,91.2d0,nloop,xmu,Q0)
      !call hoppetCachedEvolve(evolvePDF) 
      x = 1d-3
      Q = 100.0d0

      ! extract the pdf at this x,Q directly from LHAPDF
      call evolvePDF(x,Q,lhapdf)

      ! extract the pdf at this x,Q via our copy of the pdf
      call hoppetEval(x,Q,ourpdf)

      ! print it all
      write(6,*) 'x = ', x, ',   Q = ', Q
      write(6,'(a5,4a17)') 'iflv','lhapdf','hoppet'
      do iflv = -6, 6
         write(6,'(i5,4f17.7)') iflv,
     $        lhapdf(iflv),ourpdf(iflv)
      end do
      end
