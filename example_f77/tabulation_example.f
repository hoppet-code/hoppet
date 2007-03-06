      program tabulation_example
      implicit none
      !--- variables defining the grid and evolution parameters
      double precision dy, ymax, Qmin, Qmax, dlnlnQ
      integer          nloop, order
      !---------------
      external         evolvePDF
      double precision x, Q, lhapdf(-6:6), ourpdf(-6:6)
      double precision asmz, Q0
      integer          iloop,  nf, iflv

      ! initialise an LHAPDF set
      call InitPDFsetByName("cteq61.LHgrid")
      call InitPDF(0)

      ! start the dglap evolution/convolution package 
      ymax  = 14d0      ! max value of ln 1/x
      dy    = 0.1d0     ! the internal grid spacing (smaller->higher accuarcy)
                        ! 0.1 should provide at least 10^{-3} accuracy 
      Qmin  = 1.0d0     ! smallest Q value in tabulatin
      Qmax  = 1d4       ! largest Q value in tabulatin
      dlnlnQ = 0.05d0   ! tabulation spacing in dlnlnQ (<=0.05 recommended)
      nloop = 2         ! the number of loops to initialise (max=3!)
      order = -5        ! numerical interpolation order (-5 is a good choice)
      
      ! call this once at the beginning of your program
      call dglapStartExtended(ymax,dy,Qmin,Qmax,dlnlnQ,nloop,order)

      ! initialise our PDF using the LHAPDF subroutine for PDF-access
      ! at the starting scale Q0, assuming as(mz)=0.118 and evolve to 
      ! the other scales. You can replace the LHAPDF routine with your
      ! own one (with the same interface)
      !
      ! Call this routine each time you want to change 
      ! input PDFs
      asmz = 0.118d0
      Q0   = 20.0d0
      call dglapEvolve(asmz,nloop,evolvePDF,Q0)

      x = 1d-6
      Q = 10000.0d0

      ! extract the pdf at this x,Q directly from LHAPDF
      call evolvePDF(x,Q,lhapdf)

      ! extract the pdf at this x,Q via our copy of the pdf
      call dglapEval(x,Q,ourpdf)


      ! print it all
      write(6,*) 'x = ', x, ',   Q = ', Q
      write(6,'(a5,4a17)') 'iflv','lhapdf','ourpdf'
      do iflv = -6, 6
         write(6,'(i5,4f17.7)') iflv,
     $        lhapdf(iflv),ourpdf(iflv)
      end do
      end
