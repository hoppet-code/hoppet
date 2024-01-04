      program convolution_example
      implicit none
      !---------------
      external         evolvePDF
      double precision x, Q, lhapdf(-6:6), ourpdf(-6:6)
      double precision PLO_conv_ourpdf(-6:6), PNLO_conv_ourpdf(-6:6)
      double precision dy
      integer          iloop, nloop, nf, iflv

      ! initialise an LHAPDF set
      call InitPDFsetByName("cteq61.LHgrid")
      call InitPDF(0)

      ! start the hoppet evolution/convolution package 
      dy    = 0.1d0     ! the internal grid spacing (smaller->higher accuarcy)
                        ! 0.1 should provide at least 10^{-3} accuracy 
      nloop = 2         ! the number of loops to initialise (max=3!)
      call hoppetStart(dy, nloop)

      ! initialise our PDF using the LHAPDF subroutine for PDF-access
      ! (any other subroutine with same interface can be used in its place)
      call hoppetAssign(evolvePDF)

      x = 0.1d0
      Q = 30.0d0

      ! extract the pdf at this x,Q directly from LHAPDF
      call evolvePDF(x,Q,lhapdf)

      ! extract the pdf at this x,Q via our copy of the pdf
      call hoppetEval(x,Q,ourpdf)

      ! extract the convolution of the 1-loop (LO) splitting function
      ! with the currently stored pdf. The normalisation is such that to
      ! get the evolution with respect to ln(Q^2) one must multiply
      ! PLO_conv_ourpdf by alphas/(2*pi). This is the same thing that
      ! in hep-ph/0510324 is referred to as (P_0\otimes q)
      iloop = 1
      nf = 5
      call hoppetEvalSplit(x,Q,iloop,nf,PLO_conv_ourpdf)
      ! extract the convolution with the 2-loop (NLO) splitting function
      iloop = 2
      call hoppetEvalSplit(x,Q,iloop,nf,PNLO_conv_ourpdf)

      ! print it all
      write(6,*) 'x = ', x, ',   Q = ', Q
      write(6,'(a5,4a17)') 'iflv','lhapdf','ourpdf',
     $     'PLO_conv_ourpdf','PNLO_conv_ourpdf'
      do iflv = -6, 6
         write(6,'(i5,4f17.7)') iflv,
     $        lhapdf(iflv),ourpdf(iflv),
     $        PLO_conv_ourpdf(iflv),PNLO_conv_ourpdf(iflv)
      end do
      end
