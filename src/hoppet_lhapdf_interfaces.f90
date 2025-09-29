!! Module defining interfaces to LHAPDF Fortran routines. 
!!
!! Only a subset of routines is included here.
module hoppet_lhapdf_interfaces
  implicit none

  ! interface to the LHAPDF routine for accessing the PDF
  interface
    subroutine EvolvePDF(x,Q,res)
      use types; implicit none
      real(dp), intent(in)  :: x,Q
      real(dp), intent(out) :: res(*)
    end subroutine EvolvePDF
  end interface

  interface
    function alphasPDF(Q) result(as)
      use types; implicit none
      real(dp), intent(in) :: Q
      real(dp) :: as
    end function alphasPDF
  end interface

  !! Get the number of error members in the set (non-multiset version)
  interface 
    subroutine numberpdf(n)
      integer, intent(out) :: n
    end subroutine numberpdf
  end interface

  interface 
    subroutine initPDFSetByName(name)
      character(len=*), intent(in) :: name
    end subroutine initPDFSetByName
  end interface

  interface 
    subroutine initPDF(imem)
      integer, intent(in) :: imem
    end subroutine initPDF
  end interface

  interface 
    subroutine getQ2min(imem,Q2min)
      use types; implicit none
      integer, intent(in) :: imem
      real(dp), intent(out) :: Q2min
    end subroutine getQ2min
  end interface
  interface 
    subroutine getQ2max(imem,Q2max)
      use types; implicit none
      integer, intent(in) :: imem
      real(dp), intent(out) :: Q2max
    end subroutine getQ2max
  end interface

  interface 
    subroutine getxmin(imem,xmin)
      use types; implicit none
      integer, intent(in) :: imem
      real(dp), intent(out) :: xmin
    end subroutine getxmin
  end interface
  
  interface 
    subroutine getxmax(imem,xmax)
      use types; implicit none
      integer, intent(in) :: imem
      real(dp), intent(out) :: xmax
    end subroutine getxmax
  end interface

  interface 
    subroutine getorderas(order)
      integer, intent(out) :: order
    end subroutine getorderas
  end interface

  interface 
    subroutine getthreshold(nf, mquark)
      use types; implicit none
      integer, intent(in) :: nf
      real(dp), intent(out) :: mquark
    end subroutine getthreshold
  end interface

  interface 
    subroutine getnf(nf)
      integer, intent(out) :: nf
    end subroutine getnf
  end interface

  interface 
    subroutine GetPDFUncertainty(pdfset, central, err_up, err_down, err_symm)
      use types; implicit none
      real(dp), intent(in)  :: pdfset(*)
      real(dp), intent(out) :: central, err_up, err_down, err_symm
    end subroutine GetPDFUncertainty
  end interface

end module hoppet_lhapdf_interfaces
