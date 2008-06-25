**********************************************************************
*
*     File: pdfcomp.f
*
*     The subroutine PDFCOMP computes the PDFs of the single
*     partons starting from the T_i and V_i combinations
*     which are used in the evolution.
*
*     The parton densities are computing inverting eq.(15) in
*     the "Notes on Perturbative Evolution".
*
*     The combinations used in the evolution are numbered 
*     according to the following table:
*
*       1    2   3   4   5   6   7   8    9    10   11    12    13
*     Sigma  g  u_v d_v s_v c_v b_v t_v  T_3  T_8  T_15  T_24  T_35
*
*     The output PDFs are numbered according to the Les Houches
*     convention:
*
*      -6   -5   -4   -3   -2   -1   0   1   2   3   4   5   6 
*     tbar bbar cbar sbar ubar dbar  g   d   u   s   c   b   t
*
*     AG: 05.03.2007
*
*********************************************************************
*
      SUBROUTINE PDFCOMP(pdftmp,pdfout)
*
      IMPLICIT none
*
      include "../../commons/evflags.h"
*
      DOUBLE PRECISION pdftmp(13),pdfout(-6:6)
*
      IF(IVFN.EQ.1) THEN
*         
        pdfout(0) = pdftmp(2)      
*
        pdfout(1)=(10.d0*pdftmp(1) + 5.d0*pdftmp(11) + 3.d0*pdftmp(12) 
     #       - 30.d0*pdftmp(9) + 2.d0*pdftmp(13) + 10.d0*pdftmp(10)
     #       + 60.d0*pdftmp(4))/120.d0
*     
        pdfout(-1)=(10.d0*pdftmp(1) + 5.d0*pdftmp(11) + 3.d0*pdftmp(12)
     #       - 30.d0*pdftmp(9) + 2.d0*pdftmp(13) + 10.d0*pdftmp(10)
     #       - 60.d0*pdftmp(4))/120.d0
*     
        pdfout(2)=(10.d0*pdftmp(1) + 5.d0*pdftmp(11) + 3.d0*pdftmp(12) 
     #       + 30.d0*pdftmp(9) + 2.d0*pdftmp(13) + 10.d0*pdftmp(10)
     #       + 60.d0*pdftmp(3))/120.d0
*     
        pdfout(-2)=(10.d0*pdftmp(1) + 5.d0*pdftmp(11) + 3.d0*pdftmp(12)
     #       + 30.d0*pdftmp(9) + 2.d0*pdftmp(13) + 10.d0*pdftmp(10)
     #       - 60.d0*pdftmp(3))/120.d0
*     
        pdfout(3)=(10.d0*pdftmp(1) + 5.d0*pdftmp(11) + 3.d0*pdftmp(12)
     #       + 2.d0*pdftmp(13) - 20.d0*pdftmp(10) + 60.d0*pdftmp(5))
     #       /120.d0
*     
        pdfout(-3)=(10.d0*pdftmp(1) + 5.d0*pdftmp(11) + 3.d0*pdftmp(12) 
     #       + 2.d0*pdftmp(13) - 20.d0*pdftmp(10) - 60.d0*pdftmp(5))
     #       /120.d0
*     
        pdfout(4)=(10.d0*pdftmp(1) - 15.d0*pdftmp(11) + 3.d0*pdftmp(12) 
     #       + 2.d0*pdftmp(13) + 60.d0*pdftmp(6))/120.d0
*     
        pdfout(-4)=(10.d0*pdftmp(1) - 15.d0*pdftmp(11) + 3.d0*pdftmp(12)
     #       + 2.d0*pdftmp(13) - 60.*pdftmp(6))/120.d0
*     
         pdfout(5)=(5.d0*pdftmp(1) - 6.d0*pdftmp(12) + pdftmp(13) 
     #        + 30.d0*pdftmp(7))/60.d0
*     
         pdfout(-5)=(5.d0*pdftmp(1) - 6.d0*pdftmp(12) + pdftmp(13)
     #        - 30.d0*pdftmp(7))/60.d0         
*     
         pdfout(6) = (pdftmp(1) - pdftmp(13) + 6.d0*pdftmp(8))/12.d0
*     
         pdfout(-6) = (pdftmp(1) - pdftmp(13) - 6.d0*pdftmp(8))/12.d0
         
      ELSE
*
         pdfout(0) = pdftmp(2)
*     
         pdfout(1)=(3.d0*pdftmp(1) + pdftmp(11) - 6.d0*pdftmp(9) 
     #              + 2.d0*pdftmp(10)+ 12.d0*pdftmp(4))/24.d0
*     
         pdfout(-1)=(3.d0*pdftmp(1) + pdftmp(11) - 6.d0*pdftmp(9)
     #               + 2.d0*pdftmp(10) - 12.d0*pdftmp(4))/24.d0
*     
         pdfout(2)=(3.d0*pdftmp(1) + pdftmp(11) + 6.d0*pdftmp(9)
     #              + 2.d0*pdftmp(10) + 12.d0*pdftmp(3))/24.d0
*     
         pdfout(-2)=(3.d0*pdftmp(1) + pdftmp(11) + 6.d0*pdftmp(9)
     #               + 2.d0*pdftmp(10) - 12.d0*pdftmp(3))/24.d0
*     
         pdfout(3)=(3.d0*pdftmp(1) + pdftmp(11) - 4.d0*pdftmp(10) 
     #              + 12.d0*pdftmp(5))/24.d0
*     
         pdfout(-3)=(3.d0*pdftmp(1) + pdftmp(11) - 4.d0*pdftmp(10) 
     #               - 12.d0*pdftmp(5))/24.d0
*     
         pdfout(4)=(pdftmp(1) - pdftmp(11) + 4.d0*pdftmp(6))/8.d0
*     
         pdfout(-4)=(pdftmp(1) - pdftmp(11) -  4.d0*pdftmp(6))/8.d0
*     
         pdfout(5)  = 0.d0
         pdfout(-5) = 0.d0
         pdfout(6)  = 0.d0
         pdfout(-6) = 0.d0         
*
      ENDIF
*
      RETURN
      END
*     
*=========================================================================


      subroutine PDFCONVINV(PDFLH,PDFTMP)
      implicit none
***************************
*     Convert from the MLM convention for PDF ordering:

*      -6   -5   -4   -3   -2   -1   0   1   2   3   4   5   6 
*     tbar bbar cbar sbar ubar dbar  g   d   u   s   c   b   t

*     to the  convention required for the evolint.f routine:

*       1    2   3   4   5   6   7   8    9    10   11    12    13
*     Sigma  g  u_v d_v s_v c_v b_v t_v  T_3  T_8  T_15  T_24  T_35
****************************

      integer IPDF
      double precision PDFLH(-6:6),PDFTMP(13)

      PDFTMP(1)=0
      do IPDF=1,6
         PDFTMP(1)=PDFTMP(1)+(PDFLH(IPDF)+PDFLH(-IPDF))
      enddo

      PDFTMP(2)=PDFLH(0)

      PDFTMP(3)=PDFLH(2)-PDFLH(-2)
      PDFTMP(4)=PDFLH(1)-PDFLH(-1)
      PDFTMP(5)=PDFLH(3)-PDFLH(-3)
      PDFTMP(6)=PDFLH(4)-PDFLH(-4)
      PDFTMP(7)=PDFLH(5)-PDFLH(-5)
      PDFTMP(8)=PDFLH(6)-PDFLH(-6)

      PDFTMP(9)=(PDFLH(2)+PDFLH(-2)) - (PDFLH(1)+PDFLH(-1))
      PDFTMP(10)=(PDFLH(2)+PDFLH(-2)) + (PDFLH(1)+PDFLH(-1))
     #           -2.d0*(PDFLH(3)+PDFLH(-3))
      PDFTMP(11)=(PDFLH(2)+PDFLH(-2)) + (PDFLH(1)+PDFLH(-1))
     #          + (PDFLH(3)+PDFLH(-3)) - 3.d0*(PDFLH(4)+PDFLH(-4))
      PDFTMP(12)=(PDFLH(2)+PDFLH(-2)) + (PDFLH(1)+PDFLH(-1))
     #          +(PDFLH(3)+PDFLH(-3)) + (PDFLH(4)+PDFLH(-4))
     #          - 4.d0*(PDFLH(5)+PDFLH(-5))
      PDFTMP(13)=(PDFLH(2)+PDFLH(-2)) + (PDFLH(1)+PDFLH(-1))
     #          +(PDFLH(3)+PDFLH(-3)) + (PDFLH(4)+PDFLH(-4))
     #          +(PDFLH(5)+PDFLH(-5)) - 5.d0*(PDFLH(6)+PDFLH(-6))

      return
      end

*------------------------------------------------------
