! -*- mode: fortran -*-
! This is a Fortran free-form include file

#ifndef __HOPPET_PDF_TABULATE_INC__
#define __HOPPET_PDF_TABULATE_INC__

! the minimum ny we expect in a grid or subgrid in order for hard-coded
! interpolation assumptions to correct; provided as a define so that it
! can be checked at preprocessor level in "templated" code in
! pdf_tabulate_OrderNNN.F90

#define __HOPPET_tab_min_subgrid_ny__ 6
#define __HOPPET_tab_min_nQ__         5

#endif