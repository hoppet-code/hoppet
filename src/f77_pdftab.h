#ifndef __F77_PDFTAB_H__
#define __F77_PDFTAB_H__

extern "C" void dglapstart_(const double & dy,
			    const int & nloop);

extern "C" void dglapassign_(void (*pdf_subroutine)(const double & x,
						    const double & Q,
						    double * pdfres));

extern "C" void dglapeval_(const double & x,
			   const double & Q,
			   double * pdfres);


extern "C" void dglapevalsplit_(const double & x,
				const double & Q,
				const int    & iloop,
				const int    & nf,
				double * pdfres);
			    

#endif // __F77_PDFTAB_H__
