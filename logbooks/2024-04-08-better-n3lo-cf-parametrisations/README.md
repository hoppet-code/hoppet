This directory contains routines to try and extract the large-x
behaviour of the N3LO non-singlet coefficient functions, and attempts
to improve on the parametrisations in the Moch, Vogt, Vermaseren (MVV)
papers of the early 2000's.

Files of interest are:

      *) Large-x-expansions-from-MVV.nb: Contains all the large-x
      expansions of the coefficient functions (see references therein).

      *) N3LOCNS.f: Fortran implementation of the above expansions.

      *) xc*ns3*.f: The coefficient functions as provied by MVV. Small
       modifications to get the charged current piece (which is the
       problematic one) and to compute the soft part always. 

      *) print_coefficient_functions_for_fit.f: Prints the coefficient
      functions, exact, and various parametrisations to the file
      coefficient_functions_for_fit.dat.

      *) plot_exact_vs_param.gp: Gnuplot script to plot the
      coefficient functions compared to the parametrisations and
      large-x expansions.

      *) param_vs_exact.pdf: The plots.

From the above we have concluded that the parametrisations work very
well below x~0.9. The CL coefficient function in particular seems to
be working acros all x-values, whereas the C2 and C3 do not.