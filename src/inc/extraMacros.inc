! an extension to FTL_MACROS that provides stringification macros
#ifndef EXTRA_MACROS_INCLUDED
#define EXTRA_MACROS_INCLUDED

#if defined(__GFORTRAN__) || defined(NAGFOR) || defined(__flang__)
#define STRINGIFY(x) "x"
#else 
#define STRINGIFY(x) #x
#endif
#endif // EXTRA_MACROS_INCLUDED