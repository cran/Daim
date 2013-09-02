/*
 *  utils.h
 *  Daim
 *
 *  Created by Sergej Potapov on 07.06.10.
 *
 */



/* taken from pnmath Package*/
/*							*/

/* Private header file for use during compilation of Mathlib */
#ifndef MATHLIB_PRIVATE_H
#define MATHLIB_PRIVATE_H

#ifndef MATHLIB_STANDALONE
/* Mathlib in R */
# ifdef HAVE_CONFIG_H
#  include <config.h>
# endif
# if defined(HAVE_GLIBC2) && !defined(_BSD_SOURCE)
/* ensure isnan is declared */
#  define _BSD_SOURCE 1
# endif
# if defined(HAVE_GLIBC2) && !defined(_ISOC99_SOURCE)
/* ensure expm1 and log1p are declared */
#  define _ISOC99_SOURCE 1
# endif
#endif

/* need to add LDOUBLE to Rconfig.h and include Rconfig.h in pnmath.h */
#ifdef DODO
#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#else
#  define LDOUBLE double
#endif
#else
#  define LDOUBLE long double
#endif
#endif 
/* standalone */




double dmax(double *X, int n);
double d_mean(double *X, int n);
void rsort_index(double *x, int *indx, int n);
void rsort_with_x(double *x, double *indx, int n);
void rsort_xyz(double *x, double *y, double *indx, int n);
void rsort_xyzv(double *x, double *y, double *z, double *indx, int n);
void cum_sum(double *x, int size);
void my_rev_d(double *x, int *n_x);




