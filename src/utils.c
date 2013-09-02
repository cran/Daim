/*
 *  utils.c
 *  Daim
 *
 *  Created by Sergej Potapov on 07.06.10.
 *
 */

#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>


#include "utils.h"



/*	Compute the maximum of double vector
 *\param X vector 
 *\param n - length of the vector X.
*/

double dmax(double *X, int n)
{
	int i=0;
	double max=0.0;
	for (i=0; i < n; i++){
		if (X[i] > max) {
			max=X[i];
		}
	}
	return(max);
}




/*	Compute the mean of double vector
 *\param X vector 
 *\param n - length of the vector X.
 */


double d_mean(double *X, int n)
{
	int i=0;
	double max=0.0;
	for (i=0; i < n; i++){
		max +=X[i];
	}
	max /= (double) n;
	return(max);
}






static int rcmp_TW(double x, double y, Rboolean nalast)
{
    int nax = ISNAN(x), nay = ISNAN(y);
    if (nax && nay)	return 0;
    if (nax)		return nalast?1:-1;
    if (nay)		return nalast?-1:1;
    if (x < y)		return -1;
    if (x > y)		return 1;
    return 0;
}



/*	Sort the double vector x
 *\param x - vector
 *\param indx - vector of x indexes
 *\param n - length of the vector x and indx.
 */


void rsort_index(double *x, int *indx, int n)
{
    double v;
    int i, j, h, iv;
	
    for (h = 1; h <= n / 9; h = 3 * h + 1);
    for (; h > 0; h /= 3)
		for (i = h; i < n; i++) {
			v = x[i]; iv = indx[i];
			j = i;
			while (j >= h && rcmp_TW(x[j - h], v, TRUE) > 0)
			{ x[j] = x[j - h]; indx[j] = indx[j-h]; j -= h; }
			x[j] = v; indx[j] = iv;
		}
}



/*	Sort the double vector x and order vector y in sorted order of x 
 *\param x - double vector
 *\param y - double vector 
 *\param n - length of the vector x and y
 */



void rsort_with_x(double *x, double *y, int n)
{
    double v, iv;
    int i, j, h;
	
    for (h = 1; h <= n / 9; h = 3 * h + 1);
    for (; h > 0; h /= 3)
		for (i = h; i < n; i++) {
			v = x[i]; iv = y[i];
			j = i;
			while (j >= h && rcmp_TW(x[j - h], v, TRUE) > 0)
			{ x[j] = x[j - h]; y[j] = y[j-h]; j -= h; }
			x[j] = v; y[j] = iv;
		}
}



/*	Sort the double vector x and order vector y and indx in sorted order of x 
 *\param x - double vector
 *\param y - double vector 
 *\param n - length of the vector x
 */

void rsort_xyz(double *x, double *y, double *indx, int n)
{
    double v, iv, vi;
    int i, j, h;
	
    for (h = 1; h <= n / 9; h = 3 * h + 1);
    for (; h > 0; h /= 3)
		for (i = h; i < n; i++) {
			v = x[i]; iv = indx[i]; vi = y[i];
			j = i;
			while (j >= h && rcmp_TW(x[j - h], v, TRUE) > 0)
			{ x[j] = x[j - h]; indx[j] = indx[j-h]; y[j] = y[j-h]; j -= h; }
			x[j] = v; indx[j] = iv; y[j] = vi;
		}
}




/*	Sort the double vector x and order vector y, z and indx in sorted order of x 
 *\param x - double vector
 *\param y - double vector 
 *\param n - length of the vector x
 */

void rsort_xyzv(double *x, double *y, double *z, double *indx, int n)
{
    double v, iv, vi, vii;
    int i, j, h;
	
    for (h = 1; h <= n / 9; h = 3 * h + 1);
    for (; h > 0; h /= 3)
		for (i = h; i < n; i++) {
			v = x[i]; iv = indx[i]; vi = y[i]; vii = z[i];
			j = i;
			while (j >= h && rcmp_TW(x[j - h], v, TRUE) > 0)
			{ x[j] = x[j - h]; indx[j] = indx[j-h]; y[j] = y[j-h]; z[j] = z[j-h]; j -= h; }
			x[j] = v; indx[j] = iv; y[j] = vi; z[j] = vii;
		}
}






/*	calculate a cumulative sums
 *\param x - double vector
 *\param size - length of the vector x
 */

void cum_sum(double *x, int size)
{
    LDOUBLE sum = 0.;
    int i;
    for (i = 0 ; i < size ; i++) {
		sum += x[i];
		x[i] = sum;
    }
} 




/*	reverse vector elements
 *\param x - double vector
 *\param n_x - length of the vector x
 */

void my_rev_d(double *x, int *n_x)
{
	double swap = 0.0;
	int i, j=*n_x-1;
	for (i = 0; i < j; i++, j--) {
		swap = x[i];
		x[i] = x[j];
		x[j] = swap;
	}
}


