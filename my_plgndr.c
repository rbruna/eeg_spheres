/*
 *  The original plgndr function was implemented using some code from
 *  Numerical Recipes in C. However, it was not allowed to release that code
 *  as open source. The new implementation below is using some code from
 *  the GNU Scientific Library (http://www.gnu.org/software/gsl).
 *
 *  Copyright (C) 2002-2006 Robert Oostenveld
 *  Copyright (C) 2006, Thomas Hartmann
 *  Copyright (C) 2016, Ricardo Bruna
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/*  Based on FieldTrip 20160222 functions:
 *  * plgndr
 */


#include <string.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

double legendre_Pmm ( double m, double x )
{
    if ( m == 0 ) {
        return 1.0;
    } else {
        double p_mm = 1.0;
        double root_factor = sqrt ( 1.0 - x ) * sqrt ( 1.0 + x );
        double fact_coeff = 1.0;
        int i;
        for ( i = 1; i <= m; i ++ ) {
            p_mm *= -fact_coeff * root_factor;
            fact_coeff += 2.0;
        }
        return p_mm;
    }
}

double *plgndr ( int l, int m, double x ) {
    
    /* Checks the inputs. */
    if ( m < 0 ) mexErrMsgTxt ( "m must be a positive integer." );
    if ( m > l ) mexErrMsgTxt ( "l must be a positive integer equal or greater than m." );
    if ( fabs ( x ) > 1.0 ) mexErrMsgTxt ( "abs(x) cannot be greater than 1." );
    
    /* Calculates the number of iterations needed. */
    int its = l - m + 1;
    
    /* Creates an array to store the output. */
    double *p = malloc ( its * sizeof ( double ) );
    
    int ell;
    int index;
    
    /* The Legendre function is constructed iteratively. */
    for ( index = 0; index < its; index ++ ) {
        ell = index + m;
        
        /* For 1 uses the analitical construction. */
        if ( index == 0 ) {
            p [ index ] = legendre_Pmm ( m, x );
            
        /* For 2 uses a simplified version of the iterative construction. */
        } else if ( index == 1 ) {
            p [ index ] = x * ( 2 * ell - 1 ) * p [ index - 1 ];
            
        /* For values greater than 2 uses the iterative construction. */
        } else {
            p [ index ] = ( x * ( 2 * ell - 1 ) * p [ index - 1 ] - ( ell + m - 1 ) * p [ index - 2 ] ) / ( ell - m );
        }
    }
    
    /* Returns the whole array. */
    return p;
}

void mexFunction ( int nlhs, mxArray * plhs [], int nrhs, const mxArray * prhs [] ) {
    int l, m;
    double x;
    double *pd;
    
    /* Checks the inputs. */
    if ( nrhs < 3 || nrhs > 4 ) mexErrMsgTxt ( "Invalid number of arguments for PLGNDR." );
    
    /* Gets the input variables */
    l = mxGetScalar ( prhs [0] );
    m = mxGetScalar ( prhs [1] );
    x = mxGetScalar ( prhs [2] );
    
    /* Calculates the number of iterations needed. */
    int its = l - m + 1;

    double *output;
    output = plgndr ( l, m, x );
    
    /* If no fourth input returns only the last value. */
    if ( nrhs == 3 ) {
        
        /* Creates the matrix array and gets a pointer to it. */
        plhs [0] = mxCreateDoubleMatrix ( 1, 1, mxREAL );
        pd = mxGetData ( plhs [0] );
        
        /* Copies the result to the matrix array. */
        pd [0] = output [ its - 1 ];
        
    /* Otherwise returns the full array. */
    } else {
        
        /* Creates the matrix array and gets a pointer to it. */
        plhs [0] = mxCreateDoubleMatrix ( 1, its, mxREAL );
        pd = mxGetData ( plhs [0] );
        
        /* Copies the result to the matrix array. */
        memcpy ( pd, output, its * sizeof ( double ) );
    }
    
    /* Frees the memory used by the original array. */
    free ( output );
    
    return;
}
