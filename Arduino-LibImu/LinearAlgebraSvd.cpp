/* 
 * svdcomp - SVD decomposition routine. 
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a 
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is 
 * code from Numerical Recipes adapted by Luke Tierney and David Betz.
 *
 * Input to dsvd is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   m = row dimension of a
 *   n = column dimension of a
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
*/

#include <LinearAlgebraSvd.hpp>

#define MAX_ITERATION_COUNT 30   // Maximum number of iterations

//                        Internally Defined Routines 
static void Householders_Reduction_to_Bidiagonal_Form(double* A, int nrows,
    int ncols, double* U, double* V, double* diagonal, double* superdiagonal );
static int  Givens_Reduction_to_Diagonal_Form( int nrows, int ncols,
           double* U, double* V, double* diagonal, double* superdiagonal );
static void Sort_by_Decreasing_Singular_Values(int nrows, int ncols,
                                double* singular_value, double* U, double* V);


int Singular_Value_Decomposition(double* A, int nrows, int ncols, double* U, 
                      double* singular_values, double* V, double* dummy_array)
{
    Householders_Reduction_to_Bidiagonal_Form( A, nrows, ncols, U, V,
                                                singular_values, dummy_array);

    if (Givens_Reduction_to_Diagonal_Form( nrows, ncols, U, V,
                                singular_values, dummy_array ) < 0) return -1;

    Sort_by_Decreasing_Singular_Values(nrows, ncols, singular_values, U, V);
  
    return 0;
}


////////////////////////////////////////////////////////////////////////////////
// static void Householders_Reduction_to_Bidiagonal_Form(double* A, int nrows,//
//  int ncols, double* U, double* V, double* diagonal, double* superdiagonal )//
//                                                                            //
//  Description:                                                              //
//     This routine decomposes an m x n matrix A, with m >= n, into a product //
//     of the three matrices U, B, and V', i.e. A = UBV', where U is an m x n //
//     matrix whose columns are orthogonal, B is a n x n bidiagonal matrix,   //
//     and V is an n x n orthogonal matrix.  V' denotes the transpose of V.   //
//     If m < n, then the procedure may be used for the matrix A'.  The       //
//                                                                            //
//     The matrix U is the product of Householder transformations which       //
//     annihilate the subdiagonal components of A while the matrix V is       //
//     the product of Householder transformations which annihilate the        //
//     components of A to the right of the superdiagonal.                     //
//                                                                            //
//     The Householder transformation which leaves invariant the first k-1    //
//     elements of the k-th column and annihilates the all the elements below //
//     the diagonal element is P = I - (2/u'u)uu', u is an nrows-dimensional  //
//     vector the first k-1 components of which are zero and the last         //
//     components agree with the current transformed matrix below the diagonal//
//     diagonal, the remaining k-th element is the diagonal element - s, where//
//     s = (+/-)sqrt(sum of squares of the elements below the diagonal), the  //
//     sign is chosen opposite that of the diagonal element.                  //
//                                                                            //
//  Arguments:                                                                //
//     double* A                                                              //
//        On input, the pointer to the first element of the matrix            //
//        A[nrows][ncols].  The matrix A is unchanged.                        //
//     int nrows                                                              //
//        The number of rows of the matrix A.                                 //
//     int ncols                                                              //
//        The number of columns of the matrix A.                              //
//     double* U                                                              //
//        On input, a pointer to a matrix with the same number of rows and    //
//        columns as the matrix A.  On output, the matrix with mutually       //
//        orthogonal columns which is the left-most factor in the bidiagonal  //
//        decomposition of A.                                                 //
//     double* V                                                              //
//        On input, a pointer to a square matrix with the same number of rows //
//        and columns as the columns of the matrix A, i.e. V[ncols][ncols].   //
//        On output, the orthogonal matrix whose transpose is the right-most  //
//        factor in the bidiagonal decomposition of A.                        //
//     double* diagonal                                                       //
//        On input, a pointer to an array dimensioned to same as the number   //
//        of columns of the matrix A, ncols.  On output, the diagonal of the  //
//        bidiagonal matrix.                                                  //
//     double* superdiagonal                                                  //
//        On input, a pointer to an array dimensioned to same as the number   //
//        of columns of the matrix A, ncols.  On output, the superdiagonal    //
//        of the bidiagonal matrix.                                           //
//                                                                            //
//  Return Values:                                                            //
//     The function is of type void and therefore does not return a value.    //
//     The matrices U, V, and the diagonal and superdiagonal are calculated   //
//     using the addresses passed in the argument list.                       //
//                                                                            //
//  Example:                                                                  //
//     #define M                                                              //
//     #define N                                                              //
//     double A[M][N];                                                        //
//     double U[M][N];                                                        //
//     double V[N][N];                                                        //
//     double diagonal[N];                                                    //
//     double superdiagonal[N];                                               //
//                                                                            //
//     (your code to initialize the matrix A - Note this routine is not       //
//     (accessible from outside i.e. it is declared static)                   //
//                                                                            //
//     Householders_Reduction_to_Bidiagonal_Form((double*) A, nrows, ncols,   //
//                   (double*) U, (double*) V, diagonal, superdiagonal )      //
//                                                                            //
//     free(dummy_array);                                                     //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void Householders_Reduction_to_Bidiagonal_Form(double* A, int nrows,
    int ncols, double* U, double* V, double* diagonal, double* superdiagonal )
{
    int i,j,k,ip1;
    double s, s2, si, scale;
    double dum;
    double *pu, *pui, *pv, *pvi;
    double half_norm_squared;

    memcpy(U,A, sizeof(double) * nrows * ncols);
 
    diagonal[0] = 0.0;
    s = 0.0;
    scale = 0.0;
    for ( i = 0, pui = U, ip1 = 1; i < ncols; pui += ncols, i++, ip1++ ) {
        superdiagonal[i] = scale * s;
//
//                  Perform Householder transform on columns.
//
//       Calculate the normed squared of the i-th column vector starting at 
//       row i.
//
        for (j = i, pu = pui, scale = 0.0; j < nrows; j++, pu += ncols)
            scale += fabs( *(pu + i) );
       
        if (scale > 0.0) {
            for (j = i, pu = pui, s2 = 0.0; j < nrows; j++, pu += ncols) {
                *(pu + i) /= scale;
                s2 += *(pu + i) * *(pu + i);
            }
//    
//       Chose sign of s which maximizes the norm
//  
            s = ( *(pui + i) < 0.0 ) ? sqrt(s2) : -sqrt(s2);
//
//       Calculate -2/u'u
//
        half_norm_squared = *(pui + i) * s - s2;
//
//       Transform remaining columns by the Householder transform.
//
            *(pui + i) -= s;
 
            for (j = ip1; j < ncols; j++) {
                for (k = i, si = 0.0, pu = pui; k < nrows; k++, pu += ncols)
                    si += *(pu + i) * *(pu + j);
                si /= half_norm_squared;
                for (k = i, pu = pui; k < nrows; k++, pu += ncols) {
                    *(pu + j) += si * *(pu + i);
                }
            }
        }
        for (j = i, pu = pui; j < nrows; j++, pu += ncols) *(pu + i) *= scale;
        diagonal[i] = s * scale;
//       
//                  Perform Householder transform on rows.
//
//       Calculate the normed squared of the i-th row vector starting at 
//       column i.
//
        s = 0.0;
        scale = 0.0;
        if (i >= nrows || i == (ncols - 1) ) continue;
        for (j = ip1; j < ncols; j++) scale += fabs ( *(pui + j) );

        if ( scale > 0.0 ) {
            for (j = ip1, s2 = 0.0; j < ncols; j++) {
               *(pui + j) /= scale;
               s2 += *(pui + j) * *(pui + j);
            }
            s = ( *(pui + ip1) < 0.0 ) ? sqrt(s2) : -sqrt(s2);
//
//       Calculate -2/u'u
//
            half_norm_squared = *(pui + ip1) * s - s2;
//
//       Transform the rows by the Householder transform.
//
            *(pui + ip1) -= s;
            for (k = ip1; k < ncols; k++)
                superdiagonal[k] = *(pui + k) / half_norm_squared;
            if ( i < (nrows - 1) ) {
                for (j = ip1, pu = pui + ncols; j < nrows; j++, pu += ncols) {
                    for (k = ip1, si = 0.0; k < ncols; k++) 
                        si += *(pui + k) * *(pu + k);
                    for (k = ip1; k < ncols; k++) { 
                        *(pu + k) += si * superdiagonal[k];
                    }
                }
            }
            for (k = ip1; k < ncols; k++) *(pui + k) *= scale;
        }
   }

// Update V
    pui = U + ncols * (ncols - 2);
    pvi = V + ncols * (ncols - 1);
    *(pvi + ncols - 1) = 1.0;
    s = superdiagonal[ncols - 1];
    pvi -= ncols;
    for (i = ncols - 2, ip1 = ncols - 1; i >= 0; i--, pui -= ncols,
                                                      pvi -= ncols, ip1-- ) {
        if ( s != 0.0 ) {
            pv = pvi + ncols;
            for (j = ip1; j < ncols; j++, pv += ncols)
                *(pv + i) = ( *(pui + j) / *(pui + ip1) ) / s;
            for (j = ip1; j < ncols; j++) { 
                si = 0.0;
                for (k = ip1, pv = pvi + ncols; k < ncols; k++, pv += ncols)
                   si += *(pui + k) * *(pv + j);
                for (k = ip1, pv = pvi + ncols; k < ncols; k++, pv += ncols)
                   *(pv + j) += si * *(pv + i);                  
            }
        }
        pv = pvi + ncols;
        for ( j = ip1; j < ncols; j++, pv += ncols ) {
           *(pvi + j) = 0.0;
           *(pv + i) = 0.0;
        }
        *(pvi + i) = 1.0;
        s = superdiagonal[i];
    }

// Update U
    pui = U + ncols * (ncols - 1);
    for (i = ncols - 1, ip1 = ncols; i >= 0; ip1 = i, i--, pui -= ncols ) {
        s = diagonal[i];
        for ( j = ip1; j < ncols; j++) *(pui + j) = 0.0;
        if ( s != 0.0 ) {
            for (j = ip1; j < ncols; j++) { 
                si = 0.0;
                pu = pui + ncols;
                for (k = ip1; k < nrows; k++, pu += ncols)
                    si += *(pu + i) * *(pu + j);
                si = (si / *(pui + i) ) / s;
                for (k = i, pu = pui; k < nrows; k++, pu += ncols)
                    *(pu + j) += si * *(pu + i);                  
            }
            for (j = i, pu = pui; j < nrows; j++, pu += ncols){
                *(pu + i) /= s;
            }
        }
        else 
            for (j = i, pu = pui; j < nrows; j++, pu += ncols) *(pu + i) = 0.0;
        *(pui + i) += 1.0;
   }
}


////////////////////////////////////////////////////////////////////////////////
// static int Givens_Reduction_to_Diagonal_Form( int nrows, int ncols,        //
//         double* U, double* V, double* diagonal, double* superdiagonal )    //
//                                                                            //
//  Description:                                                              //
//     This routine decomposes a bidiagonal matrix given by the arrays        //
//     diagonal and superdiagonal into a product of three matrices U1, D and  //
//     V1', the matrix U1 premultiplies U and is returned in U, the matrix    //
//     V1 premultiplies V and is returned in V.  The matrix D is a diagonal   //
//     matrix and replaces the array diagonal.                                //
//                                                                            //
//     The method used to annihilate the offdiagonal elements is a variant    //
//     of the QR transformation.  The method consists of applying Givens      //
//     rotations to the right and the left of the current matrix until        //
//     the new off-diagonal elements are chased out of the matrix.            //
//                                                                            //
//     The process is an iterative process which due to roundoff errors may   //
//     not converge within a predefined number of iterations.  (This should   //
//     be unusual.)                                                           //
//                                                                            //
//  Arguments:                                                                //
//     int nrows                                                              //
//        The number of rows of the matrix U.                                 //
//     int ncols                                                              //
//        The number of columns of the matrix U.                              //
//     double* U                                                              //
//        On input, a pointer to a matrix already initialized to a matrix     //
//        with mutually orthogonal columns.   On output, the matrix with      //
//        mutually orthogonal columns.                                        //
//     double* V                                                              //
//        On input, a pointer to a square matrix with the same number of rows //
//        and columns as the columns of the matrix U, i.e. V[ncols][ncols].   //
//        The matrix V is assumed to be initialized to an orthogonal matrix.  //
//        On output, V is an orthogonal matrix.                               //
//     double* diagonal                                                       //
//        On input, a pointer to an array of dimension ncols which initially  //
//        contains the diagonal of the bidiagonal matrix.  On output, the     //
//        it contains the diagonal of the diagonal matrix.                    //
//     double* superdiagonal                                                  //
//        On input, a pointer to an array of dimension ncols which initially  //
//        the first component is zero and the successive components form the  //
//        superdiagonal of the bidiagonal matrix.                             //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The procedure failed to terminate within                  //
//                  MAX_ITERATION_COUNT iterations.                           //
//                                                                            //
//  Example:                                                                  //
//     #define M                                                              //
//     #define N                                                              //
//     double U[M][N];                                                        //
//     double V[N][N];                                                        //
//     double diagonal[N];                                                    //
//     double superdiagonal[N];                                               //
//     int err;                                                               //
//                                                                            //
//     (your code to initialize the matrices U, V, diagonal, and )            //
//     ( superdiagonal.  - Note this routine is not accessible from outside)  //
//     ( i.e. it is declared static.)                                         //
//                                                                            //
//     err = Givens_Reduction_to_Diagonal_Form( M,N,(double*)U,(double*)V,    //
//                                                 diagonal, superdiagonal ); //
//     if ( err < 0 ) printf("Failed to converge\n");                         //
//     else { ... }                                                           //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static int Givens_Reduction_to_Diagonal_Form( int nrows, int ncols,
           double* U, double* V, double* diagonal, double* superdiagonal )
{
    double epsilon;
    double c, s;
    double f,g,h;
    double x,y,z;
    double *pu, *pv;
    int i,j,k,m;
    int rotation_test;
    int iteration_count;
    for (i = 0, x = 0.0; i < ncols; i++) {
         y = fabs(diagonal[i]) + fabs(superdiagonal[i]);
         if ( x < y ) x = y;
    }

    epsilon = x * DBL_EPSILON;
    for (k = ncols - 1; k >= 0; k--) {
        iteration_count = 0;
        while(1) {
            rotation_test = 1;
            for (m = k; m >= 0; m--) { 
                if (fabs(superdiagonal[m]) <= epsilon) {rotation_test = 0; break;}
                if (fabs(diagonal[m-1]) <= epsilon) break;
        }
            if (rotation_test) {
                c = 0.0;
                s = 1.0;
                for (i = m; i <= k; i++) {  
                   f = s * superdiagonal[i];
                   superdiagonal[i] *= c;
                   if (fabs(f) <= epsilon) break;
                   g = diagonal[i];
                   h = sqrt(f*f + g*g);
                   diagonal[i] = h;
                   c = g / h;
                   s = -f / h; 
                   for (j = 0, pu = U; j < nrows; j++, pu += ncols) { 
                        y = *(pu + m - 1);
                        z = *(pu + i);
                        *(pu + m - 1 ) = y * c + z * s;
                        *(pu + i) = -y * s + z * c;
                   }
                }
         }
         z = diagonal[k];
         if (m == k ) {
            if ( z < 0.0 ) {
               diagonal[k] = -z;
               for ( j = 0, pv = V; j < ncols; j++, pv += ncols) 
                    *(pv + k) = - *(pv + k);
            }
            break;
         }
         else {
            if ( iteration_count >= MAX_ITERATION_COUNT ) return -1;
            iteration_count++;
            x = diagonal[m];
            y = diagonal[k-1];
            g = superdiagonal[k-1];
            h = superdiagonal[k];
            f = ( (y - z) * ( y + z ) + (g - h) * (g + h) )/(2.0 * h * y);
            g = sqrt( f * f + 1.0 );
            if ( f < 0.0 ) g = -g;
            f = ( (x - z) * (x + z) + h * (y / (f + g) - h) ) / x;
// Next QR Transformtion
            c = 1.0;
            s = 1.0;
            for (i = m + 1; i <= k; i++) {
                g = superdiagonal[i];
                y = diagonal[i];
                h = s * g;
                g *= c;
                z = sqrt( f * f + h * h );
                superdiagonal[i-1] = z;
                c = f / z;
                s = h / z;
                f =  x * c + g * s;
                g = -x * s + g * c;
                h = y * s;
                y *= c;
                for (j = 0, pv = V; j < ncols; j++, pv += ncols) {
                   x = *(pv + i - 1);
                   z = *(pv + i);
                   *(pv + i - 1) = x * c + z * s;
                   *(pv + i) = -x * s + z * c;
                }
                z = sqrt( f * f + h * h );
                diagonal[i - 1] = z;
                if (z != 0.0) {
                   c = f / z;
                   s = h / z;
                } 
                f = c * g + s * y;
                x = -s * g + c * y;
                for (j = 0, pu = U; j < nrows; j++, pu += ncols) {
                   y = *(pu + i - 1);
                   z = *(pu + i);
                   *(pu + i - 1) = c * y + s * z;
                   *(pu + i) = -s * y + c * z;
                }
            }
            superdiagonal[m] = 0.0;
            superdiagonal[k] = f;
            diagonal[k] = x;
         }
      } 
   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
// static void Sort_by_Decreasing_Singular_Values(int nrows, int ncols,       //
//                            double* singular_values, double* U, double* V)  //
//                                                                            //
//  Description:                                                              //
//     This routine sorts the singular values from largest to smallest        //
//     singular value and interchanges the columns of U and the columns of V  //
//     whenever a swap is made.  I.e. if the i-th singular value is swapped   //
//     with the j-th singular value, then the i-th and j-th columns of U are  //
//     interchanged and the i-th and j-th columns of V are interchanged.      //
//                                                                            //
//  Arguments:                                                                //
//     int nrows                                                              //
//        The number of rows of the matrix U.                                 //
//     int ncols                                                              //
//        The number of columns of the matrix U.                              //
//     double* singular_values                                                //
//        On input, a pointer to the array of singular values.  On output, the//
//        sorted array of singular values.                                    //
//     double* U                                                              //
//        On input, a pointer to a matrix already initialized to a matrix     //
//        with mutually orthogonal columns.  On output, the matrix with       //
//        mutually orthogonal possibly permuted columns.                      //
//     double* V                                                              //
//        On input, a pointer to a square matrix with the same number of rows //
//        and columns as the columns of the matrix U, i.e. V[ncols][ncols].   //
//        The matrix V is assumed to be initialized to an orthogonal matrix.  //
//        On output, V is an orthogonal matrix with possibly permuted columns.//
//                                                                            //
//  Return Values:                                                            //
//        The function is of type void.                                       //
//                                                                            //
//  Example:                                                                  //
//     #define M                                                              //
//     #define N                                                              //
//     double U[M][N];                                                        //
//     double V[N][N];                                                        //
//     double diagonal[N];                                                    //
//                                                                            //
//     (your code to initialize the matrices U, V, and diagonal. )            //
//     ( - Note this routine is not accessible from outside)                  //
//     ( i.e. it is declared static.)                                         //
//                                                                            //
//     Sort_by_Decreasing_Singular_Values(nrows, ncols, singular_values,      //
//                                                 (double*) U, (double*) V); //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static void Sort_by_Decreasing_Singular_Values(int nrows, int ncols,
                                double* singular_values, double* U, double* V)
{
    int i,j,max_index;
    double temp;
    double *p1, *p2;

    for (i = 0; i < ncols - 1; i++) {
        max_index = i;
        for (j = i + 1; j < ncols; j++)
            if (singular_values[j] > singular_values[max_index] ) 
                max_index = j;
        if (max_index == i) continue;
        temp = singular_values[i];
        singular_values[i] = singular_values[max_index];
        singular_values[max_index] = temp;
        p1 = U + max_index;
        p2 = U + i;
        for (j = 0; j < nrows; j++, p1 += ncols, p2 += ncols) {
            temp = *p1;
            *p1 = *p2;
            *p2 = temp;
        } 
        p1 = V + max_index;
        p2 = V + i;
        for (j = 0; j < ncols; j++, p1 += ncols, p2 += ncols) {
            temp = *p1;
            *p1 = *p2;
            *p2 = temp;
        }
   } 
}

void Singular_Value_Decomposition_Solve(double* U, double* D, double* V,  
                double tolerance, int nrows, int ncols, double *B, double* x) 
{
    int i,j,k;
    double *pu, *pv;
    double dum;

    dum = DBL_EPSILON * D[0] * (double) ncols;
    if (tolerance < dum) tolerance = dum;

    for ( i = 0, pv = V; i < ncols; i++, pv += ncols) {
        x[i] = 0.0;
        for (j = 0; j < ncols; j++)
           if (D[j] > tolerance ) {
                for (k = 0, dum = 0.0, pu = U; k < nrows; k++, pu += ncols)
                    dum += *(pu + j) * B[k];
                x[i] += dum * *(pv + j) / D[j];
           }
    } 
}

void Singular_Value_Decomposition_Inverse(double* U, double* D, double* V,  
                        double tolerance, int nrows, int ncols, double *Astar) 
{
    int i,j,k;
    double *pu, *pv, *pa;
    double dum;

    dum = DBL_EPSILON * D[0] * (double) ncols;
    if (tolerance < dum) tolerance = dum;
    for ( i = 0, pv = V, pa = Astar; i < ncols; i++, pv += ncols) 
       for ( j = 0, pu = U; j < nrows; j++, pa++) 
            for (k = 0, *pa = 0.0; k < ncols; k++, pu++)
                if (D[k] > tolerance) *pa += *(pv + k) * *pu / D[k];
}