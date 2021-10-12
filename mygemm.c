#include "mygemm.h"


/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

//Register Reuse part 1
void dgemm0(const double* A, const double* B, double* C, const int n)
{
    int i, j, k;
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            for (k=0; k<n; k++)
                C[i*n+j] += A[i*n+k] * B[k*n+j]; 
}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i=0; i<n; i++)
        for (j=0; j<n; j++) {
            register double r = C[i*n+j] ;
            for (k=0; k<n; k++)
                r += A[i*n+k] * B[k*n+j];
            C[i*n+j] = r;
    }
}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double *A, const double *B, double *C, const int n) 
{
    /* Multiply n x n matrices a and b  */
        int i, j, k;
        for (i = 0; i < n; i+=2)
            for (j = 0; j < n; j+=2)
                for (k = 0; k < n; k+=2){
                    C[i*n + j]         = A[i*n + k] * B[k*n + j] + A[i*n + k+1] * B[(k+1)*n + j] + C[i*n + j];
                    C[(i+1)*n + j]     = A[(i+1)*n + k] * B[k*n + j] + A[(i+1)*n + k+1] * B[(k+1)*n + j] + C[(i+1)*n + j];
                    C[i*n + (j+1)]     = A[i*n + k] * B[k*n + (j+1)] + A[i*n + k+1] * B[(k+1)*n + (j+1)] + C[i*n + (j+1)];
                    C[(i+1)*n + (j+1)] = A[(i+1)*n + k] * B[k*n + (j+1)] + A[(i+1)*n + k+1] * B[(k+1)*n + (j+1)] + C[(i+1)*n + (j+1)];
                }
}
    
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i = 0; i < n; i += 3) {
        for (j = 0; j < n; j += 3) {
            register int t   = i*n+j; 
            register int tt  = t+n; 
            register int ttt = tt+n; 
            register double x00 = C[t];
            register double x01 = C[t+1];
            register double x02 = C[t+2];
            register double x10 = C[tt];
            register double x11 = C[tt+1];
            register double x12 = C[tt+2];
            register double x20 = C[ttt];
            register double x21 = C[ttt+1];
            register double x22 = C[ttt+2];
            for (k = 0; k < n; k += 3) {
                register int ta   = i*n+k;
                register int tta  = ta+n;
                register int ttta = tta+n;
                register int tb   = k*n+j;
                register int ttb  = tb+n;
                register int tttb = ttb+n;

                register double R1 = A[ta]; 
                register double R2 = A[tta]; 
                register double R3 = A[ttta];
                register double R4 = B[tb]; 
                register double R5 = B[tb+1]; 
                register double R6 = B[tb+2]; 

                x00 += R1 * R4;
                x01 += R1 * R5;
                x02 += R1 * R6;
                x10 += R2 * R4;
                x11 += R2 * R5;
                x12 += R2 * R6;
                x20 += R3 * R4;
                x21 += R3 * R5;
                x22 += R3 * R6;

                R1 = A[ta+1];
                R2 = A[tta+1];
                R3 = A[ttta+1];
                R4 = B[ttb];
                R5 = B[ttb+1];
                R6 = B[ttb+2];

                x00 += R1 * R4;
                x01 += R1 * R5;
                x02 += R1 * R6;
                x10 += R2 * R4;
                x11 += R2 * R5;
                x12 += R2 * R6;
                x20 += R3 * R4;
                x21 += R3 * R5;
                x22 += R3 * R6;

                R1 = A[ta+2];
                R2 = A[tta+2];
                R3 = A[ttta+2];
                R4 = B[tttb];
                R5 = B[tttb+1];
                R6 = B[tttb+2];

                x00 += R1 * R4;
                x01 += R1 * R5;
                x02 += R1 * R6;
                x10 += R2 * R4;
                x11 += R2 * R5;
                x12 += R2 * R6;
                x20 += R3 * R4;
                x21 += R3 * R5;
                x22 += R3 * R6;
            
            }
            C[t]     = x00;
            C[t+1]   = x01;
            C[t+2]   = x02;
            C[tt]    = x10;
            C[tt+1]  = x11;
            C[tt+2]  = x12;
            C[ttt]   = x20;
            C[ttt+1] = x21;
            C[ttt+2] = x22;
        }
    }
}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i=0; i<n; i++)  {
        for (j=0; j<n; j++) {
            double sum = C[i*n+j];
            for (k=0; k<n; k++) 
                sum += A[i*n+k] * B[k*n+j];
            C[i*n+j] = sum;
        }
    } 
}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{
    /* Multiply n x n matrices a and b  */
    int i, j, k, i1, j1, k1;
    for (i = 0; i < n; i+=b)
        for (j = 0; j < n; j+=b)
            for (k = 0; k < n; k+=b)
             /* B x B mini matrix multiplications */
                for (i1 = i; i1 < i+b; i1++)
                    for (j1 = j; j1 < j+b; j1++)
                        for (k1 = k; k1 < k+b; k1++)
                            C[i1*n+j1] += A[i1*n+k1] * B[k1*n+j1];
}


void jik(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (j=0; j<n; j++) {
        for (i=0; i<n; i++) {
            double sum = 0.0;
            for (k=0; k<n; k++)
                sum += A[i*n+k] * B[k*n+j];
            C[i*n+j] = sum;
        }
    }
}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{
    /* Multiply n x n matrices a and b  */
    int i, j, k, i1, j1, k1;
    for (j = 0; j < n; j+=b)
        for (i = 0; i < n; i+=b)
            for (k = 0; k < n; k+=b)
             /* B x B mini matrix multiplications */
                for (j1 = i; j1 < j+b; j1++)
                    for (i1 = j; i1 < i+b; i1++)
                        for (k1 = k; k1 < k+b; k1++)
                            C[i1*n+j1] += C[i1*n + k1]*C[k1*n + j1];
}

void kij(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (k=0; k<n; k++) {
        for (i=0; i<n; i++) {
            double r = A[i*n+k];
            for (j=0; j<n; j++)
                C[i*n+j] += r * B[k*n+j];   
        }
    }
}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{
    /* Multiply n x n matrices a and b  */
    int i, j, k, i1, j1, k1;
    for (k = 0; k < n; k+=b)
        for (i = 0; i < n; i+=b)
            for (j = 0; j < n; j+=b)
             /* B x B mini matrix multiplications */
                for (k1 = k; k1 < k+b; k1++)
                    for (i1 = i; i1 < i+b; i1++){
                        double r = A[i1*n + k1];
                        for (j1 = j; j1 < j+b; j1++)
                            C[i1*n+j1] += r*B[k1*n + j1];
                            }                  
}


void ikj(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i=0; i<n; i++) {
        for (k=0; k<n; k++) {
            double r = A[i*n+k];
            for (j=0; j<n; j++)
                C[i*n+j] += r * B[k*n+j];
        }
    }
}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{
    /* Multiply n x n matrices a and b  */
    int i, j, k, i1, j1, k1;
    for (i = 0; i < n; i+=b)
        for (k = 0; k < n; k+=b)
            for (j = 0; j < n; j+=b)
             /* B x B mini matrix multiplications */
                for (i1 = i; i1 < i+b; i1++)
                    for (k1 = k; k1 < k+b; k1++){
                        double r = A[i1*n + k1];
                        for (j1 = j; j1 < j+b; j1++)
                            C[i1*n+j1] += r*B[k1*n + j1];
                            }      
}

void jki(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (j=0; j<n; j++) {
        for (k=0; k<n; k++) {
            double r = B[k*n+j];
            for (i=0; i<n; i++)
                C[i*n+j] += A[i*n+k] * r;
        }
    }
}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{
    /* Multiply n x n matrices a and b  */
    int i, j, k, i1, j1, k1;
    for (j = 0; j < n; j+=b)
        for (k = 0; k < n; k+=b)
            for (i = 0; i < n; i+=b)
             /* B x B mini matrix multiplications */
                for (j1 = j; j1 < j+b; j1++)
                    for (k1 = k; k1 < k+b; k1++){
                        double r = B[k1*n + j1];
                        for (i1 = i; i1 < i+b; i1++)
                            C[i1*n+j1] += A[i1*n + k1]*r;
                            } 
}

void kji(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (k=0; k<n; k++) {
        for (j=0; j<n; j++) {
            double r = B[k*n+j];
            for (i=0; i<n; i++)
                C[i*n+j] += A[i*n+k] * r;
        }
    }
}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{
    /* Multiply n x n matrices a and b  */
    int i, j, k, i1, j1, k1;
    for (k = 0; k < n; k+=b)
        for (j = 0; j < n; j+=b)
            for (k = 0; k < n; k+=b)
             /* B x B mini matrix multiplications */
                for (k1 = k; k1 < k+b; k1++)
                    for (j1 = j; j1 < j+b; j1++){
                        double r = B[k1*n + j1];
                        for (i1 = i; i1 < i+b; i1++)
                            C[i1*n+j1] += A[i1*n + k1]*r;
                            } 
}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
    /* Multiply n x n matrices a and b  */
    int i, j, k, i1, j1, k1;
    for (i = 0; i < n; i+=b)
        for (j = 0; j < n; j+=b)
            for (k = 0; k < n; k+=b)
             /* B x B mini matrix multiplications */
                for (i1 = i; i1 < i+b; i1++)
                    for (j1 = j; j1 < j+b; j1++)
                        for (k1 = k; k1 < k+b; k1++){
                            C[i*n + j]         = A[i*n + k] * B[k*n + j] + A[i*n + k+1] * B[(k+1)*n + j] + C[i*n + j];                                       
                            C[(i+1)*n + j]     = A[(i+1)*n + k] * B[k*n + j] + A[(i+1)*n + k+1] * B[(k+1)*n + j] + C[(i+1)*n + j];                    
                            C[i*n + (j+1)]     = A[i*n + k] * B[k*n + (j+1)] + A[i*n + k+1] * B[(k+1)*n + (j+1)] + C[i*n + (j+1)];                    
                            C[(i+1)*n + (j+1)] = A[(i+1)*n + k] * B[k*n + (j+1)] + A[(i+1)*n + k+1] * B[(k+1)*n + (j+1)] + C[(i+1)*n + (j+1)];
                        }
}
