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
                for (i1 = i; i1 < i+b; i++)
                    for (j1 = j; j1 < j+b; j++)
                        for (k1 = k; k1 < k+b; k++)
                            C[i1*n+j1] += A[i1*n + k1]*B[k1*n + j1];
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
                for (j1 = i; j1 < j+b; j++)
                    for (i1 = j; i1 < i+b; i++)
                        for (k1 = k; k1 < k+b; k++)
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
                for (k1 = k; k1 < k+b; k++)
                    for (i1 = i; i1 < i+b; i++){
                        double r = A[i1*n + k1];
                        for (j1 = j; j1 < j+b; j++)
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
                for (i1 = i; i1 < i+b; i++)
                    for (k1 = k; k1 < k+b; k++){
                        double r = A[i1*n + k1];
                        for (j1 = j; j1 < j+b; j++)
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
                for (j1 = j; j1 < j+b; j++)
                    for (k1 = k; k1 < k+b; k++){
                        double r = B[k1*n + j1];
                        for (i1 = i; i1 < i+b; i++)
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
                for (k1 = k; k1 < k+b; k++)
                    for (j1 = j; j1 < j+b; j++){
                        double r = B[k1*n + j1];
                        for (i1 = i; i1 < i+b; i++)
                            C[i1*n+j1] += A[i1*n + k1]*r;
                    } 
}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
    int i, j, k;
    for (i=0; i<n; i++)
        for (j=0; j<n; j++) {
            register double r=C[i*n+j];
            for (k=0; k<n; k++)
                r += A[i*n+k] * A[k*n+j];
            C[i*n+j]=r;
        }
}
