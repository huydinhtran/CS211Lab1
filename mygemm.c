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
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            for (k=0; k<n; k++)
                C[i*n+j] += A[i*n+k] * B[k*n+j]; 
}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{
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
                    C[i*n + j]         = A[i*n + k] * B[k*n + j] + A[i*n + k+1] * B[(k+1)*n + j] + C[i*n + j]
                    C[(i+1)*n + j]     = A[(i+1)*n + k] * B[k*n + j] + A[(i+1)*n + k+1] * B[(k+1)*n + j] + C[(i+1)*n + j]
                    C[i*n + (j+1)]     = A[i*n + k] * B[k*n + (j+1)] + A[i*n + k+1] * B[(k+1)*n + (j+1)] + C[i*n + (j+1)]
                    C[(i+1)*n + (j+1)] = A[(i+1)*n + k] * B[k*n + (j+1)] + A[(i+1)*n + k+1] * B[(k+1)*n + (j+1)] + C[(i+1)*n + (j+1)]
                }
}
    
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n) 
{
    for (int i = 0; i < N; i += 3) {
        for (int j = 0; j < N; j += 3) {
             register int t = i*N+j; // COLUMN 0 ROW 0
              register int tt = t+N; // COLUMN 0 ROW 1
              register int ttt = tt+N; // COLUMN 0 ROW 2
              register double rc00 = C4[t];
              register double rc01 = C4[t+1];
              register double rc02 = C4[t+2];
              register double rc10 = C4[tt];
              register double rc11 = C4[tt+1];
              register double rc12 = C4[tt+2];
              register double rc20 = C4[ttt];
              register double rc21 = C4[ttt+1];
              register double rc22 = C4[ttt+2];
              for (int k = 0; k < N; k += 3) {
                  register int ta = i*N+k;
                  register int tta = ta+N;
                  register int ttta = tta+N;
                  register int tb = k*N+j;
                  register int ttb = tb+N;
                  register int tttb = ttb+N;
                  /*
                  * C00 = A00*B00 + A01*B10 + A02*B20
                  * C01 = A00*B01 + A01*B11 + A02*B21
                  * C02 = A00*B02 + A01*B12 + A02*B22
                  *
                  * C10 = A10*B00 + A11*B10 + A12*B20;
                  * C11 = A10*B01 + A11*B11 + A12*B21;
                  * C12 = A10*B02 + A11*B12 + A12*B22;
                  *
                  * C20 = A20*B00 + A21*B10 + A22*B20;
                  * C21 = A20*B01 + A21*B11 + A22*B21;
                  * C22 = A20*B02 + A21*B12 + A22*B22;
                  */
                  register double R1 = A[ta]; // ra00
                  register double R2 = A[tta]; // r10
                  register double R3 = A[ttta]; // 20
                  register double R4 = B[tb]; // rb00
                  register double R5 = B[tb+1]; // rb01
                  register double R6 = B[tb+2]; // rb02
                  rc00 += R1 * R4;
                  rc01 += R1 * R5;
                  rc02 += R1 * R6;
                  rc10 += R2 * R4;
                  rc11 += R2 * R5;
                  rc12 += R2 * R6;
                  rc20 += R3 * R4;
                  rc21 += R3 * R5;
                  rc22 += R3 * R6;
                  R1 = A[ta+1];
                  R2 = A[tta+1];
                  R3 = A[ttta+1];
                  R4 = B[ttb];
                  R5 = B[ttb+1];
                  R6 = B[ttb+2];
                  rc00 += R1 * R4;
                  rc01 += R1 * R5;
                  rc02 += R1 * R6;
                  rc10 += R2 * R4;
                  rc11 += R2 * R5;
                  rc12 += R2 * R6;
                  rc20 += R3 * R4;
                  rc21 += R3 * R5;
                  rc22 += R3 * R6;
                  R1 = A[ta+2];
                  R2 = A[tta+2];
                  R3 = A[ttta+2];
                  R4 = B[tttb];
                  R5 = B[tttb+1];
                  R6 = B[tttb+2];
                  rc00 += R1 * R4;
                  rc01 += R1 * R5;
                  rc02 += R1 * R6;
                  rc10 += R2 * R4;
                  rc11 += R2 * R5;
                  rc12 += R2 * R6;
                  rc20 += R3 * R4;
                  rc21 += R3 * R5;
                  rc22 += R3 * R6;
                  }
                  C4[t] = rc00;
                  C4[t+1] = rc01;
                  C4[t+2] = rc02;
                  C4[tt] = rc10;
                  C4[tt+1] = rc11;
                  C4[tt+2] = rc12;
                  C4[ttt] = rc20;
                  C4[ttt+1] = rc21;
                  C4[ttt+2] = rc22;
                  }
                  }
}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n) 
{
    for (i=0; i<n; i++)  {
        for (j=0; j<n; j++) {
            sum = 0.0;
            for (k=0; k<n; k++) 
                sum += A[i][k] * B[k][j];
            C[i][j] = sum;
        }
    } 
}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{
    /* Multiply n x n matrices a and b  */
    int i, j, k;
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
    for (j=0; j<n; j++) {
        for (i=0; i<n; i++) {
            sum = 0.0;
            for (k=0; k<n; k++)
                sum += A[i][k] * B[k][j];
            C[i][j] = sum
        }
    }
}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{
    /* Multiply n x n matrices a and b  */
    int i, j, k;
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
    for (k=0; k<n; k++) {
        for (i=0; i<n; i++) {
            r = A[i][k];
            for (j=0; j<n; j++)
                C[i][j] += r * B[k][j];   
        }
    }
}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{
    /* Multiply n x n matrices a and b  */
    int i, j, k;
    for (k = 0; k < n; k+=b)
        for (i = 0; i < n; i+=b)
            for (j = 0; j < n; j+=b)
             /* B x B mini matrix multiplications */
                for (k1 = k; k1 < k+b; k++)
                    for (i1 = i; i1 < i+b; i++){
                        r = A[i1*n + k1];
                        for (j1 = j; j1 < j+b; j++)
                            C[i1*n+j1] += r*B[k1*n + j1];
                            }                  
}


void ikj(const double *A, const double *B, double *C, const int n) 
{
    for (i=0; i<n; i++) {
        for (k=0; k<n; k++) {
            r = A[i][k];
            for (j=0; j<n; j++)
                C[i][j] += r * B[k][j];
        }
    }
}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{
    /* Multiply n x n matrices a and b  */
    int i, j, k;
    for (i = 0; i < n; i+=b)
        for (k = 0; k < n; k+=b)
            for (j = 0; j < n; j+=b)
             /* B x B mini matrix multiplications */
                for (i1 = i; i1 < i+b; i++)
                    for (k1 = k; k1 < k+b; k++){
                        r = A[i1*n + k1];
                        for (j1 = j; j1 < j+b; j++)
                            C[i1*n+j1] += r*B[k1*n + j1];
                            }      
}

void jki(const double *A, const double *B, double *C, const int n) 
{
    for (j=0; j<n; j++) {
        for (k=0; k<n; k++) {
            r = B[k][j];
            for (i=0; i<n; i++)
                C[i][j] += A[i][k] * r;
        }
    }
}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{
    /* Multiply n x n matrices a and b  */
    int i, j, k;
    for (j = 0; j < n; j+=b)
        for (k = 0; k < n; k+=b)
            for (i = 0; i < n; i+=b)
             /* B x B mini matrix multiplications */
                for (j1 = j; j1 < j+b; j++)
                    for (k1 = k; k1 < k+b; k++){
                        r = B[k1*n + j1];
                        for (i1 = i; i1 < i+b; i++)
                            C[i1*n+j1] += A[i1*n + k1]*r;
                            } 
}

void kji(const double *A, const double *B, double *C, const int n) 
{
    for (k=0; k<n; k++) {
        for (j=0; j<n; j++) {
            r = B[k][j];
            for (i=0; i<n; i++)
                C[i][j] += A[i][k] * r;
        }
    }
}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{
    /* Multiply n x n matrices a and b  */
    int i, j, k;
    for (k = 0; k < n; k+=b)
        for (j = 0; j < n; j+=b)
            for (k = 0; k < n; k+=b)
             /* B x B mini matrix multiplications */
                for (k1 = k; k1 < k+b; k++)
                    for (j1 = j; j1 < j+b; j++){
                        r = B[k1*n + j1];
                        for (i1 = i; i1 < i+b; i++)
                            C[i1*n+j1] += A[i1*n + k1]*r;
                            } 
}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
    /* ijk – simple triple loop algorithm with simple single register reuse*/
    for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
        register double r=c[i*n+j];
        for (k=0; k<n; k++)
        r += a[i*n+k] * b[k*n+j];
        c[i*n+j]=r;
        }
        /* ijk – blocked version algorithm*/
        for (i = 0; i < n; i+=B)
        for (j = 0; j < n; j+=B)
        for (k = 0; k < n; k+=B)
        /* B x B mini matrix multiplications */
        for (i1 = i; i1 < i+B; i1++)
        for (j1 = j; j1 < j+B; j1++) {
            register double r=c[i1*n+j1];
            for (k1 = k; k1 < k+B; k1++)
            r += a[i1*n + k1]*b[k1*n + j1];
            c[i1*n+j1]=r;
            }
}
