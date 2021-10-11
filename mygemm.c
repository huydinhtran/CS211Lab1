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

}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{

}
//Register Reuse part 1 End

//Register Reuse part 2
void dgemm2(const double *A, const double *B, double *C, const int n) 
{

}
    
//Register Reuse part 2 End

//Register Reuse part 3
void dgemm3(const double *A, const double *B, double *C, const int n) 
{
    int i, j, k;
    for (i = 0; i < n; i += 3) {
        for (j = 0; j < n; j += 3) {
            register int t = i*n+j; // COLUMN 0 ROW 0
            register int tt = t+n; // COLUMN 0 ROW 1
            register int ttt = tt+n; // COLUMN 0 ROW 2
            register double rc00 = C[t];
            register double rc01 = C[t+1];
            register double rc02 = C[t+2];
            register double rc10 = C[tt];
            register double rc11 = C[tt+1];
            register double rc12 = C[tt+2];
            register double rc20 = C[ttt];
            register double rc21 = C[ttt+1];
            register double rc22 = C[ttt+2];
            for (k = 0; k < n; k += 3) {
                register int ta = i*n+k;
                register int tta = ta+n;
                register int ttta = tta+n;
                register int tb = k*n+j;
                register int ttb = tb+n;
                register int tttb = ttb+n;
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
            C[t] = rc00;
            C[t+1] = rc01;
            C[t+2] = rc02;
            C[tt] = rc10;
            C[tt+1] = rc11;
            C[tt+2] = rc12;
            C[ttt] = rc20;
            C[ttt+1] = rc21;
            C[ttt+2] = rc22;
        }
    }
}
//Register Reuse part 3 End

//Cache Reuse part 3
void ijk(const double *A, const double *B, double *C, const int n) 
{

}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{
    
}


void jik(const double *A, const double *B, double *C, const int n) 
{

}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void kij(const double *A, const double *B, double *C, const int n) 
{
}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{
                
}


void ikj(const double *A, const double *B, double *C, const int n) 
{

}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{
 
}

void jki(const double *A, const double *B, double *C, const int n) 
{

}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{

}

void kji(const double *A, const double *B, double *C, const int n) 
{

}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{

}
//Cache Reuse part 3 End 

//Cache Reuse part 4
void optimal(const double* A, const double* B, double *C, const int n, const int b)
{

}
