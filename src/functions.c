#include "../include/functions.h"
#include <pthread.h>
extern int threadnum;
extern int threadnumst;

int myPow(int a, int b)
{
  int tmp = a;
  if (b == 0) {
    return 1;
  } else if (b == 1) {
    return a;
  } else {
    int i = 1;

    for (i = 1; i < b; i++)
      a = a * tmp;
    return a;
  }
}
/*---------------------------------------------------------------------------*/
double wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}
/*---------------------------------------------------------------------------*/
void simpleMatrixProizvAsm(int32_t** first, int32_t** second,
                           int32_t** rezult, int size){
  int i,j,k;
  for (i = 0; i < size; i++) {
    for (k = 0; k < size; k++) {
      for (j = 0; j < size; j++) {
        asm volatile (
                      ".intel_syntax noprefix\n\t"
                      "mul edx\n\t"
                      "add %0, eax\n\t"
                      : "=m" (rezult[i][j])
                      : "a" (first[i][k]), "d"(second[k][j])
                      :
                      );
      }
    }
  }
}
/*---------------------------------------------------------------------------*/
void simpleMatrixProizv(int32_t** first, int32_t** second,
                        int32_t** rezult, int size){
  int i,j,k;
  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      for (k = 0; k < size; k++) {
        rezult[i][j] += first[i][k] * second[k][j];
      }
    }
  }
}
/*---------------------------------------------------------------------------*/
void simpleMatrixProizvCache(int32_t** first, int32_t** second,
                             int32_t** rezult, int size){
  int i,j,k;
  for (i = 0; i < size; i++) {
    for (k = 0; k < size; k++) {
      for (j = 0; j < size; j++) {
        rezult[i][j] += first[i][k] * second[k][j];
      }
    }
  }
}
/*---------------------------------------------------------------------------*/
void simpleMatrixProizvCacheOblivious(int32_t *C,  int32_t *A,  int32_t *B,
                                      int size, int rowsize)
{
  if (size == 2)
  {
    const int ind11 = 0;
    const int ind12 = 1;
    const int ind21 = rowsize;
    const int ind22 = rowsize + 1;

    C[ind11] += A[ind11] * B[ind11] + A[ind12] * B[ind21];
    C[ind12] += A[ind11] * B[ind12] + A[ind12] * B[ind22];
    C[ind21] += A[ind21] * B[ind11] + A[ind22] * B[ind21];
    C[ind22] += A[ind21] * B[ind12] + A[ind22] * B[ind22];
  } else {
    const int ind11 = 0;
    const int ind12 = size / 2;
    const int ind21 = (size / 2) * rowsize;
    const int ind22 = (size / 2) * (rowsize + 1);

    // C11 += A11 * B11
    simpleMatrixProizvCacheOblivious(C + ind11, A + ind11, B + ind11, size / 2, rowsize);
    // C11 += A12 * B21
    simpleMatrixProizvCacheOblivious(C + ind11, A + ind12, B + ind21, size / 2, rowsize);

    // C12 += A11 * B12
    simpleMatrixProizvCacheOblivious(C + ind12, A + ind11, B + ind12, size / 2, rowsize);
    // C12 += A12 * B22
    simpleMatrixProizvCacheOblivious(C + ind12, A + ind12, B + ind22, size / 2, rowsize);

    // C21 += A21 * B11
    simpleMatrixProizvCacheOblivious(C + ind21, A + ind21, B + ind11, size / 2, rowsize);
    // C21 += A22 * B21
    simpleMatrixProizvCacheOblivious(C + ind21, A + ind22, B + ind21, size / 2, rowsize);

    // C22 += A21 * B12
    simpleMatrixProizvCacheOblivious(C + ind22, A + ind21, B + ind12, size / 2, rowsize);
    // C22 += A22 * B22
    simpleMatrixProizvCacheOblivious(C + ind22, A + ind22, B + ind22, size / 2, rowsize);
  }
}
