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
    int tsize = size /2 ;
    const int ind11 = 0;
    const int ind12 = tsize;
    const int ind21 = (tsize) * rowsize;
    const int ind22 = (tsize) * (rowsize + 1);
    pthread_t tid[7];
    int* args[8][5] = {
{C + ind11, A + ind11, B + ind11, &tsize, &rowsize},//
{C + ind11, A + ind12, B + ind21, &tsize, &rowsize},//
{C + ind12, A + ind11, B + ind12, &tsize, &rowsize},//
{C + ind12, A + ind12, B + ind22, &tsize, &rowsize},//
{C + ind21, A + ind21, B + ind11, &tsize, &rowsize},//
{C + ind21, A + ind22, B + ind21, &tsize, &rowsize},
{C + ind22, A + ind21, B + ind12, &tsize, &rowsize},
{C + ind22, A + ind22, B + ind22,  &tsize, &rowsize}
}; 
    // C11 += A11 * B11
    if (threadnum < threadnumst){
    pthread_create(&tid[0],NULL,simpleMatrixProizvCacheObliviousp,&args[0]);
    threadnum++;
    } else {
    simpleMatrixProizvCacheObliviousp(C + ind11, A + ind11, B + ind11, &tsize, &rowsize);
    }
    // C11 += A12 * B21
    if (threadnum < threadnumst){
    pthread_create(&tid[1],NULL,simpleMatrixProizvCacheObliviousp,&args[1]);
    threadnum++;
    } else {
    simpleMatrixProizvCacheObliviousp(C + ind11, A + ind12, B + ind21, &tsize, &rowsize);
    }
 
    // C12 += A11 * B12
    if (threadnum < threadnumst){
    pthread_create(&tid[2],NULL,simpleMatrixProizvCacheObliviousp,&args[2]);
    threadnum++;
    } else {
    simpleMatrixProizvCacheObliviousp(C + ind12, A + ind11, B + ind12, &tsize, &rowsize);
    }
    // C12 += A12 * B22
    if (threadnum < threadnumst){
    pthread_create(&tid[3],NULL,simpleMatrixProizvCacheObliviousp,&args[3]);
    threadnum++;
    } else {
    simpleMatrixProizvCacheObliviousp(C + ind12, A + ind12, B + ind22, &tsize, &rowsize);
    }
 
    // C21 += A21 * B11
    if (threadnum < threadnumst){
    pthread_create(&tid[4],NULL,simpleMatrixProizvCacheObliviousp,&args[4]);
    threadnum++;
    } else {
    simpleMatrixProizvCacheObliviousp(C + ind21, A + ind21, B + ind11, &tsize, &rowsize);
    }

    // C21 += A22 * B21
    if (threadnum < threadnumst){
    pthread_create(&tid[5],NULL,simpleMatrixProizvCacheObliviousp,&args[5]);
    threadnum++;
    } else {
    simpleMatrixProizvCacheObliviousp(C + ind21, A + ind22, B + ind21,  &tsize, &rowsize);
    }
 
    // C22 += A21 * B12
    if (threadnum < threadnumst){
    pthread_create(&tid[6],NULL,simpleMatrixProizvCacheObliviousp,&args[6]);
    threadnum++;
    } else {
    simpleMatrixProizvCacheObliviousp(C + ind22, A + ind21, B + ind12,  &tsize, &rowsize);
    }

    // C22 += A22 * B22
    if (threadnum < threadnumst){
    pthread_create(&tid[6],NULL,simpleMatrixProizvCacheObliviousp,&args[7]);
    threadnum++;
    } else {
    simpleMatrixProizvCacheObliviousp(C + ind22, A + ind22, B + ind22,  &tsize, &rowsize);
    }
  }
}

void simpleMatrixProizvCacheObliviousp(int32_t *C,  int32_t *A,  int32_t *B,
                                      int* size, int* rowsize)
{
  if (*size == 2)
  {
    const int ind11 = 0;
    const int ind12 = 1;
    const int ind21 = *rowsize;
    const int ind22 = *rowsize + 1;

    C[ind11] += A[ind11] * B[ind11] + A[ind12] * B[ind21];
    C[ind12] += A[ind11] * B[ind12] + A[ind12] * B[ind22];
    C[ind21] += A[ind21] * B[ind11] + A[ind22] * B[ind21];
    C[ind22] += A[ind21] * B[ind12] + A[ind22] * B[ind22];
  } else {
    int tsize = *size /2 ;
    const int ind11 = 0;
    const int ind12 = tsize;
    const int ind21 = (tsize) * *rowsize;
    const int ind22 = (tsize) * (*rowsize + 1);
 
    // C11 += A11 * B11
    simpleMatrixProizvCacheObliviousp(C + ind11, A + ind11, B + ind11, &tsize, &rowsize);
    // C11 += A12 * B21
    simpleMatrixProizvCacheObliviousp(C + ind11, A + ind12, B + ind21, &tsize, &rowsize);
 
    // C12 += A11 * B12
    simpleMatrixProizvCacheObliviousp(C + ind12, A + ind11, B + ind12, &tsize, &rowsize);
    // C12 += A12 * B22
    simpleMatrixProizvCacheObliviousp(C + ind12, A + ind12, B + ind22, &tsize, &rowsize);
 
    // C21 += A21 * B11
    simpleMatrixProizvCacheObliviousp(C + ind21, A + ind21, B + ind11, &tsize, &rowsize);
    // C21 += A22 * B21
    simpleMatrixProizvCacheObliviousp(C + ind21, A + ind22, B + ind21, &tsize, &rowsize);
 
    // C22 += A21 * B12
    simpleMatrixProizvCacheObliviousp(C + ind22, A + ind21, B + ind12, &tsize, &rowsize);
    // C22 += A22 * B22
    simpleMatrixProizvCacheObliviousp(C + ind22, A + ind22, B + ind22, &tsize, &rowsize);
  }
}

