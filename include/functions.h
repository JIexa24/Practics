#ifndef FUNC_H
#define FUNC_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <inttypes.h>
#include <pthread.h>
/*---------------------------------------------------------------------------*/
typedef struct threaddata {
  int32_t *C;
  int32_t *A;
  int32_t *B;
  int size;
  int rowsize;
} dat;

double wtime();
int myPow(int a, int b);

void simpleMatrixProizv(int32_t** first, int32_t** second,
                        int32_t** reszult, int size);
void simpleMatrixProizvAsm(int32_t** first, int32_t** second,
                           int32_t** reszult, int size);
void simpleMatrixProizvCache(int32_t** first, int32_t** second,
                             int32_t** rezult, int size);
void simpleMatrixProizvCacheOblivious(int32_t *C,  int32_t *A,  int32_t *B,
                                      int size, int rowsize);
void *simpleMatrixProizvCacheObliviousp(void* ptr);

#endif
