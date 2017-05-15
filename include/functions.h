#ifndef FUNC_H
#define FUNC_H


#include <sys/time.h>
#include <inttypes.h>
/*---------------------------------------------------------------------------*/
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

#endif
