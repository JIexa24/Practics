#include "../include/functions.h"

int getrand(int32_t min, int32_t max)
{
  return (double)rand() / (RAND_MAX + 1.0) * (max - min) + min;
}
/*---------------------------------------------------------------------------*/
int main(int argc, char** argv)
{
  assert(!(argc < 2));

  int i, j;
  int32_t** matrixOne = NULL;
  int32_t** matrixTwo = NULL;
  int32_t** matrixRezult = NULL;
  int32_t* one = NULL;
  int32_t* two = NULL;
  int32_t* rezult = NULL;

  int size = atoi(argv[1]);
  int realSize = 0;
  int32_t min = -5;
  int32_t max = 5;
  for (i = 0; ; i++) {
    if ((size > myPow(2,i)) & (size < myPow(2, i + 1))) {
      realSize = myPow(2, i + 1);
      break;
    }
  }
  matrixOne = (int32_t**)malloc(sizeof(int32_t*) * realSize);
  for (i = 0; i < realSize; i++) {
    matrixOne[i] = (int32_t*)malloc(sizeof(int32_t) * realSize);
  }

  matrixTwo = (int32_t**)malloc(sizeof(int32_t*) * realSize);
  for (i = 0; i < realSize; i++) {
    matrixTwo[i] = (int32_t*)malloc(sizeof(int32_t) * realSize);
  }

  matrixRezult = (int32_t**)malloc(sizeof(int32_t*) * realSize);
  for (i = 0; i < realSize; i++) {
    matrixRezult[i] = (int32_t*)malloc(sizeof(int32_t) * realSize);
  }

  for (i = 0; i < realSize; i++) {
    for (j = 0; j < realSize; j++) {
      if ((i >= size) | (j >= size)) {
        matrixOne[i][j] = 0;
        matrixTwo[i][j] = 0;
      } else {
        matrixOne[i][j] = getrand(min,max);
        matrixTwo[i][j] = getrand(min,max);
      }
      matrixRezult[i][j] = 0;
    }
  }

  double time = wtime();
  simpleMatrixProizv(matrixOne, matrixTwo, matrixRezult, realSize);
  time = wtime() - time;
  printf("simpleMatrixProizv\t%.6lf\n" , time);

  for (i = 0; i < realSize; i++) {
    for (j = 0; j < realSize; j++) {
      matrixRezult[i][j] = 0;
    }
  }

  time = wtime();
  simpleMatrixProizvAsm(matrixOne, matrixTwo, matrixRezult, realSize);
  time = wtime() - time;
  printf("simpleMatrixProizvAsm\t%.6lf\n" , time);

  for (i = 0; i < realSize; i++) {
    for (j = 0; j < realSize; j++) {
      matrixRezult[i][j] = 0;
    }
  }

  time = wtime();
  simpleMatrixProizvCache(matrixOne, matrixTwo, matrixRezult, realSize);
  time = wtime() - time;

  printf("simpleMatrixProizvCache\t%.6lf\n" , time);

  for (i = 0; i < realSize; i++) {
    free(matrixRezult[i]);
  }

  free(matrixRezult);
  one = (int32_t*)malloc(sizeof(int32_t) * realSize * realSize);
  two = (int32_t*)malloc(sizeof(int32_t) * realSize * realSize);
  rezult = (int32_t*)malloc(sizeof(int32_t) * realSize * realSize);

  for (i = 0; i < realSize; i++) {
    for (j = 0; j < realSize; j++) {
      if ((i >= size) | (j >= size)) {
        one[i * size + j] = 0;
      } else {
        one[i * size + j] = getrand(min,max);
      }
    }
  }
  for (i = 0; i < realSize; i++) {
    for (j = 0; j < realSize; j++) {
      if ((i >= size) | (j >= size)) {
        two[i * size + j] = 0;
      } else {
        two[i * size + j] = getrand(min,max);
      }
    }
  }
  for(i=0;i < realSize; i++){
    free(matrixOne[i]);
    free(matrixTwo[i]);
  }
  free(matrixOne);
  free(matrixTwo);

  for (i = 0; i < realSize; i++) {
    for (j = 0; j < realSize; j++) {
      rezult[i * size + j] = 0;
    }
  }

  time = wtime();
  simpleMatrixProizvCacheOblivious(rezult, one, two, realSize, realSize);
  time = wtime() - time;
  printf("simpleMatrixProizvCacheOblivious\t%.6lf \n" , time);


  return 0;
}
