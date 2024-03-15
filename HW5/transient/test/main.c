#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 5
#define HALFBW 1

void matrixMult(float RHS[N+1], float A[N+1][2*HALFBW+1+1], float x[N+1], float b[N+1]);

int main(void){
  float A[N+1][2*HALFBW+1+1] = {{0, 0, 0, 0}, {0, 0, 1, 1}, {0, 1, 1, 1}, {0, 1, 1, 1}, {0, 1, 1, 1}, {0, 1, 1, 0}};
  printf("Printing A\n");
  for(int i=0; i <= N; i++){
    for(int j=0; j <= 2*HALFBW+1; j++){
      printf("%f\t", A[i][j]);
    }
    printf("\n");
  }
  float x[N+1] = {0, 1, 1, 1, 1, 1};
  float b[N+1] = {0, 0, 0, 0, 0, 0};
  float RHS[N+1];

  matrixMult(RHS, A, x, b);
  printf("Printing RHS\n");
  for(int i=0; i <= N; i++){
    printf("%f\n", RHS[i]);
  }

  return 0;
}


void matrixMult(float RHS[N+1], float A[N+1][2*HALFBW+1+1], float x[N+1], float b[N+1]){
  //Multiply banded matrix A with x

  for(int i=0; i <= N; i++){
    RHS[i] = 0.0;
  }
  int j = 0;
  for (int i=1; i <= N; i++){
    for (int jb=1; jb <= 2*HALFBW+1; jb++){
      j = i + jb - HALFBW - 1;
      if(jb > 0 && jb <= N){
        RHS[i] += A[i][jb] * x[j] + b[i];
      }
    }
  }

}