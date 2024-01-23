/*
Author: Aaryan Agrawal
Course: ENGS105, W24
Instructor: Professor Keith Paulsen

Problem 12c
Electric current flow past a circular hole in a conducting sheet
Solve PDE using FD molecule over polar grid
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solve1.h"

#define PI 3.14159

int nodemap(int i, int j, int N);

int main (void) {

int NRange[] = {5, 10};

for(int ii=0; ii < (int)sizeof(NRange)/sizeof(NRange[0]); ii++){
  int N;
  N = NRange[ii];
  printf("N = %d\n", N);
  //define N
  int n = N*N;
  //define parameters from question
  float a_param = 0.1;
  int R_param = 1;
  int k_param = 3;

  //define spacing
  float delR = (R_param-a_param)/(N-1);
  float delTheta = (PI/2)/(N-1);
  float beta = (delR/delTheta) * (delR/delTheta);

  //store theta and radius values
  float theta[N];
  float radius[N];
  for (int i=0; i<N; i++){
    theta[i] = i*delTheta;
    radius[i] = a_param + i*delR;
  }

  ///////////////////////////////////////////////////////////////////////////
  //// make b vector
  float B[n];
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      B[nodemap(i, j, N)] = 0;
    }
  }
  for(int i=N-1; i<n; i+=N){
    B[i] = cos(k_param*theta[(i+1)/N - 1]) * (delR + 1);
  }

  /*write b vector to file
  FILE *file;
  file = fopen("b.txt", "w");
  if(file == NULL){
    printf("Error opening file for B\n");
    return 1;
  }
  for(int i=0; i<n; i++){
    fprintf(file, "%f\n", B[i]);
  }
  fclose(file);*/
  
  /////////////////////////////////////////////////////////////////////////
  //// make A matrix
  float A[n][n];
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      A[i][j] = 0;
    }
  }
  float a, b, c, d;
  float rVal;
  int center, rPlus, rMinus, thetaPlus, thetaMinus;
  //iterate over all nodes, check for boundary conditions
  for (int j=0; j<N; j++){
    for(int i=0; i<N; i++){
      rVal = radius[i];

      a = - 2 - (2*beta)/(rVal*rVal);
      b = 1 + delR/(2*rVal);
      c = 1 - delR/(2*rVal);
      d = beta / (rVal*rVal);

      center = nodemap(i, j, N);
      //printf("center: %d\n", center);
      // For Banded Storage Mode 2 for A, A is n x 2N+1. Center of molecule is at index=N.
      rPlus = nodemap(i+1, j, N);
      rMinus = nodemap(i-1, j, N);
      thetaPlus = nodemap(i, j+1, N);
      thetaMinus = nodemap(i, j-1, N);

      if(center == 0){    //Type II corner
        A[center][center] = a;
        A[center][rPlus] = 2*b;
        A[center][thetaPlus] = 2*d;
        //printf("Corner\n");
      }
      else if(center >= n-N){   //Boundary I
        A[center][center] = 1;
        //printf("Boundary I\n");
      }
      else if(center > 0 && center < N-1){    //Boundary III
        A[center][center] = a;
        A[center][rPlus] = b;
        A[center][rMinus] = c;
        A[center][thetaPlus] = 2*d;
        //printf("Boundary III\n");
      }
      else if((center+1) % N == 0){   //Boundary IV
        A[center][center] = 1;
        //printf("Boundary IV\n");
      }
      else if(center% N == 0){  //Boundary II
        A[center][center] = a;
        A[center][rPlus] = 2*b;
        A[center][thetaPlus] = d;
        A[center][thetaMinus] = d;
        //printf("Boundary II\n");
      }
      else{
        A[center][center] = a;
        A[center][rPlus] = b;
        A[center][rMinus] = c;
        A[center][thetaPlus] = d;
        A[center][thetaMinus] = d;
        //printf("Interior\n");
      }
    }
  }

/*write A matrix to file*/
  FILE *file;
  file = fopen("A.txt", "w");
  if(file == NULL){
    printf("Error opening file for B\n");
    return 1;
  }
  for(int i=0; i<n; i++){
    for(int j=0; j<2*N+1; j++){
      fprintf(file, "%f\t", A[i][j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);


  /////////////////////////////////////////////////////////////////////////
  //// solve system

  float A_banded[n][2*N+1];
  for(int i=0; i<n; i++){
    for(int j=0; j<2*N+1; j++){
      A_banded[i][j] = 0;
    }
  }
  
  //// banded Storage Mode 2
  for(int i=0; i<N; i++){
    for(int j=0; j<N+i+1; j++){
      A_banded[i][N+j-i] = A[i][j];
      //printf("%d i, %d j --> banded j %d\n", i, j, N+j-i);
    }
    //printf("\n");
  }
  for(int i=N; i<n-N; i++){
    for(int j=i-N; j<N+i+1; j++){
      A_banded[i][N+j-i] = A[i][j];
      //printf("%d i, %d j --> banded j %d\n", i, j, j+i-N);
    }
    //printf("\n");
  }
  for(int i=n-N; i<n; i++){
    for(int j=i-N; j<n; j++){
      A_banded[i][N+j-i] = A[i][j];
      //printf("%d i, %d j --> banded j %d\n", i, j, j+i-N);
    }
    //printf("\n");
  }

  //write A_banded matrix to file
  file = fopen("A_banded.txt", "w");
  if(file == NULL){
    printf("Error opening file for B\n");
    return 1;
  }
  for(int i=0; i<n; i++){
    for(int j=0; j<2*N+1; j++){
      fprintf(file, "%f\t", A_banded[i][j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
  
  float A_banded_vector[n*(2*N+1)+1];
  A_banded_vector[0] = 0;
  for(int i=0; i<n; i++){
    for(int j=0; j<2*N+1; j++){
      A_banded_vector[i*(2*N+1) + j + 1] = A[i][j];
    }
  }
  /*write A_banded_vector to file
  file = fopen("A_banded_vector.txt", "w");
  if(file == NULL){
    printf("Error opening file for B\n");
    return 1;
  }
  for(int i=0; i<n*(2*N+1)+1; i++){
    fprintf(file, "%f\n", A_banded_vector[i]);
  }
  fclose(file);*/

  solve1(3, A_banded_vector, B, n, N, 2*N+1);
  /*
  //write x vector to file
  file = fopen("x.txt", "w");
  if(file == NULL){
    printf("Error opening file for B\n");
    return 1;
  }
  for(int i=0; i<n; i++){
    fprintf(file, "%f\n", x[i]);
  }
  fclose(file);
  */
}
  return 0;
}

int nodemap(int i, int j, int N){
  return (j*N + i);
}