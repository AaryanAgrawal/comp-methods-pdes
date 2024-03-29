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

int NRange[] = {120, 160};

for(int ii=0; ii < (int)sizeof(NRange)/sizeof(NRange[0]); ii++){
  FILE *file;
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
  float B[n+1];
  for(int i=0; i<n+1; i++){
    B[i] = 0;
  }
  for(int i=N; i<=n; i+=N){
    B[i] = - cos(k_param*theta[i/N - 1]);
  }

  //write b vector to file
  /*file = fopen("b.txt", "w");
  if(file == NULL){
    printf("Error opening file for B\n");
    return 1;
  }
  for(int i=0; i<n; i++){
    fprintf(file, "%.7f\n", B[i]);
  }
  fclose(file);*/
  
  /////////////////////////////////////////////////////////////////////////
  //// make A matrix
  float A[n+1][2*N+2];
  for(int i=0; i<n+1; i++){
    for(int j=0; j<2*N+2; j++){
      A[i][j] = 0;
    }
  }
  float a, b, c, d;
  float rVal;
  int rowCondition, row, center, rPlus, rMinus, thetaPlus, thetaMinus;
  //iterate over all nodes, check for boundary conditions
  for (int j=0; j<N; j++){
    for(int i=0; i<N; i++){
      rVal = radius[i];

      a = - 2 - (2*beta)/(rVal*rVal);
      b = 1 + delR/(2*rVal);
      c = 1 - delR/(2*rVal);
      d = beta / (rVal*rVal);

      //// Note: Added +1 to all indices to account for empty rows and columns
      rowCondition = nodemap(i, j, N);
      row = nodemap(i, j, N) + 1; //Added +1 to account for empty rows
      //printf("row: %d\n", row);
      // For Banded Storage Mode 2 for A, A is n x 2N+1. Center of molecule is at index=N.
      center = N + 1; //Added +1 to account for empty columns
      rPlus = (center + 1);
      rMinus = (center - 1);
      thetaPlus = (center + N);
      thetaMinus = (center - N);

      if(rowCondition == 0){    //Type II corner
        A[row][center] = a;
        A[row][rPlus] = 2;
        A[row][thetaPlus] = 2*d;
        //printf("Corner\n");
      }
      else if(rowCondition >= n-N){   //Boundary I
        A[row][center] = 1;
        //printf("Boundary I\n");
      }
      else if(rowCondition > 0 && rowCondition < N-1){    //Boundary III
        A[row][center] = a;
        A[row][rPlus] = b;
        A[row][rMinus] = c;
        A[row][thetaPlus] = 2*d;
        //printf("Boundary III\n");
      }
      else if((rowCondition+1) % N == 0){   //Boundary IV
        A[row][center] = 1;
        //printf("Boundary IV\n");
      }
      else if(rowCondition % N == 0){  //Boundary II
        A[row][center] = a;
        A[row][rPlus] = 2;
        A[row][thetaPlus] = d;
        A[row][thetaMinus] = d;
        //printf("Boundary II\n");
      }
      else{
        A[row][center] = a;
        A[row][rPlus] = b;
        A[row][rMinus] = c;
        A[row][thetaPlus] = d;
        A[row][thetaMinus] = d;
        //printf("Interior\n");
      }
    }
  }

//write A matrix to file
  /*file = fopen("A.txt", "w");
  if(file == NULL){
    printf("Error opening file for B\n");
    return 1;
  }
  for(int i=0; i<n+1; i++){
    for(int j=0; j<2*N+2; j++){
      fprintf(file, "%f\t", A[i][j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);*/


  /////////////////////////////////////////////////////////////////////////
  //// solve system
  /* Banded Storage integrated into formation of matrix A 
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

  write A_banded matrix to file
  FILE *file;
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
  fclose(file);*/
  /*
  float A_banded_vector[n*(2*N+1)+1];
  A_banded_vector[0] = 0;
  for(int i=0; i<n; i++){
    for(int j=0; j<2*N+1; j++){
      A_banded_vector[i*(2*N+1) + j + 1] = A[i][j];
    }
  }*/
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

  solve1(3, A, B, n, N, 2*N+2);
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