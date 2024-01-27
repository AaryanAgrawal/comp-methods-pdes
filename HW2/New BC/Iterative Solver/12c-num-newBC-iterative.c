/*
Author: Aaryan Agrawal
Course: ENGS105, W24
Instructor: Professor Keith Paulsen

Problem 12c - with Type III Boundary Condition - with Point Iterative Solving
Electric current flow past a circular hole in a conducting sheet
Solve PDE using FD molecule over polar grid
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592

int nodemap(int i, int j, int N);
int vector_max(double A[], int n);

int main (void) {

int NRange[] = {15};

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
  float A_param = 1;
  int D_param = -3;

  //define spacing
  double delR = (R_param-a_param)/(N-1);
  double delTheta = (PI/2)/(N-1);
  double beta = (delR/delTheta) * (delR/delTheta);

  //store theta and radius values
  double theta[N];
  double radius[N];
  for (int i=0; i<N; i++){
    theta[i] = i*delTheta;
    radius[i] = a_param + i*delR;
  }

  ///////////////////////////////////////////////////////////////////////////
  //// make b vector
  double B[n];
  for(int i=0; i<n; i++){
    B[i] = 0;
  }
  for(int i=n-N; i<n; i++){    // Boundary I
    B[i] = - 2 * D_param * beta * delTheta / radius[i-(n-N)];
  }
  for(int i=N-1; i<n; i+=N){    // Boundary IV
    B[i] = cos(k_param*theta[(i+1)/N - 1]) * (delR + 1);
  }


  //write b vector to file
  file = fopen("b.txt", "w");
  if(file == NULL){
    printf("Error opening file for B\n");
    return 1;
  }
  for(int i=0; i<n; i++){
    fprintf(file, "%f\n", B[i]);
  }
  fclose(file);

  ///////////////////////////////////////////////////////////////////////////
  //// Iterative Solve

  // Make initial guess
  double U0[n];
  for(int i=0; i<n; i++){
    U0[i] = 0;
  }
  double U1[n];
  
  double a, b, c, d;
  double rVal;
  double beta0, beta1, beta2, beta3, beta4, bij;
  int row, rPlus, rMinus, thetaPlus, thetaMinus;
  //iterate over all nodes, check for boundary conditions

  double crit = 10;
  double U0Max, U1Max;

  while(crit > 10e-5){
    for (int j=0; j<N; j++){
      for(int i=0; i<N; i++){
        rVal = radius[i];

        a = - 2 - (2*beta)/(rVal*rVal);
        b = 1 + delR/(2*rVal);
        c = 1 - delR/(2*rVal);
        d = beta / (rVal*rVal);

        //row is created to check for positions for BCs based on zero indexing
        row = nodemap(i, j, N);
        rPlus = (row + 1);
        rMinus = (row - 1);
        thetaPlus = (row + N);
        thetaMinus = (row - N);

        if(row == 0){    //Type II corner
          beta0 = - a;
          beta1 = 2*b;
          beta3 = 2*d;
          bij = B[row];
          //printf("Corner\n");
          U1[row] = (1.0/beta0) * (beta1 * U0[rPlus] + beta3 * U0[thetaPlus] - bij);
        }
        else if((row+1) % N == 0){   //Boundary IV
          beta0 = - 1;
          bij = B[row];
          //printf("Boundary IV\n");
          U1[row] = (1.0/beta0) * (- bij);
        }
        else if(row >= n-N){   //Boundary I
          beta0 = - ( a - 2*A_param*delTheta*beta/rVal );
          beta1 = b;
          beta2 = c;
          beta4 = 2*d;
          bij = B[row];
          //printf("Boundary I\n");
          U1[row] = (1.0/beta0) * (beta1 * U0[rPlus] + beta2 * U0[rMinus] + beta4 * U0[thetaMinus] - bij);
        }
        else if(row > 0 && row < N-1){    //Boundary III
          beta0 = - a;
          beta1 = b;
          beta2 = c;
          beta3 = 2*d;
          bij = B[row];
          //printf("Boundary III\n");
          U1[row] = (1.0/beta0) * (beta1 * U0[rPlus] + beta2 * U0[rMinus] + beta3 * U0[thetaPlus] - bij);
        }
        else if(row % N == 0){  //Boundary II
          beta0 = - a;
          beta1 = 2*b;
          beta3 = d;
          beta4 = d;
          bij = B[row];
          //printf("Boundary II\n");
          U1[row] = (1.0/beta0) * (beta1 * U0[rPlus] + beta3 * U0[thetaPlus] + beta4 * U0[thetaMinus] - bij);
        }
        else{
          beta0 = - a;
          beta1 = b;
          beta2 = c;
          beta3 = d;
          beta4 = d;
          bij = B[row];
          //printf("Interior\n");
          U1[row] = (1.0/beta0) * (beta1 * U0[rPlus] + beta2 * U0[rMinus] + beta3 * U0[thetaPlus] + beta4 * U0[thetaMinus] - bij);
        }
        printf("U1[%d] = %lf\t bij = %lf\n", row, U1[row], bij);
      }
    }
    
    U0Max = U0[vector_max(U0, n)];
    printf("U0Max = %lf\n", U0Max);
    U1Max = U1[vector_max(U1, n)];
    printf("U1Max = %lf\n", U1Max);
    crit = fabs(U1Max - U0Max) / fabs(U1Max);
    printf("crit = %.8e\n", crit);
    for(int i=0; i<n; i++){
      U0[i] = U1[i];
    }
  }

  //write answer to file
  file = fopen("jacobi.txt", "w");
  if(file == NULL){
    printf("Error opening file for B\n");
    return 1;
  }
  for(int i=0; i<n; i++){
    fprintf(file, "%.9lf\n", U1[i]);
  }
  fclose(file);

  ///////////////////////////////////////////////////////////////////////////
  /// Iterative Solving

  //solve1(3, A, B, n, N, 2*N+2, 991);
  
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

int vector_max(double A[], int n){
  double max = fabs(A[0]);
  int maxIndex = 0;
  for(int i=0; i<n; i++){
    if(fabs(A[i]) > max){
      max = fabs(A[i]);
      maxIndex = i;
    }
  }
  return maxIndex;
}