/*
Author: Aaryan Agrawal
Course: ENGS105, W24
Instructor: Professor Keith Paulsen

Problem 12c - with Type III Boundary Condition - with GAUSS-SEIDEL Point Iteration
Electric current flow past a circular hole in a conducting sheet
Solve PDE using FD molecule over polar grid
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592

int nodemap(int i, int j, int N);
int vector_max(double A[], int n);
float vector_rms(double A[], int n);
float vector_norm(double A[], int n);

int main (void) {

int NRange[] = {10, 20, 40};
int ARange[] = {991, 5, 1000};
for(int jj=0; jj < (int)sizeof(ARange)/sizeof(ARange[0]); jj++){
  float A_param = ARange[jj];
  if(ARange[jj] == 991){
    A_param = 0.001;
  }
  printf("A = %f\n", A_param);
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
      B[i] = - cos(k_param*theta[(i+1)/N - 1]) * (delR/2 + 1);
    }

    //write b vector to file
    /*file = fopen("b.txt", "w");
    if(file == NULL){
      printf("Error opening file for B\n");
      return 1;
    }
    for(int i=0; i<n; i++){
      fprintf(file, "%f\n", B[i]);
    }
    fclose(file);*/

    ///////////////////////////////////////////////////////////////////////////
    //// Iterative Solve

    // Make initial guess
    double U0[n], UMinus1[n];
    for(int i=0; i<n; i++){
      U0[i] = 0;
      UMinus1[i] = 0;
    }
    double U1[n];

    double a, b, c, d;
    double rVal;
    double beta0, beta1, beta2, beta3, beta4, bij;
    float spectral_radius[10000];
    double diff_vector[n], diff_vector_previous[n];
    int row, rPlus, rMinus, thetaPlus, thetaMinus;
    int iterationCnt = 0;
    float jacobiRMS[10000];
    for(int i=0; i<10000; i++){
      jacobiRMS[i] = 0.0;
      spectral_radius[i] = 0.0;
    }
    //iterate over all nodes, check for boundary conditions

    double crit = 10;

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
            beta1 = 2;
            beta3 = 2*d;
            bij = B[row];
            //printf("Corner\n");
            U1[row] = (1.0/beta0) * (beta1 * U0[rPlus] + beta3 * U0[thetaPlus] - bij);    //Jacobi
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
            U1[row] = (1.0/beta0) * (beta1 * U0[rPlus] + beta2 * U1[rMinus] + beta4 * U1[thetaMinus] - bij);
          }
          else if(row > 0 && row < N-1){    //Boundary III
            beta0 = - a;
            beta1 = b;
            beta2 = c;
            beta3 = 2*d;
            bij = B[row];
            //printf("Boundary III\n");
            U1[row] = (1.0/beta0) * (beta1 * U0[rPlus] + beta2 * U1[rMinus] + beta3 * U0[thetaPlus] - bij);
          }
          else if(row % N == 0){  //Boundary II
            beta0 = - a;
            beta1 = 2;
            beta3 = d;
            beta4 = d;
            bij = B[row];
            //printf("Boundary II\n");
            U1[row] = (1.0/beta0) * (beta1 * U0[rPlus] + beta3 * U0[thetaPlus] + beta4 * U1[thetaMinus] - bij);
          }
          else{
            beta0 = - a;
            beta1 = b;
            beta2 = c;
            beta3 = d;
            beta4 = d;
            bij = B[row];
            //printf("Interior\n");
            U1[row] = (1.0/beta0) * (beta1 * U0[rPlus] + beta2 * U1[rMinus] + beta3 * U0[thetaPlus] + beta4 * U1[thetaMinus] - bij);
          }
          //printf("U1[%d] = %lf\t bij = %lf\n", row, U1[row], bij);
        }
      }

      for(int i = 0; i<n; i++){
        diff_vector[i] = fabs(U1[i] - U0[i]);
      }
      crit = diff_vector[vector_max(diff_vector, n)];
      //printf("crit = %.8e\n", crit);

      for(int i = 0; i<n; i++){
      diff_vector_previous[i] = fabs(UMinus1[i] - U0[i]);
      }
      spectral_radius[iterationCnt] = vector_norm(diff_vector, n)/vector_norm(diff_vector_previous, n);
      iterationCnt++;

      for(int k=0; k<n; k++){
        UMinus1[k] = U0[k];
        U0[k] = U1[k];
      }
    }

    float mean_spectral_radius = 0;
    for(int i=1; i<iterationCnt; i++){
      mean_spectral_radius += spectral_radius[i];;
    }
    mean_spectral_radius = mean_spectral_radius/(iterationCnt-1);
    printf("iterationCnt = %d\n", iterationCnt);
    printf("Spectral Radius = %f\n", mean_spectral_radius);

    //write U to file
    char buf[0x100];
    char sN[0x100];
    char sA[0x100];
    itoa((int)ARange[jj], sA, 10);
    itoa(N, sN, 10);
    snprintf(buf, sizeof(buf), "N%s_A%s.txt", sN, sA);
    printf("writing to %s\n", buf);
    file = fopen(buf, "w");
    for(int fi=0; fi<n; fi++){
      fprintf(file, "%f\n", U1[fi]);
      //printf("%f\n", r[i]);
    }

    // //write RMS to file
    // snprintf(buf, sizeof(buf), "N%s_A%s_RMS_jacobi.txt", sN, sA);
    // printf("writing to %s\n", buf);
    // file = fopen(buf, "w");
    // for(int fi=0; fi<iterationCnt; fi++){
    //   fprintf(file, "%f\n", jacobiRMS[fi]);
    //   //printf("%f\n", r[i]);
    // }
    // fclose(file);

  }
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

float vector_rms(double A[], int n){
  float sum = 0;
  for(int i=0; i<n; i++){
    sum += (float)(A[i]*A[i]);
  }
  return sqrt(sum/n);
}

float vector_norm(double A[], int n){
  float sum = 0;
  for(int i=0; i<n; i++){
    sum += (float)(A[i]*A[i]);
  }
  return sqrt(sum);
}