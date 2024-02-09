/*
Author: Aaryan Agrawal
Course: ENGS105, W24
Instructor: Professor Keith Paulsen

Problem 12c - with transient behaviour - with Type I Boundary Condition - with GAUSS-SEIDEL Point Iteration
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

int main (void) {
  int NRange[] = {80};
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

    //define weight
    double weight = 1;
    double antiWeight = weight - 1;
    double r_weight = 10;

    //define spacing
    double delR = (R_param-a_param)/(N-1);
    double delTheta = (PI/2)/(N-1);
    double delTime = r_weight*delR*delR;

    //store theta and radius values
    double theta[N];
    double radius[N];
    for (int i=0; i<N; i++){
      theta[i] = i*delTheta;
      radius[i] = a_param + i*delR;
    }

    ///////////////////////////////////////////////////////////////////////////
    //// make U0 vector - Initial Conditions + Boundary Conditions Type I
    double Ut[n];   //IC at t=0, U=0
    double UtMinus[n];   //Previous time step
    for(int i=0; i<n; i++){
      Ut[i] = 0;
    }
    for(int i=N-1; i<n; i+=N){    // Boundary IV
      Ut[i] = -cos(k_param*theta[(i+1)/N - 1]);
    }

    ///////////////////////////////////////////////////////////////////////////
    //// Transient Behaviour + Iterative Solve
    
    double steadyStateCheck = 10;
    int timeCnt = 0;
    while(timeCnt*delTime < 2){

      // Make initial guess
      double U0[n];
      double U1[n];
      for(int i=0; i<n; i++){
        U0[i] = 0;
        U1[i] = 0;
      }

      double a, b, c, d, bplusc, a_anti, b_anti, c_anti, d_anti, bplusc_anti;
      double rVal;
      double beta0, beta1, beta2, beta3, beta4, bij;
      int row, rPlus, rMinus, thetaPlus, thetaMinus;
      int iterationCnt = 0;
      //iterate over all nodes, check for boundary conditions

      double crit = 10;

      while(crit > 10e-5){
        for (int j=0; j<N; j++){
          for(int i=0; i<N; i++){
            rVal = radius[i];

            a = ( - 2/(delR*delR) - 2/(delTheta*delTheta) ) * delTime * weight;
            b = ( 1/(delR*delR) + 1/(2*delR*rVal) ) * delTime * weight;
            c = ( 1/(delR*delR) - 1/(2*delR*rVal) ) * delTime * weight;
            bplusc = ( 2/(delR * delR) ) * delTime * weight;
            d = ( 1 / (delTheta*delTheta) ) * delTime * weight;

            a_anti = ( - 2/(delR*delR) - 2/(delTheta*delTheta) ) * delTime * antiWeight;
            b_anti = ( 1/(delR*delR) + 1/(2*delR*rVal) ) * delTime * antiWeight;
            c_anti = ( 1/(delR*delR) - 1/(2*delR*rVal) ) * delTime * antiWeight;
            bplusc_anti = ( 2/(delR * delR) ) * delTime * antiWeight;
            d_anti = ( 1 / (delTheta*delTheta) ) * delTime * antiWeight;

            //row is created to check for positions for BCs based on zero indexing
            row = nodemap(i, j, N);
            rPlus = (row + 1);
            rMinus = (row - 1);
            thetaPlus = (row + N);
            thetaMinus = (row - N);

            if(row == 0){    //Type II corner
              beta0 = ( a_anti - 1 );
              beta1 = bplusc_anti;
              beta3 = 2 * d_anti;
              bij = (beta0 * Ut[row] + beta1 * Ut[rPlus] + beta3 * Ut[thetaPlus]);
              //printf("Corner\n");
              beta0 = - ( a - 1 );    //Negative because of iterative solving
              beta1 = bplusc;
              beta3 = 2*d;
              U1[row] = (1.0/beta0) * (beta1 * U0[rPlus] + beta3 * U0[thetaPlus] - bij);
            }
            else if((row+1) % N == 0){   //Boundary IV
              beta0 = - 1;
              bij = Ut[row];
              //printf("Boundary IV\n");
              U1[row] = (1.0/beta0) * (- bij);
            }
            else if(row >= n-N){   //Boundary I
              beta0 = - 1;
              bij = Ut[row];
              //printf("Boundary I\n");
              U1[row] = (1.0/beta0) * (- bij);
            }
            else if(row > 0 && row < N-1){    //Boundary III
              beta0 = (a_anti-1);
              beta1 = b_anti;
              beta2 = c_anti;
              beta3 = 2*d_anti;
              bij = (beta0 * Ut[row] + beta1 * Ut[rPlus] + beta2 * Ut[rMinus] + beta3 * Ut[thetaPlus]);
              //printf("Boundary III\n");
              beta0 = - (a-1);
              beta1 = b;
              beta2 = c;
              beta3 = 2*d;
              U1[row] = (1.0/beta0) * (beta1 * U0[rPlus] + beta2 * U0[rMinus] + beta3 * U0[thetaPlus] - bij);
            }
            else if(row % N == 0){  //Boundary II
              beta0 = a_anti-1;
              beta1 = bplusc_anti;
              beta3 = d_anti;
              beta4 = d_anti;
              bij = (beta0 * Ut[row] + beta1 * Ut[rPlus] + beta3 * Ut[thetaPlus] + beta4 * Ut[thetaMinus]);
              //printf("Boundary II\n");
              beta0 = -(a-1);
              beta1 = bplusc;
              beta3 = d;
              beta4 = d;
              U1[row] = (1.0/beta0) * (beta1 * U0[rPlus] + beta3 * U0[thetaPlus] + beta4 * U0[thetaMinus] - bij);
            }
            else{
              beta0 = a_anti-1;
              beta1 = b_anti;
              beta2 = c_anti;
              beta3 = d_anti;
              beta4 = d_anti;
              bij = (beta0 * Ut[row] + beta1 * Ut[rPlus] + beta2 * Ut[rMinus] + beta3 * Ut[thetaPlus] + beta4 * Ut[thetaMinus]);
              //printf("Interior\n");
              beta0 = -(a-1);
              beta1 = b;
              beta2 = c;
              beta3 = d;
              beta4 = d;
              U1[row] = (1.0/beta0) * (beta1 * U0[rPlus] + beta2 * U0[rMinus] + beta3 * U0[thetaPlus] + beta4 * U0[thetaMinus] - bij);
            }
            //printf("U1[%d] = %lf\t bij = %lf\n", row, U1[row], bij);
          }
        }
        
        double diff_vector[n];
        for(int i = 0; i<n; i++){
          diff_vector[i] = fabs(U1[i] - U0[i]);
        }
        crit = diff_vector[vector_max(diff_vector, n)];
        //printf("crit = %.8e\n", crit);
        
        //jacobiRMS[iterationCnt] = vector_rms(diff_vector, n);
        iterationCnt++;

        for(int k=0; k<n; k++){
          U0[k] = U1[k];
        }
      }
      //printf("iterationCnt = %d\n", iterationCnt);
      
      //Save and update Ut
      for(int i=0; i<n; i++){
        UtMinus[i] = Ut[i];
        Ut[i] = U1[i];
      }

      //Check for steady state - L-infinity norm between Ut, UtMinus
      double diff_vector[n];
      for(int i = 0; i<n; i++){
        diff_vector[i] = fabs(Ut[i] - UtMinus[i]);
      }
      steadyStateCheck = diff_vector[vector_max(diff_vector, n)];
      //steadyStateCheck = 0;
      //printf("steadyStateCheck = %.8e\n", steadyStateCheck);
      timeCnt++;
    }
    printf("Total Time = %f\n", timeCnt*delTime);

    //write U to file
    char buf[0x100];
    char sWeight[0x100];
    itoa(weight*100, sWeight, 10);
    snprintf(buf, sizeof(buf), "%s.txt", sWeight);
    printf("writing to %s\n", buf);
    file = fopen(buf, "w");
    for(int fi=0; fi<n; fi++){
      fprintf(file, "%f\n", Ut[fi]);
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