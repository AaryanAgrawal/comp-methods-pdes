/*
Author: Aaryan Agrawal
Course: ENGS105, W24
Instructor: Professor Keith Paulsen

Utility Functions for time transient analysis of finite-difference molecules A & B
*/

//print out solutions to file

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solve1.h"

void applyMoleculeA(double U0[], double U1[], double U2[], int N, double r){
  //Increments a time step for Molecule A
  for(int i=1; i<N-1; i++){
    U2[i] = (1.5*r) * (U1[i-1]+U1[i+1]) + (1.0-3.0*r) * U1[i];  // l time level
    U2[i] += r * U0[i] - (0.5*r) * (U0[i+1] + U0[i-1]);  // l-1 time level
  }
  U2[0] = U1[0];
  U2[N-1] = U1[N-1];
}

void applyMoleculeB(double U0[], double U1[], double U2[], int N, double r){
  //Increments a time step for Molecule B
  double bij;
  double beta0 = - (1.0 + 5.0*r/6.0);
  double beta1 = -5.0*r/12.0;
  double beta2 = -5.0*r/12.0;
  double crit = 10.0;
  double diff_vector[N];
  double tol = 1e-6;
  int iterationCnt = 0;
  //Make Initial Guess
  double UGuess[N];
  for(int i=0; i<N; i++){
    UGuess[i] = U1[i];
  }
  //iteratively solve for each time step
  while(crit > tol){
    iterationCnt++;
    U2[0] = U1[0];
    U2[N-1] = U1[N-1];
    for(int i=1; i<N-1; i++){
      bij = (1-4.0*r/3.0) * U1[i] + (2.0*r/3.0) * (U1[i-1] + U1[i+1]);  // l time level
      bij += (-r/12.0) * (U0[i+1] + U0[i-1]) + (r/6.0) * U0[i];   // l-1 time level
      U2[i] = (1/beta0) * (beta1*U2[i-1] + beta2*UGuess[i+1] - bij);
    }
    //Calculate Vector Norm of Difference
    for(int i=0; i<N; i++){
      diff_vector[i] = fabs(U2[i] - UGuess[i]);
    }
    crit = diff_vector[vector_max(diff_vector, N)];
    //printf("crit = %f\n", crit);
    //Update UGuess
    for(int i=0; i<N; i++){
      UGuess[i] = U2[i];
    }
  }
  //printf("iterationCnt = %d\n", iterationCnt);
}

void printSolutionToFile(double U[], int N, FILE *file){
  //Prints solution to file - columns are nodes, rows are time steps
  for(int i=0; i<N; i++){
    fprintf(file, "%.15le\t", U[i]);
  }
  fprintf(file, "\n");
}

void shiftTimeLevels(double U0[], double U1[], double U2[], int N){
  //Shifts time levels for U0, U1, and U2
  for(int i=0; i<N; i++){
    U0[i] = U1[i];
    U1[i] = U2[i];
  }
}

void boundaryConditions(double U0[], double U1[], int N){
  //Applies boundary conditions to U0 and U1
  U0[0] = 0;
  U0[N-1] = 0;
  U1[0] = 0;
  U1[N-1] = 0;
}

void initialConditions(double U0[], double U1[], int N, double x[], float x0_param, float sigma_param){
  //Applies initial conditions to U0 and U1
  for(int i=0; i<N; i++){
    U0[i] = exp(-(x[i] - x0_param)*(x[i] - x0_param)/(2*sigma_param*sigma_param));
    U1[i] = U0[i];
  }
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
