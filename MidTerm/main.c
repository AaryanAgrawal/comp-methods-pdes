/*
Author: Aaryan Agrawal
Course: ENGS105, W24
Instructor: Professor Keith Paulsen

Main Program - Midterm Problem - Analyse two finite-difference molecules A & B by comparing their transient solution 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utilfunc.h"

int main(void){
  
  //Define parameters
  float D_param = 0.5;
  float L_param = 10;
  float sigma_param = 0.1;
  int x0_param = 5;

  //Define spacing
  float h = 0.1;
  float k = 0.005;

  //Calculate Number of Nodes
  int N = (int)L_param/h + 1;

  //Parameter r
  double r = D_param*k/(h*h);
  printf("r = %f\n", r);

  //Define solution files
  FILE *fileA = fopen("solutionA.txt", "w");
  FILE *fileB = fopen("solutionB.txt", "w");

  //Define x linspace
  double x[N];
  for(int i=0; i<N; i++){
    x[i] = i*h;
  }

  //Define Initial Conditions
  double U0A[N], U1A[N], U2A[N], U0B[N], U1B[N], U2B[N];
  initialConditions(U0A, U1A, N, x, x0_param, sigma_param);
  initialConditions(U0B, U1B, N, x, x0_param, sigma_param);

  //Define Boundary Conditions
  boundaryConditions(U0A, U1A, N);
  boundaryConditions(U0B, U1B, N);

  //Loop through time for Molecule A
  for(int j=0; j<=100; j++){
    applyMoleculeA(U0A, U1A, U2A, N, r);
    printSolutionToFile(U2A, N, fileA);
    shiftTimeLevels(U0A, U1A, U2A, N);
  }

  //Loop through time for Molecule B
  for(int j=0; j<=100; j++){
    applyMoleculeB(U0B, U1B, U2B, N, r);
    printSolutionToFile(U2B, N, fileB);
    shiftTimeLevels(U0B, U1B, U2B, N);
  }


  
  return 0;
}