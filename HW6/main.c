/*
Author: Aaryan Agrawal
ENGG 105: Comp Methods to PDEs I
Instructor: Prof. Paulsen

Homework 6

Solve Simplified Electrical Potential Problem using Finite Element Method
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 502    //# of Nodes
#define M 917    //# of Elements
#define HALFBW 23 // From MATLAB script
#include "solve1.h"

void scanNodes(float *xNodes, float *yNodes);
void covNoise(float xNodes[N+1], float yNodes[N+1], float sigma2, int *bcNodes, float CovB[N+1][N+1], float corLength);
void scanElems(int *elems);
int IN(int i, int j);
void ElementMatrix(int L, int *elems, float *xNodes, float *yNodes, float Ae[3][3], float be[3]);
void boundaryConditionsI(float xNodes[N+1], float yNodes[N+1], float A[N+1][2*HALFBW+2], float b[N+1]);
void gradient(int *elems, float *xNodes, float *yNodes, float *u, float grad[M+1][2]);
void matrixMult(float RHS[N+1][N+1], float A[N+1][N+1], float x[N+1][N+1]);
void printToFile(float A[N+1][N+1], char *filename);

int main(void){
  float xNodes[N+1];
  float yNodes[N+1];
  int elems[M*3+1];

  //Load Input Mesh
  scanNodes(xNodes, yNodes);
  scanElems(elems);

//////////////////////////////////////////////////////////////////////////////////////

  //Element Asssembly
  float Ae[3][3];
  float be[3];
  float b[N+1];
  float A[N+1][2*HALFBW+1+1];
  //Initialize A and b
  for (int i=0; i <= N; i++){
    b[i] = 0.0;
    for (int j=0; j <= 2*HALFBW+1; j++){
      A[i][j] = 0.0;
    }
  }
  int II, JJ, JB;
  for (int L=1; L <= M; L++) {
    ElementMatrix(L, elems, xNodes, yNodes, Ae, be);
    //print Ae be

    for (int I=1; I <= 3; I++) {

      II = elems[IN(L, I)];
      //printf("II = %d, L = %d, I = %d\n", II, L, I);
      b[II] = b[II] + be[I-1];
      for (int J=1; J <= 3; J++) {
        JJ = elems[IN(L, J)];
        //printf("II = %d, JJ = %d, L = %d, I = %d, J = %d\n", II, JJ, L, I, J);
        JB = (HALFBW + 1) + JJ - II;
        A[II][JB] = A[II][JB] + Ae[I-1][J-1];
        //printf("JB = %d, A[%d][%d] = %f\n", JB, II, JB, A[II][JB]);
      }
    }
  }


  //Line Sink
  float sinkCurrent = 1.0 / 8.0;
  b[492] += (1.0 + xNodes[492]/10.0 + yNodes[492]/5.0) * sinkCurrent;
  b[493] += (1.0 + xNodes[493]/10.0 + yNodes[493]/5.0) * sinkCurrent;

  //Boundary Conditions
  boundaryConditionsI(xNodes, yNodes, A, b);

  //Copy b to b0
  float b0[N+1];
  for(int i=0; i <= N; i++){
    b0[i] = b[i];
  }

//////////////////////////////////////////////////////////////////////////////////////

  //Print A to file
  printf("Printing A to file\n");
  FILE *fileA;
  fileA = fopen("A.dat", "w");
  for (int i=0; i <= N; i++){
    for (int j=0; j <= 2*HALFBW +1; j++){
      fprintf(fileA, "%f\t", A[i][j]);
    }
    fprintf(fileA, "\n");
  }
  fclose(fileA);
  //Print b to file
  printf("Printing b to file\n");
  FILE *file;
  file = fopen("b.dat", "w");
  for (int i=0; i <= N; i++){
    fprintf(file, "%f\n", b[i]);
  }
  fclose(file);

//////////////////////////////////////////////////////////////////////////////////////

  //Solve
  solve1(3, A, b, N, HALFBW, 2*HALFBW+2);

//////////////////////////////////////////////////////////////////////////////////////

  //Find Gradient
  float grad[M+1][2];
  gradient(elems, xNodes, yNodes, b, grad);

  //Print grad to file
  printf("Printing grad to file\n");
  file = fopen("grad.dat", "w");
  for (int i=1; i <= M; i++){
    fprintf(file, "%f\t%f\n", grad[i][0], grad[i][1]);
  }
  fclose(file);  

//////////////////////////////////////////////////////////////////////////////////////

  //Print b as u to file
  printf("Printing u to file\n");
  file = fopen("u.dat", "w");
  for (int i=0; i <= N; i++){
    fprintf(file, "%f\n", b[i]);
  }
  fclose(file);

//////////////////////////////////////////////////////////////////////////////////////

  //Load inverse matrix
  float invA[N+1][N+1];
  file = fopen("Ainv.dat", "r");
  for (int i=0; i <= N; i++){
    for (int j=0; j <= N; j++){
      fscanf(fileA, "%e", &invA[i][j]);
    }
  }
  fclose(file);
  //Transpose invA into invA_T
  float invA_T[N+1][N+1];
  for(int i=0; i <= N; i++){
    for(int j=0; j <= N; j++){
      invA_T[i][j] = invA[j][i];
    }
  }


  // ///////////// Case 1
  int sourceElem = 288;
  float var = 0.5;
  float sigma2 = var * var;
  float bNoise[N+1];
  float CovB[N+1][N+1];
  float CovU_temp[N+1][N+1];
  float CovU[N+1][N+1];
  for(int i=0; i<=N; i++){
    bNoise[i] = 0.0;
  }

  for(int j=0; j<3; j++){
    bNoise[elems[IN(sourceElem, j+1)]] = b[elems[IN(sourceElem, j+1)]] * var;
  }

  for(int i = 0; i <= N; i++){
    for(int j = 0; j <= N; j++){
      CovB[i][j] = bNoise[i] * bNoise[j];
    }
  }

  matrixMult(CovU_temp, invA, CovB); 
  matrixMult(CovU, CovU_temp, invA_T);

  //Print CovU, CovB to file
  printf("Printing CovU to file\n");
  printToFile(CovU, "CovU_I.dat");
  printf("Printing CovB to file\n");
  printToFile(CovB, "CovB_I.dat");

  // ///////////// Case 2
  int sinkNodes[] = {492, 493};
  var = 0.4;
  sigma2 = var * var;

  for(int i=0; i<=N; i++){
    bNoise[i] = 0.0;
  }

  for(int j=0; j<2; j++){
    bNoise[sinkNodes[j]] = b[sinkNodes[j]] * var;
  }

  for(int i = 0; i <= N; i++){
    for(int j = 0; j <= N; j++){
      CovB[i][j] = bNoise[i] * bNoise[j];
    }
  }

  matrixMult(CovU_temp, invA, CovB); 
  matrixMult(CovU, CovU_temp, invA_T);

  //Print CovU, CovB to file
  printf("Printing CovU to file\n");
  printToFile(CovU, "CovU_II.dat");
  printf("Printing CovB to file\n");
  printToFile(CovB, "CovB_II.dat");

  // ///////////// Case 3
  var = 0.5;
  sigma2 = var * var;
  FILE *fp;
  fp = fopen("hw44.dnd.dat", "r");
  if (fp == NULL){
    printf("Error: file not found\n");
    exit(1);
  }
  int bcNodes[N];
  float temp;
  int bcNodesCnt = 0;
  while(fscanf(fp, "%d %f", &bcNodes[bcNodesCnt], &temp) != EOF){
    bcNodesCnt++;
  }
  fclose(fp);

  for(int i=0; i<=N; i++){
    bNoise[i] = 0.0;
  }

  for(int j=0; j<33; j++){
    bNoise[bcNodes[j]] = var;
  }

  //For Correlation Length = 0
  for(int i = 0; i <= N; i++){
    for(int j = 0; j <= N; j++){
      CovB[i][j] = bNoise[i] * bNoise[j];
    }
  }
  matrixMult(CovU_temp, invA, CovB); 
  matrixMult(CovU, CovU_temp, invA_T);
  //Print CovU, CovB to file
  printf("Printing CovU to file\n");
  printToFile(CovU, "CovU_IIIA.dat");
  printf("Printing CovB to file\n");
  printToFile(CovB, "CovB_IIIA.dat");

  //For Correlation Length = 2, 5
  float corLength[] = {2.0/6.5, 5.0/6.5};

  covNoise(xNodes, yNodes, sigma2, bcNodes, CovB, corLength[0]);
  matrixMult(CovU_temp, invA, CovB);
  matrixMult(CovU, CovU_temp, invA_T);
  //Print CovU, CovB to file
  printf("Printing CovU to file\n");
  printToFile(CovU, "CovU_IIIB.dat");
  printf("Printing CovB to file\n");
  printToFile(CovB, "CovB_IIIB.dat");

  covNoise(xNodes, yNodes, sigma2, bcNodes, CovB, corLength[1]);
  matrixMult(CovU_temp, invA, CovB);
  matrixMult(CovU, CovU_temp, invA_T);
  //Print CovU, CovB to file
  printf("Printing CovU to file\n");
  printToFile(CovU, "CovU_IIIC.dat");
  printf("Printing CovB to file\n");
  printToFile(CovB, "CovB_IIIC.dat");

}

void covNoise(float xNodes[N+1], float yNodes[N+1], float sigma2, int *bcNodes, float CovB[N+1][N+1], float corLength){

  float r;
  int I, J;
  //Set CovB to 0
  for(int i=0; i <= N; i++){
    for(int j=0; j <= N; j++){
      CovB[i][j] = 0.0;
    }
  }
  //i and j loops from 0 to 32
  for(int i=0; i < 33; i++){
    I = bcNodes[i];
    for(int j=0; j < 33; j++){
      J = bcNodes[j];
      r = sqrt((xNodes[I] - xNodes[J])*(xNodes[I] - xNodes[J]) + (yNodes[I] - yNodes[J])*(yNodes[I] - yNodes[J]));
      CovB[I][J] = sigma2 * (1 + r/corLength) * exp(-r/corLength);
    }
  }
  return sigma2 * exp(-r);
}

void printToFile(float A[N+1][N+1], char *filename){
  FILE *file;
  file = fopen(filename, "w");
  for (int i=0; i <= N; i++){
    for (int j=0; j <= N; j++){
      fprintf(file, "%e\t", A[i][j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);
}

void gradient(int *elems, float *xNodes, float *yNodes, float *u, float grad[M+1][2]){
  float dudx, dudy;
  float dely[3], delx[3];  //Zero Indexed
  float area = 0.0;

  //For each element
  for(int L=1; L<=M; L++){

    //Find delx, dely and area
    area = 0.0;
    for (int i=0; i<3; i++) {
      dely[i] = yNodes[elems[IN(L, (i+1)%3+1)]] - yNodes[elems[IN(L, (i+2)%3+1)]];
      delx[i] = xNodes[elems[IN(L, (i+1)%3+1)]] - xNodes[elems[IN(L, (i+2)%3+1)]];
    }
    for (int i=0; i<3; i++) {
      area += 0.5 * (xNodes[elems[IN(L, i+1)]] * dely[i]);
    }

    //Find gradient to be placed at centroid
    grad[L][0] = 0.0;
    grad[L][1] = 0.0;
    for(int i=0; i<3; i++){
      dudx = dely[i]/(2*area);
      dudy = -delx[i]/(2*area);
      grad[L][0] += u[elems[IN(L, i+1)]] * dudx;
      //printf("grad[%d][0] = %f, area=%f, node j=%d \n", L, grad[L][0], area, elems[IN(L, i+1)]);
      grad[L][1] += u[elems[IN(L, i+1)]] * dudy;
    }
  }
}

void boundaryConditionsI(float xNodes[N+1], float yNodes[N+1], float A[N+1][2*HALFBW+2], float b[N+1]){

  FILE *fp;
  fp = fopen("hw44.dnd.dat", "r");
  if (fp == NULL){
    printf("Error: file not found\n");
    exit(1);
  }
  int bcNodes[N];
  float temp;
  int i = 0;
  while(fscanf(fp, "%d %f", &bcNodes[i], &temp) != EOF){
    i++;
  }
  fclose(fp);

  for(int j = 0; j<i; j++){
    // printf("bcNodes[%d] = %d\n", j, bcNodes[j]);
    for(int k = 0; k<2*HALFBW+2; k++){
      A[bcNodes[j]][k] = 0.0;
    }
    A[bcNodes[j]][HALFBW+1] = 1.0;
    b[bcNodes[j]] = 0.0;
  
  }
}

void scanNodes(float *xNodes, float *yNodes) {
  FILE *file;
  file = fopen("hw44.nod.dat", "r");
  if (file == NULL){
    printf("Error: file not found\n");
    exit(1);
  }
  int temp;
  for(int i = 1; i <= N; i++){
    fscanf(file, "%d %f %f", &temp, &xNodes[i], &yNodes[i]);
  }
  fclose(file);
}

void scanElems(int *elems) {
  FILE *file;
  file = fopen("hw44.ele.dat", "r");
  if (file == NULL){
    printf("Error: file not found\n");
    exit(1);
  }
  int temp;
  for (int i=1; i <= M; i++){
    fscanf(file, "%d %d %d %d %d", &temp, &elems[IN(i, 1)], &elems[IN(i, 2)], &elems[IN(i, 3)], &temp);
  }
  fclose(file);
}

int IN(int i, int j) {
  // Takes in L (Element Number) and I (Node Number). Outputs Global Node Number
  return (i-1)*3 + j;
}

void ElementMatrix(int L, int *elems, float *xNodes, float *yNodes, float Ae[3][3], float be[3]) {
  float dely[3], delx[3];  //Zero Indexed
  float area = 0.0;
  float t_coeff = 0.0;
  for (int i=0; i<3; i++) {
    dely[i] = yNodes[elems[IN(L, (i+1)%3+1)]] - yNodes[elems[IN(L, (i+2)%3+1)]];
    delx[i] = xNodes[elems[IN(L, (i+1)%3+1)]] - xNodes[elems[IN(L, (i+2)%3+1)]];
  }
  for (int i=0; i<3; i++) {
    area += 0.5 * (xNodes[elems[IN(L, i+1)]] * dely[i]);
    // if(L==1){
    // printf("xNodes = %f, elems = %d, L = %d, i+1 = %d, dely = %f\n", xNodes[elems[IN(L, i+1)]], elems[IN(L, i+1)], L, i+1, dely[i]);
    // }
  }
  //printf("area = %f\n", area);
  
  for(int j=0; j<3; j++){
    t_coeff += (1.0 + xNodes[elems[IN(L, j)]]/10.0 + yNodes[elems[IN(L, j)]]/5.0);
  }
  t_coeff = t_coeff / 3.0;

  //Ae Matrix
  for (int i=0; i < 3; i++){
    for (int j=0; j < 3; j++){
      Ae[i][j] = t_coeff * (- dely[i]*dely[j] - delx[i]*delx[j]) / (4 * area);
    }
    //printf("\n");
  }

  //be Vector w/ Type II Boundary Conditions
  for (int j=0; j<3; j++){
    be[j] = 0.0;
  }

  float sourceCurrent = 1.0;
  int elemPoint = 288; //From MATLAB script
  float x0 = -0.206052;
  float y0 = 0.56379;
  float a, b, c, p;
  if(L==288){
    for(int j=0; j<3; j++){
      a = ( xNodes[elems[IN(L, (j+1)%3+1)]] * yNodes[elems[IN(L, (j+2)%3+1)]] - xNodes[elems[IN(L, (j+2)%3+1)]] * yNodes[elems[IN(L, (j+1)%3+1)]] ) / (2.0*area);
      b = dely[j]/(2.0*area);
      c = -delx[j]/(2.0*area);
      p = a + b*x0 + c*y0;
      be[j] += - p * sourceCurrent;
    }
  }
}

void matrixMult(float RHS[N+1][N+1], float A[N+1][N+1], float x[N+1][N+1]){
  //Multiply banded matrix A with x


  for(int i=0; i <= N; i++){    //rows
    for(int j=0; j <= N; j++){    //columns
      RHS[i][j] = 0.0;

      for(int k=0; k <= N; k++){
        RHS[i][j] += A[i][k] * x[k][j];
      }
    }
  }

}