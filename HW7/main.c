/*
Author: Aaryan Agrawal
ENGG 105: Comp Methods to PDEs I
Instructor: Prof. Paulsen

Homework 7

Solve Simplified Electrical Potential Problem using Finite Element Method and Parameter Minimization
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 502    //# of Nodes
#define M 917    //# of Elements
#define HALFBW 23 // From MATLAB script
#include "solve1.h"
#define N_POINTS 12

void scanNodes(float *xNodes, float *yNodes);
void covNoise(float xNodes[N+1], float yNodes[N+1], float sigma2, int *bcNodes, float CovB[N+1][N+1], float corLength);
void scanElems(int *elems);
int IN(int i, int j);
void ElementMatrix(int L, int *elems, float *xNodes, float *yNodes, float Ae[3][3], float be[3]);
void boundaryConditionsI(float xNodes[N+1], float yNodes[N+1], float A[N+1][2*HALFBW+2], float b[N+1]);
void gradient(int *elems, float *xNodes, float *yNodes, float *u, float grad[M+1][2]);
void matrixMult(float RHS[N+1][N+1], float A[N+1][N+1], float x[N+1][N+1]);
void printToFile(float A[N+1][N+1], char *filename);
void samplingMatrix(float xNodes[N+1], float yNodes[N+1], int elems[3*M+1], float samplingElem[N_POINTS], float xSam[N_POINTS], float ySam[N_POINTS], float u[N+1], float S[N_POINTS+1][N+1], float measU[N_POINTS]);
void matrixVectorMult(float RHS[N+1], float A[N+1][N+1], float x[N+1]);

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
  float Afull[N+1][N+1];
  //Initialize A and b
  for (int i=0; i <= N; i++){
    b[i] = 0.0;
    for (int j=0; j <= 2*HALFBW+1; j++){
      A[i][j] = 0.0;
    }
    for (int j=0; j <= N; j++){
      Afull[i][j] = 0.0;
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
        Afull[II][JJ] = Afull[II][JJ] + Ae[I-1][J-1];
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
  //Transpose Afull into Afull_T
  float Afull_T[N+1][N+1];
  for(int i=0; i <= N; i++){
    for(int j=0; j <= N; j++){
      Afull_T[i][j] = Afull[j][i];
    }
  }

//////////////////////////////////////////////////////////////////////////////////////

  /////////////// Type I Nodes
  float var = 0.1;
  float sigma2 = var * var;
  float bNoise[N+1];
  float CovB[N+1][N+1];
  float CovU[N+1][N+1];
  float CovU_temp[N+1][N+1];
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

  //For Correlation Length = 0.3
  float corLength = 0.3;

  covNoise(xNodes, yNodes, sigma2, bcNodes, CovB, corLength);
  matrixMult(CovU_temp, invA, CovB);
  matrixMult(CovU, CovU_temp, invA_T);

  ///////////////// Interior Nodes
  float covBII[N+1][N+1];
  float covUII[N+1][N+1];
  float covUII_temp[N+1][N+1];
  var = 1.0;
  sigma2 = var * var;
  corLength = 0.4;
  covNoise(xNodes, yNodes, sigma2, , covBII, corLength);

//////////////////////////////////////////////////////////////////////////////////////

  ////Find Sampling Matrix S based on sampling points
  //Input sampling point coordinates
    float xSam[N_POINTS] = {-0.6, -0.3, 0.0, 0.3, 0.6, -0.3, 0.0, 0.3, 0.6, 0.0, 0.3, 0.6};
    float ySam[N_POINTS] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.4, 0.6, 0.6, 0.6};
  //Load 12 sampling points elements from hw44.samplingElem.dat
  int samplingElem[12];
  file = fopen("hw44.samplingElem.dat", "r");
  for (int i=0; i < 12; i++){
    fscanf(file, "%d", &samplingElem[i]);
  }
  fclose(file);
  //Find S
  float S[N_POINTS][N+1];
  float measU[N_POINTS];
  samplingMatrix(xNodes, yNodes, elems, samplingElem, xSam, ySam, b, S, measU);
  //Find transpose of S
  float S_T[N+1][N_POINTS];
  for(int i=0; i < N_POINTS; i++){
    for(int j=0; j < N; j++){
      S_T[j][i] = S[i][j];
    }
  }

  ////Generate CovDelta and WDelta
  float CovDelta[N_POINTS][N_POINTS];
  float WDelta[N_POINTS][N_POINTS];
  float sigma2 = 0.05*0.05;

  //Due to independence and zero mean properties of noise
  for(int i=0; i < N_POINTS; i++){
    for(int j=0; j < N_POINTS; j++){
      if(i==j){
        CovDelta[i][j] = sigma2;
        WDelta[i][j] = 1.0 / sigma2;
      }
      else{
        CovDelta[i][j] = 0.0;
        WDelta[i][j] = 0.0;
      }
    }
  }

  //////////////Representor Approach
  ////Solving for delta and U in Model-Data Misfit
  float delta[N_POINTS];
  float lambda[N_POINTS];
  float temp[N+1][N_POINTS];
  float temp2[N_POINTS][N_POINTS];
  float bi[N_POINTS];
  float ui[N_POINTS];

  float URep[N_POINTS][N_POINTS];
  float RRep[N_POINTS][N_POINTS];

  //initialize delta to 0.0
  for (int i=0; i < N_POINTS; i++){
    delta[i] = 0.0;
  }

  for (int i=0; i < N_POINTS; i++){
    //Unit Misfit
    delta[i] = 1.0;

    //Find lamba - Equation 14.66 from Textbook
    //multiply S_T with WDelta
    for (int j=0; j < N; j++){
      for (int k=0; k < N_POINTS; k++){
        temp[j][k] = 2.0 * S_T[j][k];
      }
    }
    for(int j=0; j < N_POINTS; j++){
      for(int k=0; k < N_POINTS; k++){
        temp[j][k] = 2.0 * temp[j][k];
      }
    }
    matrixVectorMult(temp2, temp, delta);
    matrixVectorMult(lambda, invA_T, temp2);

    //Find bi - Equation 14.67 from Textbook
    matrixVectorMult(bi, S, lambda);
    for(int j=0; j < N_POINTS; j++){
      bi[j] = 0.5 * bi[j];
    }

    //Find ui - Equation 14.68 from Textbook
    matrixVectorMult(ui, invA, bi);

  }

}

void matrixVectorMult(float RHS[N+1], float A[N+1][N+1], float x[N+1]){
  //Multiply banded matrix A with x
  for(int i=0; i <= N; i++){    //rows
    RHS[i] = 0.0;
    for(int j=0; j <= N; j++){    //columns
      RHS[i] += A[i][j] * x[j];
    }
  }
}

void samplingMatrix(float xNodes[N+1], float yNodes[N+1], int elems[3*M+1], float samplingElem[N_POINTS], float xSam[N_POINTS], float ySam[N_POINTS], float u[N+1], float S[N_POINTS+1][N+1], float measU[N_POINTS]){
  float x0, y0;
  float a, b, c;
  float phi[3];
  //Set S to 0.0
  for(int i=0; i <= N_POINTS; i++){
    for(int j=0; j <= N; j++){
      S[i][j] = 0.0;
    }
  }
  //Find S
  float dely[3], delx[3];  //Zero Indexed
  float area = 0.0;

  //Loop over N_POINTS
  for(int k=0; k < N_POINTS; k++){
    //Find delx, dely and area
    for (int i=0; i<3; i++) {
      dely[i] = yNodes[elems[IN(samplingElem[k], (i+1)%3+1)]] - yNodes[elems[IN(samplingElem[k], (i+2)%3+1)]];
      delx[i] = xNodes[elems[IN(samplingElem[k], (i+1)%3+1)]] - xNodes[elems[IN(samplingElem[k], (i+2)%3+1)]];
    }
    for (int i=0; i<3; i++) {
      area += 0.5 * (xNodes[elems[IN(samplingElem[0], i+1)]] * dely[i]);
    }

    //Find weighing functions phi
    x0 = xSam[k];
    y0 = ySam[k];

    for(int i=0; i<3; i++){
      a = ( xNodes[elems[IN(samplingElem[k], (i+1)%3+1)]] * yNodes[elems[IN(samplingElem[k], (i+2)%3+1)]] - xNodes[elems[IN(samplingElem[k], (i+2)%3+1)]] * yNodes[elems[IN(samplingElem[k], (i+1)%3+1)]] ) / (2.0*area);
      b = dely[i]/(2.0*area);
      c = -delx[i]/(2.0*area);
      phi[i] = a + b*x0 + c*y0;
    }

    for(int jj=0; jj<3; jj++){
      S[k][elems[IN(samplingElem[k], jj)]] = phi[jj];
      measU[k] = phi[jj] * u[elems[IN(samplingElem[k], jj)]];
    }

  }
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

void covNoiseInt(float xNodes[N+1], float yNodes[N+1], float sigma2, float CovB[N+1][N+1], float corLength){

  float r;
  int I, J;
  //Set CovB to 0
  for(int i=0; i <= N; i++){
    for(int j=0; j <= N; j++){
      CovB[i][j] = 0.0;
    }
  }
  //i and j loops from 0 to N
  for(int i=1; i <= N; i++){
    I = i;
    for(int j=1; j <= N; j++){
      J = j;
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