/*
Author: Aaryan Agrawal
ENGG 105: Comp Methods to PDEs I
Instructor: Prof. Paulsen

Solve Bioheat Equation over cross-section of patient using Finite Element Method

1. Load Input Mesh
2. Load Material Properties
2. Element Assembly
  1. Element Matrix Ae - Determine if element is triangle or quadrilateral
  2. Forcing Vector be - Using heat generation file
3. Global Assembly
4. Boundary Conditions
5. Solve


1. Add boundary conditions to both matrices
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "solve1.h"
#define N 969
#define M 1089
#define HALFBW 48
#define GAUSS 0.577350269


void applyBCI(int nodeBC[N], float A[N+1][2*HALFBW+1+1], float b[N+1]);
void matrixMult(float RHS[N+1], float A[N+1][2*HALFBW+1+1], float x[N+1], float b[N+1]);
void applyBC(int nodeBC[N], int neighNodeBC[N][2], float coeffBC[N][2], float A[N+1][2*HALFBW+1+1], float AA[N+1][2*HALFBW+1+1], float b[N+1], float xNodes[N+1], float yNodes[N+1], int mat[M+1], float matK[9], float theta, float antiTheta, float delt);
void scanNodes(float *xNodes, float *yNodes);
void scanElemHeat(float elemheat[M+1]);
void scanElems(int *elems, int *mat);
void inputMatProp(float matK[9], float matC[9], float matM[9]);
int IN(int i, int j);
void quadElementMatrix(float *xL, float *yL, float Ae[4][4], float be[4], float Km, float Fm, float Gm, float Me[4][4], float Cm);
void elementBasis(float   Z, float E, float P[4], float *DJ, float *xL, float *yL, float *dPX, float *dPY);
int checkShape(int L, int *elems);
void getElemCoords(int L, int *elems, float *xNodes, float *yNodes, float xL[4], float yL[4]);
float getGauss(int m, int i);
void triElementMatrix(int L, int *elems, float *xNodes, float *yNodes, float Ae[4][4], float be[4], float Km, float Fm, float Gm, float Me[4][4], float Cm);
float sideLength(float x1, float y1, float x2, float y2);

int main(void){

  float xNodes[N+1];
  float yNodes[N+1];
  int elems[M*4+1];
  int mat[M+1];

  // Load Input Mesh
  scanNodes(xNodes, yNodes);  // x, y coordinates for each node
  scanElems(elems, mat);  //elems - incidence list, mat - material list for each element

  // Load Material Properties
  float matK[8+1];  //Material Parameter K for each material # 3 -> 8
  float matM[8+1];  //Material Parameter M
  float matC[8+1];  //Material Parameter C
  //initialize matK and matM
  for(int i=0; i<=8; i++){
    matK[i] = 0.0;
    matM[i] = 0.0;
    matC[i] = 0.0;
  }
  inputMatProp(matK, matC, matM);
  // printf("matK\n");
  // for(int i=0; i<=8; i++){
  //   printf("%d\n", matK[i]);
  // }
  // printf("matM\n");
  // for(int i=0; i<=8; i++){
  //   printf("%d\n", matM[i]);
  // }


////////////////////////////////////////////////////////////////////////////////////

  //Load Heat Generation
  float elemHeat[M+1];
  scanElemHeat(elemHeat);
  
  //Load Boundary Conditions
  int nodeBC[N];
  int neighNodeBC[N][2];
  float coeffBC[N][2];
  //initialize nodeBC, neighNodeBC, coeffBC
  for(int i=0; i<N; i++){
    nodeBC[i] = 0;
    for(int j=0; j<2; j++){
      neighNodeBC[i][j] = 0;
      coeffBC[i][j] = 0.0;
    }
  }
  scanBC(nodeBC, neighNodeBC, coeffBC);

////////////////////////////////////////////////////////////////////////////////////

//Time transient
  float delt = 5.0;
  float theta = 0.7;
  float antiTheta = 1 - theta;

////////////////////////////////////////////////////////////////////////////////////

  float Ae[4][4];
  float Me[4][4];
  float be[4];
  float b[N+1];
  float A[N+1][2*HALFBW+1+1];
  float AA[N+1][2*HALFBW+1+1];
  //Initialize A and b
  for (int i=0; i <= N; i++){
    b[i] = 0.0;
    for (int j=0; j <= 2*HALFBW+1; j++){
      A[i][j] = 0.0;
      AA[i][j] = 0.0;
    }
  }

  float xL[4];
  float yL[4];
  int shape;
  float Km, Fm, Gm, Cm;
  int II, JJ, JB;
  int matNum;

  // For each element
  for(int L=1; L <= M; L++){

    for(int i=0; i<4; i++){
      be[i] = 0.0;
      for(int j=0; j<4; j++){
        Ae[i][j] = 0.0;
        Me[i][j] = 0.0;
      }
    }

    getElemCoords(L, elems, xNodes, yNodes, xL, yL);

    matNum = mat[L];
    Km = matK[matNum];
    Fm = - matM[matNum];
    Gm = - elemHeat[L];
    Cm = - matC[matNum];
    //printf("Km = %f, Fm = %f, Gm = %f, Cm = %f\n", Km, Fm, Gm, Cm);

    shape = checkShape(L, elems);
    
    if(shape==1){      //Quadrilateral
      quadElementMatrix(xL, yL, Ae, be, Km, Fm, Gm, Me, Cm);
    }
    else{      //Triangle 
      triElementMatrix(L, elems, xNodes, yNodes, Ae, be, Km, Fm, Gm, Me, Cm);
    }

    //print Ae and Me
    if(shape==3){
      printf("Ae\n");
      for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
          printf("%f\t", Ae[i][j]);
        }
        printf("\n");
      }
      printf("Me\n");
      for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
          printf("%f\t", Me[i][j]);
        }
        printf("\n");
      }
    }

    //Global Assembly
    for (int I=1; I <= 4; I++) {
      II = elems[IN(L, I)];
      //printf("II = %d, L = %d, I = %d\n", II, L, I);
      b[II] = b[II] + delt*be[I-1];
      for (int J=1; J <= 4; J++) {
        JJ = elems[IN(L, J)];
        //printf("II = %d, JJ = %d, L = %d, I = %d, J = %d\n", II, JJ, L, I, J);
        JB = (HALFBW + 1) + JJ - II;
        A[II][JB] = A[II][JB] - antiTheta*delt*Ae[I-1][J-1] + Me[I-1][J-1];
        AA[II][JB] = AA[II][JB] + theta*delt*Ae[I-1][J-1] + Me[I-1][J-1];

        //printf("JB = %d, A[%d][%d] = %f\n", JB, II, JB, A[II][JB]);
      }
    }
  }

////////////////////////////////////////////////////////////////////////////////////

  //Boundary Conditions
  applyBC(nodeBC, neighNodeBC, coeffBC, A, AA, b, xNodes, yNodes, mat, matK, theta, antiTheta, delt);

////////////////////////////////////////////////////////////////////////////////////

  FILE *fp;

  //Initial Conditions
  float x[N+1];
  for (int i=0; i <= N; i++){
    x[i] = 0.0;
  }

  //k --> k+1 Equation RHS [B*u + R]
  float RHS[N+1];
  matrixMult(RHS, A, x, b);
  //applyBCI(nodeBC, AA, RHS);


  ////Start time stepping
  //Copy RHS
  float tempRHS[N+1];
  for (int i=0; i <= N; i++){
    tempRHS[i] = RHS[i];
  }
  //Copy A
  float tempA[N+1][2*HALFBW+1+1];
  for (int i=0; i <= N; i++){
    for (int j=0; j <= 2*HALFBW+1; j++){
      tempA[i][j] = A[i][j];
    }
  }

  solve1(1, AA, tempRHS, N, HALFBW, 2*HALFBW+2);  //AA destroyed for LU Decomp

  FILE *ftime;
  ftime = fopen("transient.dat", "w");
  if(ftime == NULL){
    printf("Error opening transient.dat\n");
    exit(1);
  }

  //For each time step
  for(int tCnt = 0; tCnt < 2500; tCnt++){
  
    for (int i=0; i <= N; i++){
      fprintf(ftime, "%f\t", x[i]);
    }
    fprintf(ftime, "\n");

    matrixMult(RHS, A, x, b);
    //applyBCI(nodeBC, tempA, RHS);

    //Advance time step
    solve1(2, AA, RHS, N, HALFBW, 2*HALFBW+2);

    for (int i=0; i <= N; i++){
      x[i] = RHS[i];
    }
  }

  fclose(ftime);
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

void applyBCI(int nodeBC[N], float A[N+1][2*HALFBW+1+1], float b[N+1]){
  //Apply Type I BC
  int i = 0;
  while(nodeBC[i] != 0){
    A[nodeBC[i]][HALFBW + 1] = 1.0;
    for(int j=0; j <= 2*HALFBW+1; j++){
      if(j != HALFBW + 1){
        A[nodeBC[i]][j] = 0.0;
      }
    }
    b[nodeBC[i]] = 0.0;
    i++;
  }
}

void applyBC(int nodeBC[N], int neighNodeBC[N][2], float coeffBC[N][2], float A[N+1][2*HALFBW+1+1], float AA[N+1][2*HALFBW+1+1], float b[N+1], float xNodes[N+1], float yNodes[N+1], int mat[M+1], float matK[9], float theta, float antiTheta, float delt){

  //Scan file with boundary element numbers (created using MATLAB)
  int bElem[2*M];
  int bElemCnt=0;
  FILE *file;
  file = fopen("bElem.dat", "r");
  if(file==NULL){
    printf("Error opening bElem.dat\n");
    exit(1);
  }
  while(fscanf(file, "%d, %d", &bElem[2*bElemCnt], &bElem[2*bElemCnt+1]) != EOF){
    bElemCnt++;
  }
  fclose(file);

  // for(int i=0; i<bElemCnt; i++){
  //   printf("%d %d\n", bElem[i*2], bElem[i*2+1]);
  // }

  float aI, aII, cI, cII, LI, LII, KI, KII;
  int B, C, D;
  //Apply Type III BC
  int i = 0;
  float h, Ta;
  h = coeffBC[i][0];
  Ta = coeffBC[i][1];
  //printf("h = %f, Ta = %f\n", h, Ta);

  while(nodeBC[i] != 0){

    B = nodeBC[i];
    C = neighNodeBC[i][0];
    D = neighNodeBC[i][1];

    LI = sideLength(xNodes[B], yNodes[B], xNodes[C], yNodes[C]);
    LII = sideLength(xNodes[B], yNodes[B], xNodes[D], yNodes[D]);
    //printf("LI = %f, LII = %f\n", LI, LII);

    KI = matK[mat[bElem[i*2]]];
    //printf("KI = %f\n", KI);
    KII = matK[mat[bElem[i*2+1]]];
    //printf("KII = %f\n", KII);
    aI = -h/KI;
    aII = -h/KII;
    cI = h * Ta / KI;
    cII = h * Ta / KII;

    A[B][(HALFBW + 1) + C - B] += - antiTheta * delt * (aI * LI / 6);
    A[B][HALFBW + 1] += - antiTheta * delt * (LI*aI + LII*aII) / 3;
    A[B][(HALFBW + 1) + D - B] += - antiTheta * delt * (aII * LII / 6);

    AA[B][(HALFBW + 1) + C - B] += theta * delt * (aI * LI / 6);
    AA[B][HALFBW + 1] += theta * delt * (LI*aI + LII*aII) / 3;
    AA[B][(HALFBW + 1) + D - B] += theta * delt * (aII * LII / 6);

    b[B] += delt * ( - (cI * LI / 2) - (cII * LII / 2) );

    i++;
  }
}

float sideLength(float x1, float y1, float x2, float y2){
  return sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}

void scanBC(int nodeBC[N], int neighNodeBC[N][2], float coeffBC[N][2]){
  FILE *file;
  file = fopen("bpeltr4.dat", "r");
  if(file == NULL){
    printf("Error: file not found\n");
    exit(1);
  }
  int temp;
  int i = 0;
  while( fscanf(file, "%d,%d,%d,%d,%d, %f,%f", &temp, &nodeBC[i], &temp, &neighNodeBC[i][0], &neighNodeBC[i][1], &coeffBC[i][0], &coeffBC[i][1]) != EOF){
    i++;
  }
  fclose(file);
}

void scanElemHeat(float elemheat[M+1]){
  FILE *file;
  file = fopen("ppelt4.dat", "r");
  if(file == NULL){
    printf("Error: file not found\n");
    exit(1);
  }
  int temp;
  for(int i=1; i <= M; i++){
    fscanf(file, "%d,  %d,    %f", &temp, &temp, &elemheat[i]);
  }
  fclose(file);
}

void triElementMatrix(int L, int *elems, float *xNodes, float *yNodes, float Ae[4][4], float be[4], float Km, float Fm, float Gm, float Me[4][4], float Cm){
  float dely[3], delx[3];  //Zero Indexed
  float area = 0.0;
  for (int i=0; i<3; i++) {
    dely[i] = yNodes[elems[IN(L, (i+1)%3+1)]] - yNodes[elems[IN(L, (i+2)%3+1)]];
    delx[i] = xNodes[elems[IN(L, (i+1)%3+1)]] - xNodes[elems[IN(L, (i+2)%3+1)]];
  }
  for (int i=0; i<3; i++) {
    area += 0.5 * (xNodes[elems[IN(L, i+1)]] * dely[i]);
    //printf("xNodes = %f, elems = %d, L = %d, i+1 = %d, dely = %f\n", xNodes[elems[IN(L, i+1)]], elems[IN(L, i+1)], L, i+1, dely[i]);
  }
  //printf("area = %f\n", area);
  //printf("Km = %f, Fm = %f, Gm = %f\n", Km, Fm, Gm);
  //Ae Matrix
  for (int i=0; i < 3; i++){
    for (int j=0; j < 3; j++){
      if (i == j) {
        Ae[i][j] = Km*(- dely[i]*dely[i] - delx[i]*delx[i]) / (4 * area) + (Fm*area / 6);
        Me[i][j] = Cm*area / 6;
      }
      else {
        Ae[i][j] = Km*(- dely[i]*dely[j] - delx[i]*delx[j]) / (4 * area) + (Fm*area / 12);
        Me[i][j] = Cm*area / 12;
      }
      //printf("Ae[%d][%d] = %f\t", i, j, Ae[i][j]);
    }
    //printf("\n");
  }

  //be Vector w/ Type II Boundary Conditions
  for (int j=0; j<3; j++){
    be[j] = Gm*area / 3;
  }

}\

void quadElementMatrix(float *xL, float *yL, float Ae[4][4], float be[4], float Km, float Fm, float Gm, float Me[4][4], float Cm){
  float Z, E;
  float P[4];
  float DJ;
  float dPX[4];
  float dPY[4];
  
  for(int m=1; m<=4; m++){
    Z = getGauss(m, 1);
    E = getGauss(m, 2);

    elementBasis(Z, E, P, &DJ, xL, yL, dPX, dPY);

    //Gauss Point Matrix                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    for(int i=0; i<4; i++){
      for(int j=0; j<4; j++){
        Ae[i][j] = Ae[i][j] + ( -Km*(dPX[i]*dPX[j] + dPY[i]*dPY[j]) + (Fm*P[i]*P[j]) ) * DJ;
        Me[i][j] = Me[i][j] + (Cm*P[i]*P[j]) * DJ;
        //printf("Add = %f\n", (-Km*(dPX[i]*dPX[j] + dPY[i]*dPY[j]) + (Fm*P[i]*P[j]) ) * DJ);
      }
      be[i] = be[i] + DJ*Gm*P[i];
    }
  }
  // printf("Ae\n");
  // for(int i=0; i<4; i++){
  //   for(int j=0; j<4; j++){
  //     printf("%f\t", Ae[i][j]);
  //   }
  //   printf("\n");
  // }
  // printf("be\n");
  // for(int i=0; i<4; i++){
  //   printf("%f\n", be[i]);
  // }

}

void elementBasis(float Z, float E, float P[4], float *DJ, float *xL, float *yL, float *dPX, float *dPY){

  float dPZ[4];
  float dPE[4];
  float dXZ = 0.0;
  float dXE = 0.0;
  float dYZ = 0.0;
  float dYE = 0.0;

  //Bases
  P[0] = 0.25*(1-Z)*(1-E);
  P[1] = 0.25*(1+Z)*(1-E);
  P[2] = 0.25*(1+Z)*(1+E);
  P[3] = 0.25*(1-Z)*(1+E);
  dPZ[0] = -0.25*(1-E);
  dPZ[1] = 0.25*(1-E);
  dPZ[2] = 0.25*(1+E);
  dPZ[3] = -0.25*(1+E);
  dPE[0] = -0.25*(1-Z);
  dPE[1] = -0.25*(1+Z);
  dPE[2] = 0.25*(1+Z);
  dPE[3] = 0.25*(1-Z);

  //Jacobian
  for(int i=0; i<4; i++){
    dXZ = dXZ + dPZ[i]*xL[i];
    dXE = dXE + dPE[i]*xL[i];
    dYZ = dYZ + dPZ[i]*yL[i];
    dYE = dYE + dPE[i]*yL[i];
  }
  *DJ = dXZ * dYE - dXE * dYZ;

  //Derivatives
  for(int i=0; i<4; i++){
    dPX[i] = (1/ *DJ) * (dYE * dPZ[i] - dYZ * dPE[i]);
    dPY[i] = (1/ *DJ) * (dXZ * dPE[i] - dXE * dPZ[i]);
  }
  // printf("DJ = %f\n", *DJ);
  // printf("dPX = %f, %f, %f, %f\n", dPX[0], dPX[1], dPX[2], dPX[3]);
  // printf("dPY = %f, %f, %f, %f\n", dPY[0], dPY[1], dPY[2], dPY[3]);

}

int checkShape(int L, int *elems){
  if(elems[IN(L, 4)] == elems[IN(L, 3)]){
    return 2; //Triangle
  }
  else{
    return 1; //Quadrilateral
  }
}

void getElemCoords(int L, int *elems, float *xNodes, float *yNodes, float xL[4], float yL[4]){
  for(int i = 1; i <= 4; i++){
    //printf("i=%d, elems[IN(L, i)] = %d\n", i, elems[IN(L, i)]);
    xL[i-1] = xNodes[elems[IN(L, i)]];
    yL[i-1] = yNodes[elems[IN(L, i)]];
  }
}

float getGauss(int m, int i){
  switch(m){
    case 1:
      if(i==1){
        return - GAUSS;
      }
      else{
        return - GAUSS;
      }
      break;
    case 2:
      if(i==1){
        return GAUSS;
      }
      else{
        return - GAUSS;
      }
      break;
    case 3:
      if(i==1){
        return GAUSS;
      }
      else{
        return GAUSS;
      }
      break;
    case 4:
      if(i==1){
        return - GAUSS;
      }
      else{
        return GAUSS;
      }
      break;
    default:
      printf("Error: Invalid Gauss Point\n");
      exit(1);
  }
}

void inputMatProp(float matK[9], float matC[9], float matM[9]){
  FILE *file;
  file = fopen("matprop.dat", "r");
  if (file == NULL){
    printf("Error: file not found\n");
    exit(1);
  }
  int temp;
  for(int i = 3; i <= 8; i++){
    fscanf(file, "%d, %f, %f, %f", &temp, &matK[i], &matC[i], &matM[i]);
  }
  fclose(file);
}

void scanElems(int *elems, int *mat) {
  FILE *file;
  file = fopen("epeltr4.dat", "r");
  if (file == NULL){
    printf("Error: file not found\n");
    exit(1);
  }
  int temp;
  for (int i=1; i <= M; i++){
    fscanf(file, "%d,%d,%d,%d,%d,%d", &temp, &elems[IN(i, 1)], &elems[IN(i, 2)], &elems[IN(i, 3)], &elems[(IN(i, 4))], &mat[i]);
  }
  fclose(file);
}

int IN(int i, int j){
  // elems[IN(i, j)] returns Global Node #
  return (i-1)*4 + j;
}

void scanNodes(float *xNodes, float *yNodes) {
  FILE *file;
  file = fopen("npeltr4.dat", "r");
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