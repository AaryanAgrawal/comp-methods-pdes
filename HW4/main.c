/*
Author: Aaryan Agrawal
ENGG 105: Comp Methods to PDEs I
Instructor: Prof. Paulsen

Solve Simplified Flow Problem using Finite Element Method
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 28    //# of Nodes
#define M 38    //# of Elements
#define HALFBW 6
#include "solve1.h"

void scanNodes(float *xNodes, float *yNodes);
void scanElems(int *elems);
int IN(int i, int j);
void ElementMatrix(int L, int *elems, float *xNodes, float *yNodes, float Ae[3][3], float be[3]);
void boundaryConditionsI(float xNodes[N], float yNodes[N], float A[N+1][2*HALFBW+2], float b[N+1]);
void gradient(int *elems, float *xNodes, float *yNodes, float *u, float grad[M+1][2]);

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
    for (int I=1; I <= 3; I++) {
      II = elems[IN(L, I)];
      //printf("II = %d, L = %d, I = %d\n", II, L, I);
      b[II] = b[II] + be[I-1];
      for (int J=1; J <= 3; J++) {
        JJ = elems[IN(L, J)];
        printf("II = %d, JJ = %d, L = %d, I = %d, J = %d\n", II, JJ, L, I, J);
        JB = (HALFBW + 1) + JJ - II;
        A[II][JB] = A[II][JB] + Ae[I-1][J-1];
        //printf("JB = %d, A[%d][%d] = %f\n", JB, II, JB, A[II][JB]);
      }
    }
  }
  //Boundary Conditions
  boundaryConditionsI(xNodes, yNodes, A, b);

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
  file = fopen("grad_III.dat", "w");
  for (int i=1; i <= M; i++){
    fprintf(file, "%f\t%f\n", grad[i][0], grad[i][1]);
  }
  fclose(file);

//////////////////////////////////////////////////////////////////////////////////////

  //Print b as u to file
  printf("Printing u to file\n");
  file = fopen("u_III.dat", "w");
  for (int i=0; i <= N; i++){
    fprintf(file, "%f\n", b[i]);
  }
  fclose(file);

  return 0;
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
      dudx = -dely[i]/(2*area);
      dudy = -delx[i]/(2*area);
      grad[L][0] += u[elems[IN(L, i+1)]] * dudy;
      printf("grad[%d][0] = %f, area=%f, node j=%d \n", L, grad[L][0], area, elems[IN(L, i+1)]);
      grad[L][1] += u[elems[IN(L, i+1)]] * dudx;
    }
  }
}

void boundaryConditionsI(float xNodes[N], float yNodes[N], float A[N+1][2*HALFBW+2], float b[N+1]){
  //Upper Edge
  int nodesBUpper[] = {5, 10, 15, 20, 24, 27, 28, 25, 21};
  for (int i=0; i<9; i++){
    int I = nodesBUpper[i];
    b[I] = 1.0;
    for (int j=0; j<=2*HALFBW+1; j++){
      if(j==HALFBW+1){
        A[I][j] = 1.0;
      }
      else{
        A[I][j] = 0.0;
      }
    }
  }

  //Lower Edge
  int nodesBLower[4] = {1, 6, 11, 16};
  for (int i=0; i<4; i++){
    int I = nodesBLower[i];
    b[I] = 0.0;
    for (int j=0; j<=2*HALFBW+1; j++){
      if(j==HALFBW+1){
        A[I][j] = 1.0;
      }
      else{
        A[I][j] = 0.0;
      }
    }
  }

  //Boundary Conditions Type I
  //Left Edge
  float Length;
  Length = yNodes[5] - yNodes[1];
  for (int i=1; i<=5; i++){
    b[i] = yNodes[i] / Length;
    for (int j=1; j<=2*HALFBW+1; j++){
      if(j==HALFBW+1){
        A[i][j] = 1.0;
      }
      else{
        A[i][j] = 0.0;
      }
    }
  }
}

void scanNodes(float *xNodes, float *yNodes) {
  FILE *file;
  file = fopen("hw5.nodes.dat", "r");
  if (file == NULL){
    printf("Error: file not found\n");
    exit(1);
  }
  int temp;
  for(int i = 1; i <= 28; i++){
    fscanf(file, "%d %f %f", &temp, &xNodes[i], &yNodes[i]);
  }
  fclose(file);
}

void scanElems(int *elems) {
  FILE *file;
  file = fopen("hw5.elems.dat", "r");
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
  for (int i=0; i<3; i++) {
    dely[i] = yNodes[elems[IN(L, (i+1)%3+1)]] - yNodes[elems[IN(L, (i+2)%3+1)]];
    delx[i] = xNodes[elems[IN(L, (i+1)%3+1)]] - xNodes[elems[IN(L, (i+2)%3+1)]];
  }
  for (int i=0; i<3; i++) {
    area += 0.5 * (xNodes[elems[IN(L, i+1)]] * dely[i]);
    //printf("xNodes = %f, elems = %d, L = %d, i+1 = %d, dely = %f\n", xNodes[elems[IN(L, i+1)]], elems[IN(L, i+1)], L, i+1, dely[i]);
  }
  //printf("area = %f\n", area);
  
  //Ae Matrix
  for (int i=0; i < 3; i++){
    for (int j=0; j < 3; j++){
      if (i == j) {
        Ae[i][j] = (- dely[i]*dely[i] - delx[i]*delx[i]) / (4 * area);
      }
      else {
        Ae[i][j] = (- dely[i]*dely[j] - delx[i]*delx[j]) / (4 * area);
      }
      printf("Ae[%d][%d] = %f\t", i, j, Ae[i][j]);
    }
    printf("\n");
  }

  //be Vector w/ Type II Boundary Conditions
  for (int j=0; j<3; j++){
    be[j] = 0.0;
  }

}