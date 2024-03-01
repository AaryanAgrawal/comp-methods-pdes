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
#define HALFBW 2
#include "solve1.h"

int test(void){
  float xNodes[N+1];
  float yNodes[N+1];
  int elems[M*3+1];

  //Load Input Mesh
  scanNodes(xNodes, yNodes);
  scanElems(elems);

  //Element Asssembly
  float Ae[3][3];
  float be[3];
  float b[N+1];
  float A[N+1][N+1];
  //Initialize A and b
  for (int i=0; i <= N; i++){
    b[i] = 0.0;
    for (int j=0; j <= N; j++){
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
        A[II][JJ] = A[II][JJ] + Ae[I-1][J-1];
        //printf("JB = %d, A[%d][%d] = %f\n", JB, II, JB, A[II][JB]);
      }
    }
  }


  //Boundary Conditions Type I
  //Left Edge
  float Length;
  Length = yNodes[5] - yNodes[1];
  for (int i=1; i<=5; i++){
    b[i] = yNodes[i] / Length;
    for (int j=1; j<=N; j++){
      if(j==i){
        A[i][j] = 1.0;
      }
      else{
        A[i][j] = 0.0;
      }
    }
  }


  //Upper Edge
  int nodesBUpper[6] = {5, 10, 15, 20, 24, 27};
  for (int i=0; i<6; i++){
    int I = nodesBUpper[i];
    b[I] = 1.0;
    for (int j=0; j<=N; j++){
      if(j==I){
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
    for (int j=0; j<=N; j++){
      if(j==I){
        A[I][j] = 1.0;
      }
      else{
        A[I][j] = 0.0;
      }
    }
  }


  //Print A to file
  printf("Printing A to file\n");
  FILE *fileA;
  fileA = fopen("Atest.dat", "w");
  for (int i=0; i <= N; i++){
    for (int j=0; j <= N; j++){
      fprintf(fileA, "%f\t", A[i][j]);
    }
    fprintf(fileA, "\n");
  }
  fclose(fileA);
  //Print b to file
  printf("Printing b to file\n");
  FILE *file;
  file = fopen("btest.dat", "w");
  for (int i=0; i <= N; i++){
    fprintf(file, "%f\n", b[i]);
  }
  fclose(file);

  return 0;
}
