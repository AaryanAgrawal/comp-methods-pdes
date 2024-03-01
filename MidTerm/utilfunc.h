#ifndef UTILFUNC_H
#define UTILFUNC_H

void applyMoleculeA(float U0[], float U1[], float U2[], int N, float r);
void applyMoleculeB(float U0[], float U1[], float U2[], int N, float r);
void printSolutionToFile(float U[], int N, FILE *file);
void shiftTimeLevels(float U0[], float U1[], float U2[], int N);
void boundaryConditions(float U0[], float U1[], int N);
void initialConditions(float U0[], float U1[], int N, float x[], float x0_param, float sigma_param);
int vector_max(float A[], int n);

#endif // UTILFUNC_H