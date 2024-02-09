#ifndef UTILFUNC_H
#define UTILFUNC_H

void applyMoleculeA(double U0[], double U1[], double U2[], int N, double r);
void applyMoleculeB(double U0[], double U1[], double U2[], int N, double r);
void printSolutionToFile(double U[], int N, FILE *file);
void shiftTimeLevels(double U0[], double U1[], double U2[], int N);
void boundaryConditions(double U0[], double U1[], int N);
void initialConditions(double U0[], double U1[], int N, double x[], float x0_param, float sigma_param);
int vector_max(double A[], int n);

#endif // UTILFUNC_H