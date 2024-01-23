#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A_ROWS 1200
#define A_COLS 6

int main() {
  FILE *file;
  double matA[A_ROWS][A_COLS];
  int i, j;

  file = fopen("hw0.dat", "r");
  if (file == NULL) {
    printf("Failed to open the file.\n");
    return 1;
  }

  // Read the matrix elements
  for (i = 0; i < A_ROWS; i++) {
    for (j = 0; j < A_COLS; j++) {
      fscanf(file, "%lf", &matA[i][j]);
    }
  }
  // Close the file
  fclose(file);

  //find max value of matA
  double maxA = matA[0][0];
  for(i=0; i<A_ROWS; i++) {
    for(j=0; j<A_COLS; j++) {
      if(matA[i][j] > maxA) {
        maxA = matA[i][j];
      }
    }
  }
  printf("Maximum of A = %e\n", maxA);

  //find min value of matA
  double minA = matA[0][0];
  for(i=0; i<A_ROWS; i++) {
    for(j=0; j<A_COLS; j++) {
      if(matA[i][j] < minA) {
        minA = matA[i][j];
      }
    }
  }
  printf("Minimum of A = %e\n", minA);

  //find mean squared of matA
  double mean_squaredA = 0;
  for(i=0; i<A_ROWS; i++) {
    for(j=0; j<A_COLS; j++) {
      mean_squaredA += matA[i][j] * matA[i][j];
    }
  }
  mean_squaredA /= (A_ROWS * A_COLS);
  printf("Mean Squared of A = %e\n", mean_squaredA);

  //transpose matA
  double matA_T[A_COLS][A_ROWS];
  for(i=0; i<A_ROWS; i++) {
    for(j=0; j<A_COLS; j++) {
      matA_T[j][i] = matA[i][j];
    }
  }

  ////find B = A^T * A
  double matB[A_COLS][A_COLS];
  for(i=0; i<A_COLS; i++) {
    for(j=0; j<A_COLS; j++) {
      matB[i][j] = 0;
      for(int k=0; k<A_ROWS; k++) {
        matB[i][j] += matA_T[i][k] * matA[k][j];
      }
    }
  }
  //write be to an ascii file named B.dat
  file = fopen("B.dat", "w");
  if (file == NULL) {
    printf("Failed to open the file.\n");
    return 1;
  }
  for(i=0; i<A_COLS; i++) {
    for(j=0; j<A_COLS; j++) {
      fprintf(file, "%e\t", matB[i][j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);

  //read B.dat
  double matB_read[A_COLS][A_COLS];
  file = fopen("B.dat", "r");
  if (file == NULL) {
    printf("Failed to open the file.\n");
    return 1;
  }
  for(i=0; i<A_COLS; i++) {
    for(j=0; j<A_COLS; j++) {
      fscanf(file, "%lf", &matB_read[i][j]);
    }
  }

  //Diagonal of B
  double diagonalB[A_COLS];
  for(i=0; i<A_COLS; i++) {
    for(j=0; j<A_COLS; j++) {
      fscanf(file, "%lf", &matB_read[i][j]);
      if(i == j) {
        diagonalB[i] = matB_read[i][j];
      }
    }
  }
  fclose(file);
  
  //write diagonal to diagonal.dat
  file = fopen("diagonalB.dat", "w");
  if (file == NULL) {
    printf("Failed to open the file.\n");
    return 1;
  }
  for(i=0; i<A_COLS; i++) {
    fprintf(file, "%e\n", diagonalB[i]);
  }
  fclose(file);

  //Symmetry check on B
  double Bsym[A_COLS][A_COLS];
  double Bsym_mean = 0;
  double Bsym_max = 0;
  double Bsym_rms = 0;

  for(i=0; i<A_COLS; i++){
    for(j=0; j<A_COLS; j++){
      Bsym[i][j] = abs(matB_read[i][j] - matB_read[j][i]);
      Bsym_mean += Bsym[i][j];
      if(Bsym[i][j] > Bsym_max) {
        Bsym_max = Bsym[i][j];
      }
      Bsym_rms += Bsym[i][j] * Bsym[i][j];
    }
  }

  Bsym_mean /= (A_COLS * A_COLS);
  Bsym_rms /= (A_COLS * A_COLS);
  Bsym_rms = sqrt(Bsym_rms);
  printf("\n----------------------------------------------------------------------------------------------------------------\n");
  printf("Mean of Symmetry Check for B = %e\n", Bsym_mean);
  printf("Max of Symmetry Check for B = %e\n", Bsym_max);
  printf("RMS of Symmetry Check for B = %e\n", Bsym_rms);


  //Symmetry Check on B -- Normalized
  double Bsym_norm[A_COLS][A_COLS];
  double Bsym_norm_mean = 0;
  double Bsym_norm_max = 0;
  double Bsym_norm_rms = 0;

  for(i=0; i<A_COLS; i++){
    for(j=0; j<A_COLS; j++){
      Bsym_norm[i][j] = abs( (matB[i][j] - matB[j][i])*2 / (matB[i][j] + matB[j][i]) );
      Bsym_norm_mean += Bsym_norm[i][j];
      if(Bsym_norm[i][j] > Bsym_norm_max) {
        Bsym_norm_max = Bsym_norm[i][j];
      }
      Bsym_norm_rms += Bsym_norm[i][j] * Bsym_norm[i][j];
    }
  }

  Bsym_norm_mean /= (A_COLS * A_COLS);
  Bsym_norm_rms /= (A_COLS * A_COLS);
  Bsym_norm_rms = sqrt(Bsym_norm_rms);
  printf("\n----------------------------------------------------------------------------------------------------------------\n");
  printf("Mean of Symmetry Check for B -- Normalized = %e\n", Bsym_norm_mean);
  printf("Max of Symmetry Check for B -- Normalized = %e\n", Bsym_norm_max);
  printf("RMS of Symmetry Check for B -- Normalized = %e\n", Bsym_norm_rms);


  //print 3rd row of B
  printf("\n----------------------------------------------------------------------------------------------------------------\n");
  printf("3rd row of B =\n");
  for(i=0; i<A_COLS; i++) {
    printf("%e\t", matB[2][i]);
  }
  printf("\n");
  //save to 3rdrowB.dat
  file = fopen("3rdrowB.dat", "w");
  if (file == NULL) {
    printf("Failed to open the file.\n");
    return 1;
  }
  for(i=0; i<A_COLS; i++) {
    fprintf(file, "%e\t", matB[2][i]);
  }
  fclose(file);

  //print 2nd subdiagonal of B
  printf("\n----------------------------------------------------------------------------------------------------------------\n");
  printf("2nd subdiagonal of B =\n");
  for(i=0; i<A_COLS-1; i++) {
    printf("%e\t", matB[i+1][i]);
  }
  printf("\n");
  //save to 2ndsubdiagonalB.dat
  file = fopen("2ndsubdiagonalB.dat", "w");
  if (file == NULL) {
    printf("Failed to open the file.\n");
    return 1;
  }
  for(i=0; i<A_COLS-1; i++) {
    fprintf(file, "%e\t", matB[i+1][i]);
  }
  fclose(file);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////find C = A * A^T
  //malloc matC
  double **matC = (double **)malloc(A_ROWS * sizeof(double *));
  if (matC == NULL) {
    perror("malloc 1 failed");
    return 1;
  }
  for(i=0; i<A_ROWS; i++) {
    matC[i] = malloc(A_ROWS * sizeof(double));
    if (matC[i] == NULL) {
      perror("malloc 2 failed");
    }
  }

  //find C = A * A^T
  for(i=0; i<A_ROWS; i++) {
    for(j=0; j<A_ROWS; j++) {
      matC[i][j] = 0;
      for(int k=0; k<A_COLS; k++) {
        matC[i][j] += matA[i][k] * matA_T[k][j];
      }
    }
  }

  //write C to C.dat
  file = fopen("C.dat", "w");
  if (file == NULL) {
    printf("Failed to open the file.\n");
    return 1;
  }
  for(i=0; i<A_ROWS; i++) {
    for(j=0; j<A_ROWS; j++) {
      fprintf(file, "%e\t", matC[i][j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);

  //read C.dat
  double **matC_read = (double **)malloc(A_ROWS * sizeof(double *));
  if (matC_read == NULL) {
    perror("malloc 1 failed");
    return 1;
  }
  for(i=0; i<A_ROWS; i++) {
    matC_read[i] = malloc(A_ROWS * sizeof(double));
    if (matC_read[i] == NULL) {
      perror("malloc 2 failed");
    }
  }
  
  file = fopen("C.dat", "r");
  if (file == NULL) {
    printf("Failed to open the file.\n");
    return 1;
  }
  for(i=0; i<A_ROWS; i++) {
    for(j=0; j<A_ROWS; j++) {
      fscanf(file, "%lf", &matC_read[i][j]);
    }
  }
  
  //Diagonal of C
  double diagonalC[A_ROWS];
  for(i=0; i<A_ROWS; i++) {
    for(j=0; j<A_ROWS; j++) {
      fscanf(file, "%lf", &matC_read[i][j]);
      if(i == j) {
        diagonalC[i] = matC_read[i][j];
      }
    }
  }
  fclose(file);

  //write diagonal to diagonalC.dat
  file = fopen("diagonalC.dat", "w");
  if (file == NULL) {
    printf("Failed to open the file.\n");
    return 1;
  }
  for(i=0; i<A_ROWS; i++) {
    fprintf(file, "%e\n", diagonalC[i]);
  }
  fclose(file);

  //Symmetry check on C
    double **Csym = (double **)malloc(A_ROWS * sizeof(double *));
  if (Csym == NULL) {
    perror("malloc 1 failed");
    return 1;
  }
  for(i=0; i<A_ROWS; i++) {
    Csym[i] = malloc(A_ROWS * sizeof(double));
    if (Csym[i] == NULL) {
      perror("malloc 2 failed");
    }
  }
  double Csym_mean = 0;
  double Csym_max = 0;
  double Csym_rms = 0;

  for(i=0; i<A_ROWS; i++){
    for(j=0; j<A_ROWS; j++){
      Csym[i][j] = abs(matC_read[i][j] - matC_read[j][i]);
      Csym_mean += Csym[i][j];
      if(Csym[i][j] > Csym_max) {
        Csym_max = Csym[i][j];
      }
      Csym_rms += Csym[i][j] * Csym[i][j];
    }
  }

  Csym_mean /= (A_ROWS * A_ROWS);
  Csym_rms /= (A_ROWS * A_ROWS);
  Csym_rms = sqrt(Csym_rms);
  printf("\n----------------------------------------------------------------------------------------------------------------\n");
  printf("Mean of Symmetry Check for C = %e\n", Csym_mean);
  printf("Max of Symmetry Check for C = %e\n", Csym_max);
  printf("RMS of Symmetry Check for C = %e\n", Csym_rms);

  //Symmetry Check on C -- Normalized
  double Csym_norm_mean = 0;
  double Csym_norm_max = 0;
  double Csym_norm_rms = 0;

  for(i=0; i<A_ROWS; i++){
    for(j=0; j<A_ROWS; j++){
      Csym[i][j] = abs( (matC[i][j] - matC[j][i])*2 / (matC[i][j] + matC[j][i]) );
      Csym_norm_mean += Csym[i][j];
      if(Csym[i][j] > Csym_norm_max) {
        Csym_norm_max = Csym[i][j];
      }
      Csym_norm_rms += Csym[i][j] * Csym[i][j];
    }
  }

  Csym_norm_mean /= (A_ROWS * A_ROWS);
  Csym_norm_rms /= (A_ROWS * A_ROWS);
  Csym_norm_rms = sqrt(Csym_norm_rms);

  printf("\n----------------------------------------------------------------------------------------------------------------\n");
  printf("Mean of Symmetry Check for C -- Normalized = %e\n", Csym_norm_mean);
  printf("Max of Symmetry Check for C -- Normalized = %e\n", Csym_norm_max);
  printf("RMS of Symmetry Check for C -- Normalized = %e\n", Csym_norm_rms);

  //print 3rd row of C
  printf("\n----------------------------------------------------------------------------------------------------------------\n");
  printf("3rd row of C =\n");
  for(i=0; i<A_ROWS; i++) {
    printf("%e\t", matC[2][i]);
  }
  printf("\n");
  //save to 3rdrowC.dat
  file = fopen("3rdrowC.dat", "w");
  if (file == NULL) {
    printf("Failed to open the file.\n");
    return 1;
  }
  for(i=0; i<A_ROWS; i++) {
    fprintf(file, "%e\t", matC[2][i]);
  }
  fclose(file);

  //print 2nd subdiagonal of C
  printf("\n----------------------------------------------------------------------------------------------------------------\n");
  printf("2nd subdiagonal of C =\n");
  for(i=0; i<A_ROWS-1; i++) {
    printf("%e\t", matC[i+1][i]);
  }
  printf("\n");
  //save to 2ndsubdiagonalC.dat
  file = fopen("2ndsubdiagonalC.dat", "w");
  if (file == NULL) {
    printf("Failed to open the file.\n");
    return 1;
  }
  for(i=0; i<A_ROWS-1; i++) {
    fprintf(file, "%e\t", matC[i+1][i]);
  }

  //free matC
  for(i=0; i<A_ROWS; i++) {
    free(matC[i]);
  }
  free(matC);
  //free Csym
  for(i=0; i<A_ROWS; i++) {
    free(Csym[i]);
  }
  free(Csym);
  //free matC_read
  for(i=0; i<A_ROWS; i++) {
    free(matC_read[i]);
  }
  free(matC_read);

  return 0;
}
