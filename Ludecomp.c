#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LUdecomp.h"

double **createMatrix(int N) {
  double **M = (double **) malloc(N*sizeof(double*));
  for (int i = 0; i < N; i++)
    M[i] = (double*) malloc(N*sizeof(double));
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      M[i][j] = (i == j) ? 1.0 : 0.0;
    return M;
  }

void LUdestroy(LUdecomp* LU){
    double **M = LU->LU;
    int N = LU->N;
    for (int i = 0; i < N; i++)
      free(M[i]);
    free(M);
    free(LU->mutate);
    free(LU);
  }

LUdecomp *LUdecompose(int N, const double **A) {
  LUdecomp *LU = (LUdecomp*) malloc(sizeof(LUdecomp));
  LU->N = N;
  LU->LU = createMatrix(N);
  LU->mutate = (short *) malloc(N*sizeof(short));
  LU->d = +1;

  // Clone A into LU
  double **A_ = LU->LU;
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      A_[i][j] = A[i][j];

    for (int i = 0; i < N; i++)
      LU->mutate[i] = (short) i;

  //
  //Partial pevoting
  //

  // compute aij <--- Bij on and above diagonal
    for(int j=0; j < N; j++){
      for (int i=0; i <= j; i++){
        double sum =0.0;
        for(int k = 0; k < i; k++) sum += A_[i][k] * A_[k][j];
        A_[i][j] = A_[i][j] - sum; 
      }

      // p:initial pivot value, n: initial pivot row
      double  p = fabs(A_[j][j]);
      int n = j;

      // compute aij <--- Aij on and bellow diagonal
      for(int i=j+1; i < N; i++){
        double sum = 0.0;
        for (int k = 0; k < j; k++) sum += A_[i][k] * A_[k][j];
        
        A_[i][j] = A_[i][j] - sum;

        //if better pivot found, then record the new pivot
        if(fabs(A_[i][j]) > p){
          p = fabs(A_[i][j]);
          n = i;
        }

        if(p == 0){
          return(0);
        }
      }

      if(n != j){
        short swap = LU->mutate[j];
        LU->mutate[j] = LU->mutate[n];
        LU->mutate[n] = swap;

        double *temp = A_[j];
        A_[j] = A_[n];
        A_[n] = temp;

        LU->d = -(LU->d);
      }
      
      // perform divisoin bellow the diagonal
      for(int i= j+1; i < N; i++)
        A_[i][j] = A_[i][j] / A_[j][j];
    }
    return LU;
}

void LUsolve(LUdecomp *decomp, const double *b, double *x){
  int N = decomp->N;
  // forward subatitution
  double y[N];
  y[0] = b[(decomp->mutate)[0]];
  for(int i=1; i < N; i++){
    double sum = 0.0;
    for(int j=0; j < i; j++) sum+= (decomp->LU)[i][j] * y[j];
    y[i] = b[(decomp->mutate)[i]] - sum;
  }

  // backward solving 
  x[N-1] = y[N-1]/ (decomp->LU)[N-1][N-1];
  for(int i = N-2; i >= 0; i--){
    double sum = 0.0;
    for(int j = i+1; j < N; j++) sum+= (decomp->LU)[i][j] * x[j];
    x[i] = (y[i] - sum)/(decomp->LU)[i][i];
  }
}

