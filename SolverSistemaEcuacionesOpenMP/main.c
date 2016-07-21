#include <stdio.h>
#include <stdlib.h>
#include "Matrix.h"
#include "MyRand.h"
#include "MonteCarlo.h"
int main()
{

    srand(time(NULL));
    int N=400;

    Matrix *A = NewMatrix(N,N);
    Vector *b = NewVector(N);
    Vector *X = NewVector(N);
    double Inf=1,Sup=1e3;

    omp_set_num_threads(4);
    RandomMatrixDiagonalDominante(A,b, Inf, Sup);


    double ParametroRelajacion=1;
    double Tolerancia = 1e-3;
    MontecarloAlmostOptimal(A,b,X, ParametroRelajacion, Tolerancia);
   /* printf("Matrix A:\n");
    PrintMatrix(A);
    printf("\nVector b:\n");
    PrintVector(b);
    printf("\nSolución X:\n");
    PrintVector(X);
*/
    FreeMatrix(A);
    FreeVector(b);
    FreeVector(X);
    return 0;
}
