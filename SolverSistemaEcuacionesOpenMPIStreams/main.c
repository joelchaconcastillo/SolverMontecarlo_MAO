/**
    Autor: Joel Chacon castillo
    Fecha: 02/Junio/2016
    Descripcion: Este fichero realiza la llamada
                 de la función realiza el proceso de
                 montecarlo.
**/
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "Matrix.h"
#include "MyRand.h"
#include "MonteCarlo.h"
int main(int argc, char **argv)
{

    int Rank, NumeroProcesos;
    /***Inicializar el ambiente de MPI*/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &NumeroProcesos);
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);


    srand(time(NULL));
  int N=350;

if(Rank==ROOT)
{
    Matrix *A = NewMatrix(N,N);
    Vector *b = NewVector(N);
    Vector *X = NewVector(N);
    double Inf=1,Sup=1e3;

//    omp_set_num_threads(4);
    /***
        Genera de forma aleatoria una matriz que diagonalmente dominante.
    **/
    RandomMatrixDiagonalDominante(A,b, Inf, Sup);


    double ParametroRelajacion=1;
    double Tolerancia = 1e-3;


    /**
        El nodo
    **/
    MontecarloAlmostOptimal(A,b,X, ParametroRelajacion, Tolerancia);
    /*printf("Matrix A:\n");
    PrintMatrix(A);
    printf("\nVector b:\n");
    PrintVector(b);
    printf("\nSolución X:\n");
    PrintVector(X);*/

    FreeMatrix(A);
    FreeVector(b);
    FreeVector(X);
}else
{

    Matrix *L  = (Matrix *) malloc(sizeof(Matrix));
    Vector *LSum = (Vector *) malloc(sizeof(Vector));
    Vector *F = (Vector *) malloc(sizeof(Vector));
    Vector *Streams = (Vector *) malloc(sizeof(Vector));
    Matrix *P = (Matrix *) malloc(sizeof(Matrix));
    double Tolerancia;
    int packsize, SizeStream;
    /***

      Recibir la información....

    ***/
      DesenpaquetarMatrix(L, LSum, F, P, &Tolerancia,Streams, &SizeStream);
      Vector *X = NewVector(LSum->Size);

      int IdThread = Rank-1;

      /**
        Implementar la herramienta de HybridTaus
        como generador de números aleatorios..
      **/
        ptrHybridTaus HyT = (ptrHybridTaus)malloc(sizeof(HybridTaus));
        HyT->z1=2*Rank;
        HyT->z2=3*Rank;
        HyT->z3=4*Rank;
        HyT->z4=5*Rank;

      for(int j = (IdThread)  * SizeStream ; j <  Streams->data[IdThread]  ; j++)
      {
          double Xm=0;
          for(int i = 0 ; i  < L->n; i++)
          {
              Xm+=GetTrajectory(L, LSum, F,P, j, Tolerancia, HyT);
          }
          X->data[j] = Xm/L->n;
      }
      MPI_Send(X->data, X->Size, MPI_DOUBLE, ROOT,TAG_RESULT, MPI_COMM_WORLD);

      free(HyT);

}


 MPI_Finalize();
    return 0;
}
