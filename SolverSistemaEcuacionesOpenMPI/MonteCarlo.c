/**
    Autor: Joel Chacon castillo
    Fecha: 02/Junio/2016
    Descripcion: Este fichero distribuye con
                 el modelo maestro esclavo las trayectorias de
                 la cadena de markov, cada cadena pertenece a su respectiva variable
**/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "omp.h"
#include "Matrix.h"
#include "MyRand.h"
#include "MonteCarlo.h"
double getNorma(double *X, int Size)
{
    double Norma=0;
    for(int i = 0; i < Size; i++)
    Norma+=pow(X[i],2);
    return sqrt(Norma);
}
void ComputeL(Matrix *A, Matrix *L, double ParametroRelajacion)
{
    #pragma omp parallel for
    for(int i = 0; i < L->n; i++)
    {
        for(int j = 0; j < L->m; j++)
        {
            if(i==j)
                L->data[i*L->m + j] = 1.0-ParametroRelajacion;
            else
                    L->data[i*L->m + j] = -ParametroRelajacion*( A->data[i*A->m+j] / A->data[i*A->m+i]  );
        }
    }
}
void ComputeLSUM(Matrix *L, Vector *LSum)
{
    #pragma omp parallel for
    for(int i = 0; i < L->n; i++)
    {
        for(int j = 0; j < L->m; j++)
        {
            LSum->data[i]+=fabs(L->data[i*L->m + j]);
        }
    }
}
double GetSumaFila(Matrix *P, int IndexFila, int IndexLimiteFila)
{
    double Suma=0;
    for(int j = 0; j < IndexLimiteFila; j++) Suma+=P->data[IndexFila*P->m + j];
    return Suma;
}
double GetTrajectory(Matrix *L, Vector *LSum, Vector *F,Matrix *P,  int Index, double Tolerancia, ptrHybridTaus HyT )
{

    /**Asignar valores iniciales*/
    double X=0, W=1;
    X+=W*F->data[Index];
    double Epsilon =  Tolerancia;

    while(fabs(W) > Epsilon)
    {
        double Eta=  drand(HyT);// ( (double)rand() )/RAND_MAX;

        int j = 0;
        while(Eta > GetSumaFila(P, Index, ++j)) ;
        j--;
        W= W*(Signo(L->data[Index*L->m + j]) )*LSum->data[Index];
        X= X+W*F->data[j];
        Index=j;
    }
    return X;
}
void GetMatrixTransition(Matrix *L, Vector *LSum, Matrix *P)
{
    for(int i = 0; i < P->n; i++)
    {
        for(int j = 0; j < P->m; j++)
        {
            P->data[i*P->m + j] = fabs(L->data[i*P->m + j]) / LSum->data[i];
        }
    }
}
void MontecarloAlmostOptimal(Matrix *A, Vector *b, Vector *X, double ParametroRelajacion, double Tolerancia)
{

    /** Calcular la matriz L*/

    Matrix *L = NewMatrix(A->n, A->m);

    ComputeL(A,L, ParametroRelajacion);

    Vector *F = NewVector(b->Size);
    for(int i = 0; i < F->Size; i++) F->data[i] = ParametroRelajacion*(b->data[i]/A->data[i*A->m + i]);


    /**Calcular el vector suma que corresponde a cada fila  */
    Vector * LSum= NewVector(b->Size);
    ComputeLSUM(L, LSum);
    /**Generar la matriz de probabilidades Matriz de transición P*/
    Matrix *P=NewMatrix(A->n, A->m);
    GetMatrixTransition(L, LSum, P);


    /***Calcular la trayectoria*/

    int NThreads, Rank;
    MPI_Comm_size(MPI_COMM_WORLD, &NThreads);
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

    /**
      No contemplar al proceso maestro...
    */
    //NThreads--;
    /**Primero obtener el tamaño de bloques que procesará cada hilo*
    Vector *Streams = NewVector(NThreads-1);
    int SizeStream = L->m/NThreads;
    int Cont=0;

    for(int i = 0; i < NThreads; i++)
    {
        Streams->data[i] = (i+1)*SizeStream;
        Cont +=SizeStream;
    }
     //Agregar la diferencia faltante que corresponde al último hilo
        Streams->data[Streams->Size-1]+=L->m-Cont;
  /**
        Distribuir el cálculo de cada variable en cada proceso.
  **/
  int packsize = GetPackSize(L, LSum, F, P);
  unsigned char *package =  EmpaquetarMatrix(L, LSum, F, P, Tolerancia, packsize);
  int FLAG =1;

  /**
  Enviar la información de la matriz a procesar...
  **/

  for(int i = 0; i  < NThreads; i++)
  {
    if(i== ROOT) continue;
    MPI_Send(&packsize, 1, MPI_INT, i,INFORMATIONMATRIX, MPI_COMM_WORLD);
    MPI_Send(package, packsize, MPI_UNSIGNED_CHAR, i,INFORMATIONMATRIX, MPI_COMM_WORLD);
  }

  /**
    Enviar el índice de las variables que se van a procesar...
  **/

  for(int i = 0; i  < F->Size; i++)
  {
    int TemporalRank = i%(NThreads-1) + 1;
    MPI_Send(&FLAG, 1, MPI_INT, TemporalRank,TAG_FLAG, MPI_COMM_WORLD);
    MPI_Send(&i, 1, MPI_INT, TemporalRank,TAG_INDEX_VARIABLE, MPI_COMM_WORLD);
  }


  for(int i = 0; i  < F->Size; i++)
  {
    int TemporalRank = i%(NThreads-1) + 1;
    MPI_Recv( &(X->data[i] ), 1, MPI_DOUBLE, TemporalRank, TAG_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  /**
    Indicar que ya se calcularon todas las variables
  **/

  FLAG = 0;
  for(int i = 0; i < NThreads;i++)
  {
      if(i==ROOT) continue;
      MPI_Send(&FLAG, 1, MPI_INT, i,TAG_FLAG, MPI_COMM_WORLD);
  }

    FreeMatrix(L);
    FreeMatrix(P);
    FreeVector(LSum);
    FreeVector(F);
}
int GetPackSize(Matrix *L, Vector *LSum, Vector *F, Matrix *P)
{
  /**
      Calcular el tamaño de la información que se va a enviar..
  **/
    int increment=0, packsize=0;
    MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &increment);
    packsize+=increment;
    MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &increment);
    packsize+=increment;
    MPI_Pack_size(L->n*L->m, MPI_DOUBLE, MPI_COMM_WORLD, &increment);
    packsize+=increment;
    MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &increment);
    packsize+=increment;
    MPI_Pack_size(LSum->Size, MPI_DOUBLE, MPI_COMM_WORLD, &increment);
    packsize+=increment;
    MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &increment);
    packsize+=increment;
    MPI_Pack_size(F->Size, MPI_DOUBLE, MPI_COMM_WORLD, &increment);
    packsize+=increment;
    MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &increment);
    packsize+=increment;
    MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &increment);
    packsize+=increment;
    MPI_Pack_size(P->n*P->m, MPI_DOUBLE, MPI_COMM_WORLD, &increment);
    packsize+=increment;
    MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &increment);
    packsize+=increment;
    return packsize;

}
unsigned char * EmpaquetarMatrix(Matrix *L, Vector *LSum, Vector *F, Matrix *P, double Tolerancia, int packsize)
{


  int position =0;
  /**
    Reservar memoria para el paquete...
  **/
  unsigned char *package = (unsigned char *) malloc(packsize);

  /**
      Agregar la información al paquete
  **/

    MPI_Pack(&(L->n), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
    MPI_Pack(&(L->m), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
    MPI_Pack(L->data, L->n*L->m, MPI_DOUBLE, package, packsize, &position, MPI_COMM_WORLD );
    MPI_Pack(&(LSum->Size), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
    MPI_Pack(LSum->data, LSum->Size, MPI_DOUBLE, package, packsize, &position, MPI_COMM_WORLD );
    MPI_Pack(&(F->Size), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
    MPI_Pack(F->data, F->Size, MPI_DOUBLE, package, packsize, &position, MPI_COMM_WORLD );
    MPI_Pack(&(P->n), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
    MPI_Pack(&(P->m), 1, MPI_INT, package, packsize, &position, MPI_COMM_WORLD );
    MPI_Pack(P->data, P->n*P->m, MPI_DOUBLE, package, packsize, &position, MPI_COMM_WORLD );
    MPI_Pack(&(Tolerancia), 1, MPI_DOUBLE, package, packsize, &position, MPI_COMM_WORLD );
    return package;

}
void DesenpaquetarMatrix(Matrix *L, Vector *LSum, Vector *F, Matrix *P, double *Tolerancia)
{
    int position=0;
		int increment, packsize=0;
		MPI_Recv(&packsize, 1, MPI_INT, ROOT, INFORMATIONMATRIX, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    unsigned char *package = (unsigned char *) malloc(packsize);
    MPI_Recv(package, packsize, MPI_UNSIGNED_CHAR, ROOT, INFORMATIONMATRIX, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Unpack(package, packsize, &position, &(L->n), 1, MPI_INT, MPI_COMM_WORLD );
    MPI_Unpack(package, packsize, &position, &(L->m), 1, MPI_INT, MPI_COMM_WORLD );
    L->data = (double *) malloc( L->n*L->m * sizeof(double));
    MPI_Unpack(package, packsize, &position, L->data, L->n*L->m, MPI_DOUBLE, MPI_COMM_WORLD );

    MPI_Unpack(package, packsize, &position, &(LSum->Size), 1, MPI_INT, MPI_COMM_WORLD );
    LSum->data = (double *) malloc( LSum->Size * sizeof(double));
    MPI_Unpack(package, packsize, &position, LSum->data, LSum->Size, MPI_DOUBLE, MPI_COMM_WORLD );

    MPI_Unpack(package, packsize, &position, &(F->Size), 1, MPI_INT, MPI_COMM_WORLD );
    F->data = (double *) malloc( F->Size * sizeof(double));
    MPI_Unpack(package, packsize, &position, F->data, F->Size, MPI_DOUBLE, MPI_COMM_WORLD );

    MPI_Unpack(package, packsize, &position, &(P->n), 1, MPI_INT, MPI_COMM_WORLD );
    MPI_Unpack(package, packsize, &position, &(P->m), 1, MPI_INT, MPI_COMM_WORLD );
    P->data = (double *) malloc( P->n*P->m * sizeof(double));
    MPI_Unpack(package, packsize, &position, P->data, P->n*P->m, MPI_DOUBLE, MPI_COMM_WORLD );

    MPI_Unpack(package, packsize, &position, Tolerancia, 1, MPI_DOUBLE, MPI_COMM_WORLD );

}
