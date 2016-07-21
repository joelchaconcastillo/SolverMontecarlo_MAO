#include <stdio.h>
#include <stdlib.h>
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
   // #pragma omp parallel for
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
   // #pragma omp parallel for
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


    int NThreads = omp_get_max_threads();

    /**Primero obtener el tamaño de bloques que procesará cada hilo*/
    Vector *Streams = NewVector(NThreads);
    int SizeStream = L->m/NThreads;
    int Cont=0;

    for(int i = 0; i < NThreads; i++)
    {
        Streams->data[i] = (i+1)*SizeStream;
        Cont +=SizeStream;
    }
     //Agregar la diferencia faltante que corresponde al último hilo
        Streams->data[Streams->Size-1]+=L->m-Cont;


    #pragma omp parallel default(none) shared(Tolerancia, L, P, X, F, LSum, NThreads, Streams, SizeStream)
    {

        int IdThread = omp_get_thread_num();
        for(int j = IdThread  * SizeStream ; j <  Streams->data[IdThread]  ; j++)
        {
            ptrHybridTaus HyT = (ptrHybridTaus)malloc(sizeof(HybridTaus));
            HyT->z1=2*omp_get_thread_num();
            HyT->z2=3*omp_get_thread_num();
            HyT->z3=4*omp_get_thread_num();
            HyT->z4=5*omp_get_thread_num();
            double Xm=0;
            for(int i = 0 ; i  < L->n; i++)
            {
                Xm+=GetTrajectory(L, LSum, F,P, j, Tolerancia, HyT);
            }
            X->data[j] = Xm/L->n;
            free(HyT);
        }
    }



    FreeMatrix(L);
    FreeMatrix(P);
    FreeVector(LSum);
    FreeVector(Streams);
    FreeVector(F);
}
