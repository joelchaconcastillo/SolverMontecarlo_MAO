#include <stdio.h>
#include <stdlib.h>
#include "Matrix.h"
Matrix * NewMatrix(int n, int m)
{
    Matrix *M = (Matrix *) malloc( sizeof(Matrix));
    M->data = ( double *) malloc( n*m*sizeof(double));
    M->n=n;
    M->m=m;
    for(int i = 0; i < n*m; i++) M->data[i] = 0;
    return M;
}
Vector * NewVector(int n)
{
    Vector * V = ( Vector *) malloc( sizeof(Vector));
    V->Size=n;
    V->data = (double *) malloc( n*sizeof(double));
    for(int i=0; i < n; i++) V->data[i] = 0;
    return V;
}
void ResetVector(Vector *X){
    for(int i = 0 ; i < X->Size; i++) X->data[i]=0;
}
void cpyRow(Matrix * A, Matrix *B, int Index    )
{
    for(int i = 0; i < A->m; i++)
    A->data[Index*A->m + i ] = B->data[Index*B->m + i ];
}
void cpyVector(Vector *A, Vector *B)
{
    for(int i = 0; i < A->Size; i++)
    A->data[i] = B->data[i];
    A->Size = B->Size;
}
void cpyVectorMatrix(Vector *A, Matrix *B, int Index)
{
    for(int i = 0; i < A->Size; i++)
    B->data[ Index*B->m + i] = A->data[i];
}
void cpyMatrixVector(Vector *A, Matrix *B, int Index)
{
    for(int i = 0; i < A->Size; i++)
        A->data[i] = B->data[ Index*A->Size + i];

}
void ResetMatrix(Matrix *M){
    for(int i = 0; i < M->m * M->n; i++)
        M->data[i]=0;
}
void PrintMatrix(Matrix * M)
{
    for(int i = 0; i < M->n; i++)
    {
        for(int j = 0; j < M->m; j++)
        {
            printf("%f ", M->data[i*M->m + j]);
        }
        printf("\n");
    }
}
void PrintVector(Vector *V)
{
    for(int i = 0; i < V->Size; i++)
    printf("%f\n", V->data[i]);
    printf("\n");
}
void FreeMatrix(Matrix *M)
{
    free(M->data);
    free(M);
}
void FreeVector(Vector *V)
{
    free(V->data);
    free(V);
}
