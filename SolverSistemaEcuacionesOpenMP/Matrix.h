#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

typedef struct Matrix
{
	int n;
	int m;
    double *data;
}Matrix;
typedef struct Vector
{
	int Size;
	double *data;
}Vector;
typedef struct
{
	Vector *Elite;
	Vector *Media;
	Vector *Varianza;
	int Flag;
}Packet;
Matrix * NewMatrix(int n, int m);
Vector * NewVector(int n);
void ResetVector(Vector *X);
void cpyRow(Matrix * A, Matrix *B, int Index);
void cpyVector(Vector *A, Vector *B);
void cpyVectorMatrix(Vector *A, Matrix *B, int Index);
void cpyMatrixVector(Vector *A, Matrix *B, int Index);
void PrintMatrix(Matrix * M);
void PrintVector(Vector *V);
void FreeMatrix(Matrix *M);
void FreeVector(Vector *V);


#endif // MATRIX_H_INCLUDED
