#ifndef MONTECARLO_H_INCLUDED
#define MONTECARLO_H_INCLUDED
#define Signo(x)(x < 0) ? -1.0 : 1.0
#define ROOT 0
#define TAG_FLAG 10
#define TAG_INDEX_VARIABLE 11
#define TAG_RESULT 12
#define INFORMATIONMATRIX 1000
double getNorma(double *X, int Size);
void ComputeL(Matrix *A, Matrix *L, double ParametroRelajacion);
void ComputeLSUM(Matrix *L, Vector *LSum);
double GetSumaFila(Matrix *P, int IndexFila, int IndexLimiteFila);
double GetTrajectory(Matrix *L, Vector *LSum, Vector *F,Matrix *P,  int Index, double Tolerancia, ptrHybridTaus HyT );
void GetMatrixTransition(Matrix *L, Vector *LSum, Matrix *P);
void MontecarloAlmostOptimal(Matrix *A, Vector *b, Vector *X, double ParametroRelajacion, double Tolerancia);
int GetPackSize(Matrix *L, Vector *LSum, Vector *F, Matrix *P, Vector *Streams, int SizeStream);
unsigned char * EmpaquetarMatrix(Matrix *L, Vector *LSum, Vector *F, Matrix *P, double Tolerancia, int packsize, Vector *Streams,  int SizeStream);
void DesenpaquetarMatrix(Matrix *L, Vector *LSum, Vector *F, Matrix *P, double *Tolerancia, Vector *Streams,  int *SizeStream);
#endif // MONTECARLO_H_INCLUDED
