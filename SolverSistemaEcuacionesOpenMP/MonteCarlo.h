#ifndef MONTECARLO_H_INCLUDED
#define MONTECARLO_H_INCLUDED
#define Signo(x)(x < 0) ? -1.0 : 1.0
double getNorma(double *X, int Size);
void ComputeL(Matrix *A, Matrix *L, double ParametroRelajacion);
void ComputeLSUM(Matrix *L, Vector *LSum);
double GetSumaFila(Matrix *P, int IndexFila, int IndexLimiteFila);
double GetTrajectory(Matrix *L, Vector *LSum, Vector *F,Matrix *P,  int Index, double Tolerancia, ptrHybridTaus HyT );
void GetMatrixTransition(Matrix *L, Vector *LSum, Matrix *P);
void MontecarloAlmostOptimal(Matrix *A, Vector *b, Vector *X, double ParametroRelajacion, double Tolerancia);


#endif // MONTECARLO_H_INCLUDED
