#ifndef MYRAND_H_INCLUDED
#define MYRAND_H_INCLUDED
#define SIZEMATRIX 500
#define MaxRandom 1000000
typedef struct InHybridTaus
{
    unsigned z1, z2, z3, z4;
} *ptrHybridTaus, HybridTaus;

unsigned TausStep(unsigned *z, int S1, int S2, int S3, unsigned M);
unsigned LCGStep(unsigned *z, unsigned A, unsigned C);
double drand(ptrHybridTaus pseed);
/*
Marsaglia polar Method
*/

double HybridRandom();

double box_muller(double m, double s);
double box_mullerTaus(ptrHybridTaus HyT, double m, double s);
double GenerarAleatorio(double min, double max);
void RandomMatrixDiagonalDominante(Matrix *A, Vector *b, double Max, double Min);
#endif // MYRAND_H_INCLUDED
