#include <stdlib.h>
#include <stdio.h>
#ifndef _OPENMP
#define omp_get_num_threads() 1
#define omp_get_thread_num() 1
#define omp_get_max_threads() 1
#define omp_get_wtime() 1
#endif // _OPENMP

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP
#include "Matrix.h"
#include "MyRand.h"
unsigned TausStep(unsigned *z, int S1, int S2, int S3, unsigned M)
{
    unsigned b=((((*z) << S1) ^ (*z)) >> S2);
    return (*z) = ((((*z) & M) << S3) ^ b);
}
unsigned LCGStep(unsigned *z, unsigned A, unsigned C)
{
    return (*z)=(A*(*z)+C);
}
double drand(ptrHybridTaus pseed)
{
    // Combined period is lcm(p1,p2,p3,p4)~ 2^121
    return 2.3283064365387e-10 * ( // Periods
    TausStep(&(pseed->z1), 13, 19, 12, 4294967294UL) ^ // p1=2^31-1
    TausStep(&(pseed->z2), 2, 25, 4, 4294967288UL) ^ // p2=2^30-1
    TausStep(&(pseed->z3), 3, 11, 17, 4294967280UL) ^ // p3=2^28-1
    LCGStep(&(pseed->z4), 1664525, 1013904223UL) // p4=2^32
);
}
double HybridRandom()
{

    srand(time(0));
    //printf("%d", rand());
    ptrHybridTaus HyT = (ptrHybridTaus)malloc(sizeof(HybridTaus));
    HyT->z1=2*omp_get_thread_num();
    HyT->z2=3*omp_get_thread_num();
    HyT->z3=4*omp_get_thread_num();
    HyT->z4=5*omp_get_thread_num();
    //Cargar números psudoaleatorios

   // for(register int i=0; i < MaxRandom; i++)
        return drand(HyT);
}
/*
Marsaglia polar Method
*/
double box_muller(double m, double s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * ((double)rand()/(double)RAND_MAX) - 1.0;
			x2 = 2.0 * ((double)rand()/(double)RAND_MAX) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}
double box_mullerTaus(ptrHybridTaus HyT, double m, double s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * drand(HyT) - 1.0;
			x2 = 2.0 * drand(HyT) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return( m + y1 * s );
}
double GenerarAleatorio(double min, double max)
{
    double f =(double)rand() / RAND_MAX;
    return min + f * (max - min);
}
void RandomMatrixDiagonalDominante(Matrix *A,Vector *b, double Min, double Max)
{
    Vector *SumaFila= NewVector(A->n);
    for(int i = 0; i < A->n; i++)
    {
        for(int j = 0; j < A->m; j++)
        {
            if(i==j) continue;
            A->data[i*A->m + j] = GenerarAleatorio(Min,Max);
            SumaFila->data[i]+=fabs(A->data[i*A->m + j]);
        }
        b->data[i] = GenerarAleatorio(Min,Max);
    }
    for(int i = 0; i < A->n; i++)
    {
        A->data[i*A->m+i]=30*SumaFila->data[i];
    }
    FreeVector(SumaFila);
}
