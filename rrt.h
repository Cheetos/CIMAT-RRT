#ifndef RRT_H
#define RRT_H

#include <time.h>
#include "rand.h"
#include "matrix.h"
#include "numerical_methods.h"
#include "global.h"

typedef struct stIndividual
{
	double X[DIM];
	double fitness;
	double penalty;
	int v;
}Individual;

extern Individual P[];
extern Individual S[];
extern double M[];
extern double lowerBound;
extern double upperBound;

void Init(double, double);
void Verify_Boundaries(Individual &);
void Build_RRT();
void Extend(Individual &, int);
void Build_RRT();
void Obtain_Sample(Individual *, Individual *, int);
void Update_Mean(Individual *, int);
void Update_Covariance(Individual *, int);

Individual Strongest(Individual, Individual);
Individual Best_Individual(Individual *, int);

double Get_Covariance(Individual *, int, int, int);
double Distance(Individual &, Individual &);
double Evaluation_Function(Individual &);
double Penalty(Individual &);

#endif
