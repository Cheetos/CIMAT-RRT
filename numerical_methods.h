#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H

#include "matrix.h"

CMatrix<double> Cholesky(CMatrix<double> &);
bool Solve_by_Jacobi(CMatrix<double>, CMatrix<double> &, double, int);
bool Jacobi_Rotation(CMatrix<double> &, CMatrix<double> &, CMatrix<double> &, double, int);
double Jacobi_Infinite_Norm(CMatrix<double> &, int &, int &);
void Swap(double &, double &);

#endif
 
