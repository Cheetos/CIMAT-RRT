#include "numerical_methods.h"

CMatrix<double> Cholesky(CMatrix<double> &A)
{
	int n;
	double sum;
	CMatrix<double> L = A;

	if(L.rows() != L.columns())
		return L;

	n = L.columns();

	for(int i=0; i<n; i++)
	{
		sum = 0.0;
		for(int k=0; k<i; k++)
			sum = sum + L[i][k]*L[i][k];

		L[i][i] = sqrt(A[i][i] - sum);

		for(int j=i+1; j<n; j++)
		{
			sum = 0.0;
			for(int k=0; k<i; k++)
				sum = sum + L[i][k]*L[j][k];

			L[j][i] = (A[i][j] - sum)/L[i][i];
			L[i][j] = 0.0;
		}
	}

	return L;
}

bool Solve_by_Jacobi(CMatrix<double> A, CMatrix<double> &B, double tol, int maxIterations) {

    CMatrix<double> U(A.rows(), A.columns());
    CMatrix<double> Ut(A.columns(), A.rows());
    CMatrix<double> S(A.rows(), A.columns());
	CMatrix<double> At = A.Transpose();

	A = At*A;
	B = At*B;

    if(!Jacobi_Rotation(A,S,U,tol,maxIterations))
        return false;

    B = U*B;
    for(int i=0; i<B.rows(); i++)
    {
        if(fabs(S[i][i]) < tol)
            return 0;

        B[i][0] = B[i][0]/S[i][i];
    }

    Ut = U.Transpose();
    B = Ut*B;

    return true;
}

bool Jacobi_Rotation(CMatrix<double> &A, CMatrix<double> &D, CMatrix<double> &Q, double tol, int maxIterations) {
    int r, c, n, nIterations = 0;
    double max, S, C, T, k, num, den;

    CMatrix<double> Qtemp(A.rows(), A.columns());
    n = A.rows();

    D = A;
    Q.Identity();

    do
    {
        max = Jacobi_Infinite_Norm(A,r,c);
        if(max < tol)
            return true;

       k = (A[r][r] - A[c][c])/(2.0*A[r][c]);

       if(k < tol)
            T = -1.0*(k + sqrt(1.0+(k*k)));
       else
            T = sqrt(1.0+(k*k)) - k;

       C = 1/sqrt(T*T+1.0);
       S = T*C;

       D = A;
       Qtemp = Q;

       D[c][c] = A[c][c]*C*C + A[r][r]*S*S - 2.0*A[r][c]*S*C;
       D[r][r] = A[c][c]*S*S + A[r][r]*C*C + 2.0*A[r][c]*S*C;
       D[r][c] = 0.0;
       D[c][r] = 0.0;

       for(int i=0; i<n; i++)
       {
            Qtemp[c][i] = Q[c][i]*C - Q[r][i]*S;
            Qtemp[r][i] = Q[r][i]*C + Q[c][i]*S;
       }

       for(int i=0; i<n; i++)
       {
            if(i != c && i != r)
            {
                D[c][i] = A[c][i]*C - A[r][i]*S;
                D[r][i] = A[c][i]*S + A[r][i]*C;
                D[i][c] = D[c][i];
                D[i][r] = D[r][i];
            }
       }

       A = D;
       Q = Qtemp;
    }
    while(++nIterations < maxIterations);

    return true;
}

double Jacobi_Infinite_Norm(CMatrix<double> &A, int &r, int &c) {
    double max = 0.0;

    r = 1;
    c = 0;

    for(int i=0; i<A.rows(); i++)
    {
        for(int j=0; j<i; j++)
        {
            if(fabs(A[i][j]) > max)
            {
                max = fabs(A[i][j]);
                r = i;
                c = j;
            }
        }
    }

    return max;
}

void Swap(double &a, double &b)
{
	double temp = a;
	a = b;
	b = temp;
}