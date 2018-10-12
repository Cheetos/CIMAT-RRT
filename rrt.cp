#include "rrt.h"

Individual P[N_INDIVIDUALS];
Individual S[SAMPLE];
double M[DIM];

CMatrix<double> COV(DIM,DIM);

double lowerBound;
double upperBound;

void Init(double l, double u)
{
	lowerBound = l;
	upperBound = u;

	ini_rand(time(NULL));

	for(int i=0; i<N_INDIVIDUALS; i++)
	{
		for(int j=0; j<DIM; j++)
			P[i].X[j] = lowerBound + (upperBound - lowerBound)*random();

		Verify_Boundaries(P[i]);

		if(i > 0)
			Extend(P[i],i);
		
		P[i].fitness = Evaluation_Function(P[i]);
		P[i].penalty = Penalty(P[i]);
	}

	Build_RRT();
}

void Extend(Individual &A, int n)
{
	int pos;
	double d, k;
	double m;

	Individual B;

	pos = 0;
	k = Distance(A,P[0]);

	for(int i=1; i<n; i++)
	{
		if(P[i].v != -1)
		{
			d = Distance(A,P[i]);
			if(d < k)
			{
				k = d;
				pos = i;
			}
		}
	}


	m = 0.0;
	for(int i=0; i<DIM; i++)
	{
		B.X[i] = A.X[i] - P[pos].X[i];
		m = m + B.X[i]*B.X[i];
	}

	m = sqrt(m);

	for(int i=0; i<DIM; i++)
		B.X[i] = B.X[i]/m;

	for(int i=0; i<DIM; i++)
	{
		B.X[i] = DELTA*B.X[i];
		A.X[i] = P[pos].X[i] + B.X[i];
	}

	A.v = pos;
}

void Build_RRT()
{
	CMatrix<double> Temp(N_INDIVIDUALS,DIM);
	CMatrix<double> L(DIM,DIM);
	CMatrix<double> I(N_INDIVIDUALS,DIM);

	Individual best;

	best = Strongest(best,Best_Individual(P,N_INDIVIDUALS));
	
	for(int g=0; g<200; g++)
	{
		Obtain_Sample(S,P,SAMPLE);

		Update_Mean(S,SAMPLE);
		Update_Covariance(S,SAMPLE);

    	//Check if COV has 0's in the diagonal
		for(int i=0; i<DIM; i++)
		{
			if(fabs(COV[i][i]) < 0.0001)
			{
				for(int j=0; j<DIM; j++)
				{
					COV[j][i] = 0.0;
					COV[i][j] = 0.0;
				}

				COV[i][i] = 1.0;
			}
		}

		L = Cholesky(COV);

		for(int i=0; i<N_INDIVIDUALS; i++)
			for(int j=0; j<DIM; j++)
				Temp[i][j] = rndNormal();

		I = Temp*L.Transpose();

		for(int i=0; i<N_INDIVIDUALS; i++)
		{
			for(int j=0; j<DIM; j++)
					P[i].X[j] = I[i][j] + M[j];

			Verify_Boundaries(P[i]);

			if(i > 0)
				Extend(P[i],i);

			P[i].fitness = Evaluation_Function(P[i]);
			P[i].penalty = Penalty(P[i]);
		}

		//best = Strongest(best,Best_Individual(P,N_INDIVIDUALS));
		//P[0] = best;
	}

	best = Best_Individual(P,N_INDIVIDUALS);

	for(int i=0; i<DIM; i++)
		cout << "X(" << i << ") = " << best.X[i] << endl;

	cout << "F = " << best.fitness << endl;
}

void Obtain_Sample(Individual *dst, Individual *src, int n)
{
	int a, b;

	for(int i=0; i<n; i++)
	{
		a = (int)(random()*N_INDIVIDUALS);
		b = (int)(random()*N_INDIVIDUALS);

		while(b == a)
			b = (int)(random()*n);

		dst[i] = Strongest(src[a],src[b]);
	}
}

Individual Strongest(Individual a, Individual b)
{
	if(a.penalty < EPS && b.penalty < EPS)
	{
		if(a.fitness < b.fitness)
			return a;
		else
			return b;
	}
	else
	{
		if(a.penalty < b.penalty)
			return a;
		else
			return b;
	}
}

Individual Best_Individual(Individual *A, int n)
{
	Individual best = A[0];

	for(int i=1; i<n; i++)
		best = Strongest(best,A[i]);

	return best;
}

void Update_Mean(Individual *A, int n)
{
	for(int i=0; i<DIM; i++)
	{
		M[i] = 0.0;
		for(int j=0; j<n; j++)
			M[i] = M[i] + A[j].X[i];

		M[i] = M[i]/n;
	}
}

void Update_Covariance(Individual *A, int n)
{
	for(int i=0; i<DIM; i++)
	{
		for(int j=i; j<DIM; j++)
		{
			COV[i][j] = Get_Covariance(A,n,i,j);
			COV[j][i] = COV[i][j];
		}
	}
}

double Get_Covariance(Individual *A, int n, int a, int b)
{
	double cov = 0.0;

	for(int i=0; i<n; i++)
		cov = cov + (A[i].X[a] - M[a])*(A[i].X[b] - M[b]);

	cov = cov/n;

	return cov;
}

void Verify_Boundaries(Individual &A)
{
	for(int i=0; i<DIM; i++)
	{
		if(A.X[i] < lowerBound || A.X[i] > upperBound)
        	A.X[i] = lowerBound + (upperBound - lowerBound)*random();
	}
}

double Distance(Individual &A, Individual &B)
{
	double d = 0.0;

	for(int i=0; i<DIM; i++)
		d = d + (A.X[i] - B.X[i])*(A.X[i] - B.X[i]);

	return d;
}

double Evaluation_Function(Individual &A)
{
	double f;
	double f1;

	//f = A.X[0];
	//f = (A.X[0] - 10.0)*(A.X[0] - 10.0) + (A.X[1] - 200.0)*(A.X[1] - 200.0);

	/***************************************
	*				Sphere				   *
	****************************************/

	f = 0.0;
	for(int i=0; i<DIM; i++)
		f = f + A.X[i]*A.X[i];

	return f;
}

double Penalty(Individual &A)
{
	double g[10];
	double s = 0.0;

	/*g[0] = 10000.0 - (A.X[0] - 200.0)*(A.X[0] - 200.0) - (A.X[1] - 0.0)*(A.X[1] - 0.0);
	g[1] = 4900.0 - (A.X[0] - 250.0)*(A.X[0] - 250.0) - (A.X[1] - 250.0)*(A.X[1] - 250.0);
	g[2] = 10000.0 - (A.X[0] - 0.0)*(A.X[0] - 0.0) - (A.X[1] - 400.0)*(A.X[1] - 400.0);

	for(int i=0; i<3; i++)
		if(g[i] > 0.0)
        	s = s + g[i];*/

	return s;
}
