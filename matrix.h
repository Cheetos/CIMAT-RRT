#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

template <class T>
class CArray 
{
    private:
        T *data;
        int size;
    public:
        CArray() 
		{
            data = NULL;
            size = 0;
        }
        
        CArray(int n) 
		{
            data = new T[n];
            size = n;
        }
        
        T &operator [] (int n)
		{
            if(n < 0 || n >= size)
                throw 1;
            return data[n];
        }
};

template <class T>
class CMatrix 
{
    private:
        int nRows;
        int nColumns;
        CArray<T> *data;
    public:
        /******************************************************
        *                   Constructors                      *
        *******************************************************/

        CMatrix() 
		{
            data = NULL;
            nRows = 0;
            nColumns = 0;
        }
        
        CMatrix(int r, int c) 
		{
            nRows = r;
            nColumns = c;  
            data = new CArray<T> [r];
            for(int i=0; i<r; i++)
                data[i] = CArray<T>(c);
        }

        /******************************************************
        *               Constructor by Copy                   *
        *******************************************************/

        CMatrix(const CMatrix<T> &b) 
		{
            if(b.nRows > 0 && b.nColumns > 0) 
			{
                nRows = b.nRows;
                nColumns = b.nColumns;
                data = new CArray<T> [b.nRows];
                for(int i=0; i<b.nRows; i++)
                    data[i] = CArray<T>(b.nColumns);

                for(int i=0; i<nRows; i++)
                    for(int j=0; j<nColumns; j++)
                        data[i][j] = b.data[i][j];
            }
            else 
			{
                data = NULL;
                nRows = 0;
                nColumns = 0;
            }
        }

        /******************************************************
        *               Destructor Definition                 *
        *******************************************************/

        ~CMatrix() 
		{
            if(nRows > 0)
                delete [] data;
        }

        /******************************************************
        *           Operators Definition []                   *
        *******************************************************/

        CArray<T> &operator [] (int n) 
		{
            if(n < 0 || n >= nRows)
                throw 1;
        
            return data[n];
        }

        /******************************************************
        *           Operators Definition =                *
        *******************************************************/

        CMatrix<T> &operator = (const CMatrix<T> &b) 
		{
            if(nRows > 0)
                delete [] data;

            if(b.nRows > 0 && b.nColumns > 0) 
			{
                nRows = b.nRows;
                nColumns = b.nColumns;
                data = new CArray<T> [nRows];
                for(int i=0; i<nRows; i++)
                    data[i] = CArray<T>(nColumns);

                for(int i=0; i<nRows; i++)
                    for(int j=0; j<nColumns; j++)
                        data[i][j] = b.data[i][j];
            }
            else 
			{
                data = NULL;
                nRows = 0;
                nColumns = 0;
            }
                    
            return *this;
        }
        
        CMatrix<T> &operator = (const T val) 
		{
            for(int i=0; i<nRows; i++)
                for(int j=0; j<nColumns; j++)
                    data[i][j] = val;
            return *this;
        }

        /******************************************************
        *           Operators Definition +,-,*                *
        *******************************************************/

        CMatrix<T> operator + (const CMatrix<T> &b) 
		{
            if(nRows != b.nRows || nColumns != b.nColumns)
                return *this;

            CMatrix<T> temp(nRows,nColumns);

            for(int i=0; i<nRows; i++)
                for(int j=0; j<nColumns; j++)
                    temp.data[i][j] = data[i][j] + b.data[i][j];

            return temp;
        }

        CMatrix<T> operator - (const CMatrix<T> &b)
		{
            if(nRows != b.nRows || nColumns != b.nColumns)
                return *this;
                
            CMatrix<T> temp(nRows,nColumns);
            
            for(int i=0; i<nRows; i++)
                for(int j=0; j<nColumns; j++)
                    temp.data[i][j] = data[i][j] - b.data[i][j];
            
            return temp;
        }
        
        CMatrix<T> operator * (const CMatrix<T> &b) 
		{
            if(nColumns != b.nRows)
                return *this;
            
            CMatrix<T> temp(nRows,b.nColumns);
            temp = 0;
            
            for(int i=0; i<nRows; i++)
                for(int j=0; j<b.nColumns; j++)
                    for(int k=0; k<nColumns; k++)
                        temp.data[i][j] = temp.data[i][j] + data[i][k] * b.data[k][j];
            
            return temp;
        }

        CMatrix<T> operator * (const T val) 
		{
            CMatrix<T> temp(nRows,nColumns);

            for(int i=0; i<nRows; i++)
                for(int j=0; j<nColumns; j++)
                        temp.data[i][j] = data[i][j]*val;
            
            return temp;
        }

        /******************************************************
        *           Operators Definition +=,-=                *    
        *******************************************************/

        CMatrix<T> &operator += (const CMatrix<T> &b) 
		{
            if(nRows != b.nRows || nColumns != b.nColumns)
                return *this;

            for(int i=0; i<nRows; i++)
                for(int j=0; j<nColumns; j++)
                    data[i][j] = data[i][j] + b.data[i][j];

            return *this;
        }

        CMatrix<T> &operator -= (const CMatrix<T> &b) 
		{
            if(nRows != b.nRows || nColumns != b.nColumns)
                return *this;

            for(int i=0; i<nRows; i++)
                for(int j=0; j<nColumns; j++)
                    data[i][j] = data[i][j] - b.data[i][j];

            return *this;
        }
        
        /******************************************************
        *           Operators Definition <<                   *
        *******************************************************/

        friend ostream &operator << (ostream &os, const CMatrix<T> &b) 
		{
            for(int i=0; i<b.nRows; i++)
			{
                for(int j=0; j<b.nColumns; j++)
                    os << setiosflags(ios::fixed) << setw(10) << setprecision(8) << b.data[i][j] << " ";
                os << endl;
            }
            os << endl;
            
            return os;
        }

        /******************************************************
        *               Methods Definition                     *
        *******************************************************/

		void Resize(int r, int c)
		{
			if(nRows > 0)
				delete [] data;

			nRows = r;
            nColumns = c;  
            data = new CArray<T> [r];
            for(int i=0; i<r; i++)
                data[i] = CArray<T>(c);
		}

        CMatrix<T> Transpose() 
		{
            CMatrix<T> temp(nColumns,nRows);

            for(int i=0; i<nRows; i++)
                for(int j=0; j<nColumns; j++)
                    temp.data[j][i] = data[i][j];

            return temp;
        }

		CMatrix<T> Inverse()
		{
			int r;
			T max, p;

			CMatrix<T> B = *this;
			CMatrix<T> A(nRows,2*nColumns);

			if(nRows != nColumns)
				return B;

			A = 0;
			for(int i=0; i<nRows; i++)
			{
				A.data[i][i+nColumns] = 1;
				for(int j=0; j<nColumns; j++)
					A.data[i][j] = data[i][j];
			}

			for(int i=0; i<nRows; i++)
			{
				max = A.data[i][i];
				r = i;

				for(int j=i; j<nRows; j++)
				{
					if(fabs(A.data[j][i]) > max)
					{
						max = fabs(A.data[j][i]);
						r = j;
					}
				}

				A.Swap_Rows(i,r);
				for(int j=i; j<2*nColumns; j++)
					A.data[i][j] = A.data[i][j]/max;

				for(int j=0; j<nRows; j++)
				{
					p = A.data[j][i];
					if(i != j)
					{
						for(int k=i; k<2*nColumns; k++)
						{
							A.data[j][k] = A.data[j][k] - p*A.data[i][k]/A.data[i][i];
							if(fabs(A.data[j][k]) < 0.000000001)
								A.data[j][k] = 0.0;
						}
					}
				}
			}

			for(int i=0; i<nRows; i++)
				for(int j=0; j<nColumns; j++)
					B.data[i][j] = A.data[i][j+nColumns];

			return B;
		}

        void Identity()
		{
            if(nRows == nColumns) 
			{
                for(int i=0; i<nRows; i++)
                    for(int j=0; j<nColumns; j++)
                        if(i == j)
                            data[i][j] = 1;
                        else
                            data[i][j] = 0;
            }
        }

        double Infinite_Norm(int &r, int &c) 
		{
            double max = 0.0;
            for(int i=0; i<nRows; i++)
                for(int j=0; j<nColumns; j++)
                    if(fabs(double(data[i][j])) > max)
					{
                        max = fabs(double(data[i][j]));
                        r = i;
                        c = j;
                    }

            return max;
        }

        double Norm2() 
		{
            double max = 0.0;
            for(int i=0; i<nRows; i++)
                for(int j=0; j<nColumns; j++)
                    max = max + double(data[i][j]) * double(data[i][j]);

            return sqrt(max);
        }
        
        void Delete_Row(int r)
        {
			if(r <= nRows && nRows > 0)
			{
				CMatrix<T> temp = *this;
				
				delete [] data;
				
				data = new CArray<T> [nRows-1];
				for(int i=0; i<nRows-1; i++)
					data[i] = CArray<T>(nColumns);
				
				int k = 0;
				for(int i=0; i<nRows; i++)
				{
					if(i == r) continue;
					
					for(int j=0; j<nColumns; j++)
						data[k][j] = temp.data[i][j];
						
					k++;
				}
	
				nRows--;
			}
		}
		
		void Delete_Column(int c)
        {
			if(c <= nColumns && nColumns > 0)
			{
				CMatrix<T> temp = *this;
				
				delete [] data;
				
				data = new CArray<T> [nRows];
				for(int i=0; i<nRows; i++)
					data[i] = CArray<T>(nColumns-1);
				
				
				for(int i=0; i<nRows; i++)
				{	
					int k = 0;
					for(int j=0; j<nColumns; j++)
					{
						if(j == c)	continue;
						data[i][k++] = temp.data[i][j];
					}
				}
	
				nColumns--;
			}
		}
		
		void Swap_Rows(int r1, int r2)
		{
			T temp;

			if(r1 >= 0 && r1 < nRows && r2 >= 0 && r2 < nRows)
			{
				for(int i=0; i<nColumns; i++)
				{
					temp = data[r1][i];
					data[r1][i] = data[r2][i];
					data[r2][i] = temp;
				}
			}
		}
		
        int columns() {return nColumns;}
        int rows() {return nRows;}
};

#endif
