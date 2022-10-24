
#include "ReadMatrix.h"

double norma(double* A, double*b);

int main()
{
    double* A;
    double* b;
    string filename;
    int n = -1; //размер матрицы
    int m = -1; //размер вывода
    int k = -1; // формула
    cin>>n>>m>>k;

    if(k == 0)
    {
        if(!getline(cin, filename))
        {
            cout<<"-2"<<endl;
            return 0;
        }
    }

    A = (double*) malloc(sizeof(double)*n*n);
    b = (double*) malloc(sizeof(double)*n);
    
    ReadMatrix(A,n,k, filename);
    for(int i = 0; i < n; i++)
    {
        for( int j = 0; j <= (n+1)/2; j++)
        {
            b[i] += A[(n*i) + (2*k + 1)];
        }
        
    }

    return 0;
}

double norma(double* A, double*b)
{
    for( int i = 0; i < N; i++)
    {
        for(int j =  0; j< N; j++)
        {
            x[i]+=A[i*N + j];
        }
    }
    int max = x[0];
    for( int i = 1; i < N; i++)
    {
        if(max <= x[i])
        {
            max = x[i];
        }
    }
    return max;
}
