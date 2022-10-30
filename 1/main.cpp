
#include "ReadMatrix.h"
#include "Solve.h"
#include <cmath>
#include <ctime>

double norma_nevyaski(double* A, double*b, double* x, int N);
double norma_pogreshnosty(double* x,double* x_real, int N);


int main()
{
    double* A;
    double* b;
    double* x;
    double* x_real;
    
    string filename;
    int n = -1; //размер матрицы
    int m = -1; //размер вывода
    int k = -1; // формула
    cin>>n>>m>>k;
    if(k == 0)
    {
        cin>>filename;
    }

    A = (double*) malloc(sizeof(double)*n*n);
    b = (double*) malloc(sizeof(double)*n);
    x = (double*) malloc(sizeof(double)*n);
    x_real = (double*) malloc(sizeof(double)*n);

    
    if(ReadMatrix(A,n,k, filename)!=0)
    {
        free(A);
        free(b);
        free(x);
        free(x_real);
        return 0;
    }

    for(int i = 0; i < n; i++)
    {
        for( int j = 0; j <= (n+1)/2; j++)
        {
            b[i] += A[(n*i) + (2*k + 1)];
        }
        x[i] = 0;
        x_real[i] = (i+1)%2;
    }
    PrintMatrix(A, n, n, m);
    PrintMatrix(b, 1, n, m);
    PrintMatrix(x, 1, n, m);

    int start = clock();
    int res = Solve(n, A, b, x);
    int end = clock(); 
    int time = (end - start)/CLOCKS_PER_SEC;// время работы  в секундах
    cout<<"время работы: "<<time<<endl;

    norma_nevyaski(A,b,x,n);
    norma_pogreshnosty(x,x_real, n);

    PrintMatrix(b, 1, n, m);
    PrintMatrix(x, 1, n, m);
    PrintMatrix(A, n, n, m);

    free(A);
    free(b);
    free(x);
    free(x_real);
    return 0;
}

double norma_nevyaski(double* A, double*b, double* x, int N)
{   //  ||Ax-b|| / ||b||
    double sum = 0; //ищем максимум по суммам модулей коэф. строк.
    double b_max = 0;
    double max  = 0;
    for( int i = 0; i < N; i++) // ||Ax-b|| и ||b||
    {
        for(int j =  0; j< N; j++)
        {
            sum+=A[i*N + j]* x[i];
        }
        sum - b[i];
        if(( abs(b[i]) - b_max) > EPS) // ||b||
        {
            b_max = abs(b[i]);
        }
        sum = abs(sum);
        if( (sum - max) > EPS)
        {
            max = sum;
        }
        sum = 0;
    }
    max/= b_max; //норма
    cout<< "норма невязки: "<<scientific<<max<< endl; 
    return max;
}

double norma_pogreshnosty(double* x,double* x_real, int N)
{
    double norma  = 0;
    for( int i = 0; i < N; i++) // ||x - x_real||
    {
        if(( abs(x[i] - x_real[i])- norma) > EPS)
        {
            norma = abs(x[i] - x_real[i]);
        }
    }
    cout<< "норма погрешности: "<<scientific<< norma<< endl; 
    return norma;
}