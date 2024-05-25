#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <float.h>
#include <math.h>
#include <iomanip>
#include <algorithm>

using namespace std;

// вариант 1000

void Generate_d(double* d, double* x, int N);
void Generate_x(double* MATRIX, int N);
void Generate_f(double* f, double* x, int N);
double func(double x);
double b(double x);
void Progonka(int N, double* p, double* q, double* d, double* y, double*x, double*f);
void fourier(double* C, int N, double*y, double*x, double*f);
void Get_Coef(  double *C, int N, double*f, double*x);
double lambda_n(int n, int N);
double dot_f_phi(int N, int j, double*f, double*x);

//splot '1.txt' using 1:2:3 with points, '1.txt' using 1:2:4 with points

int main()
{
    double *y,*x, *p, *q, *d, *f,*C;
    int N = -1; // число узлов
    int constant = -1;

    if(!(scanf("%d %d",&N, &constant)==1))
    {
        //ошибка чтения
        return -1;
    }

    if(N<=0)
    {
        //ошибка 
        cout<<N<<endl;
        cout<<" N >= 0"<<endl;
        return -1;
    }
    
    y = (double*) malloc(sizeof(double)*(N+1));
    x = (double*) malloc(sizeof(double)*(N+1));
    p = (double*) malloc(sizeof(double)*(N+1));
    q = (double*) malloc(sizeof(double)*(N+1));
    d = (double*) malloc(sizeof(double)*(N+1));
    f = (double*) malloc(sizeof(double)*(N+1));
    C = (double*) malloc(sizeof(double)*(N+1));

    Generate_x(x, N);
    Generate_f(f, x, N);

    if(constant == 1)
    {
        Generate_d(d, x, N);
        Progonka(N, p, q, d, y, x,f);
    }
    else if(constant == 0)
    {
        fourier(C,N,y,x,f);
    }
    Write(N,y);
    

    free(x);
    free(y);
    free(p);
    free(q);
    free(d);
    free(f);
    free(C);
    return 0;
}

int Write(int N, double* y)
{
    
}

void Generate_d(double* d, double* x, int N)
{
    for(int k = 0; k<=N; k++)
    {
        d[k] = 2 + b(x[k])*(1./(double)(N*N)) ;
    }
    return;
}

void Generate_f(double* f, double* x, int N)
{
    for(int k = 0; k<=N; k++)
    {
        f[k] = func(x[k])*(1./(double)(N*N)) ;
    }
    return;
}

void Progonka(int N, double* p, double* q, double* d, double* y, double*x, double*f)
{
    p[0] =(double)(1./d[0]);
    q[0] = (double)(f[0]/d[0]);
    for(int k  = 0; k < N; k++)
    {
        p[k+1] = 1./(d[k]-p[k]);
        q[k+1] = (f[k] +q[k])/(d[k] - p[k]);
    }

    //обратно

    y[N] = (f[N] + q[N])/(d[N] - p[N]);
    for(int k = N - 1; k >= 0;k--)
    {
        y[k] = p[k+1]*y[k+1] + q[k+1];
    }
    return;
}

double func(double x)
{
    return (x*x);
}

double b(double x)
{
    return 2*x + 1;
}

void Generate_x(double* MATRIX, int N)
{
    for(int i = 0; i < N; i++)
    {
        MATRIX[i] = i*(1/(double)N) ;
    }
    MATRIX[N] = 1.0;
    return;
}



void fourier(double* C, int N, double*y, double*x, double*f)
{
    double res = 0.0;
    Get_Coef(C, N, f, x);
    for(int k = 0; k <= N; k++)
    {
        for(int i = 1; i < N; i++)
        {
            res += C[i] * cos(M_PI*(i+0.5)*x[k]);
        }
        y[k] = res;
    }
    return;
}

void Get_Coef(  double *C, int N, double*f, double*x)
{
    for(int i = 0; i < N; i++)
    {
        C[i] = dot_f_phi(N,i, f, x) *2.0 / lambda_n(i,N);
    }
    return;
}

double lambda_n(int n, int N)
{
    double h = (1/(double)N);
    return 2*N*N*(1 - cos(M_PI*h*(n + 0.5)));
}

double dot_f_phi(int N, int j, double*f, double*x)
{
    double res = 0.0;
    double h = 1/(double)N;
    res += f[0] * cos(0.) /2.0;
    for(int i = 1; i < N; i++)
    {
       res += f[i] * cos(M_PI*(j+0.5)*(x[i]));
    }
    return res*h;
}



