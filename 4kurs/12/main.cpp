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
//plot '1.txt' using 1:2 with linespoints, '1.txt' using 1:3 with linespoints


void Generate_d(double* d, double* x, int N);
void Generate_x(double* MATRIX, int N);
void Generate_f(double* f, double* x, int N);
double func(double x);
double ans(double x);
double b(double x);
int Write(int N, double* y);
void Progonka(int N, double* alpha, double* beta, double* c, double* y, double*f);
void fourier(double* C, int N, double*y, double*x, double*f);
void Get_Coef(  double *C, int N, double*f, double*x);
double lambda_n(int n, int N);
double dot_f_phi(int N, int j, double*f, double*x);

int main(int argc, char* argv[])
{
    double *y,*x, *p, *q, *d, *f,*C;
    int N = -1; // число узлов
    int constant = -1;// 1 --- прогонка, 0--- Фурье

    if(!((argc == 3)&&
    (sscanf(argv[1],"%d",&N)==1)&&
    (sscanf(argv[2],"%d",&constant)==1)))
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
        Progonka(N, p, q, d, y, f);
        cout<<"Прогонка"<<endl;
    }
    if(constant == 0)
    {
        fourier(C,N,y,x,f);
        cout<<"Фурье"<<endl;
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
    ofstream out;
    out.open("1.txt");
    if(out.is_open())
    {
        out<<setprecision(15)<<fixed;
        for(int i = 0; i <=N; i++)
        {
            out<<setw(5)<<i/(double)N<<setw(25)<<y[i]<<setw(25)<<ans(i/(double)N)<<setw(25)<<abs(y[i] - ans(i/(double)N))<<endl;
        }
        out. close();
        return 1;
    }
    out.close();
    return -1;
}

void Generate_d(double* d, double* x, int N)
{
    for(int k = 1; k<N; k++)
    {
        d[k] = 2.*(N*N) + b(x[k]);
        //cout<<d[k]<< " "<<x[k]<<endl;
    }
    d[0] = N + (1./(double)(2*N))*b(x[0]);//1*N
    d[N] = 1;
    return;
}

void Generate_f(double* f, double* x, int N)
{
    for(int k = 0; k<=N; k++)
    {
        f[k] = func(x[k]);
    }
    return;
}

void Progonka(int N, double* alpha, double* beta, double* c, double* y, double*f)
{
    double a =(double) (N*N);
    double b = (double)(N*N);
    double b0 = N;
    double aN = 0.;
    f[0] = f[0]*0.5/(double)N;//0
    f[N] = 0;
    alpha[0] = 0;
    beta[0] = 0;

    alpha[1] =(b0/(double)c[0]);
    beta[1] = (f[0]/(double)c[0]);

    for(int i  = 1; i < N; i++)
    {
        alpha[i+1] = b/(double)(c[i] - a*alpha[i]);
        beta[i+1] = (f[i] + a*beta[i])/(double)(c[i] - a* alpha[i]);
        //cout<<alpha[i]<<" "<<beta[i]<<" "<<i<<endl;
    }

    //обратно

    y[N] = (f[N] + aN* beta[N])/(double)(c[N] - aN* alpha[N]);
    //cout<<y[N]<<" "<<endl;
    for(int i = N - 1; i > -1;i--)
    {
        y[i] = alpha[i+1]*y[i+1] + beta[i+1];
        //cout<<y[i]<<" "<<endl;
    }
    return;
}

double func(double x)
{
    //return 2 + b(x)*(-x*x +1);
    return cos(5*M_PI*x*0.5) *(25*M_PI*M_PI*0.25  + b(x));
}

double b(double x)
{

    x = x;
    return x;
    //return 1;
}
double ans(double x)
{
    //return -x*x +1;
    return cos(5*M_PI*x*0.5);
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
    Get_Coef(C, N, f, x);
    for(int k = 0; k <= N; k++)
    {
        y[k] = 0;
        for(int i = 0; i < N; i++)
        {
            y[k] += C[i] * cos(M_PI*(i+0.5)*x[k]);
        }
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
    return 2*N*N*(1 - cos(M_PI*h*(n + 0.5))) + b(h);
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



