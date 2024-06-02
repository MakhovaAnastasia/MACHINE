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

//plot '1.txt' using 1:2 with linespoints, '1.txt' using 1:3 with linespoints


int f(double x, double* y, double* z, int M);
int ans(double x,double* z, int M);
int NU(double*z, int M);
void K_1(double x, double* y, double* k1, double h, double M);
void K_2(double x, double* y,  double* k1,double* k2, double h, double M, double*tmp);
void K_3(double x, double* y, double* k2, double* k3, double h, double M, double*tmp);
void K_4(double x, double* y, double* k3, double* k4, double h, double M, double*tmp);
void y_n(double*y, double* k1, double* k2, double* k3, double* k4, int M);
void E(double*z, double* k1, double* k2, double* k3, double* k4, int M);
void Generate_x(double* MATRIX, double h, int N);
double normVV(double*A, double*B, int M);
double normV(double*A,  int M);



int main(int argc, char* argv[])
{
    double **y, *F, *alpha, *beta, *c;

    double h = -0.1;//timestep
    double t = -0.1;//timestep
    int N = 0;//num of timesteps
    if(!((argc == 2)&&
    (sscanf(argv[1],"%d",&N)==1)))
    {
        //ошибка чтения
        return -1;
    }

    if(N<=0)
    {
        //ошибка 
        cout<<N<<endl;
        cout<<" N>= 0"<<endl;
        return -1;
    }
    h = ((1)/(double)N);
    t = h = ((1)/(double)N);

    y = (double**) malloc(sizeof(double*)*(N+1));
    for(int i = 0; i<= N; i++)
    {
        y[i] = (double*) malloc(sizeof(double)*(N+1));
    }

    F = (double*) malloc(sizeof(double)*(N+1));
    alpha = (double*) malloc(sizeof(double)*(N+1));
    beta = (double*) malloc(sizeof(double)*(N+1));
    c = (double*) malloc(sizeof(double)*(N+1));

    Yavno(y,h,N,t);
    ofstream out;
    out.open("1.txt");
    if(out.is_open())
    {
        out<<setprecision(15)<<fixed;
        for(int m = 0; m <=N; m++)
        {
            for(int n = 0; n<=N; n++)
            {
                out<<setw(25)<< y[n][m];
            }
            out<<endl;
        }
    }
    out.close();

    Progonka(N,alpha,beta,c,y,F,h,t);
    ofstream out2;
    out2.open("2.txt");
    if(out2.is_open())
    {
        out2<<setprecision(15)<<fixed;
        for(int m = 0; m <=N; m++)
        {
            for(int n = 0; n<=N; n++)
            {
                out2<<setw(25)<< y[n][m];
            }
            out2<<endl;
        }
    }
    out2.close();



    for(int i = 0; i<=N; i++)
    {
        free(y[i]);
    }
    free(y);
    free(F);
    free(c);
    free(alpha);
    free(beta);

    return 0;
}

double u_0(double x)
{
    return 1-x;
}

void Yavno(double** y, double h, int N, double t)
{
    //n==0
    for(int m = 0; m<= N; m++)
    {
        y[0][m] = u_0(h*m);
    }
    //n >=1
    for(int  n = 0; n<N; n++)
    {
        //left
        y[n+1][0] = y[n][0] - t*(f(n*t,0.) - b(0.)*y[n][0] - (2./(double)(h*h))*(y[n][1] - y[n][0]));

        //right
        y[n][N] = 0;

        //other
        for(int m = 1; m < N; m++)
        {
            y[n+1][m] = t*f(n*t, m*h) + y[n][m] *(1 - t*b(m*h) - 2.*t/(h*h)) + (t/(h*h))*(y[n][m+1] + y[n][m-1]);
        }
        
    }
}

double b(double x)
{
    return x;
}

double f(double x, double y)
{

    return 1;
}

int ans(double x,double* z, int M)
{
   
    return 1;
}

void Generate_d(double* d, int N, double t, double h)
{
    for(int k = 1; k<N; k++)
    {
        d[k] = 1 - 2.*(N*N)*t;
    }
    d[0] = -N + 0.5*h*b(0) + 1/t;
    d[N] = 1;
    return;
}

void Generate_f(double* F, double** y, int N, int n, double h, double t)
{
    for(int m = 0; m<=N; m++)
    {
        F[m] = t*f(n*t, m*h) + y[n][m]*(1 - b(m*h)*t);
    }
    F[0] = (F[0] + y[n][0]/t)*0.5*h;
    F[N] = 0;
    return;
}

void Progonka(int N, double* alpha, double* beta, double* c, double** y, double*f, double h, double t)
{
    double a = t*(double) (N*N);
    double b = t*(double)(N*N);
    double b0 = -N;
    double aN = 0.;

    alpha[0] = 0;
    beta[0] = 0;

    for(int m = 0; m<= N; m++)
    {
        y[0][m] = u_0(m*h);
    }
    for(int n = 0; n<N; n++)
    {
        Generate_d(c,N, t, h);
        Generate_f(f,y,N,n,h, t);
        alpha[1] =(b0/(double)c[0]);
        beta[1] = (f[0]/(double)c[0]);

        for(int i  = 1; i < N; i++)
        {
            alpha[i+1] = b/(double)(c[i] - a*alpha[i]);
            beta[i+1] = (f[i] + a*beta[i])/(double)(c[i] - a* alpha[i]);
        //cout<<alpha[i]<<" "<<beta[i]<<" "<<i<<endl;
        }

        //обратно

        y[n+1][N] = (f[N] + aN* beta[N])/(double)(c[N] - aN* alpha[N]);
        //cout<<y[N]<<" "<<endl;
        for(int i = N - 1; i > -1;i--)
        {
            y[n+1][i] = alpha[i+1]*y[n+1][i+1] + beta[i+1];
            //cout<<y[i]<<" "<<endl;
        }
    }
    return;
}


void Generate_x(double* MATRIX, double h, int N)
{
    for(int i = 0; i <= N; i++)
    {
        MATRIX[i] = i*h ;
    }
    return;
}

double normVV(double*A, double*B, int M)
{
    double max = -1;
    for(int i = 0; i < M; i++)
    {
        if(abs(A[i] - B[i]) > max)
        {
            max = abs(A[i] - B[i]);
        }
    }
    return max;
}
double normV(double*A,  int M)
{
    double max = -1;
    for(int i = 0; i<M; i++)
    {
        if(abs(A[i] ) > max)
        {
            max = abs(A[i]);
        }
    }
    return max;
}



