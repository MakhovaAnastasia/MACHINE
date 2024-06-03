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

//splot '1.txt' using 1:2:3 with points, '1.txt' using 1:2:4 with points
//splot '2.txt' using 1:2:3 with points, '2.txt' using 1:2:4 with points

double u_0(double x);
void Yavno(double** y, double h, int N, double t);
double b(double x);
double f(double x, double y);
double ans(double x,double y);
void Generate_d(double* d, int N, double t, double h);
void Generate_f(double* F, double** y, int N, int n, double h, double t);
void Progonka(int N, double* alpha, double* beta, double* c, double** y, double*f, double h, double t);



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
    t = ((1)/(double)N);

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
        for(int n = 0; n<=N; n++)
        {
            for(int m = 0; m<=N; m++)
            {
                //out<<setw(25)<<m*h<<setw(25)<<y[1][m]<<setw(25)<<ans(1*t,m*h)<<endl;
                out<<setw(25)<<m*h<<setw(25)<<n*t<<setw(35)<<y[n][m]<<setw(35)<<ans(n*t,m*h)<<endl;
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
        for(int n = 0; n<=N; n++)
        {
            for(int m = 0; m<=N; m++)
            {
                out2<<setw(25)<<m*h<<setw(25)<<n*t<<setw(35)<<y[n][m]<<setw(35)<<ans(n*t,m*h)<<endl;
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
        y[n+1][0] = t*h*0.5*f(n*t,0.) + y[n][0] *(-b(0.)*t*h*0.5 +1 - t/(h)) + y[n][1]*t/(h);

        //right
        y[n+1][N] = 0;

        //other
        for(int m = 1; m < N; m++)
        {
            y[n+1][m] = t*f(n*t, m*h) + y[n][m] *(1 - t*b(m*h) - 2.*t/(h*h)) + (t/(h*h))*(y[n][m+1] + y[n][m-1]);
        }
        
    }
}

double u_0(double x)
{
    double k = 1;
    return cos(M_PI*(k+0.5)*x);
}

double b(double x)
{
    return 2*x;
}

double f(double x, double y)
{
    return b(y)*ans(x,y);
}

double ans(double x,double y)
{
    double k = 1;
    return cos(M_PI*(k+0.5)*y)*exp(-M_PI*M_PI*(k+0.5)*(k+0.5)*x);
}

void Generate_d(double* d, int N, double t, double h)
{
    for(int k = 1; k<N; k++)
    {
        d[k] = t*b(k*h) + 2.*(double)(N*N)*t + 1;
        //cout<<d[k]<<" "<<k<<endl;
    }
    d[0] = -(double)N + 0.5*h*b(0.) + 0.5*h/t;
    d[N] = 1;
    return;
}

void Generate_f(double* F, double** y, int N, int n, double h, double t)
{
    for(int m = 0; m<N; m++)
    {
        F[m] = t*f((n+1)*t, m*h) + y[n][m];
    }
    F[0] = -(f((n+1)*t,0.) + y[n][0]/t)*0.5*h;
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

    //nu
    for(int m = 0; m<= N; m++)
    {
        y[0][m] = u_0(m*h);
    }

    Generate_d(c,N, t, h);
    for(int n = 0; n<N; n++)
    {
       
        Generate_f(f,y,N,n,h, t);
        alpha[1] =(b0/(double)c[0]);
        beta[1] = (f[0]/(double)c[0]);
        
        for(int i  = 1; i < N; i++)
        {
            alpha[i+1] = b/(double)(c[i] - a*alpha[i]);
            beta[i+1] = (f[i] + a*beta[i])/(double)(c[i] - a* alpha[i]);
        }

        //обратно

        y[n+1][N] = (f[N] + aN* beta[N])/(double)(c[N] - aN* alpha[N]);
    
        for(int i = N - 1; i > -1;i--)
        {
            y[n+1][i] = alpha[i+1]*y[n+1][i+1] + beta[i+1];
        
        }
    }
    return;
}





