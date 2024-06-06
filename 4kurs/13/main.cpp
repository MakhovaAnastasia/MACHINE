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
void Yavno(double** y, double h, int N, double t, int T);
double b(double x);
double f(double t, double x);
double ans(double t,double x);
void Generate_d(double* d, int N, double t, double h);
void Generate_f(double* F, double** y, int N, int n, double h, double t);
void Progonka(int N, double* alpha, double* beta, double* c, double** y, double*f, double h, double t, int T);



int main(int argc, char* argv[])
{
    double **y, *F, *alpha, *beta, *c;

    double h = -0.1;//step
    double t = -0.1;//timestep
    int N = 0;//num of steps
    int T = 0;
    if(!((argc == 3)&&
    (sscanf(argv[1],"%d",&N)==1)&&
    (sscanf(argv[2],"%d",&T)==1)))
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
    h = ((1.)/(double)N);
    t = ((1.)/(double)T);

    y = (double**) malloc(sizeof(double*)*(T+1));
    for(int i = 0; i<= T; i++)
    {
        y[i] = (double*) malloc(sizeof(double)*(N+1));
    }

    F = (double*) malloc(sizeof(double)*(N+1));
    alpha = (double*) malloc(sizeof(double)*(N+1));
    beta = (double*) malloc(sizeof(double)*(N+1));
    c = (double*) malloc(sizeof(double)*(N+1));

     double max  = -1;
    if(t<(h*h)*0.5)
    {
        cout<<"Явная схема устойчива"<<endl;
        Yavno(y,h,N,t, T);

        ofstream out;
        out.open("1.txt");
        if(out.is_open())
        {
            out<<setprecision(15)<<fixed;
            for(int n = 0; n<=T; n++)
            {
                for(int m = 0; m<=N; m++)
                {
                    //out<<setw(25)<<m*h<<setw(25)<<y[1][m]<<setw(25)<<ans(1*t,m*h)<<endl;
                    out<<setw(25)<<m*h<<setw(25)<<n*t<<setw(35)<<y[n][m]<<setw(35)<<ans(n*t,m*h)<<endl;
                    if(fabs(y[n][m] - ans(n*t,m*h)) > max)
                    {
                        max = fabs(y[n][m] - ans(n*t,m*h));
                    }
                }
                out<<endl;
            }
        }
        out.close();

    }
    else {
        cout<<"Явная схема не устойчива"<<endl;
    }
    cout<<"Явно: "<<max<<endl;
    max  = -1;
    Progonka(N,alpha,beta,c,y,F,h,t, T);
    ofstream out2;
    out2.open("2.txt");
    if(out2.is_open())
    {
        out2<<setprecision(15)<<fixed;
        for(int n = 0; n<=T; n++)
        {
            for(int m = 0; m<=N; m++)
            {
                out2<<setw(25)<<m*h<<setw(25)<<n*t<<setw(35)<<y[n][m]<<setw(35)<<ans(n*t,m*h)<<endl;
                if(fabs(y[n][m] - ans(n*t,m*h)) > max)
                {
                    max = fabs(y[n][m] - ans(n*t,m*h));
                }
            }
            out2<<endl;
        }
    }
    out2.close();

    cout<<"Неявно: "<<max<<endl;

    for(int i = 0; i<=T; i++)
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

void Yavno(double** y, double h, int N, double t, int T)
{
    //n==0
    for(int m = 0; m<= N; m++)
    {
        y[0][m] = u_0(h*m);
    }
    //n >=1
    for(int  n = 0; n<T; n++)
    {
        //left
        y[n+1][0] =t*f(n*t,0.) + (y[n][1] -y[n][0])*2.*t/(h*h) + y[n][0]*(1- b(0.)*t);

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
    //double k = 2;
    //return cos(M_PI*(k+0.5)*x);
    return (x-1.)*exp(x);
}

double b(double x)
{
    x = x;
    return 1;
    //return x;
}

double f(double t, double x)
{
    //return b(y)*ans(x,y);
    return (-((x-1.)*exp(x))*exp(-t) - exp(x)*(x +1.)*exp(-t)) +1.*((x-1.)*exp(x))*exp(-t);
}

double ans(double t,double x)
{
    //double k = 2;
    //return cos(M_PI*(k+0.5)*y)*exp(-M_PI*M_PI*(k+0.5)*(k+0.5)*x);
    return ((x-1.)*exp(x))*exp(-t);
}

void Generate_d(double* d, int N, double t, double h)
{
    for(int k = 1; k<N; k++)
    {
        d[k] = b(k*h) + 2./(h*h)+ 1/t;
        //cout<<d[k]<<" "<<k<<endl;
    }
    d[0] = 2./(h*h) + 1./t +b(0.);
    d[N] = 1.;
    return;
}

void Generate_f(double* F, double** y, int N, int n, double h, double t)
{
    for(int m = 0; m<N; m++)
    {
        F[m] = f((n+1)*t, m*h) + y[n][m]/t;
    }
    F[0] = f((n+1)*t,0. ) + y[n][0]/t;
    F[N] = 0;
    return;
}

void Progonka(int N, double* alpha, double* beta, double* c, double** y, double*f, double h, double t, int T)
{
    double a = 1./(h*h);
    double b = 1./(h*h);
    double b0 = 2./(h*h);
    double aN = 0.;


    alpha[0] = 0;
    beta[0] = 0;

    //nu
    for(int m = 0; m<= N; m++)
    {
        y[0][m] = u_0(m*h);
    }

    Generate_d(c,N, t, h);
    for(int n = 0; n<T; n++)
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
        for(int i = N - 1; i > -1;i--)
        {
            y[n+1][i] = alpha[i+1]*y[n+1][i+1] + beta[i+1];

        }
        double nevyaz = 0.;
        for(int i = 0; i<=N; i++)
        {
            if(i== 0)
            {
                nevyaz += y[n+1][0]*c[0] -b0*y[n+1][1] - f[0];
                continue;
            }
            if(i == N)
            {
                nevyaz += y[n+1][N]*c[N] - aN*y[n+1][N-1] - f[N];
                continue;
            }
            nevyaz += y[n+1][i]*c[i] -b*y[n+1][i+1] - a*y[n+1][i-1] - f[i];


        }
        //cout<< "Невязка = "<<nevyaz<<endl;
    }

    return;
}







