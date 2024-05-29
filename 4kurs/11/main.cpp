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
    double *y, *x, *k1, *k2, *k3, *k4, *res, *tmp;
    double h = -0.1;//timestep
    int N = 0;//num of timesteps
    int M = 2;//num of equations
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
    h = ((M_PI)/(double)N);

    y = (double*) malloc(sizeof(double)*(N+1));
    x = (double*) malloc(sizeof(double)*(N+1));
    k1 = (double*) malloc(sizeof(double)*(N+1));
    k2 = (double*) malloc(sizeof(double)*(N+1));
    k3 = (double*) malloc(sizeof(double)*(N+1));
    k4 = (double*) malloc(sizeof(double)*(N+1));
    res = (double*) malloc(sizeof(double)*(N+1));
    tmp = (double*) malloc(sizeof(double)*(N+1));

    Generate_x(x, h, N);

    NU(y,M);

    ofstream out;
    out.open("1.txt");
    if(out.is_open())
    {
        out<<setprecision(15)<<fixed;
        for(int i = 1; i <=N; i++)
        {
            K_1(x[i],y,k1,h,M);
            K_2(x[i],y,k1,k2,h,M, tmp);
            K_3(x[i],y,k2,k3,h,M, tmp);
            K_4(x[i],y,k3,k4,h,M, tmp);

            y_n(y,k1,k2,k3,k4,M);
            ans(x[i], res, M);

            out<<setw(5)<<x[i];
            for(int j = 0; j<M; j++)
            {
                out<<setw(25)<<y[j];
            }   
            for(int j = 0; j<M; j++)
            {
                out<<setw(25)<<res[j];
            }     
            out<<setw(25)<<normVV(y, res,M);
            E(res,k1,k2,k3,k4,M);
            out<<setw(25)<<normV(res,M)<<endl;
        }
    }
    out.close();

    free(x);
    free(y);
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(res);
    free(tmp);
    return 0;
}


int f(double x, double* y, double* z, int M)
{
    if(M!= 2)
    {return -1;}
    //z[0] = -2*x;

    y[0] = y[0];//просто, чтобы использовать y, если он не нужен
    x = x;

    //z[0] = y[0]; //exp
    //z[0] = 5*x*x*x*x +1;
    //z[1] = 4*x*x*x + 2*x -1;
    //z[0] = 5*x*x*x*x;
    //z[1] = 4*x*x*x;
    z[0] = y[1];
    z[1] = -y[0];
    return 1;
}

int ans(double x,double* z, int M)
{
    if(M!= 2)
    {return -1;}
    //z[0] = -x*x + 1;
    //z[0] = exp(x);

    //z[0] = x*x*x*x*x + x;
    //z[1] = x*x*x*x + x*x - x;

    //z[0] = x*x*x*x*x;
    //z[1] = x*x*x*x +1;

    z[0] = sin(x);
    z[1] = cos(x);

    //z[0] = exp(x);
    //z[1] = exp(x);
    return 1;
}

int NU(double*z, int M)
{
    if(M!= 2)
    {return -1;}
    z[0] = 0;
    z[1] = 1;
    return 1;
}
void K_1(double x, double* y, double* k1, double h, double M)
{
    f(x, y,k1, M);
    for(int i = 0; i<M; i++)
    {
        k1[i] *= h;
    }
    return;
}
void K_2(double x, double* y,  double* k1,double* k2, double h, double M, double* tmp)
{
    for(int i = 0; i<M; i++)
    {
        tmp[i] = y[i]+0.5*k1[i];
    }
    f(x+0.5*h, tmp, k2, M);
    for(int i = 0; i<M; i++)
    {
        k2[i] *= h;
    }
    return;
}

void K_3(double x, double* y, double* k2, double* k3, double h, double M, double* tmp)
{
    for(int i = 0; i<M; i++)
    {
        tmp[i] = y[i]+0.5*k2[i];
    }
    f(x+0.5*h, tmp, k3, M);
    for(int i = 0; i<M; i++)
    {
        k3[i] *= h;
    }
    return;
}

void K_4(double x, double* y, double* k3, double* k4, double h, double M, double* tmp)
{
    for(int i = 0; i<M; i++)
    {
        tmp[i] = y[i]+k3[i];
    }
    f(x+h, tmp, k4, M);
    for(int i = 0; i<M; i++)
    {
        k4[i] *= h;
    }
    return;
}

void y_n(double*y, double* k1, double* k2, double* k3, double* k4, int M)
{
    for(int i = 0; i< M; i++)
    {
        y[i] = y[i] + (k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.;
    }
    return;
}

void E(double*z, double* k1, double* k2, double* k3, double* k4, int M)
{
    for(int i = 0; i< M; i++)
    {
        z[i] = 2.*(k1[i] - k2[i] - k3[i] + k4[i])/3.;
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



