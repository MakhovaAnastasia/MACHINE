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

long double f(long double t,long double x, long double y); //функция
long double ans(long double t,long double x, long double y);
int Generate_x(long double* MATRIX, int N); //создаем узлы сетки
int Get_Coef(long double *C, long double *C_2, int N, long double t); //получаем коэфициенты с_m
long double dot_f_phi(int N, int i,  int j, long double t); //скалярное произведение f и phi_m
long double dot_c_phi(long double *C, int  N, int i,  int j); //скалярное произведение C_1 и phi_m
long double fourier(long double* C, int N, long double x, long double y, long double t); //значение ряда фурье в точке x
int Write(long double* X,long double* Y, int N, long double* C,long double* C_2, int T); //выведем в файл X, f(X), fourier(X), |f(X) - fourier(X)|

//splot '1.txt' using 1:2:3 with points, '1.txt' using 1:2:4 with points

/*
reset
set term gif animate delay 20
set output "animate.gif"
FILE(i) = sprintf("%d.txt",i)
set xrange [0:1]
set yrange [0:1]
set zrange [-1:1]
do for [i = 0:9]{
    splot FILE(i) using 1:2:3 with points, FILE(i) using 1:2:4 with points
}
*/

int main(int argc, char* argv[])
{
    long double *X,*Y, *C_1, *C_2;
    int N = -1; // число узлов
    int T = -1;
    
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
        cout<<" N >= 0"<<endl;
        return -1;
    }
    C_1 = (long double*) malloc(sizeof(long double)*(N*N));
    C_2 = (long double*) malloc(sizeof(long double)*(N*N));
    X = (long double*) malloc(sizeof(long double)*(N+1));
    Y = (long double*) malloc(sizeof(long double)*(N+1));

    Generate_x(X, N);
    Generate_x(Y, N);
    Write(X,Y,N,C_1,C_2, T);


    free(X);
    free(Y);
    free(C_1);
    free(C_2);

    return 0;
}

long double f(long double t,long double x, long double y)
{
    
    return M_PI*M_PI*(2.5)*(2.5)*ans(t,x,y);
}
long double ans(long double t,long double x, long double y)
{
    
    return cos(5*M_PI*x/2.0)* cos(5*M_PI*y/2.0)*exp(-M_PI*M_PI*(2.5)*(2.5)*t*0.001);
    //return cos(M_PI*x*(0.5))*cos(M_PI*0.5*y)*exp(-M_PI*M_PI*0.25*t);
}

int Generate_x(long double* MATRIX, int N)
{
    for(int i = 0; i < N; i++)
    {
        MATRIX[i] = i*(1/(long double)N) ;
    }
    MATRIX[N] = 1.0;
    return 1;
}

long double fourier(long double* C, int N, long double x, long double y, long double t)
{
    long double res = 0.0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j< N; j++)
        {
            res += C[i*N + j] * cos(M_PI*(i+0.5)*x)* cos(M_PI*(j+0.5)*y);
        }
    }
    return res;
}

int Get_Coef(long double *C, long double *C_2, int N, long double t)
{
    for(int j = 0; j < N; j++)
    {   
        for(int i = 0; i < N; i++)
        {
            C[N*i+j] = dot_f_phi(N,i,j,t)*2.0;
        }
    }

    for(int j = 0; j < N; j++)
    {   
        for(int i = 0; i < N; i++)
        {
            C_2[N*i+j] = dot_c_phi(C,N,i,j)*2.0;
        }
    }

    return 1;
}


long double dot_f_phi(int N, int i,  int j, long double t)
{
    long double res = 0.0;
    long double h = 1/(long double)N;
    res += f(t, 0., j*h) * cos(0.) /2.0;
    for(int k = 1; k < N; k++)
    {
        res += f(t,k*h, j*h) * cos(M_PI*(i+0.5)*(k*h));
    } 
    return res*h;
}

long double dot_c_phi(long double *C, int  N, int i,  int j)
{
    long double res = 0.0;
    long double h = 1/(long double)N;
    res += C[i*N] * cos(0.) /2.0;
    for(int k = 1; k < N; k++)
    {
        res += C[i*N+k] * cos(M_PI*(j+0.5)*(k*h));
    } 
    return res*h;
}



int Write(long double* X,long double* Y, int N, long double* C,long double* C_2, int T)
{
    for(int k = 0; k<T; k++)
    {
        Get_Coef(C,C_2,N,k);
        ofstream out;
        out.open(to_string(k)+".txt");
        if(out.is_open())
        {
            long double h = 1/(long double)N;
            h/=3.0;
            out<<setprecision(15)<<fixed;
            for(int i = 0; i <= N; i++)
            {
                for(int j = 0; j<= N; j++)
                {
                    for(int l = 0; l<= 2; l++)
                    {
                        //if(max < fabs(f(k,X[i]+l*h, Y[j]+l*h) - fourier(C,N,X[i]+l*h,Y[j]+l*h)))
                        //{
                        //    max = fabs(f(k,X[i]+l*h, Y[j]+l*h) - fourier(C,N,X[i]+l*h,Y[j]+l*h));
                        //}
                        out<<setw(25)<<X[i]+l*h<<setw(25)<<Y[j]+l*h<<setw(25)<<f(k,X[i]+l*h, Y[j]+l*h)<<setw(25)<<fourier(C_2,N,X[i]+l*h,Y[j]+l*h,k)<<setw(25)<<fabs(f(k,X[i]+l*h,Y[j]+l*h) - fourier(C_2,N,X[i]+l*h,Y[j]+l*h, k))<<endl;
                    }

               }
            }
            out.close();
        }
    /*
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                cout<<C[i*N+j]<<" ";
            }
            cout<<endl;
        }
    */
        
    }
    return 0;
}

