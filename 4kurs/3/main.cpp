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
int find_p(long double* X,long double*Y,  long double* C);
long double f(long double x, long double y); //функция
int Generate_x(long double* MATRIX, int N); //создаем узлы сетки
int Get_Coef(long double *C, int N); //получаем коэфициенты с_m
long double dot_f_phi(int N, int i,  int j); //скалярное произведение f и phi_m
long double fourier(long double* C, int N, long double x, long double y); //значение ряда фурье в точке x
int Write(long double* X,long double* Y, int N, long double* C); //выведем в файл X, f(X), fourier(X), |f(X) - fourier(X)|

//plot '1.txt' using 1:2 with linespoints, '1.txt' using 1:3 with lines

int main()
{
    long double *X,*Y, *C;
    int N = -1; // число узлов

    C = (long double*) malloc(sizeof(long double)*(300*300));
    X = (long double*) malloc(sizeof(long double)*(301));
    Y = (long double*) malloc(sizeof(long double)*(301));

    if(!(scanf("%d",&N)==1))
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
    
    Generate_x(X, N);
    Generate_x(Y, N);
    Get_Coef(C, N);
    Write(X,Y,N,C);

    //find_p(X,Y,C);

    free(X);
    free(Y);
    free(C);

    return 0;
}

long double f(long double x, long double y)
{
    //return  (-x +1)* (-y +1);
    //return -(x*x)+1;
    //return cos(5*M_PI*x/2.0)* cos(5*M_PI*y/2.0);
    return cos(M_PI*x*(1-0.5))* cos(M_PI*y*(1-0.5));

    //if((x<=0.50001)&&(x>= 0.499998) &&(y<=0.50001)&&(y>= 0.499998))
    //{
    //    return 1;
    //}
    //return 0;
}

int Generate_x(long double* MATRIX, int N)
{
    for(int i = 0; i < N; i++)
    {
        MATRIX[i] = i*(1/(long double)N) ;
    }
    return 1;
}

long double fourier(long double* C, int N, long double x, long double y)
{
    long double res = 0.0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j< N; j++)
        {
            res += C[i*N + j] * cos(M_PI*(i-0.5)*x)* cos(M_PI*(j-0.5)*y);
        }
    }
    return res;
}

int Get_Coef(long double *C, int N)
{
    long double h = 1/(long double)N;
    for(int i = 1; i < N; i++)
    {   
        C[N*i+0] = f(i*h, 0.);
        for(int j = 1; j < N; j++)
        {
            C[N*i+j] = dot_f_phi(N,i,j)*2.0;//;  /dot_phi_phi(N,i)/2.0
            C[N*i+0]-= C[N*i+j];
        }
    }
    return 1;
}

long double dot_f_phi(int N, int i,  int j)
{
    long double res = 0.0;
    long double h = 1/(long double)N;
    res += f(0., 0.) * cos(0.)* cos(0.)/2.0;
    for(int m = 0; m < N; m++)
    {
        for(int k = 0; k < N; k++)
        {
            res += f(i*h, j*h) * cos(M_PI*(m-0.5)*(i*h))* cos(M_PI*(k-0.5)*(j*h));
            if((k == 0) || (m == 0 ))
            {
                res /= 2.0;
            }
        } 
    }
    return res*h*h;
}



int Write(long double* X,long double* Y, int N, long double* C)
{
    ofstream out;
    out.open("1.txt");
    long double max = 0;
    if(out.is_open())
    {
        out<<setprecision(15)<<fixed;
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j< N; j++)
            {
                if(max < fabs(f(X[i], Y[j]) - fourier(C,N,X[i],Y[j])))
                {
                    max = fabs(f(X[i], Y[j]) - fourier(C,N,X[i],Y[j]));
                }
                cout<<C[i*N+j]<<endl;
                out<<setw(25)<<X[i]<<setw(25)<<Y[j]<<setw(25)<<f(X[i], Y[j])<<setw(25)<<fourier(C,N,X[i],Y[j])<<setw(25)<<fabs(f(X[i],Y[j]) - fourier(C,N,X[i],Y[j]))<<endl;
            }
        }
        out. close();
        return max;
    }
    out.close();
    return -1;
}

int find_p(long double* X,long double*Y,  long double* C)
{
    int NUM[6] = {10,50,100,150,200,300};
    long double a = 0.;
    long double b = 0.;
    long double c = 0.;
    long double d = 0.;
    ofstream out;
    out.open("2.txt");
    if(out.is_open())
    {
        for(int j = 0; j < 6; j++)
        {
            Generate_x(X, NUM[j]);
            Generate_x(Y, NUM[j]);
            Get_Coef(C, NUM[j]);
            long double max = 0;
            out<<setprecision(15)<<fixed;

            for(int i = 0; i < NUM[j]; i++)
            {
                for(int k = 0; k < NUM[j]; k++)
                {
                    if(max < fabs(f(X[i],Y[k]) - fourier(C, NUM[j], X[i],Y[k])))
                    {
                        max = fabs(f(X[i],Y[k]) - fourier(C, NUM[j], X[i],Y[k]));
                    }
                }
            }

            if(j==1)
            {
                a = log((long double)NUM[j]);
                b = log(1/max);
            }
            if(j==2)
            {
                c = log((long double)NUM[j]);
                d = log(1/max);
            }
            out<<setw(25)<<NUM[j]<<setw(25)<<log((long double)NUM[j])<<setw(25)<<log(1/max)<<endl;
        }
        cout<<"p = "<<(long double)(d - b)/(long double)(c - a)<<endl;
        out.close();
        return 0;
    }
    out. close();
    return -1;
}

