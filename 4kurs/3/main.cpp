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
int find_p(long double* X,long double*Y,  long double* C, long double* C_2);
long double f(long double x, long double y); //функция
int Generate_x(long double* MATRIX, int N); //создаем узлы сетки
int Get_Coef(long double *C, long double *C_2, int N); //получаем коэфициенты с_m
long double dot_f_phi(int N, int i,  int j); //скалярное произведение f и phi_m
long double dot_c_phi(long double *C, int  N, int i,  int j); //скалярное произведение C_1 и phi_m
long double fourier(long double* C, int N, long double x, long double y); //значение ряда фурье в точке x
int Write(long double* X,long double* Y, int N, long double* C); //выведем в файл X, f(X), fourier(X), |f(X) - fourier(X)|

//splot '1.txt' using 1:2:3 with points, '1.txt' using 1:2:4 with points

int main()
{
    long double *X,*Y, *C_1, *C_2;
    int N = -1; // число узлов

    C_1 = (long double*) malloc(sizeof(long double)*(300*300));
    C_2 = (long double*) malloc(sizeof(long double)*(300*300));
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
    Get_Coef(C_1,C_2, N);
    Write(X,Y,N,C_2);

    find_p(X,Y,C_1, C_2);

    free(X);
    free(Y);
    free(C_1);
    free(C_2);

    return 0;
}

long double f(long double x, long double y)
{
    //return  (-x +1)* (-y +1);
    return (-(x*x)+1)*(-(y*y)+1);
    //return cos(5*M_PI*x/2.0)* cos(5*M_PI*y/2.0);
    //return cos(M_PI*x*(1+0.5))*cos(M_PI*y*(1+0.5));

    //if((x<=0.5)&&(x>= 0.4) &&(y<=0.5)&&(y>= 0.4))
    //{
    //    return 1;
    //}
    return 0;
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

long double fourier(long double* C, int N, long double x, long double y)
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

int Get_Coef(long double *C, long double *C_2, int N)
{
    for(int j = 0; j < N; j++)
    {   
        for(int i = 0; i < N; i++)
        {
            C[N*i+j] = dot_f_phi(N,i,j)*2.0;
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


long double dot_f_phi(int N, int i,  int j)
{
    long double res = 0.0;
    long double h = 1/(long double)N;
    res += f(0., j*h) * cos(0.) /2.0;
    for(int k = 1; k < N; k++)
    {
        res += f(k*h, j*h) * cos(M_PI*(i+0.5)*(k*h));
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



int Write(long double* X,long double* Y, int N, long double* C)
{
    ofstream out;
    out.open("1.txt");
    long double max = 0;
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
                    if(max < fabs(f(X[i]+l*h, Y[j]+l*h) - fourier(C,N,X[i]+l*h,Y[j]+l*h)))
                    {
                        max = fabs(f(X[i]+l*h, Y[j]+l*h) - fourier(C,N,X[i]+l*h,Y[j]+l*h));
                    }
                    out<<setw(25)<<X[i]+l*h<<setw(25)<<Y[j]+l*h<<setw(25)<<f(X[i]+l*h, Y[j]+l*h)<<setw(25)<<fourier(C,N,X[i]+l*h,Y[j]+l*h)<<setw(25)<<fabs(f(X[i]+l*h,Y[j]+l*h) - fourier(C,N,X[i]+l*h,Y[j]+l*h))<<endl;
                }

            }
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
        out. close();
        return max;
    }
    out.close();
    return -1;
}

int find_p(long double* X,long double*Y,  long double* C, long double* C_2)
{
    int NUM[4] = {10,20,50,70};
    long double a = 0.;
    long double b = 0.;
    long double c = 0.;
    long double d = 0.;
    ofstream out;
    out.open("2.txt");
    if(out.is_open())
    {
        for(int j = 0; j < 4; j++)
        {
            Generate_x(X, NUM[j]);
            Generate_x(Y, NUM[j]);
            Get_Coef(C,C_2, NUM[j]);
            long double max = 0;
            out<<setprecision(15)<<fixed;

            long double h = 1/(long double)NUM[j];
            h/=3.0;

            for(int i = 0; i <= NUM[j]; i++)
            {
                for(int k = 0; k <=NUM[j]; k++)
                {
                    for(int l = 0; l<= 2; l++)
                    {
                        if(max < fabs(f(X[i]+l*h,Y[k]+l*h) - fourier(C_2, NUM[j], X[i]+l*h,Y[k]+l*h)))
                        {
                            max = fabs(f(X[i]+l*h,Y[k]+l*h) - fourier(C_2, NUM[j], X[i]+l*h,Y[k]+l*h));
                        }
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

