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
int find_p();
long double f(long double x); //функция
int Generate_x(long double* X, int N); //создаем узлы сетки
int Get_Coef(long double *C, int N); //получаем коэфициенты с_m
long double dot_f_phi(int N, int j); //скалярное произведение f и phi_m
long double dot_phi_phi(int N, int j); // скалярное произведение phi_m и phi_m
long double fourier(long double* C, int N, long double x); //значение ряда фурье в точке x
int Write(long double* X, int N, long double* C); //выведем в файл X, f(X), fourier(X), |f(X) - fourier(X)|

//plot '1.txt' using 1:2 with linespoints, '1.txt' using 1:3 with lines

int main(int argc, char* argv[])
{
    long double *X, *C;
    int N = -1; // число узлов

    if(!(sscanf(argv[1],"%d",&N)==1))
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
    
    C = (long double*) malloc(sizeof(long double)*(N));
    X = (long double*) malloc(sizeof(long double)*(N+1));

    Generate_x(X, N);
    Get_Coef(C, N);
    Write(X, N, C);

    free(X);
    free(C);

    return 0;
}


int Generate_x(long double* X, int N)
{
    for(int i = 0; i < N; i++)
    {
        X[i] = i*(1/(long double)N) ;
    }
    return 1;
}

int Get_Coef(long double *C, int N)
{
    C[0] = 0;
    for(int i = 1; i < N; i++)
    {
        C[i] = dot_f_phi(N,i)*2.0;//;  /dot_phi_phi(N,i)/2.0
    }
    return 1;
}

long double dot_f_phi(int N, int j)
{
    long double res = 0.0;
    long double h = 1/(long double)N;
    res += f(0.) * cos(0.) /2.0;
    for(int i = 1; i < N; i++)
    {
       res += f(i*h) * cos(M_PI*(j-0.5)*(i*h));
    }
    return res*h;
}

long double dot_phi_phi(int N, int j)
{
    long double res = 0.0;
    long double h = 1/(long double)N;
    res += cos(0.) * cos(0.) /2.0;
    for(int i = 1; i < N; i++)
    {
        res += cos(M_PI*(j-0.5)*(i*h )) * cos(M_PI*(j-0.5)*(i*h));
    }
    return res*h;
}

long double fourier(long double* C, int N, long double x)
{
    long double res = 0.0;
    for(int i = 1; i < N; i++)
    {
        res += C[i] * cos(M_PI*(i-0.5)*x);
    }
    return res;
}


long double f(long double x)
{
    //return  -x +1;
    //return -(x*x)+1;
    //return cos(5*M_PI*x/2.0);
    return cos(M_PI*x*(1-0.5));
}

int Write(long double* X, int N, long double* C)
{
    ofstream out;
    out.open("1.txt");
    long double max = 0;
    if(out.is_open())
    {
        out<<setprecision(15)<<fixed;
        for(int i = 0; i < N; i++)
        {
            if(max < fabs(f(X[i]) - fourier(C,N,X[i])))
            {
                max = fabs(f(X[i]) - fourier(C,N,X[i]));
            }
            cout<<C[i]<<endl;
            out<<setw(25)<<X[i]<<setw(25)<<f(X[i])<<setw(25)<<fourier(C,N,X[i])<<setw(25)<<fabs(f(X[i]) - fourier(C,N,X[i]))<<endl;
        }
        cout<<setw(25)<<N<<setw(25)<<log((long double)N)<<setw(25)<<log(1/max)<<endl;
        out. close();
        return max;
    }
    out.close();
    return -1;
}

int find_p()
{
    ofstream out;
    out.open("2.txt");
    if(out.is_open())
    {
        out<<setprecision(15)<<fixed;
        for(int i = 0; i < 5; i++)
        {
            out<<0.00499905/0.00249988<<endl;
        }
        out.close();
        return 0;
    }
    out. close();
    return -1;
}

