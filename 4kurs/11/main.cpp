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
void Progonka(int N, double* p, double* q, double* d, double* y, double*x, double*f);
void fourier(double* C, int N, double*y, double*x, double*f);
void Get_Coef(  double *C, int N, double*f, double*x);
double lambda_n(int n, int N);
double dot_f_phi(int N, int j, double*f, double*x);

//splot '1.txt' using 1:2:3 with points, '1.txt' using 1:2:4 with points

int main(int argc, char* argv[])
{
    double *y,*x, *k1, *k2, *k3, *k4;
    int N = -1; // число узлов
    int constant = -1;// 1 --- прогонка, 0--- Фурье

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
        cout<<" N >= 0"<<endl;
        return -1;
    }
    
    y = (double*) malloc(sizeof(double)*(N+1));
    x = (double*) malloc(sizeof(double)*(N+1));
    k1 = (double*) malloc(sizeof(double)*(N+1));
    k2 = (double*) malloc(sizeof(double)*(N+1));
    k3 = (double*) malloc(sizeof(double)*(N+1));
    k4 = (double*) malloc(sizeof(double)*(N+1));
    Generate_x(x, N);
    Generate_f(f, x, N);

    if(constant == 1)
    {
        Generate_d(d, x, N);
        Progonka(N, p, q, d, y, x,f);
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
    free(k1);
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
            out<<setw(5)<<x[i]<<setw(25)<<y[i]<<setw(25)<<ans(i/(double)N)<<setw(25)<<abs(y[i] - ans(i/(double)N))<<endl;
        }
        out. close();
        return 1;
    }
    out.close();
    return -1;
}

double f(double x, double y)
{
    return 1;
}

void K1(int N,double*x, double* y, double* k1)
{
    for(int i = 0; i < N; i++)
    {
        k1[i] = f(x[i], y[i]);
    }
    return;
}





