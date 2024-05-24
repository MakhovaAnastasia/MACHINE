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


int no(double A, int n, double* y, int num, double EPS);
int no1(double A, int n, double* y,double EPS);
int no2(double A, int n, double* y, double EPS);
int no3(double A, int n, double* y, double EPS);
int no4(double A, int n, double* y);
int no5(double A, int n, double* y,double EPS);
int no6(double A, int n, double* y);
double En(double A, int n, double* y, double EPSe);

int Write(double* y);

int main(void)
{
    int A[3] = {1,10,1000};
    int N[4] = {1, 2, 3, 6};
    int m[6] = {1,1,2,2,2,1};
    double EPS = 1e-20;
    double EPSe = -700;
    double* y;
    y = (double*) malloc(sizeof(double) * 1000000);
    cout<<setw(15)<<"# "<<setw(15)<<"E1 "<<setw(15)<<"E2 "<<setw(15)<<"E3 "<<setw(15)<<"E6 "<<setw(15)<<"m "<<setw(15)<<"A "<<endl;
    for(int i = 1;i <=6; i++)
    {
        cout<<"-------------------------------------------------------------------------------------------------------------"<<endl;
        for(int k = 0; k < 3; k++)//Ak
        {
            cout<<setw(15)<<i;
            for(int j = 0; j < 4; j++)//Ej
            {
                no(A[k], N[j], y, i, EPS);
                cout<<setw(15)<<En(A[k], N[j], y, EPSe);

            }
            cout<<setw(15)<<m[i-1]<<setw(15)<<A[k]<<endl;
        }
    }
    Write(y);
    free(y);
    return 0;
}

double En(double A, int n, double* y, double EPSe)
{
    double h = pow(10, -n);
    int N = (int)(1/h);
    double max = -1;
    for(int k = 0; k< N; k++)
    {
           if(-A*h*k < EPSe)
           {
               return -k;
           }
        if(fabs(y[k] - exp(-A*h*k)) > max)
        {
            max = abs(y[k] - exp(-A*h*k));
        }
    }
    return max;
}

int Write(double* y)
{
    ofstream out;
    out.open("1.txt");
    if(out.is_open())
    {
        out<<setprecision(15)<<fixed;
        no4(10, 4, y);
        for(int i = 0; i <pow(10, 4); i++)
        {
            out<<i<<" "<<y[i]<<endl;
        }
        out. close();
        return 1;
    }
    out.close();
    return -1;
}

int no(double A, int n, double* y, int num, double EPS)
{
    switch (num)
    {
    case 1:
        no1(A, n, y, EPS);
        break;
    case 2:
        no2(A, n, y, EPS);
        break;
    case 3:
        no3(A, n, y, EPS);
        break;
    case 4:
        no4(A, n, y);
        break;
    case 5:
        no5(A, n, y, EPS);
        break;
    case 6:
        no6(A, n, y);
        break;
    default:
        break;
    }
    return 1;
}

int no1(double A, int n, double* y, double EPS)
{
    double h = pow(10, -n);
    int N = (int)(1/h);
    y[0] = 1.;
    for(int k = 1; k < N; k++)
    {
        y[k] = y[k - 1]*(1 - A*h);
        if(fabs(y[k])< EPS)
        {
            y[k] = 0.;
        }
    }
    return 1;
}

int no2(double A, int n, double* y, double EPS)
{
    double h = pow(10, -n);
    int N = (int)(1/h);
    y[0] = 1.;
    for(int k = 1; k < N; k++)
    {
        y[k] = y[k - 1] / (1 + A*h);
        if(fabs(y[k])< EPS)
        {
            y[k] = 0.;
        }
    }
    return 1;
}

int no3(double A, int n, double* y, double EPS)
{
    double h = pow(10, -n);
    int N = (int)(1/h);
    y[0] = 1.;
    for(int k = 1; k < N; k++)
    {
        y[k] = y[k - 1] * (2 - A * h) / (2 + A * h);
        if(fabs(y[k])< EPS)
        {
            y[k] = 0.;
        }
    }
    return 1;
}

int no4(double A, int n, double* y)
{
    double h = pow(10, -n);
    int N = (int)(1/h);
    y[0] = 1.;
    y[1] = 1 - A*h;
    for(int k = 2; k < N; k++)
    {
        if(fabs(y[k-1])> 1e+10)
        {
            y[k] = 1e+10;
        }
        else {
            y[k] = y[k - 2] - 2*A*h*y[k-1];
        }

    }
    return 1;
}

int no5(double A, int n, double* y,double EPS )
{
    double h = pow(10, -n);
    int N = (int)(1/h);
    y[0] = 1.;
    y[1] = 1 - A*h;
    for(int k = 2; k < N; k++)
    {
        if((y[k-1] < EPS) && (y[k-2] < EPS))
        {
            y[k] = 0;
            continue;
        }
        else {
            if((y[k-1] < EPS))
            {
                y[k] = (- 0.5*y[k - 2])/ (1.5 + A*h);
                continue;
            }
            if((y[k-2] < EPS))
            {
                y[k] = (2*y[k - 1])/ (1.5 + A*h);
                continue;
            }
        }
        y[k] = (2*y[k - 1] - 0.5*y[k - 2])/ (1.5 + A*h);
    }
    return 1;
}

int no6(double A, int n, double* y)
{
    double h = pow(10, -n);
    int N = (int)(1/h);
    y[0] = 1.;
    y[1] = 1 - A*h;
    for(int k = 2; k < N; k++)
    {
        if((fabs(y[k-1])> 1e+10)||(fabs(y[k-2])> 1e+10))
        {
            y[k] = 1e+10;
        }
        else{
             y[k] = (2*A*h - 3)*y[k - 2] + 4*y[k - 1];
        }

    }
    return 1;
}

