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

int no(double A, int n, double* y, int num);
int no1(double A, int n, double* y);
int no2(double A, int n, double* y);
int no3(double A, int n, double* y);
int no4(double A, int n, double* y);
int no5(double A, int n, double* y);
int no6(double A, int n, double* y);
double En(double A, int n, double* y);

//plot '2.txt' using 1:2 with linespoints, '2.txt' using 1:3 with lines

int main(void)
{
    int A[3] = {1,10,1000};
    int N[4] = {1, 2, 3, 6};
    double* y;

    cout<<setw(15)<<"# "<<setw(15)<<"E1 "<<setw(15)<<"E2 "<<setw(15)<<"E3 "<<setw(15)<<"E6 "<<setw(15)<<"m "<<setw(15)<<"A "<<endl;
    for(int i = 1;i <=6; i++)
    {
        cout<<"-------------------------------------------------------------------------------------------------------------"<<endl;
        for(int k = 0; k < 3; k++)//Ak
        {
            cout<<setw(15)<<i;
            for(int j = 0; j < 4; j++)//Ej
            {
                double h = pow(10, -N[j]);
                int col = (int)(1/h);
                y = (double*) malloc(sizeof(double) * col);
                no(A[k], N[j], y, i);
                cout<<setw(15)<<En(A[k], N[j], y);
            }
            cout<<setw(15)<<"m "<<setw(15)<<A[k]<<endl;
        }
    }
    
    return 0;
}

double En(double A, int n, double* y)
{
    double h = pow(10, -n);
    int N = (int)(1/h);
    double max = -1;
    for(int k = 0; k< N; k++)
    {
        if(abs(y[k] - exp(-A*h*k)) > max)
        {
            max = abs(y[k] - exp(-A*h*k));
        }
    }
    return max;
}

int no(double A, int n, double* y, int num)
{
    switch (num)
    {
    case 1:
        no1(A, n, y);
        break;
    case 2:
        no2(A, n, y);
        break;
    case 3:
        no3(A, n, y);
        break;
    case 4:
        no4(A, n, y);
        break;
    case 5:
        no5(A, n, y);
        break;
    case 6:
        no6(A, n, y);
        break;
    default:
        break;
    }
    return 1;
}

int no1(double A, int n, double* y)
{
    double h = pow(10, -n);
    int N = (int)(1/h);
    y[0] = 1.;
    for(int k = 1; k < N; k++)
    {
        y[k] = y[k - 1]*(1 - A*h);
    }
    return 1;
}

int no2(double A, int n, double* y)
{
    double h = pow(10, -n);
    int N = (int)(1/h);
    y[0] = 1.;
    for(int k = 1; k < N; k++)
    {
        y[k] = y[k - 1] / (1 + A*h);
    }
    return 1;
}

int no3(double A, int n, double* y)
{
    double h = pow(10, -n);
    int N = (int)(1/h);
    y[0] = 1.;
    for(int k = 1; k < N; k++)
    {
        y[k] = y[k - 1] * (2 - A * h) / (2 + A * h);
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
        y[k] = y[k - 2] - 2*A*h*y[k-1];
    }
    return 1;
}

int no5(double A, int n, double* y)
{
    double h = pow(10, -n);
    int N = (int)(1/h);
    y[0] = 1.;
    y[1] = 1 - A*h;
    for(int k = 2; k < N; k++)
    {
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
        y[k] = (2*A*h - 3)*y[k - 2] + 4*y[k - 1];
    }
    return 1;
}

