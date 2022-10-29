#include "Solve.h"

int Solve(int n,double* A,double* b,double* x)
{
    //умножаем слева на Tij
    for( int i = n-2; i > 1; i--)
    {
        for(int j = n-1; j > i-1; j--)
        {
            TA(i,j, A, n); //Tij*A
            Tb(i,j,b, n); //Tij*b
        }
    }
    // обратный ход метода Гаусса
    x[n-1] = b[n-1];
    for( int i = n-2; i > 0; i--)
    {
        x[n-1] = b[n-1];
        for(int j = i+1; j < n; j++)
        {
            x[n-1] -= A[i*n +j]* x[j];
        }
    }
    return 1;
}

int TA(int i, int j, double* A, int n)
{
    double x = 0;
    double y = A[i*n + j];
    double cos_phi = 0;
    double sin_phi = 0;
    for(int k = i; k < j; k++)
    {
        x+= (A[i*n+k])* (A[i*n+k]);
    }
    x = sqrt(x);
    cos_phi = x / (sqrt(x*x +y*y));
    sin_phi = -y / (sqrt(x*x +y*y));
    //умножение
    for(int k = 0; k < n; k++) //столбцы А
    {
        double xi = A[k*n + i];
        double xj = A[k*n + j];
        for(int l = 0; l < n; l++) //строки
        {
            if((k==i)&&(l==j)) //A[j;i] = 0
            {
                A[l*n + k] = 0;
                continue;
            }
            if(l == i)
            {
                A[l*n + k] = xi*cos_phi - xj*sin_phi;
                continue;
            }
            if(l == j)
            {
                A[l*n + k] = xi*sin_phi + xj*cos_phi;
                continue;
            }
        }
    }
    return 1;
}

int Tb(int i, int j, double* b, int n)
{
    double x = 0;
    double y = b[i];
    double cos_phi = 0;
    double sin_phi = 0;
    for(int k = i; k < j; k++)
    {
        x+= (b[k])* (b[k]);
    }
    x = sqrt(x);
    cos_phi = x / (sqrt(x*x +y*y));
    sin_phi = -y / (sqrt(x*x +y*y));
    //умножение
    double xi = b[i];
    double xj = b[j];
    b[i] = xi*cos_phi - xj*sin_phi;
    b[j] = xi*sin_phi + xj*cos_phi;
    //b[i] = sqrt(x*x+b[j]*b[j]); // b[j] еще не изменен тк i>j
    //b[j] = 0;
    return 1;
}

