#include "Solve.h"
#include "ReadMatrix.h"
//problem- 2

int Solve(int n,double* A,double* x, double EPS)
{
    for(int i = 0;i<n; i++)
    {
        //s[k]
        for(int j = k+2; j<n;j++)
        {
            s+=(abs(A[j*n+k]))*(abs(A[j*n+k]));
        }
        s = sqrt(s);
        for(int j = k; j<n;j++)
        {
            x[j] = A[(n*j) +k];
            if(j==k)
            {
                x[j]-=sqrt(s+abs(A[A[(k+1)*n+k]]))
            }
        }
    }

    //

    NearlyTriangle(n, A, EPS);
    QRRotate(n,A,x,EPS);
    return 0;
}

int TA(int i, int j, double* A, double*b, int n)
{
    //определим угол поворота
    double x = A[i*n + i];
    double y = A[j*n + i];//[i*n + j];
    double cos_phi = x / (sqrt(x*x + y*y));
    double sin_phi =  -y / (sqrt(x*x +y*y));
    double xi = 0;
    double xj = 0;
    //умножение TA
    for(int k = i; k < n; k++) //столбцы А
    {
        xi = A[i*n+k];
        xj = A[j*n + k];
        //A[i*n + k] = xi*cos_phi - xj*sin_phi;
        //A[j*n + k] = xi*sin_phi + xj*cos_phi;
        for(int l = i; l < n; l++) //строки
        {
            if((k==i)&&(l==j)) //A[j;i] = 0
            {
                A[l*n + k] = 0;
                continue;
            }
            if((k==i)&&(l==i)) //A[i;i] = ||a(i, j)||
            {
                A[l*n + k] = sqrt(x*x + y*y);
                continue;
            }
            if(l == i) //A[i,k] = A[i,k]*cos- A[j,k]*sin
            {
                A[l*n + k] = xi*cos_phi - xj*sin_phi;
                continue;
            }
            if(l == j) //A[i,k] = A[i,k]*sin + A[j,k]*cos
            {
                A[l*n + k] = xi*sin_phi + xj*cos_phi;
                continue;
            }
            //иначе: A[l,k] =A[l,k]
        }

    }

    //Tb
    xi = b[i];
    xj = b[j];
    b[i] = xi*cos_phi - xj*sin_phi;
    b[j] = xi*sin_phi + xj*cos_phi;

    return 0;
}

