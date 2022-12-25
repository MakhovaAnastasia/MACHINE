#include "Solve.h"
#include "ReadMatrix.h"
//problem- 2

int Solve(int n,double* A,double* x, double EPS)
{
    for(int k = 0;k < n-1; k++) //n-2 steps
    {
        double s= 0.;
        double norm_x = 0;
        double b_x = 0;
        double norm_a = 0;
        for(int j = k+2; j < n; j++)
        {
            s+=(abs(A[j*n+k]))*(abs(A[j*n+k]));
        }
        s = sqrt(s);
        //x
        for(int j = k+1; j < n; j++)
        {
            x[j] = A[(n*j) + k];
            if(j==k+1)
            {
                norm_a = sqrt(s + abs(A[(k+1)*n+k])*abs(A[(k+1)*n+k]));
                x[j]-=norm_a;
                norm_x = sqrt((x[j]*x[j])+s);
                if(norm_x < EPS)    return -1;
            }
            x[j]/= norm_x;
        }
        //U(x)AU(x)
        for(int j = k; j<n; j++)
        {
            double scalar = 0;
            //(||a||,0 0 0 ... 0 )
            for(int i = k+1; i < n; i++)
            {
                if((i==k+1) && (j == k))
                {
                    A[i*n+j] = norm_a;
                    continue;
                }  
                if(j == k)
                {
                    A[i*n+j] = 0;
                    continue;
                }  
                scalar+=A[i*n+j]*x[i];
            }
            //U(x)A по столбцам
            for(int i = k+1; i < n; i++)
            {   
                if(j != k)
                {
                    A[i*n+j] = A[i*n+j]-2*scalar*x[i];
                }
            }
        }
        //AU(x)~~ по лемме 11 по строкам
        for(int i = 0; i < n; i++)
        {   
            double scalar = 0;
            for(int j = k+1; j<n; j++)
            {
                scalar+=A[i*n+j]*x[j];
            }
            for(int j = k+1; j<n; j++)
            {
                A[i*n+j] = A[i*n+j]-2*scalar*x[j];
            }
        }
    }
    QRRotate(n,A,x,EPS);
    return 0;
}

int TA(int i, int j, double* A, double*b, int n,double EPS)
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

