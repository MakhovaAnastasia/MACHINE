#include "Solve.h"
#include "ReadMatrix.h"

int Solve(int n,double* A,double* b,double* x)
{
    //умножаем слева на Tij
    for( int i = 0; i < n-1; i++)
    {
        for(int j = i+1; j < n; j++)
        {
            int res  = TA(i,j, A,b, n); //Tij*A Tij*b
            if(res == -1)
            {
                return -1;
            }
        }
    }
    // обратный ход метода Гаусса
    for( int i = n-1; i >0; i--)
    {
        x[i] = b[i];
        for(int j = i+1; j < n; j++)
        {
            x[i] -= A[i*n +j]* x[j];
        }
        if(abs(A[i*n + i]) < EPS) //на диагонали  нашли 0
        {
            //cout<<"нет точного ответа. x["<<i<<"] = 1 "<<endl;
            x[i] = 1;
        }
        else{
            x[i]/= A[i*n +i];
        }
    }
    return 0;
}

int TA(int i, int j, double* A, double*b, int n)
{
    //определим угол поворота
    double x = A[i*n + i];
    double y = A[j*n + i];//[i*n + j];
    double root = sqrt(x*x + y*y);
    if(abs(root) < EPS)
    {
        return -1;
    }
    double cos_phi = x / root;
    double sin_phi =  -y / root;
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

