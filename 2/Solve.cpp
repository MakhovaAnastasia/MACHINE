#include "Solve.h"
#include "ReadMatrix.h"
//problem- 2

int Solve(int n,double* A,double* x, double EPS, double* Q_cos, double* Q_sin)
{
    //почти треугольный вид
    for(int k = 0;k < n-1; k++) //n-2 steps
    {
        double s= 0.;
        double norm_x = 0;
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
                if(norm_x < EPS)
                {
                      cout<<"norm";
                    return -1;
                }
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
            cout<<"!!";
    int k = n-1;
    do{
                cout<<"!!";
        double s = A[k*n +k];
        double nA  = 0;
        //QR со сдвигом: A -sI = QR
        //R = ПT(i,i+1)(A-sI)
        for(int j = 0; j< k; j++)
        {
            A[j*n+j] -=s;
        }
        for(int i = 0; i< n-1; i++)
        {
            int res=TA(i, i+1, A,k,EPS, Q_cos, Q_sin);
            if(res == -1)
            {
                cout<<"do";
                return -1;
            }
        }
        //A = RQ +sI = RПT(i,i+1) +sI
        for(int i = n-2; i>=0; i++)
        {
            //T(i,j)^t = T(j,i)
            int res=TA(i+1, i, A,k,EPS, Q_cos, Q_sin);
            if(res == -1)
            {
                  cout<<"do2";
                return -1;
            }
        }
        for(int j = 0; j< n; j++)
        {
            A[j*n+j] +=s;
        }
        for(int i = 0; i< n; i++)
        {
            double sum = 0;
            for( int j = 0; j <n; j++)
            {
                sum += abs(A[n*i +j]);
            }
            if((i == 0)||(nA < sum))
            {
                nA = sum;
            }
        }
        if(abs(A[k*n +(k-1)])< EPS*nA)
        {
           cout<<x[k]<<" "<< s<<endl;
                   x[k] = s;//новое собственное значение
            k--;
            continue;
        }
    }while(k > 1);
    //у матрицы 2х2 легко найти собственные значения:
    double D = (A[0*n +0] +A[1*n +1])*(A[0*n +0] +A[1*n +1]) - 4*((A[0*n +0]*A[1*n +1]) + (A[0*n +1] +A[1*n +0]));
    if(D < -EPS)//нет решения
    {
                cout<<"!!";
        return -1;
    }
    if(abs(D) < EPS)//D = 0
    {
        cout<<"!!";
        x[1] = (A[0*n +0] +A[1*n +1])/2.0;
        x[0] = (A[0*n +0] +A[1*n +1])/2.0;
        return 0;
    }
            cout<<"!!";
        x[1] = ((A[0*n +0] +A[1*n +1]) + sqrt(D))/2.0;
        x[0] = ((A[0*n +0] +A[1*n +1]) - sqrt(D))/2.0;
    return 0;
}

int TA(int i, int j, double* A,int n,double EPS,double* Q_cos, double* Q_sin)
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
    Q_cos[i] = cos_phi;
    Q_sin[i] = sin_phi;
    double xi = 0;
    double xj = 0;
    //умножение TA

    for(int k = i; k < n; k++) //столбцы А
    {
        xi = A[i*n+k];
        xj = A[j*n + k];
        if(k==i)
        {
            A[i*n + i] = root;
            if(abs(A[i*n + i]) < EPS)
            {
                A[i*n + i] = 0;
            }
            A[j*n + i] = 0;
        }
        else
        {
            A[i*n + k] = xi*cos_phi - xj*sin_phi;
            if(abs(A[i*n + k]) < EPS)
            {
                A[i*n + k] = 0;
            }

            A[j*n + k] = xi*sin_phi + xj*cos_phi;
            if(abs(A[j*n + k]) < EPS)
            {
                A[j*n + k] = 0;
            }
        }
    }
    return 0;
}



