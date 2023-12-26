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

double Richardson(double*X, double* True_value, double* A, double *F, double tau, double q, int N, int mIter);

void dot(const double* A, const double* X, double* Temp, int N);
double norm_h(const double* F, double* Temp, int N);
double lambda_n(int n, int N, double p);


//plot '2.txt' using 1:2 with linespoints, '2.txt' using 1:3 with lines

int main(int argc, char* argv[])
{
    double *True_value, *X, *A, *F;
    int N = -1; // число узлов
    double p = 0.;
    
    int mIter = 100;
    //double eps = 0.00000000001;


    if(!((argc == 3)&&
    (sscanf(argv[1],"%lf",&p)==1)&&
    (sscanf(argv[2],"%d",&N)==1)))
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

    if(p<0)
    {
        //ошибка 
        cout<<p<<endl;
        cout<<" p >= 0"<<endl;
        return -1;
    }

    double a = -(N)*(N);//
    double b  = (2*(N)*(N) +p);//
    double m = lambda_n(1,N,p); //p;
    double M = lambda_n(N-1,N,p);//(4*(N+1)*(N+1) +p);
    double tau = 2/(double)(m+M);
    double q = (M-m)/(double)(M+m);


    //cout<<q<<endl;


    True_value = (  double*) malloc(sizeof(  double)*(N+1));
    X = (  double*) malloc(sizeof(  double)*(N+1));
    A = (  double*) malloc(sizeof(  double)*((N+1)*(N+1)));
    F = (  double*) malloc(sizeof(  double)*(N+1));

    //True_value
    for(int i = 0; i < N+1; i++)
    {
        if((i == 0)||(i==N))
        //if(i%2 ==0)
        {
            True_value[i] = 0.;
        }
        else{
            True_value[i] = 1.;
        }
    }

    for(int i  = 0; i< N+1; i++)
    {
        for(int j = 0; j< N+1; j++)
        {
            if(i == j)
            {
                A[(N+1)*i+j] = b;

                cout<<setw(5)<<A[(N+1)*i+j]<<" ";
                continue;
            }
            if((j -1 == i) ||(i-1 == j))
            {
                A[(N+1)*i+j] = a;

                cout<<setw(5)<<A[(N+1)*i+j]<<" ";
                continue;
            }
            else{
                A[(N+1)*i+j] = 0;
            }
           cout<<setw(5)<<A[(N+1)*i+j]<<" ";
        }
        cout<<endl;
    }
    //X
    for(int i = 0; i< N+1; i++)
    {
        if((i== 0)|| (i==N))
        {
            X[i] = 0;
        }
        else{
            X[i] = 0.5;
        }
    }

    //F
    for(int i = 0; i < N+1; i++)
    {
        F[i] = 0;
        for(int j = 1; j<N;j++)
        {
            F[i]+=A[(N+1)*i+j];
        }

        //cout<<"         "<<F[i]<<endl;
    }

    /* 2 */
    cout<<"невязка (Ричардсон): "<<Richardson(X,True_value, A, F, tau, q, N, mIter)<<endl;


    free(True_value);
    free(X);
    free(A);
    free(F);

    return 0;
}

double Richardson(double*X, double* True_value, double* A, double *F, double tau, double q, int N, int mIter)
{
    ofstream out;
    out.open("4.txt");
    //double norm_0 = norm_h(True_value, X, N);
    double nevyaz = 0;
    double qk = q;
    double* Mass;
    Mass = (double*) malloc(sizeof(double)*(N+1));

    if(out.is_open())
    {
        out<<setprecision(15)<<fixed;
        for(int k = 1; k <= mIter ; k++)
        {
            dot(A,X, Mass, N);//Ax

            for(int i = 0; i < N+1; i++)
            {
                //x^(k) = (I-tA)x^(k-1) +tF
                X[i] = X[i]- tau*Mass[i] + tau*F[i];

            }

            dot(A,X,Mass,N);
            // ||True - X||, ||F-Ax||
            out<<setw(25)<<k<<setw(25)<<norm_h(True_value, X, N)<<setw(25)<<norm_h(F,Mass,N)<<endl;

           //out<<setw(25)<<k<<setw(25)<<norm_h(True_value, X, N)<<setw(25)<<qk*norm_0<<endl;
           qk*=q;

        }


        //nevyazka ||F-Ax||
        dot(A, X, Mass, N);
        nevyaz = norm_h(F, Mass, N);

        free(Mass);
        out. close();
        return nevyaz;
    }
    out.close();
    return -1;
    
}


double norm_h(const double* F, double* Temp, int N)
{
    double res = 0;
    for(int i = 0; i< N+1; i++)
    {
        res += (F[i] - Temp[i])*(F[i] - Temp[i])/(double) N;
    }

    return sqrt(res);
}

void dot(const double* A, const double* X, double* Temp, int N)
{
    for(int i = 0; i< N+1; i++)
    {
        Temp[i] = 0;
        for(int j = 0; j<N+1; j++)
        {
            Temp[i] += A[(N+1)*i +j] * X[j];
        }
    }
    return;
}

double lambda_n(int n, int N, double p)
{
    return p - (2*N*N*(cos(M_PI*(n) / (double)N) - 1));
}


