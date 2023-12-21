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
double f(int i,   double p, int N); //функция
double lambda_n(int n, int N, double p);
int Get_Coef(  double *C, int N, double p);//получаем коэфициенты с_m
double dot_f_phi(int N, int j,   double p); //скалярное произведение f и phi_m
double y_k(  double* C, int N, int k); //значение ряда y_k
double Write(int N,   double* C,   double* True_val); //выведем в файл i, true_val, y, |true_val - y|
double Err(int N,   double* C,   double p);
double Richardson(double*X, double* X_new, double* True_value, double* A, double *F, double tau, double q, int N, double eps, int mIter);
//plot '1.txt' using 1:2 with linespoints, '1.txt' using 1:3 with lines

int main(int argc, char* argv[])
{
    double *C, *True_value, *X, *A, *F, *X_new;
    int N = -1; // число узлов
    double p = 0.;
    
    int mIter = 50;
    double eps = 0.00000000001;


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

    double a = -N*N;
    double b  = (2*N*N +p);
    double m = min(b - fabs(2*a), b - fabs(a));
    double M = (-2)*a + b;
    double tau = 2/(double)(m+M);
    double q = (M-m)/(double)(M+m);

        //cout<<M<<endl<<m<<endl<<q<<endl<<tau<<endl;

    C = (  double*) malloc(sizeof(  double)*(N));
    True_value = (  double*) malloc(sizeof(  double)*(N+1));
    X = (  double*) malloc(sizeof(  double)*(N+1));
    X_new = (  double*) malloc(sizeof(  double)*(N+1));
    A = (  double*) malloc(sizeof(  double)*((N+1)*(N+1)));
    F = (  double*) malloc(sizeof(  double)*(N+1));

    //True_value
    for(int i = 0; i < N+1; i++)
    {
        if(i%2 == 0)
        {
            True_value[i] = 0.;
        }
        else{
            True_value[i] = 1.;
        }
    }

    //F,A
    for(int i  = 0; i< N+1; i++)
    {
        if((i ==0) ||(i== N))
        {
            F[i] = a;
        }
        if(i%2 == 0)
        {
            F[i] = 2*a;
        }
        else{
            F[i] = b;
        }

        for(int j = 0; j< N+1; j++)
        {
            if(i == j)
            {
                A[N*i+j] = b;
                cout<<A[N*i+j]<<" ";
                continue;
            }
            if((j -1 == i) ||(i-1 == j))
            {
                A[N*i+j] = a;
                cout<<A[N*i+j]<<" ";
                continue;
            }
            else{
                A[N*i+j] = 0;
            }
            cout<<A[N*i+j]<<" ";
        }
        cout<<endl;
    }
    //X, X_new
    for(int i = 0; i< N+1; i++)
    {
        if((i== 0)|| (i==N))
        {
            X[i] = 0;
        }
        else{
            X[i] = 1;
        }
        X_new[i] = 0;
    }

    /* 1 */
    Get_Coef(C, N, p);
    Write(N, C, True_value);
    cout<<"невязка: "<<Err(N,C,p)<<endl;

    /* 2 */
    cout<<"невязка (Ричардсон): "<<Richardson(X,X_new,True_value, A, F, tau, q, N, eps, mIter)<<endl;

    free(True_value);
    free(C);
    free(X);
    free(X_new);
    free(A);
    free(F);

    return 0;
}

double Richardson(double*X, double* X_new, double* True_value, double* A, double *F, double tau, double q, int N, double eps, int mIter)
{
    ofstream out;
    out.open("2.txt");
    double sum = 0;
    double norm_0 = 0;
    double norm_k = 0;
    double nevyaz = 0;
    if(out.is_open())
    {
        out<<setprecision(15)<<fixed;
        for(int k = 1; k < mIter ; k++)
        {
            //x^(k+1) = (I-tA)x^k +tF
            for(int i = 0; i < N+1; i++)
            {
                // I-tA * X
                for(int j = 0; j < N+1; j++)
                {
                    if(j == i)
                    {
                        sum+=(1 - tau*A[N*i +j])*X[j];
                    }
                    else{
                        sum+=(- tau*A[N*i +j])*X[j];
                    }
                }
                // +tF
                X_new[i] = sum + tau*F[i];
                sum = 0;

                //norm
                if(k==1)
                {
                    norm_0 += (True_value[i] - X[i])*(True_value[i] - X[i])/(double)N;
                }
                norm_k += (True_value[i] - X_new[i])*(True_value[i] - X_new[i])/(double)N;


            }
            out<<setw(25)<<k<<setw(25)<<sqrt(norm_k)<<setw(25)<<pow(q,k)*sqrt(norm_0)<<endl;
            //swap
            for(int l = 0; l<N+1; l++)
            {
                X[l] = X_new[l];
            }
            norm_k = 0;
        }

        //nevyazka
        for(int i = 0; i < N+1;i++)
        {
            double res = 0;
            for(int j = 0;j < N+1; j++)
            {
                res+= A[N*i+j]*X[j];
            }

            nevyaz += (F[i] - res)*(F[i] - res)/(double)N;
        }
        out. close();
        return sqrt(nevyaz);
    }
    out.close();
    return -1;
    
}


double f(int i,   double p, int N)
{
    if(i == 0)
    {
        return (-1)*N*N;
    }
    if(i%2 == 0)
    {
        return (-2)*N*N;
    }
    else{
        return (2*N*N) + p;
    }
}

double lambda_n(int n, int N, double p)
{
    return p - (2*N*N*(cos(M_PI*n / (double)N) - 1));
}

int Get_Coef(  double *C, int N, double p)
{
    for(int i = 0; i < N; i++)
    {
        C[i] = dot_f_phi(N,i, p) *2.0 / lambda_n(i,N,p);//;  /dot_phi_phi(N,i)
    }
    return 1;
}

double dot_f_phi(int N, int j,   double p)
{
      double res = 0.0;
      double h = 1/(  double)N;
    for(int i = 0; i < N; i++)
    {
       res += f(i, p, N) * sin(M_PI*(j)*(i*h));
    }
    return res*h;//
}

double y_k(  double* C, int N, int k)
{
      double res = 0.0;
      double h = 1/(  double)N;
    for(int i = 1; i < N; i++)
    {
        res += C[i] * sin(M_PI*(i)*k*h);
    }
    return res;
}



  double Write(int N,   double* C,   double* True_value)
{
    ofstream out;
    out.open("1.txt");
      double max = 0;

    if(out.is_open())
    {
        out<<setprecision(15)<<fixed;
        for(int i = 0; i < N; i++)
        {
            if(max < fabs(True_value[i] - y_k(C,N,i)))
            {
                max = fabs(True_value[i] - y_k(C,N,i));
            }
            out<<setw(25)<<True_value[i]<<setw(25)<<y_k(C,N,i)<<setw(25)<<fabs(True_value[i] - y_k(C,N,i))<<endl;
        }
        out. close();
        return max;
    }
    out.close();
    return -1;
}

  double Err(int N,   double* C,   double p)
{
    double max = 0;
    double a = (-1)*(N*N);
    double b = (2*N*N) + p;

    
    for(int i = 1; i < N; i++)
    {
        max += (a*y_k(C,N,i-1) + b*y_k(C,N,i) + a*y_k(C,N,i+1) - f(i,p,N))*(a*y_k(C,N,i-1) + b*y_k(C,N,i) + a*y_k(C,N,i+1) - f(i,p,N)) /(double)N;
    }

    return sqrt(max);
}



