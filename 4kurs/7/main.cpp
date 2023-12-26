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

double find_q(double* A,int N);
double BSolver(double*X, const double*A, const double*B, const double*F, double tau, int N, double p, int mIter,  double* Y, double* Temp, double q);
double norm_h(const double* F, double* Temp, int N);
void dot(const double* A,const double* X, double* Temp, int N);
double lambda_n(int n, int N, double p);
double dot_f_phi(double* Temp, int N, int j);
void y_k(double* Y, double* Temp, int N, double p);

//plot '3.txt' using 1:2 with linespoints, '3.txt' using 1:3 with linespoints


int main(int argc, char* argv[])
{
    double *X, *A,*B, *F;
    double* Y, *Temp;
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

    double a = -(N)*(N)/M_PI/M_PI;//B
    double b  = (2*(N/M_PI)*(N/M_PI) + p);//B
    //double m = lambda_n(1,N,p); 
    //double M = lambda_n(N-1,N,p);
    double tau = 0.8;//2/(double)(m+M);
    //double q = (M-m)/(double)(M+m);


    X = new double[(N+1)];//(  double*) malloc(sizeof(  double)*(N+1));
    A = new double[(N+1)*(N+1)];//(  double*) malloc(sizeof(  double)*((N+1)*(N+1)));
    B = new double[(N+1)*(N+1)];//(  double*) malloc(sizeof(  double)*((N+1)*(N+1)));
    F = new double[(N+1)];//(  double*) malloc(sizeof(  double)*(N+1));
    Y = new double[(N+1)];//(  double*) malloc(sizeof(  double)*(N+1));
    Temp = new double[(N+1)];//(  double*) malloc(sizeof(  double)*(N+1));

    //A,B
    for(int i  = 0; i< N+1; i++)
    {
        for(int j = 0; j< N+1; j++)
        {
            if(i == j)
            {
                B[(N+1)*i+j] = b;//2*(N*N/M_PI/M_PI) + 1 +sin(M_PI*M_PI*(i)/(double) N)*sin(M_PI*M_PI*(i)/(double) N);
                A[(N+1)*i+j] = 2*(N*N/M_PI/M_PI) + 1 +sin(M_PI*M_PI*(i)/(double) N)*sin(M_PI*M_PI*(i)/(double) N);
                cout<<setw(5)<<A[(N+1)*i+j]<<" ";
                continue;
            }
            if((j -1 == i) ||(i-1 == j))
            {
                B[(N+1)*i+j] = a;
                A[(N+1)*i+j] = a;
                cout<<setw(5)<<A[(N+1)*i+j]<<" ";
                continue;
            }
            else{
                B[(N+1)*i+j] = 0;
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
        Y[i] = 0;
        Temp[i] = 0;

       // cout<<"         "<<F[i]<<endl;
    }

    double q =  find_q(A,N);
    /* 3 */
    cout<<"Невязка: "<<BSolver(X,A,B,F,tau,N,p,mIter, Y,Temp, q)<<endl;

    delete [] X;//free(X);
    delete [] A;//free(A);
    delete [] B;//free(B);
    delete [] F;//free(F);
    delete [] Temp;//free(Temp);
    delete [] Y;//free(Y);


    return 0;
}

double find_q(double* A,int N)
{
    double sum;
    double max = -1;
    for(int i = 0; i< N+1; i++)
    {   
        sum = 0;
        for(int j = 0; j< N+1; j++)
        {
            if(i!= j)
            {
                sum+= fabs(A[(N+1)*i+j]);
            }
        }

        if(sum/A[(N+1)*i+i]> max)
        {
            max = sum/fabs(A[(N+1)*i+i]);
        }
    }
    return max;
}

double BSolver(double*X,const double*A,const  double*B,const double*F, double tau, int N, double p, int mIter, double* Y, double* Temp,double q)
{   
    ofstream out;
    out.open("33.txt");
    out<<setprecision(15)<<fixed;
    double qk = B[(N+1)*2+2];
    qk= q;

    dot(A,X,Temp, N);
    double norm_0  = norm_h(F,Temp,N);
    for(int k = 1; k < mIter; k++)
    {
        dot(A,X,Temp,N);//Temp := A*x^k-1

        for(int i = 0; i<N+1; i++)
        {
            Temp[i] = F[i] - Temp[i];//Temp:= F - AX^k-1 (=  By^k)
        }
        
        y_k(Y,Temp,N,p);// Y методом фурье из By^k = Temp

        for(int i = 0; i<N+1; i++)
        {
            X[i] = X[i] + tau*Y[i];
        }
        
        dot(A,X,Temp, N);// Temp:= Ax
        out<<setw(25)<<k<<setw(25)<<norm_h(F,Temp,N)<<setw(25)<<qk*norm_0<<endl;
        qk*=q;
    }
    out. close();
    dot(A,X,Temp, N);// Temp:= Ax
    return norm_h(F,Temp,N);
}

double norm_h(const double* F, double* Temp, int N)
{
    double res = 0;
    for(int i = 0; i< N+1; i++)
    {
        if(res < fabs(F[i] - Temp[i]))
        {
            res = fabs(F[i] - Temp[i]);
        }
        //res += (F[i] - Temp[i])*(F[i] - Temp[i])/(double) N;
    }
    return res;
    //return sqrt(res);
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

double dot_f_phi(double* Temp, int N, int j)
{
    double res = 0.0;
    double h = 1/(  double)N;
    for(int i = 1; i < N; i++)
    {
       res += Temp[i] * sin(M_PI*(j)*((i)*h));
    }
    return res*h;//
}

void y_k(double* Y, double* Temp, int N, double p)
{
    double h = 1/(  double)N;
    for(int i = 0; i < N+1; i++)
    {
        Y[i] =0;
        for(int j = 1; j<N; j++)
        {
            Y[i] += dot_f_phi(Temp, N, j) *2.0 * sin(M_PI*(j)*(i)*h)/ lambda_n(j,N,p);
        }
    }
    return;
}





