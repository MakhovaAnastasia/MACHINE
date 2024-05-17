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

double lambda_n(int n, int N);
double y_k_n(int k,int n, int N);
int check_solution(int N, double eps);


int main(int argc, char* argv[])
{
    int N = -1; // число узлов
    double eps = 1e-5;
    int err = -1;
    double** A;
    double max = -1.;
    double max2 = -1.;
    int maxn = -1;
    double pogr = 0.;
    double res = 0;


    if(!((argc == 2)&&
    (sscanf(argv[1],"%d",&N)==1)))
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
    A =(double**) malloc(N* sizeof(double*));
    for(int i = 0; i<N; i++)
    {
        A[i] = (double*) malloc(N*sizeof(double));
    }

    err = check_solution(N, eps);

    if(err  == 0)
    {
        cout<<"True"<<endl;
    }
    else{
        cout<<"False"<<endl;
    }

    //проверим скалярное произведение
    if(N<= 20)
    {
        double    h = (1/(double)N);

        cout<<setprecision(4)<<fixed;

        for(int i = 0; i < N; i++) //k = 1...N-1
        {
            for(int j = 0; j < N; j++)
            {
                res = 0;
                for(int n = 0; n < N; n++ )
                {
                    res+=y_k_n(i,n,N)*y_k_n(j,n,N);
                }
                A[i][j] = res*h;
                cout<<setw(15)<<scientific<<A[i][j];
                if(max < (res*h))
                {
                    max = res*h;
                }
            }
            cout<<endl;
        }
        cout<<"||(y,y)_h||_h  =  "<<max<<endl;
    }

    //--------------------------------------------
    //самый несобственный вектор: (A-lambdaE)y V 0

    for(int n = 0; n < N; n++)
    {
        pogr = 0;
        for(int i = 0; i < N; i++) //k = 1...N-1
        {
            res = 0;
            for(int j = 0; j < N; j++)
            {

                if(i==j)
                {
                    res +=(A[i][i] - lambda_n(n, N))*A[n][j];
                }
                else {
                     res +=A[i][j]*A[n][j];
                }
            }
            pogr += abs(res);
        }
        if(pogr > max2)
        {
            max2 = pogr;
            maxn = n;
        }
    }
    cout<<endl<<"самый 'несобственный' вектор: "<<maxn<<endl;
    for(int i = 0;i< N; i++)
    {
        cout<<A[maxn][i]<<" ";
    }
    cout<<endl;
    //-----------------------------------------------

    for(int i = 0; i<N; i++)
    {
        free(A[i]);
    }
    free(A);
    return 0;
}


double lambda_n(int n, int N)
{
    double h = (1/(double)N);
    return 2*N*N*(1 - cos(M_PI*h*(n + 0.5)));
}

double y_k_n(int k,int n, int N)
{
    double h = (1/(double)N);
    return  cos(M_PI*k*h*(n + 0.5)) ;
}

int check_solution(int N, double eps)
{
    double h = (1/(double)N);
    double res = 100;

    for(int n = 0; n < N; n++)
    {
        for(int k = 1; k < N; k++)
        {
            res = ((y_k_n(k+1, n, N)  - 2*y_k_n(k, n, N) + y_k_n(k-1, n, N))/(double)(h*h)) + lambda_n(n,N) *y_k_n(k, n, N);


            if(fabs(res) > eps)
            {
                cout<<"k = "<<k<<endl<<"n = "<<n<<endl<<"res = "<<res<<endl;
                return -1;
            }
        }
    }
    return 0;
}




