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

//plot '2.txt' using 1:2 with linespoints, '2.txt' using 1:3 with lines

int main(int argc, char* argv[])
{
    int N = -1; // число узлов
    double eps = 1e-10;
    int res = -1;


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

    res = check_solution(N, eps);

    if(res  == 0)
    {
        cout<<"True"<<endl;
    }
    else{
        cout<<"False"<<endl;
    }

    return 0;
}


double lambda_n(int n, int N)
{
    double h = (1/(double)N);
    return 2/((double) h*h)*(1 -cos((2*M_PI*n + M_PI)/((double)(2*N -1))));
}

double y_k_n(int k,int n, int N)
{
    return  cos((M_PI*k*(2*n - 1))/ ((double)(2*N - 1))) ;
}

int check_solution(int N, double eps)
{
    double h = (1/(double)N);
    double res = 100;

    for(int n = 0; n < N; n++)
    {
        for(int k = 1; k < N; k++)
        {
            res = ((y_k_n(k+1, n, N)  - 2*y_k_n(k, n, N) + y_k_n(k-1, n, N))/(double)(2*h*h)) + lambda_n(n,N) *y_k_n(k, n, N);

            if(fabs(res) > eps)
            {
                cout<<"k = "<<k<<endl<<"n = "<<n<<endl;
                return -1;
            }
        }
    }
    return 0;
}




