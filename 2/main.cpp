
#include "ReadMatrix.h"
#include "Solve.h"
#include <cmath>
#include <ctime>
#include <cstdlib>

double n_1(double* A, double* x, int N);
double n_2(double* A, double* x, int N);


int main(int argc, char* argv[])
{
    double* A ;
    double* x;
    double* Q_cos;
    double* Q_sin;
    
    string filename;

    int n = -1; //размер матрицы
    int m = -1; //размер вывода
    double EPS = 0.0; //точность нахождения с.з.
    int k = -1; // формула

    if(!(((argc == 5)||(argc == 6))&&
    (sscanf(argv[1],"%d",&n)==1)&&
    (sscanf(argv[2],"%d",&m)==1)&&
    (sscanf(argv[3],"%le",&EPS)==1)&&
    (sscanf(argv[4],"%d",&k)==1)))
    {
        //ошибка чтения
        return -1;
    }
    if((n<=0)||(m<=0)||(n<m))
    {
        cout<<"0<m<=n!!!"<<endl;
        return -1;
    }
    if((k == 0)&&(argc == 6))
    {
        filename = argv[5];
    }

    A = (double*) malloc(sizeof(double)*n*n);
    x = (double*) malloc(sizeof(double)*n);
    Q_cos = (double*) malloc(sizeof(double)*(n-1));
    Q_sin = (double*) malloc(sizeof(double)*(n-1));
    
    if(ReadMatrix(A,n,k, filename)!=0)
    {
        cout<<"ошибка чтения"<<endl;
        free(A);
        free(x);
        free(Q_cos);
        free(Q_sin);
        return 0;
    }
    for(int i = 0; i<n ; i++)
    {
        x[i] = 0;
        if(i != n-1)
        {
            Q_cos[i] = 0;
            Q_sin[i] = 0;
        }

    }

    cout<<"A--------"<<endl;
    PrintMatrix(A, n, n, m);

    int start = clock();
    int res = Solve(n, A, x, EPS, Q_cos, Q_sin);
    if(res == -1)
    {
        cout<<"делим на ноль. пришлось выйти"<<endl;
    }
    int end = clock(); 
    int time = (end - start)/(CLOCKS_PER_SEC/100);// время работы  в секундах
    cout<<"время работы(сотые доли сек.): "<<time<<endl;

    //printf("невязка в первом инварианте:  %10.3e\n",n_1(A,x,n));
    //printf("невязка во втором инварианте: %10.3e\n",n_2(A,x,n));


    cout<<"x(собственные значения)--------"<<endl;
    PrintMatrix(x, 1, n, m);
    cout<<"A--------"<<endl;
    PrintMatrix(A, n, n, m);

    free(A);
    free(x);
    free(Q_cos);
    free(Q_sin);

    return 0;
}

double n_1(double* A, double* x, int N)
{
    double sum = 0.;
    for(int i  = 0; i< N; i++)
    {
        sum+=(A[N*i+i] - x[i]);
    }
    return abs(sum);
}

double n_2(double* A, double* x, int N)
{
    double length = 0.;
    double sz = 0.;
    for(int i = 0; i< N*N; i++)
    {
        length+=(A[i]*A[i]);
        if(i<N)
        {
            sz+= (x[i]*x[i]);
        }
    }
    length = sqrt(length);

    return abs(length - sz);
}
