
#include "ReadMatrix.h"
#include "Solve.h"
#include <cmath>
#include <ctime>
#include <cstdlib>

double norma_nevyaski(double* A, double*b, double* x, int N);
double norma_pogreshnosty(double* x,double* x_real, int N);


int main(int argc, char* argv[])
{
    double* A;
    double* b;
    double* x;
    double* x_real;
    
    string filename;

    int n = -1; //размер матрицы
    int m = -1; //размер вывода
    int k = -1; // формула

    if(!(((argc == 4)||(argc == 5))&&
    (sscanf(argv[1],"%d",&n)==1)&&
    (sscanf(argv[2],"%d",&m)==1)&&
    (sscanf(argv[3],"%d",&k)==1)))
    {
        //ошибка чтения
        return -1;
    }
    if((n<=0)||(m<=0)||(n<m))
    {
        cout<<n<<m<<k<<endl;
        cout<<"0<m<=n!!!"<<endl;
        return -1;
    }
    if((k == 0)&&(argc == 5))
    {
        filename = argv[4];
    }

    A = (double*) malloc(sizeof(double)*n*n);
    if (A==NULL) exit (1);
    b = (double*) malloc(sizeof(double)*n);
    x = (double*) malloc(sizeof(double)*n);
    x_real = (double*) malloc(sizeof(double)*n);

    
    if(ReadMatrix(A,n,k, filename)!=0)
    {
        cout<<"ошибка чтения"<<endl;
        free(A);
        free(b);
        free(x);
        free(x_real);
        return 0;
    }

    for(int i = 0; i < n; i++)
    {
        b[i] = 0;
        for( int j = 0; j <= ((n+1)/2)-1; j++)
        {
            b[i] += A[(n*i) + (2*j)];
        }

        x[i] = 0;
        x_real[i] = (i+1)%2;
    }
    cout<<"A--------"<<endl;
    PrintMatrix(A, n, n, m);
    cout<<"b--------"<<endl;
    PrintMatrix(b, 1, n, m);
    cout<<"x--------"<<endl;
    PrintMatrix(x, 1, n, m);

    int start = clock();
    int res = Solve(n, A, b, x);
    if(res == -1)
    {
        cout<<"делим на ноль. пришлось выйти"<<endl;
    }
    int end = clock(); 
    int time = (end - start)/(CLOCKS_PER_SEC/100);// время работы  в секундах
    cout<<"время работы(сотые доли сек.): "<<time<<endl;
if(res == 0)
{
    printf("норма невязки:  %10.3e\n",norma_nevyaski(A,b,x,n));
    printf("норма погрешности: %10.3e\n",norma_pogreshnosty(x,x_real, n));

    cout<<"b--------"<<endl;
    PrintMatrix(b, 1, n, m);
    cout<<"x--------"<<endl;
    PrintMatrix(x, 1, n, m);
    cout<<"A--------"<<endl;
    PrintMatrix(A, n, n, m);
}
    free(A);
    free(b);
    free(x);
    free(x_real);
    return 0;
}

double norma_nevyaski(double* A, double*b, double* x, int N)
{   //  ||Ax-b|| / ||b||
    double sum = 0; //ищем максимум по суммам модулей коэф. строк.
    double b_max = 0;
    double max  = 0;
    for( int i = 0; i < N; i++) // ||Ax-b|| и ||b||
    {
        for(int j =  0; j< N; j++)
        {
            sum+=A[i*N + j]* x[j];
        }
        sum -= b[i];
        if(( abs(b[i]) - b_max) > EPS) // ||b||
        {
            b_max = abs(b[i]);
        }
        sum = abs(sum);

        if( (sum - max) > EPS)
        {
            max = sum;
        }
        sum = 0;
    }
    if(abs(b_max) > EPS)
    {
        max/= b_max; //норма
    }
    else {
        cout<< "||b|| = 0"<<endl;
        return -1;
    }
    //cout<< "норма невязки: "<<scientific<<max<< endl; 
    return max;
}

double norma_pogreshnosty(double* x,double* x_real, int N)
{
    double norma  = 0;
    for( int i = 0; i < N; i++) // ||x - x_real||
    {
        if(( abs(x[i] - x_real[i])- norma) > EPS)
        {
            norma = abs(x[i] - x_real[i]);
        }
    }
    //cout<< "норма погрешности: "<<scientific<< norma<< endl; 
    return norma;
}
