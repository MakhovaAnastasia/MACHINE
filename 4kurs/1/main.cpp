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

long double f(long double x); //функция
int first_table(long double* X, long double a, long double b, int n, int k); //создаем первый файл
int second_table(long double*X,long double*XX,long double*P,
                long double* A,long double* c,
                long double*L, int n); //создаем второй файл
int Write1(long double* X, int n); //вывод в первый файл
int Write2(long double* XX, long double* P, long double* L, int n);//вывод во второй
int Generate(long double* X, int a, int b, int n, int k); //создадим узлы x_i
int Ln(long double* XX, long double* L,int n); //считаем полином в форме Лагранжа
int Pn(long double* XX, long double* P,long double* A,long double* c,int n); //считаем полином в алгебраической форме
long double PHI(long double* XX,int i, int j, int n); //функция Ф

int main(int argc, char* argv[])
{
    long double* X, *XX, *P, *L, *A, *c;
    long double a = 0; //отрезок [a,b]
    long double b = 0; 
    int n = -1; // число узлов
    int k = -1; // тип узлов

    if(!((argc == 5)&&
    (sscanf(argv[1],"%Lf",&a)==1)&&
    (sscanf(argv[2],"%Lf",&b)==1)&&
    (sscanf(argv[3],"%d",&n)==1)&&
    (sscanf(argv[4],"%d",&k)==1)))
    {
        //ошибка чтения
        return -1;
    }

    if(b < a) //если концы в другом порядке, то поменяем их
    {
        long double buf = a;
        a = b;
        b = buf;
    }
    if((n<=0)||(k>=4)||(k<=0))
    {
        //ошибка 
        cout<<n<<' '<<k<<endl;
        cout<<"k = 1, 2, 3  !!!  n >= 0"<<endl;
        return -1;
    }
    
    A = (long double*) malloc(sizeof(long double)*n*n);
    if (A==NULL) exit (1);
    c = (long double*) malloc(sizeof(long double)*n);

    X = (long double*) malloc(sizeof(long double)*n);
    XX = (long double*) malloc(sizeof(long double)*3*n);
    P = (long double*) malloc(sizeof(long double)*3*n);
    L = (long double*) malloc(sizeof(long double)*3*n);

    // 1-ая таблица
    first_table(X, a, b, n, k);

    // 2-ая таблица
    second_table(X, XX, P, A, c, L, n);

    free(X);
    free(XX);
    free(P);
    free(L);

    return 0;
}

int first_table(long double*X, long double a, long double b, int n, int k)
{
     //создаем узлы
    if(Generate(X,a,b,n,k)!= k)
    {
        cout<<"error in Generator"<<endl;
        return -1;
    }

    //запишем в файл результат
    
    Write1(X, n);

    return 0;
}

int Generate(long double* X, int a, int b, int n, int k)
{
    if(k==1) // равноотстоящие точки
    {
        for(int i = 0; i < n; i++)
        {
            X[i] = a + (((b-a)*i) / (long double)(n-1));
        }
        return 1;
    }
    if(k==2) // чебышевские
    {
        for(int i = 0; i < n; i++)
        {
            X[n-1-i] = (a+b)/2.0 + (b-a)* cosl((2.*i+1.)* M_PI / (2.0*n) )/2.0;
        }
        return 2;
    }
    if(k==3) // случайные
    {
        X[0] = a;
        X[n-1] = b;
        for(int i = 1; i < n-1; i++)
        {
            X[i] =((long double)(rand()) /(long double)(RAND_MAX))*(b-a)+ a;
        }
        sort(X,X+n);
        return 3;
    }
    return -1;
}



int second_table(long double*X,long double*XX,long double*P,
                long double* A,long double* c,
                long double*L, int n)
{
     //создаем доп узлы
    for(int i = 0; i < n-1; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            XX[(3*i)+j] = X[i] + (((X[i+1]-X[i])*j) / (long double)(3));
        }
        
    }
    XX[3*(n-1)] = X[n-1];

    //P
    Pn(XX, P, A, c, n);
    //L
    Ln(XX,L,n);

    //запишем в файл результат
    
    Write2(XX, P, L, n);

    return 0;
}

int Pn(long double* XX, long double* P,long double* A,long double* c,int n)
{
    //матрица Ван дер Монда
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            if(j == 0)
            {
                A[(n*i)+j] = 1;
            }
            else
            {
                A[(n*i)+j] = A[(n*i)+j-1] * XX[3*i];
            }
        }
        c[i] = f(XX[i*3]);
    }

    //прямой ход метода Гаусса
    for(int j = 0; j<n;j++)
    {
        for(int i = j+1;i<n;i++)
        {
            long double temp = A[n*i+j]/A[j*n+j];
            
            for(int k = j; k< n;k++)
            {
                A[i*n+k] -= A[j*n+k]*temp;
            }
            c[i] -= temp *c[j];
        }
        
    }

    // обратный ход метода Гаусса
    for( int i = n-1; i >= 0; i--)
    {
        P[i] = c[i];
        for(int j = i+1; j < n; j++)
        {
            P[i] -= A[i*n +j]* P[j];
        }
        P[i]/= A[i*n +i];
    }

    //считаем значения в точках
    for(int j = 0; j<n;j++)
    {
        A[j] = P[j];
    }

    for(int i = 0; i<= 3*(n-1);i++)
    {
        long double power = 1.;
        for(int j = 0; j < n; j++)
        {
            P[i]= P[i]+(A[j]*power);
            power= power* XX[i];
        }
    }

    return 0;
}

int Ln(long double* XX, long double* L,int n)
{
    long double sum = 0;

    for(int i = 0; i <= 3*(n-1); i++)
    {
        for(int j = 0; j < n ; j++)
        {
            sum += f(XX[j*3]) * PHI(XX, i, j, n);
        }
        L[i] = sum;
        sum  = 0;
    }
    return  0;
}

long double PHI(long double* XX,int i, int j, int n)
{
    long double prod = 1;
    for(int k = 0; k < n ; k++)
    {
        if(k!=j)
        {
            prod *= (XX[i]-XX[3*k])/((XX[3*j]-XX[3*k]));
        }
    }
    return prod;
}


long double f(long double x)
{
    //return -x*x*x + 4*x*x+ x-8;

    //return abs(x);
    return 1/(1+25*x*x);
}

int Write1(long double* X, int n)
{
    ofstream out;
    out.open("1.txt");
    if(out.is_open())
    {
        out<<setprecision(15)<<fixed;
        for(int i = 0; i < n; i++)
        {
            out<<setw(25)<<X[i]<<setw(25)<<f(X[i])<<endl;
        }
        out. close();
        return 0;
    }
    out. close();
    return -1;
}

int Write2(long double* XX, long double* P, long double* L, int n)
{
    ofstream out;
    out.open("2.txt");
    if(out.is_open())
    {
        out<<setprecision(15)<<fixed;
        for(int i = 0; i <3*n-2; i++)
        {
            out<<setw(25)<<XX[i]<<setw(25)<<L[i]<<setw(25)<<P[i]<<setw(25)<<f(XX[i])<<endl;
        }
        out. close();
        return 0;
    }
    out. close();
    return -1;
}
