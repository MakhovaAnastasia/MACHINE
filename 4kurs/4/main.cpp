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

#define eps 0.00000001

int test();
double f(double x);
double xn(int n, double x);
double xn_int(int n, double a, double b);
double Integrate_Simpson( double a,double b, double(*f)(double));
double Integrate_Gauss( double a,double b, double(*f)(double));
double Integrate_Simpson_N( double a,double b, double(*f)(double), int N);
double Integrate_Gauss_N( double a,double b, double(*f)(double), int N);
double Simpson_err( double a,double b, double(*f)(double), int N);
double Gauss_err( double a,double b, double(*f)(double), int N);
double deriv_max( double a,double b, double(*f)(double), int k);
double deriv(double x, double(*f)(double), int k);

double x0(double x) {return 1;}
double x1(double x) {return x;}
double x2(double x) {return pow(x, 2);}
double x3(double x) {return pow(x, 3);}
double x5(double x) {return pow(x, 5);}
double x9(double x) {return pow(x, 9);}



int main()
{
    double a = 0; 
    double b = 0;

    if(!(scanf("%lf",&a)==1)&&
        !(scanf("%lf",&b)==1))
    {
        //ошибка чтения
        return -1;
    }

    if(b < a)
    {
        double temp = a;
        a = b;
        b = temp;
    }
    test();


    return 0;
}

int test()
{
    int NUM[6] = {0,1,2,3,5,9};
    double (*func[6])(double) = {x0, x1, x2, x3, x5, x9};
    int a = 1;
    int b = 1.1;
    ofstream out;
    out.open("1.txt");
    if(out.is_open())
    {
        out<<setprecision(15)<<fixed;
        out<<setw(5)<<"n"<<setw(20)<<"Simpson"<<setw(20)<<"Simpson(formula)"<<setw(20)<<"Gauss"<<setw(20)<<"Gauss(formula)"<<endl;
        for(int j = 0; j < 6; j++)
        {
            double true_val = xn_int(NUM[j], a, b);
            out<<setw(5)<<NUM[j]<<setw(20)<<fabs(Integrate_Simpson(a, b, func[j]) - true_val)
            <<setw(20)<<Simpson_err(a, b, func[j], 1)
            <<setw(20)<<fabs(Integrate_Gauss(a, b, func[j]) - true_val)
            <<setw(20)<<Gauss_err(a, b, func[j], 1)<<endl;

        }
        out. close();
        return 0;
    }
    out.close();
    return -1;

}

double f(double x)
{
    return x;
}


double xn_int(int n, double a, double b)
{
    return ((pow(b,n+1)/(n+1)) - (pow(a,n+1)/(n+1)));
}

double Integrate_Simpson( double a,double b, double(*f)(double))
{
    return ((b-a)/6.)*(f(a)+ 4*f((a+b)/2.) +f(b));
}

double Integrate_Gauss( double a,double b, double(*f)(double))
{
    double x1 = (a+b)/2.;
    double x0 = (a+b)/2. - (b-a)/2. *sqrt(3./5.);
    double x2 = (a+b)/2. + (b-a)/2. *sqrt(3./5.);

    return ((b-a)/18.)*(5*f(x0)+ 8*f(x1) + 5*f(x2));
}

double Integrate_Simpson_N( double a,double b, double(*f)(double), int N)
{
    double h = (b - a)/(double) N;
    double res = 0.;
    for(int i = 0; i < N; i++)
    {
        res += Integrate_Simpson(a+(h*i), a+(h*(i+1)), f);
    }
    return res;
}

double Integrate_Gauss_N( double a,double b, double(*f)(double), int N)
{
    double h = (b - a)/(double) N;
    double res = 0.;
    for(int i = 0; i < N; i++)
    {
        res += Integrate_Gauss(a+(h*i), a+(h*(i+1)), f);
    }
    return res;
}

double Simpson_err( double a,double b, double(*f)(double), int N)
{
    return deriv_max(a,b,f,4)*pow(b-a, 5)/ (2880 * pow(N,4));
}

double Gauss_err( double a,double b, double(*f)(double), int N)
{
    return deriv_max(a,b,f,6)*pow(b-a, 7)/ (2016000* pow(N,6));
}

double deriv_max( double a,double b, double(*f)(double), int k)
{
    int L = 50;
    double h = (b - a)/(double) L;
    double max = -1;
    double res = -1;
    for(int i =  0; i < L; i++)
    {
        res = deriv(a+ h*i, f, k);
        if(fabs(res) > max)
        {
            max = fabs(res);
        }
    }
    return max;
}

double deriv(double x, double(*f)(double), int k)
{
    if(k == 1)
    {
        return (f(x+eps) - f(x - eps))/(2*eps);
    }

    return (deriv(x+eps, f, k-1) - deriv(x - eps, f, k-1))/(2*eps);
}