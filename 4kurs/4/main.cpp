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

#define eps 0.01

int test1();
int test2(int N);
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

double x0(double x) {return pow(x, 0);}
double x1(double x) {return x;}
double x2(double x) {return pow(x, 2);}
double x3(double x) {return pow(x, 3);}
double x5(double x) {return pow(x, 5);}
double x9(double x) {return pow(x, 9);}

double costest(double x) { return cos(100*x);}
//double exptest(double x) {return pow((double)M_E, (-1000)*x);}
double sqrttest(double x) {return 1./(sqrt(1 - x*x));}



int main(int argc, char* argv[])
{
    double a = 0; 
    double b = 0;
    int N = -1;


    if(!((argc == 4)&&
    (sscanf(argv[1],"%lf",&a)==1)&&
    (sscanf(argv[2],"%lf",&b)==1)&&
    (sscanf(argv[3],"%d",&N)==1)))
    {
        //ошибка чтения
        return -1;
    }

    if(N <= 0)
    {
        cout<<a<<" "<<b<<" "<<N;
        cout<<"ошибка. изменено на N = 1"<<endl;
        N = 1;
    }

    if(b < a)
    {
        double temp = a;
        a = b;
        b = temp;
    }

/* Задание 1 */
    if(N==1)
    {
        cout<<setprecision(15)<<fixed;
        cout<<setw(25)<<" "<<setw(25)<<"Simpson"<<setw(25)<<"Gauss"<<endl;
        cout<<setw(25)<<"answer:"<<setw(25)<<Integrate_Simpson(a, b, f)<<setw(25)<<Integrate_Gauss(a, b, f)<<endl;
        cout<<setw(25)<<"err:"<<setw(25)<<Simpson_err(a, b, f, 1)<<setw(25)<<Gauss_err(a, b, f, 1)<<endl;
        test1();
    }
/* Задание 2 */
    else{
        cout<<setprecision(15)<<fixed;
        cout<<setw(25)<<" "<<setw(25)<<"Simpson"<<setw(25)<<"Gauss"<<endl;
        cout<<setw(25)<<"answer:"<<setw(25)<<Integrate_Simpson_N(a, b, f, N)<<setw(25)<<Integrate_Gauss_N(a, b, f, N)<<endl;
        cout<<setw(25)<<"err:"<<setw(25)<<Simpson_err(a, b, f, N)<<setw(25)<<Gauss_err(a, b, f, N)<<endl;
        test2(N);
    }

    return 0;
}

int test1()
{
    int NUM[6] = {0,1,2,3,5,9};
    double (*func[6])(double) = {x0, x1, x2, x3, x5, x9};
    double a = 1;
    double b = 1.1;
    ofstream out;
    out.open("1.txt");
    if(out.is_open())
    {
        out<<setprecision(15)<<fixed;
        out<<setw(5)<<"n"<<setw(25)<<"Simpson"<<setw(25)<<"Simpson(formula)"<<setw(25)<<"Gauss"<<setw(25)<<"Gauss(formula)"<<endl;
        for(int j = 0; j < 6; j++)
        {
            double true_val = xn_int(NUM[j], a, b);
            out<<setw(5)<<NUM[j]<<setw(25)<<fabs(Integrate_Simpson(a, b, func[j])- true_val)
            <<setw(25)<<Simpson_err(a, b, func[j], 1)
            <<setw(25)<<fabs(Integrate_Gauss(a, b, func[j])- true_val)
            <<setw(25)<<Gauss_err(a, b, func[j], 1)<<endl;

        }
        out. close();
        return 0;
    }
    out.close();
    return -1;

}

 int test2(int N)
 {
     double (*func[2])(double) = {costest,  sqrttest};//exptest,
     double true_val[3] = {0,  M_PI}; //0.001,
     double a[3] = {0,  -1}; //0,
     double b[3] = {M_PI, 1};//1,
     ofstream out;
     out.open("2.txt");
     if(out.is_open())
     {
         out<<setprecision(15)<<fixed;
         out<<setw(5)<<"func"<<setw(25)<<"Simpson_N"<<setw(25)<<"Simpson_N(formula)"<<setw(25)<<"Gauss_N"<<setw(25)<<"Gauss_N(formula)"<<endl;
         for(int j = 0; j < 3; j++)
         {
             out<<setw(5)<<j<<setw(25)<<fabs(Integrate_Simpson_N(a[j], b[j], func[j], N) - true_val[j])
             <<setw(25)<<Simpson_err(a[j], b[j], func[j], N)
             <<setw(25)<<fabs(Integrate_Gauss_N(a[j], b[j], func[j], N) - true_val[j])
             <<setw(25)<<Gauss_err(a[j], b[j], func[j], N)<<endl;

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
    return ((pow(b,n+1)/(double)(n+1)) - (pow(a,n+1)/(double)(n+1)));
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
    return deriv_max(a,b,f,4)*pow(b-a, 5) / ((double)2880 * pow((double)N,4));
}

double Gauss_err( double a,double b, double(*f)(double), int N)
{
    return deriv_max(a,b,f,6)*pow(b-a, 7) / ((double)2516000* pow((double)N,6));
}

double deriv_max( double a,double b, double(*f)(double), int k)
{
    int L = 100;
    double h = (b - a- 2*(k+1)*eps)/(double) L;
    double max = -1;
    double res = -1;
    for(int i =  0; i < L; i++)
    {
        res = deriv(a+(k+1)*eps+ h*i, f, k);
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
        return (f(x+eps) - f(x-eps))/(double)(2.*eps);
    }

    return (deriv(x+eps, f, k-1) - deriv(x - eps, f, k-1))/(double)(2.*eps);
}
