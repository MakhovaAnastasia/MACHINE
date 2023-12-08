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

#define eps 0.001

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
int find_p(void);

double x0(double x) {return pow(x, 0);}
double x1(double x) {return x;}
double x2(double x) {return pow(x, 2);}
double x3(double x) {return pow(x, 3);}
double x5(double x) {return pow(x, 5);}
double x9(double x) {return pow(x, 9);}

double costest(double x) { return cos(100*x);}
double exptest(double x) {return pow((double)M_E, (-100)*x);}
double exptest2(double x) {return pow(M_E, x);}
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
        find_p();
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
     double (*func[3])(double) = {costest, exptest, sqrttest};//
     double true_val[3] = {0, 0.001, M_PI}; //
     double a[3] = {0,0, -1+eps}; //
     double b[3] = {M_PI, 1, 1-eps};//
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
    double res = (deriv(x+eps, f, k-1) - deriv(x - eps, f, k-1))/(double)(2.*eps);
    if(res < eps)
    {
        res =0;
    }
    return res;
}


int find_p()
{

    int NUM[4] = {10,20,40,80};

    long double a = 0.;
    long double b = 0.;
    long double c = 0.;
    long double d = 0.;


    long double a_g = 0.;
    long double b_g = 0.;
    long double c_g = 0.;
    long double d_g = 0.;

    int flag_s = 0;
    int flag_g = 0;


    //exp2
    double true_val = (double)M_E - 1.0;
    double aa = 0;
    double bb = 1;


    ofstream out;
    out.open("3.txt");
    if(out.is_open())
    {
        for(int j = 0; j < 4; j++)
        {

            double max_s = fabs(Integrate_Simpson_N(aa, bb, exptest2, NUM[j]) - true_val);
            if(max_s <  0.00000000000001)
            {
                flag_s = 1;
            }

            if(flag_s ==0)
            {
                if(j==1)
                {
                    a = log((long double)NUM[j]);
                    b = log(1/max_s);
                }
                if(j==2)
                {
                    c = log((long double)NUM[j]);
                    d = log(1/max_s);
                }
            }

       //out<<setw(25)<<NUM[j]<<setw(25)<<log((long double)NUM[j])<<setw(25)<<log(1/max_s)<<setw(25)<<log(1/max_g)<<endl;
        }

        if(flag_s == 0)
        {

             if(fabs(c - a) > 0.00000000000001)
             {
                 cout<<"pSimpson = "<<(long double)(d - b)/(long double)(c - a)<<endl;

             }
             else{ flag_s = 1;}

         }
         if(flag_s == 1){
             cout<<"pSimpson нельзя вычислить"<<endl;
        }

        for(int j = 0; j < 4; j++)
        {

            double max_s = fabs(Integrate_Simpson_N(aa, bb, exptest2, NUM[j]) - true_val);
            double max_g = fabs(Integrate_Gauss_N(aa, bb, exptest2, NUM[j]) - true_val);

            if(max_g <  0.00000000000001)
            {
                flag_g = 1;
            }
            else{
                flag_g = 0;
            }

            if(flag_g == 0)
            {
                if(j==0)
                {
                    a_g = log((long double)NUM[j]);
                    b_g= log(1/max_g);
                }
                if(j==1)
                {
                    c_g = log((long double)NUM[j]);
                    d_g = log(1/max_g);
                }
                out<<setw(25)<<NUM[j]<<setw(25)<<log((long double)NUM[j])<<setw(25)<<log(1/max_s)<<setw(25)<<log(1/max_g)<<endl;
            }
            cout<<NUM[j]<<" "<<max_s<<" "<< max_g<<endl;

        }
        if(flag_g == 1)
        {
            if(fabs(c_g - a_g) >  0.00000000000001)
            {
                cout<<"pGauss = "<<(long double)(d_g - b_g)/(long double)(c_g - a_g)<<endl;

            }
            else{ flag_g = 1;}

        }
        if(flag_g == 1){
           // cout<<"pGauss нельзя вычислить"<<endl;
        }

        out.close();
        return 0;
    }
    out. close();
    return -1;
}
