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

double f(double x1, double x2);
int Triangular(double Lx,double Ly,int Nx,int Ny,double** V,int** T,int** vR,int** gR);
double Integrate(int N, double **V,int** T);
int find_p(double **V,int** T,int** vR,int** gR);

int main(int argc, char* argv[])
{
    double Lx = 1;
    double Ly = 1;
    int Nx1 = 100;
    int Ny1 = 100;
    int Nx = 100;
    int Ny = 100;

    if(!((argc == 5)&&
    (sscanf(argv[1],"%lf",&Lx)==1)&&
    (sscanf(argv[2],"%lf",&Ly)==1)&&
    (sscanf(argv[3],"%d",&Nx)==1)&&
    (sscanf(argv[4],"%d",&Ny)==1)))
    {
        //ошибка чтения
        return -1;
    }
    if(Lx <= 0)
    {
        cout<<"ошибка. изменено на abs(Lx)"<<endl;
        Lx = fabs(Lx);
    }
    if(Ly <= 0)
    {
        cout<<"ошибка. изменено на abs(Ly)"<<endl;
        Ly = fabs(Ly);
    }
    if(Nx <= 0)
    {
        cout<<"ошибка. изменено на Ny = 1"<<endl;
        Ny = 1;
    }

    if(Nx <= 0)
    {
        cout<<"ошибка. изменено на Nx = 1"<<endl;
        Nx = 1;
    }


    double **V;
    int **T, **vR, **gR;

    V = (double**) malloc(sizeof(double*)*((Nx1 + 1)*(Ny1 +1))); //вершины
    for(int i = 0; i< (Nx1 + 1)*(Ny1 +1); i++)
    {
        V[i] = (double*) malloc(sizeof(double)*(2));
    }

    T = (int**) malloc(sizeof(int*)*((Nx1)*(Ny1)*2)); //треугольники
    for(int i = 0; i< ((Nx1)*(Ny1)*2); i++)
    {
        T[i] = (int*) malloc(sizeof(int)*(3));
    }

    vR = (int**) malloc(sizeof(int*)*((Nx1*Ny1*4+Ny1) - ((Nx1 +Ny1)*2))); //ребра внутренние
    for(int i = 0; i< (Nx1*Ny1*4+Ny1) - ((Nx1 +Ny1)*2); i++)
    {
        vR[i] = (int*) malloc(sizeof(int)*(2));
    }

    gR = (int**) malloc(sizeof(int*)*((Nx1 +Ny1)*2)); //ребра граничные
    for(int i = 0; i< (Nx1 +Ny1)*2; i++)
    {
        gR[i] = (int*) malloc(sizeof(int)*(2));
    }



    find_p(V,T,vR,gR);


    /*   Триангуляция */
    Triangular(Lx, Ly, Nx, Ny,V,T,vR,gR);

    //Интеграл на [0,1]x[0, 1]

    double Sn = Integrate(Nx,V, T);
    cout<<"Sn = "<<Sn<<endl;
    cout<<"R = "<<scientific<<fabs(Sn  - (13/12.))<<endl;



    for(int i = 0; i< (Nx1 + 1)*(Ny1 +1); i++)
    {
        free(V[i]);
    }
    free(V);

    for(int i = 0; i< ((Nx1)*(Ny1)*2); i++)
    {
        free(T[i]);
    }
    free(T);
     
    for(int i = 0; i< (Nx1*Ny1*4+Ny1) - ((Nx1 +Ny1)*2); i++)
    {
        free(vR[i]);
    }
    free(vR);
    
    for(int i = 0; i< (Nx1 +Ny1)*2; i++)
    {
        free(gR[i]);
    }
    free(gR);

    return 0;
}


int Triangular(double Lx,double Ly,int Nx,int Ny,double **V,int** T,int** vR,int** gR)
{
    for(int i = 0; i < Ny+1;i++)
    {
        for(int j = 0; j < Nx+1; j++)
        {
            V[(i*(Nx+1))+j][0] = j*(Lx/(double)Nx);
            V[(i*(Nx+1))+j][1] = Ly - i*(Ly/(double)Ny);
        }
    }
    int k = 0;
    int n = 0;
    int m = 0;

    for(int i = 0; i < Ny;i++)
    {
        for(int j = 0; j < Nx; j++)
        {
            k = i*(Nx)*2 +2*j;
            //upper triangle
            T[k][0]= i*(Nx+1)+j;

            T[k][1]= i*(Nx+1)+(j+1);

            T[k][2]= (i+1)*(Nx+1)+j;

            if((i == 0))
            {
                if(j == 0)
                {
                    gR[n][0] = T[k][0];
                    gR[n][1] = T[k][2];
                    n ++;
                }
                else{
                    vR[m][0] = T[k][0];
                    vR[m][1] = T[k][2];
                    m++;
                }

                gR[n][0] = T[k][0];
                gR[n][1] = T[k][1];
                n++;
            }
            else
            {
                if(j == 0)
                {
                    gR[n][0] = T[k][0];
                    gR[n][1] = T[k][2];
                    n ++;
                }
                else{
                    vR[m][0] = T[k][0];
                    vR[m][1] = T[k][2];
                    m++;
                }

                vR[m][0] = T[k][0];
                vR[m][1] = T[k][1];
                m++;
            }

            vR[m][0] = T[k][1];
            vR[m][1] = T[k][2];
            m++;
            k++;

            //lower triangle

            T[k][0]= i*(Nx+1)+(j+1);

            T[k][1]= (i+1)*(Nx+1)+j;

            T[k][2]= (i+1)*(Nx+1)+(j+1);

            if((i == Ny -1))
            {
                if(j == Nx - 1)
                {
                    gR[n][0] = T[k][0];
                    gR[n][1] = T[k][2];
                    n ++;
                }

                gR[n][0] = T[k][1];
                gR[n][1] = T[k][2];
                n++;
            }
            else
            {
                if(j == Nx - 1)
                {
                    gR[n][0] = T[k][0];
                    gR[n][1] = T[k][2];
                    n ++;
                }

                //vR[m][0] = T[k][1];
                //vR[m][1] = T[k][2];
                //m++;
            }
            
            //vR[m][0] = T[k][0];
            //vR[m][1] = T[k][1];
            //m++;
            k++;
        }
    }

    ofstream out;
    out.open("1.txt");
    if(out.is_open())
    {
        out<<(Nx + 1)*(Ny +1)<<endl;//число вершин
        out<<(Nx*Ny*2)<<endl; //число треугольников
        out<<m<<endl;
        out<<n<<endl;
        out<<endl;
        for(int i = 0; i < (Nx + 1)*(Ny +1); i++)
        {
            out<<i<<" : < "<< V[i][0]<<", "<<V[i][1]<<" >"<<endl;
        }
        out<<endl;
        for(int i = 0; i < (Nx*Ny*2); i++)
        {
            out<<i<<" : < "<< T[i][0]<<", "<<T[i][1]<<", "<<T[i][2]<<" >"<<endl;
        }
        out<<endl;
        for(int i = 0; i < m; i++)
        {
            out<<i<<" : < "<< vR[i][0]<<", "<<vR[i][1]<<" >"<<endl;
        }
        out<<endl;
        for(int i = 0; i < n; i++)
        {
            out<<i<<" : < "<< gR[i][0]<<", "<<gR[i][1]<<" >"<<endl;
        }


        out. close();
        return 0;
    }
    out.close();
    return -1;
}


double f(double x1, double x2)
{
    //return (pow(x1,4) + pow(x1,2)*pow(x2,2)+ pow(x2,4));
    return x1+pow(x1,2) +x1*x2;
}

double Integrate(int N, double **V,int** T)
{
    double sum = 0;
    double A1,A2,B1,B2,C1,C2;
    double S = (1/(double)N)*(1/(double)N)*0.5;
    for(int i = 0; i <(N*N*2);i++)
    {
        A1 = fabs(V[T[i][0]][0] + V[T[i][1]][0])/2.;
        A2 = fabs(V[T[i][0]][1] + V[T[i][1]][1])/2.;
        B1 = fabs(V[T[i][1]][0] + V[T[i][2]][0])/2.;
        B2 = fabs(V[T[i][1]][1] + V[T[i][2]][1])/2.;
        C1 = fabs(V[T[i][2]][0] + V[T[i][0]][0])/2.;
        C2 = fabs(V[T[i][2]][1] + V[T[i][0]][1])/2.;
        sum+=(1/3.)*S*(f(A1,A2)+f(B1,B2)+f(C1,C2));
    }
    return sum;
}

int find_p(double **V,int** T,int** vR,int** gR)
{
    double true_val = 13/12.;//23/45.;

    int NUM[4] = {10,20,40,80};

    long double a = 0.;
    long double b = 0.;
    long double c = 0.;
    long double d = 0.;
    int flag = 0;
    ofstream out;
    out.open("3.txt");
    if(out.is_open())
    {
        for(int j = 0; j < 4; j++)
        {
            Triangular(1,1,NUM[j],NUM[j],V,T,vR,gR);
            double max = fabs(Integrate(NUM[j], V, T) - true_val);
            if(max <  0.00000000000001)
            {
                flag= 1;
            }

            if(flag ==0)
            {
                if(j==1)
                {
                    a = log((long double)NUM[j]);
                    b = log(1/max);
                }
                if(j==2)
                {
                    c = log((long double)NUM[j]);
                    d = log(1/max);
                }
                out<<setw(25)<<NUM[j]<<setw(25)<<log((long double)NUM[j])<<setw(25)<<log(1/max)<<endl;
            }

        }    
            if(flag == 0)
            {
                if(fabs(c - a) > 0.00000000000001)
                {
                    cout<<"p = "<<(long double)(d - b)/(long double)(c - a)<<endl;
                }
                else{ flag = 1;}

            }
            if(flag == 1){
            cout<<"p нельзя вычислить" <<endl;
            }
        out.close();
        return 0;
    }
    out. close();
    return -1;
}
