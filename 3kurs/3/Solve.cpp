#include "Solve.h"
#include "ReadMatrix.h"
#include "synchronize.h"

int Solve(int n,double* A,double* b,double* x, double EPS, int thread_num,int total_threads, double* Cos_p, double* Sin_p)
{
    //int err_code = 0;
    //умножаем слева на Tij
    //double cos_phi= 0.0;
    //double sin_phi = 0.0;
    double xi = 0;
    double xj = 0;
    int ncols = n/total_threads;
    if ((n%total_threads)!=0)
        ncols++;

    int first_row ;//= ((n)*thread_num)/total_threads -1;
    first_row = ncols*thread_num;
        if(first_row == -1) first_row = 0;
    int  last_row ;//= ((n)*(thread_num +1))/(total_threads) -1;
    last_row = ncols*(thread_num+1);
    if(last_row > n)
        last_row = n;
    printf("%d - %d- %d\n",first_row, last_row, thread_num);
    double xx = 0.0;
    double y= 0.0;
    double root= 0.0;
    //synchronize(total_threads);
    for( int i = 0; i < n-1; i++)
    {
        //-----------
        if((i>= first_row) &&(i< last_row))
        {
            for(int j = i+1; j < n; j++)
            {
                xx = A[i*n + i];
                y = A[j*n + i];//[i*n + j];
                root = sqrt(xx*xx+ y*y);
                if(abs(root) < EPS)
                {
                    Cos_p[i] = -100;
                     break;
                }
                Cos_p[j] = xx / root;
                Sin_p[j] =  -y / root;
                xi = A[i*n+i];
                xj = A[j*n + i];
                    A[i*n + i] = root;
                    //if(abs(A[i*n + i]) < EPS)
                    //{
                    //    A[i*n + i] = 0;
                    //}
                    A[j*n + i] = 0;

            }
          }
        synchronize(total_threads);
        if(abs(Cos_p[i]) > 2)
        {
            return -1;
        }

            for(int j = i+1; j < n; j++)
            {
                if(thread_num == 0 ){
                    //Tb
                    xi = b[i];
                    xj = b[j];
                    b[i] = xi*Cos_p[j] - xj*Sin_p[j];
                   b[j] = xi*Sin_p[j] + xj*Cos_p[j];
                    if(abs(b[i])< EPS)   b[i] = 0.0;
                    if(abs(b[j])< EPS)   b[j] = 0.0;
                }
                //умножение TA
                for(int k = first_row; k < last_row; k++) //столбцы А
                {
                    xi = A[i*n+k];
                    xj = A[j*n + k];
                    if(k==i)
                    {
                        continue;
                    }
                    else
                    {
                        A[i*n + k] = xi*Cos_p[j] - xj*Sin_p[j];
                        if(abs(A[i*n + k]) < EPS)
                        {
                            A[i*n + k] = 0;
                        }
                        A[j*n + k] = xi*Sin_p[j] + xj*Cos_p[j];
                        if(abs(A[j*n + k]) < EPS)
                        {
                            A[j*n + k] = 0;
                        }
                    }
                }
            }
                         synchronize(total_threads);
        /*
        for(int j = i+1; j < n; j++)
        {
            //Tij*A Tij*b
            //определим угол поворота
            //synchronize(total_threads);
            xx = A[i*n + i];
            y = A[j*n + i];//[i*n + j];
            root = sqrt(xx*xx+ y*y);
            if(abs(root) < EPS)
            {
                return -1;
            }
            cos_phi = xx / root;
            sin_phi =  -y / root;
            if(thread_num == 0 ){
                //Tb
                xi = b[i];
                xj = b[j];
                b[i] = xi*cos_phi - xj*sin_phi;
               b[j] = xi*sin_phi + xj*cos_phi;
                if(abs(b[i])< EPS)   b[i] = 0.0;
                if(abs(b[j])< EPS)   b[j] = 0.0;
            }
            synchronize(total_threads);
            //умножение TA
            for(int k = first_row; k < last_row; k++) //столбцы А
            {
                xi = A[i*n+k];
                xj = A[j*n + k];
                if(k==i)
                {
                    A[i*n + k] = root;
                    //if(abs(A[i*n + i]) < EPS)
                    //{
                    //    A[i*n + i] = 0;
                    //}
                    A[j*n + k] = 0;
                }
                else
                {
                    A[i*n + k] = xi*cos_phi - xj*sin_phi;
                    if(abs(A[i*n + k]) < EPS)
                    {
                        A[i*n + k] = 0;
                    }
                    A[j*n + k] = xi*sin_phi + xj*cos_phi;
                    if(abs(A[j*n + k]) < EPS)
                    {
                        A[j*n + k] = 0;
                    }
                }
            } 

            synchronize(total_threads);
        }*/
    }
    synchronize(total_threads);
    // обратный ход метода Гаусса
    if((thread_num == 0))//&&(err_code != -1))
    {
        for( int i = n-1; i >= 0; i--)
        {
            x[i] = b[i];
            for(int j = i+1; j < n; j++)
            {
                x[i] -= A[i*n +j]* x[j];
            }
            if(abs(A[i*n + i]) < EPS) //на диагонали  нашли 0. матрица вырождена
            {
                //cout<<"нет точного ответа. x["<<i<<"] = 1 "<<endl;
                //err_code = -1;
                return -1;
            }
            else{
                x[i]/= A[i*n +i];
                if((abs(x[i]))<EPS)
                {
                    x[i] = 0;
                }
                //cout<<"x["<<i<<"] = "<<x[i]<<endl;
            }
        }
    }

    return 0;
}

