#include "Solve.h"
#include "ReadMatrix.h"
#include "synchronize.h"

void Solve(int n,double* A,double* b,double* x, double EPS, int thread_num,int total_threads, int* err)
{
    //умножаем слева на Tij
    for( int i = 0; i < n-1; i++)
    {
        for(int j = i+1; j < n; j++)
        {
            //Tij*A Tij*b
        if(*err != -1)
        {
            //определим угол поворота
            double x = A[i*n + i];
            double y = A[j*n + i];//[i*n + j];
            double root = sqrt(x*x + y*y);
            if(abs(root) < EPS)
            {
                *err = -1;
            }
            else{
                double cos_phi = x / root;
                double sin_phi =  -y / root;
                double xi = 0;
                double xj = 0;
                if(thread_num == 0 ){
                    //Tb
                    xi = b[i];
                    xj = b[j];

                    b[i] = xi*cos_phi - xj*sin_phi;
                    b[j] = xi*sin_phi + xj*cos_phi;
                    if(abs(b[i])< EPS)   b[i] = 0.0;
                    if(abs(b[j])< EPS)   b[j] = 0.0;
                    //cout<<"b["<<i<<"] = "<<b[i]<<" ("<<i<<j<<endl;
                }
                synchronize(total_threads);
                //умножение TA
                int first_row = ((n+1)*thread_num)/total_threads -1;
                if(first_row == -1) first_row = 0;
                int last_row = ((n+1)*(thread_num +1))/(total_threads) -1;
                //if((i == 0) &&(j == 1))
                //{
                //    cout<<first_row<<" - "<<last_row<<"  для "<<thread_num<<endl;
                //}
                for(int k = first_row; k < last_row; k++) //столбцы А
                {
                    xi = A[i*n+k];
                    xj = A[j*n + k];
                    if(k==i)
                    {
                        A[i*n + i] = root;
                        if(abs(A[i*n + i]) < EPS)
                        {
                            A[i*n + i] = 0;
                        }
                        A[j*n + i] = 0;
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

                }
            }
        
            synchronize(total_threads);
        }
    }
    synchronize(total_threads);
    // обратный ход метода Гаусса
    if((thread_num == 0)&&(*err != -1))
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
                *err = -1;
                break;
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

    return;
}

