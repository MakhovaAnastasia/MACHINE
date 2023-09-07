#include "Solve.h"
#include "ReadMatrix.h"
//problem- 2

int Solve(const int n,double* A,double* x, double EPS, double* Q_cos, double* Q_sin)
{
    //норма
    double nA  = 0;
    for(int i = 0; i< n; i++)
    {
        double sum = 0;
        for( int j = 0; j < n; j++)
        {
            sum += abs(A[i*n +j]);
        }
        if((i == 0)||(nA < sum))
        {
            nA = sum;
        }
    }
    //почти треугольный вид
    for(int k = 0;k < n-2; k++) //n-2 steps
    {
        double s= 0.;
        double norm_x = 0;
        double norm_a = 0;
        for(int j = k+2; j < n; j++)
        {
            s+=(abs(A[j*n+k]))*(abs(A[j*n+k]));
        }
        norm_a = sqrt(s + (abs(A[(k+1)*n+k])*abs(A[(k+1)*n+k])));
        //x
        for(int j = k+1; j < n; j++)
        {
            x[j] = A[(n*j) + k];
            if(j==k+1)//x[0]
            {
                
                x[j]-=norm_a;
                norm_x = sqrt((abs(x[j])*abs(x[j]))+s);
                if(norm_x < EPS)
                {
                     for(int i = k+1; i < n; i++)
                     {
                         x[i] = 0;
                     }
                     break;
                }
            }
            x[j]/= norm_x;
        }
        //U(x)AU(x)
        for(int j = k; j<n; j++)//столбец
        {
            double scalar = 0;

            for(int i = k+1; i < n; i++)//строка
            {
                            //(||a||,0 0 0 ... 0 )
                if((i==k+1) && (j == k))
                {
                    A[i*n+j] = norm_a;
                    continue;
                }  
                if((i!=k+1) &&(j == k))
                {
                    A[i*n+j] = 0;
                    continue;
                }
                if(j!= k) 
                {
                    scalar+=A[i*n+j]*x[i];
                }
            }
            //U(x)A по столбцам
            for(int i = k+1; i < n; i++)//строки
            {   
                if(j != k)
                {
                    A[i*n+j] = A[i*n+j]-2*scalar*x[i];
                    if(abs(A[i*n+j])< EPS)
                    {
                        A[i*n+j] = 0;
                    }
                }
            }
        }
            //AU(x)~~ по лемме 11 по строкам
            for(int i = 0; i < n; i++)//строка
            {   
                double scalar = 0;
                for(int j = k+1; j<n; j++)
                {
                    scalar+=A[i*n+j]*x[j];
                }
                for(int j = k+1; j<n; j++)
                {
                    A[i*n+j] = A[i*n+j]-2*scalar*x[j];
                    if(abs(A[i*n+j])< EPS)
                    {
                        A[i*n+j] = 0;
                    }
                }
            }
                //PrintMatrix(A, n, n, n);
    }

//QR-rotate
    int k = n-1;
    int steps = 0;
    if(n> 2){
    do{
        //double s = 0;//A[k*n +k];
        //QR со сдвигом: A -sI = QR
        //R = ПT(i,i+1)(A-sI)
            //cout<<" s ="<< s<<" k = "<<k<<endl;
        //сдвиг
        //for(int j = 0; j<= k; j++)
        //{
        //    A[j*n+j] -=s;
        //}
            //PrintMatrix(A, n, n, n);

        for(int i = 0; i< k; i++)
        {
            int res=TA(i, i+1, A,n,k+1,EPS, Q_cos, Q_sin);
                //   cout<<"!!"<<endl;
                //PrintMatrix(A, n, n, n);
            if(res == -1)
            {
                //cout<<"do";
                return -1;
            }
        }
        //A = RQ +sI = RПT(i,i+1) +sI
        //for(int i = k-1; i>=0; i--)
        for(int i = 0; i< k; i++)
        {
            int res=T2(i, i+1, A,n,EPS, Q_cos, Q_sin);
                //cout<<"??"<<Q_cos[i]<<" "<<Q_sin[i]<<endl;
                //PrintMatrix(A, n, n, n);
            if(res == -1)
            {
                    // cout<<"do2";
                return -1;
            }
        }

        //обратный сдвиг
        //for(int j = 0; j<= k; j++)
        //{
        //    A[j*n+j] +=s;
        //}

                //cout<<abs(A[k*n +(k-1)]) <<" "<<EPS*nA<<endl;
                //PrintMatrix(A, n, n, n);
                //cout<<"x["<<k<<"]="<<A[k*n+k]<<endl;
                //cout<<"nA = "<<nA<<" A="<<A[k*n +(k-1)]<<endl;
        if(abs(A[k*n +(k-1)])< EPS*nA)
        {
           
                   x[k] = A[k*n+k];//s;//новое собственное значение
                        //PrintMatrix(A, n, n, n);
                        //cout<<x[k]<<" !!!!"<< s<<" "<<nA<<" "<<k<<endl;
            k--;
            steps = 0;
        }
        steps++;
    }while((k > 1)&&(steps < 1000));
    if((k>1)&&(steps >= 1000))
    {
        return -1;//не сходится быстро
    }
    }
    if(n>1)
    {
        //у матрицы 2х2 легко найти собственные значения:
        double D = ((A[0*n +0] -A[1*n +1])*(A[0*n +0] -A[1*n +1])) + (4*A[0*n +1] *A[1*n +0]);
        cout<<"D = "<<D<<endl;
        if(D < -EPS)//нет решения
        {

            return n-2;
        }
        if(abs(D) < EPS)//D = 0
        {
            x[1] = (A[0*n +0] +A[1*n +1])/2.0;
            x[0] = (A[0*n +0] +A[1*n +1])/2.0;
            return 0;
        }
            x[1] = ((A[0*n +0] +A[1*n +1]) + sqrt(D))/2.0;
            x[0] = ((A[0*n +0] +A[1*n +1]) - sqrt(D))/2.0;
        return 0;
    }
    else{//n = 1
        x[0] = A[0*n +0];
        return 0;
    }
    return -1;
}

int TA(int i, int j, double* A,int n,int m, double EPS,double* Q_cos, double* Q_sin)
{
    //определим угол поворота
    double x = A[i*n + i];
    double y = A[j*n + i];
    double root = sqrt(x*x + y*y);
    double cos_phi = 0;
    double sin_phi = 0;
    if(abs(root) < EPS)
    {
        root = 0;
        cos_phi = (x > 0 ? 1. : -1.);
    }
    else{
        cos_phi = x / root;
        sin_phi =  -y / root;
    }
    if(abs(cos_phi)<EPS)
    {
        cos_phi = 0;
    }
        if(abs(sin_phi)<EPS)
    {
        sin_phi = 0;
    }

    Q_cos[i] = cos_phi;
    Q_sin[i] = sin_phi;
        //cout<<i<<" x ="<<x<<" "<<j<<"y = "<<y<<" sin = "<<sin_phi<<" cos = "<<cos_phi<<" r =  "<<root<<endl;
    double xi = 0;
    double xj = 0;
    //умножение TA

    for(int k = i; k < m; k++) //столбцы А
    {
        xi = A[i*n+k];
        xj = A[j*n + k];
        if(k==i)
        {
            A[i*n + i] = root;
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
    return 0;
}

int T2(int i, int j, double* A,int n, double EPS,double* Q_cos, double* Q_sin)
{
    double cos_phi = Q_cos[i];
    double sin_phi =  Q_sin[i];
    double xi = 0;
    double xj = 0;
    //умножение АТ

    for(int k = 0; k < i+2; k++) //строки А
    {
        xi = A[k*n+i];
        xj = A[k*n + j];

            A[k*n + i] = xi*cos_phi - xj*sin_phi;
            if(abs(A[k*n + i]) < EPS)
            {
                A[k*n + i] = 0;
            }

            A[k*n + j] = xi*sin_phi + xj*cos_phi;
            if(abs(A[k*n + j]) < EPS)
            {
                A[k*n + j] = 0;
            }
    }
    return 0;
}




