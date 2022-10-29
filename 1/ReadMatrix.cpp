#include "ReadMatrix.h"

int ReadMatrix(double*A, int N, int K, string FileName)
{
    if(K == 0)
    {
        return Read_from_file(A,N, FileName);
    }
    else{
        return Read_by_func(A,N,K);
    }
    return -1;
}

int Read_by_func(double* A, int N, int K)
{
    switch(K)
    {
        case(1):
            for(int i = 0; i < N; i++)
            {
                for(int j = 0; j < N; j++)
                {
                    A[i+(N*j)] = N- max(i,j)+1;
                }
            } 
            return K;      
            break;

        case(2):
            for(int i = 0; i < N; i++)
            {
                for(int j = 0; j < N; j++)
                {
                    A[i+(N*j)] = max(i,j);
                }
            } 
            return K;
            break;

        case(3):
            for(int i = 0; i < N; i++)
            {
                for(int j = 0; j < N; j++)
                {
                    A[i+(N*j)] = abs(i-j);
                }
            }
            return K;
            break;

        case(4):
            for(int i = 0; i < N; i++)
            {
               for(int j = 0; j < N; j++)
                {
                    A[i+(N*j)] = 1/(i+j-1);
                }
            }
            return K;
            break;

        default:    //error
            return -1;
            break;

    };
    return -1;
}

int Read_from_file(double*A, int N, string FileName)
{
    double new_number;
    int count_num = 0;
    ifstream in(FileName);
    if(in.is_open())
    {
        while((!in.eof())|| (count_num >= N*N))
        {
            in>>new_number;
            A[count_num] = new_number;
            count_num++;
        }
        if (count_num < N*N)
        {
            //меньшее число элементов
            return -11;
        }
        
    }
    else{
        return -10; //ошибка при открытии файла
    }
    in.close();
    if(in.fail())
    {
        return -12; //ошибка при закрытии
    }
    return 0;
}

int PrintMatrix(double* M, int l, int n, int m)
{
    if(l>m)
    {
        l = m;
    }
    if(n>m)
    {
        n = m;
    }
    for(int i = 0; i < l; i++)
    {
        for(int j = 0; j < n; j++)
        {
            cout<<" %10.3e"<<M[i*l+n]; 
        }
        cout<<endl;
    } 
    cout<<endl;
}