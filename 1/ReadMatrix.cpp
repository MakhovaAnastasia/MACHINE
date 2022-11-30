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
                    A[(N*i)+j] = N- max(i+1,j+1)+1;
                }
            } 
            return 0;      
            break;

        case(2):
            for(int i = 0; i < N; i++)
            {
                for(int j = 0; j < N; j++)
                {
                    A[(N*i)+j] = max(i+1,j+1);
                }
            } 
            return 0;
            break;

        case(3):
            for(int i = 0; i < N; i++)
            {
                for(int j = 0; j < N; j++)
                {
                    A[(N*i)+j] = abs(i-j);
                }
            }
            return 0;
            break;

        case(4):
            for(int i = 0; i < N; i++)
            {
               for(int j = 0; j < N; j++)
                {
                    A[(N*i)+j] = 1/(double)(i+1+j+1-1);

                }
            }
            return 0;
            break;

        default:    //error
            return -1;
            break;

    };
    return -1;
}

int Read_from_file(double*A, int N, string FileName)
{
    string new_number{};
    int count_num = 0;
    ifstream in(FileName);
    if(in.is_open())
    {
        while((!in.eof())&&(count_num <= N*N))
        {
            in>>new_number;
            if(!isValid(new_number))
            return -13;
            A[count_num] = stoi(new_number);
            count_num++;
        }
        if (count_num < N*N)
        {
            cout<<"недостаточно данных"<<endl;
            //меньшее число элементов
            return -11;
        }
        
    }
    else{
        cout<<"ошибка при открытии файла"<<endl;
        return -10; //ошибка при открытии файла
    }
    in.close();
    if(in.fail())
    {
        cout<<"ошибка при закрытии файла"<<endl;
        return -12; //ошибка при закрытии
    }
    return 0;
}

bool isValid(string input)
{
    int points = 0;
    for(int i = 0; i < input.length();i++ )
    {
        if((input[i]!=0)||(input[i]!=1)||(input[i]!=2)||(input[i]!=3)
        ||(input[i]!=4)||(input[i]!=5)||(input[i]!=6)||(input[i]!=7)||(input[i]!=8)
        ||(input[i]!=9))
        {
            if(input[i]=='.')
            {
                if(points!=0)
                {
                    return false;
                }
                points++;
            }
            return false;
        }
    }
    return true;
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
            cout<<" "<<scientific<<M[i*l+j];
        }
        cout<<endl;
    } 
    cout<<endl;
    return 1;
}