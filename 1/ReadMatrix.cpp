#include "ReadMatrix.h"

int ReadMatrix(double*A, int N, int K, string FileName)
{
    if(K == 0)
    {
        return Read_from_file(A,N, FileName);
    }
    if((K==1)||(K==2)||(K==3)||(K==4))
    {
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                A[(N*i)+j] = f(K,N,i,j);
            }
        }
        return 0;
    }

    return -1;
}


double f(int K,int N,int i,int j)
{
    double r =0;
    switch(K)
    {
        case(1):
            r = N- max(i+1,j+1)+1;      
            break;

        case(2):
            r = max(i+1,j+1);
            break;

        case(3):
            r = abs(i-j);
            break;

        case(4):
            r = 1/(double)(i+1+j+1-1);
            break;

        default:    //error
            return -1;
            break;

    };
    return r;
}


int Read_from_file(double*A, int N, string FileName)
{
    string new_number{};
    int count_num = 0;
    ifstream in(FileName);
    if(in.is_open())
    {
        while((!in.eof())&&(count_num < N*N))
        {
            in>>new_number;
            if(!isValid(new_number))
            {   
                return -13;
            }
            A[count_num] = stod(new_number);
            cout<<A[count_num]<<" "<<new_number;
            count_num++;
        }
        if (count_num < N*N)
        {
            //cout<<"недостаточно данных"<<endl;
            //меньшее число элементов
            return -11;
        }
    }
    else{
        //cout<<"ошибка при открытии файла"<<endl;
        return -10; //ошибка при открытии файла
    }
    in.close();
    return 0;
}

bool isValid(string input)
{
    int points = 0;
    for(int i = 0; i < (int)input.length();i++ )
    {
        if((input[i]!='0')&&(input[i]!='1')&&(input[i]!='2')&&(input[i]!='3')
        &&(input[i]!='4')&&(input[i]!='5')&&(input[i]!='6')&&(input[i]!='7')&&(input[i]!='8')
        &&(input[i]!='9')&&(input[i]!='-'))
        {
            if(input[i]=='.')
            {
                if(points!=0)
                {
                    return false;
                }
                points++;
            }
            else{
                return false;
            }
        }
    }
    return true;
}

