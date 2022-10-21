#include "ReadMatrix.h"

int ReadMatrix(double*A, int N, int K, string FileName)
{
    if(K == 0)
    {
        Read_from_file(A,N, FileName);
    }
    else{
        Read_by_func(A,N,K);
    }
    return 0;
}

void Read_by_func(double* A, int N, int K)
{
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < N; j++)
        {
            A[i+(N*j)] = f(K,N, i, j);
        }
    }
}

double f(int k, int n, int i, int j)
{
    switch(k)
    {
        case(1):
            return n- max(i,j)+1;
            break;
        case(2):
            return max(i,j);
            break;
        case(3):
            return abs(i-j);
            break;
        case(4):
            return 1/(i+j-1);
            break;
        default:    //error
            return -1;
            break;
        
    };
}