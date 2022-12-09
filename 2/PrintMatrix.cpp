#include "ReadMatrix.h"

int PrintMatrix(double* M, int l, int n, int m)
{
    int L = l;
    int N = n;
    if(l>m)
    {
        L = m;
    }
    if(n>m)
    {
        N = m;
    }
    for(int i = 0; i < L; i++)
    {
        for(int j = 0; j < N; j++)
        {
            printf(" %10.3e",M[i*n+j]);
            //cout<<" "<<scientific<<M[i*n+j];
        }
        cout<<endl;
    } 
    cout<<endl;
    return 1;
}
