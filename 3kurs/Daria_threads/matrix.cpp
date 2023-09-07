#include "matrix.hpp"
#include <stdio.h>
#include <math.h>

double f(int k, int n, int i, int j)
{
    switch (k)
    {
        case 1:
            return n-fmax(i, j)+1;
        case 2:
            return fmax(i, j);
        case 3:
            return fabs(i-j);
        case 4:
            return 1.0/(i+j+1);
        default:
            return -1;
    }
}

int enter_matrix(double* a, int n, int k, FILE* fin)
{
    int i, j;
    
    if (fin)
    {
        for (i = 0; i < n; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                if (fscanf(fin, "%lf", &a[i*n+j]) != 1)
                    return -1;
            }
        }
    }
    else
    {
        for (i = 0; i < n; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                a[i*n+j] = f(k, n, i, j);
            }
        }
    }
    
    return 0;
}

void print_matrix(double* a, int n, int m)
{
    int i, j;
    
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            printf("%10.3e ", a[i*n+j]);
        }
        printf("\n");
    }
}

double error_norm(double* a, double* a_inv, int n)
{
    int i, j, k;
    double temp, sum = 0.0, max = 0.0;
        
    for (i = 0; i < n; ++i)
    {
        sum = 0.0;
        
        for (j = 0; j < n; ++j)
        {
            temp = 0.0;
            
            for (k = 0; k < n; ++k)
                temp += a[i*n+k] * a_inv[k*n+j];
            
            if (i == j)
                temp -= 1.0;
            
            sum += fabs(temp);
        }
        
        if (sum > max)
            max = sum;
    }
    
    return max;
}

double matrix_norm(double* a, int n)
{
    int i, j;
    double sum = 0.0, max = 0.0;
        
    for (i = 0; i < n; ++i)
    {
        sum = 0.0;
        
        for (j = 0; j < n; ++j)
            sum += fabs(a[i*n+j]);
        
        if (sum > max)
            max = sum;
    }
    
    return max;
}
