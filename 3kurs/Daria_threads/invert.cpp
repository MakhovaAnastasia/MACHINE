#include "invert.hpp"
#include "matrix.hpp"
#include <math.h>

void invert(double *a, double *a_inv, double *x, int n, int thread_num, int threads_count, int *continue_flag,
            int *return_flag)
{
    int i, j, k;
    int begin_col;
    int last_col;
    double norm = 1;
    
    double eps = 1e-20;
        
    if (thread_num == 0)
    {
        norm = matrix_norm(a, n);
        
        for (i = 0; i < n; ++i)
        {
            x[i] = 0;
            
            for (j = 0; j < n; ++j)
            {
                if (i == j)
                    a_inv[i*n+j] = 1;
                else
                    a_inv[i*n+j] = 0;
                
                a[i*n+j] /= norm; 
            }
        }
    }
    
    
    for (i = 0; i < n; ++i)
    {
        synchronize(threads_count);
        double norm_a1 = 0;

        if (thread_num == 0)
        {
            double s = 0.0;
            for (j = i+1; j < n; ++j)
                s += a[j*n+i] * a[j*n+i];//(12)

            norm_a1 = sqrt(a[i*n+i]*a[i*n+i] + s);//(13)

            if (norm_a1 < eps)
            {
                *return_flag = 0;
            }
            else
            {
                x[i] = a[i*n+i] - norm_a1;
                
                for (j = i+1; j < n; ++j)
                    x[j] = a[j*n+i]; //(14)
                
                double norm_x = sqrt(x[i]*x[i] + s);//(15)

                if (norm_x < eps)
                {
                    *continue_flag = 0;
                }
                else
                {
                    *continue_flag = 1;
                    
                    for (j = i; j < n; ++j)
                        x[j] /= norm_x; //(16)
                }
            }
        }
        
        synchronize(threads_count);
        
        if (!*return_flag)
            return;
        
        if(!*continue_flag)
            continue;
 
        begin_col = (n-i) * thread_num;
        begin_col = begin_col/threads_count + i;
        last_col = (n-i) * (thread_num+1);
        last_col = last_col/threads_count + i;
                    
        for (k = begin_col; k < last_col; ++k) //лемма 10-11
        {
            double sum = 0.0;
            for (j = i; j < n; ++j)
                sum += x[j] * a[j*n+k];

            sum *= 2.0;
            for (j = i; j < n; ++j)
                a[j*n+k] -= sum * x[j];
        }

        synchronize(threads_count);

        begin_col = n * thread_num/threads_count;
        last_col = n * (thread_num+1)/threads_count;
        
        for (k = begin_col; k < last_col; ++k)//лемма 10-11
        {
            double sum = 0.0;
            for (j = i; j < n; ++j)
                sum += x[j] * a_inv[j*n+k];

            sum *= 2.0;
            for (j = i; j < n; ++j)
                a_inv[j*n+k] -= sum * x[j];
        }
    }
    
    synchronize(threads_count);

    begin_col = n * thread_num/threads_count;
    last_col = n * (thread_num+1)/threads_count;
    
    for (k = begin_col; k < last_col; ++k)
    {
        for (i = n-1; i >= 0; --i)
        {
            double temp = a_inv[i*n+k];

            for (j = i+1; j < n; ++j)
                temp -= a[i*n+j] * a_inv[j*n+k];

            a_inv[i*n+k] = temp/a[i*n+i];
        }
    }//Обратный Гаусс
    
    synchronize(threads_count);

    if (thread_num == 0)
    {
        for (i = 0; i < n; ++i)
            for (j = 0; j < n; ++j)
                a_inv[i*n+j] /= norm;
    }
    
    synchronize(threads_count);
}

void synchronize(int total_threads)
{
    static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
    static int threads_in = 0;
    static int threads_out = 0;
    
    pthread_mutex_lock(&mutex);
    
    threads_in++;
    if (threads_in >= total_threads)
    {
        threads_out = 0;
        pthread_cond_broadcast(&condvar_in);
    } else
        while (threads_in < total_threads)
            pthread_cond_wait(&condvar_in,&mutex);
    
    threads_out++;
    if (threads_out >= total_threads)
    {
        threads_in = 0;
        pthread_cond_broadcast(&condvar_out);
    } else
        while (threads_out < total_threads)
            pthread_cond_wait(&condvar_out,&mutex);
    
    pthread_mutex_unlock(&mutex);
}
