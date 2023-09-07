#include "solution.hpp"
#include "input.hpp"
#include <iostream>
#include <math.h>
#include <vector>
#include <future>

const double eps = 1e-16;

using namespace std;

void UToA(double * A, double * x, int i, int tem, int n, int start, int finish) {
    
    double scalyar = 0;

    for (int ii = start; ii < finish; ii++) {
        for (int j = 0; j < tem; j++) {
            scalyar += A[(i + j) * n + ii + i] * x[j];
        }
        for (int j = 0; j < tem; j++) {
            A[(i + j) * n + ii + i] -= 2 * x[j] * scalyar;
        }
        scalyar = 0;
    }
}

double FindNorm(double * X, int n) {
    double sum_X = 0;
    for (int i = 0; i < n; i++)
        sum_X += X[i] * X[i];
    return sqrt(sum_X);
}

double MatrixNorm(double * A, int n) {
    double sum = 0, max = 0;
    for (int i = 0; i < n; i++) {
        sum = 0.0;
        for (int j = 0; j < n; j++)
            sum += abs(A[i * n + j]);
        if (sum > max)
            max = sum;
    }
    return max;
}

void FindReflectionMatrix(int n, double * A_i, double * x, double * TMP) {
    ///find vector x
    double norm_TMP = FindNorm(A_i, n);
    
    for (int i = 0; i < n; i++) {
        TMP[i] = A_i[i];
    }
    TMP[0] -= norm_TMP;
    norm_TMP = FindNorm(TMP, n);
    for (int i = 0; i < n; i++) {
        if (abs(norm_TMP) > eps) {
            x[i] = TMP[i] / norm_TMP;
        }
    }
}

void SolveMatrix(double * A, double * B, double * X, int n, int thread_num, int thread_count, int k,
                 double * x, double * A_i, double * TMP, int * error) {
    
    double scalyar = 0;
    int tr;

    ///if the matrix A is bad, we need to multiply matrix on norm
    double normA = 1;
    
    if (thread_num == 0) {
        normA = MatrixNorm(A, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                A[i * n + j] /= normA;
            B[i] /= normA;
        }
    }
    
    for (int i = 0, tem = n; i < n - 1; i++, tem--) {
        
        synchronize(thread_count);
        
        if (thread_num == 0) {
            for (int k = 0; k < tem; k++) {
                A_i[k] = A[(k + i) * n + i];
            }
            FindReflectionMatrix(tem, A_i, x, TMP);
        }
        
        ///multiply matrix U(x) and A
        synchronize(thread_count);
        
        tr = tem / thread_count;
        int start = thread_num * tr;
        int finish;
        
        if (thread_num == thread_count - 1){
            finish = tem;
        } else {
            finish = thread_num * tr + tr;
        }
        
        UToA(A, x, i, tem, n, start, finish);
        
        synchronize(thread_count);
        
        ///multiply matrix U(x) and B
        if (thread_num == 0) {
            for (int j = 0; j < tem; j++) {
                scalyar += B[i + j] * x[j];
            }
            for (int j = 0; j < tem; j++) {
                B[i + j] -= 2 * x[j] * scalyar;
            }
            scalyar = 0;
        }
    }
    
    synchronize(thread_count);
    
    ///reverse of the Gauss algorithm
    if (thread_num == 0) {
        double sum = 0;
        for (int i = n - 1; i >= 0; i--) {
            sum = B[i];
            for (int j = i + 1; j < n; j++)
                sum -= X[j] * A[i * n + j];
            if (abs((A[i * n + i])) < eps) {
                * error = 1;
            }
            else
                X[i] = sum / A[i * n + i];
        }
    }
    synchronize(thread_count);
}

void synchronize(int total_threads) {
    
    static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
    static int threads_in = 0;
    static int threads_out = 0;
    
    pthread_mutex_lock(&mutex);
    
    threads_in++;
    if (threads_in >= total_threads) {
        threads_out = 0;
        pthread_cond_broadcast(&condvar_in);
    } else
        while (threads_in < total_threads)
            pthread_cond_wait(&condvar_in,&mutex);
    
    threads_out++;
    if (threads_out >= total_threads) {
        threads_in = 0;
        pthread_cond_broadcast(&condvar_out);
    } else
        while (threads_out < total_threads)
            pthread_cond_wait(&condvar_out,&mutex);
    
    pthread_mutex_unlock(&mutex);
}
