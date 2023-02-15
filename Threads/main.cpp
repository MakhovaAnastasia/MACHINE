#include "input.hpp"
#include "discrepancy.hpp"
#include "solution.hpp"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <chrono>
#include <pthread.h>
#include "sys/time.h"
#include "stdexcept"

using namespace std;
using namespace std::chrono;

long double get_time() {
    struct timeval took_time;
    gettimeofday(&took_time, 0);
    return took_time.tv_sec + took_time.tv_usec / 1000000.0;
}

typedef struct {
    double * A;
    double * B;
    double * X;
    int n;
    int thread_num;
    int total_threads;
    int k;
    double * x;
    double * A_i;
    double * TMP;
    int * error_catch;
    long double time_thread;
} ARGS;

void * SolveMatrix(void * arguments) {
    
    ARGS * arg = (ARGS*)arguments;
    long double took_time;
    
    synchronize(arg -> total_threads);
    took_time = get_time();

    SolveMatrix(arg -> A, arg -> B, arg -> X, arg -> n, arg -> thread_num, arg -> total_threads, arg -> k,
                arg -> x, arg -> A_i, arg -> TMP, arg -> error_catch);
    
    synchronize(arg -> total_threads);
    arg -> time_thread = get_time() - took_time;

    return NULL;
}

int main(int argc, const char * argv[]) {
    
    int n;///matrix dimension
    int m;///submatrix dimension to print
    int k;///formula number
    int thread_count;///number of cores
    long double time;
    
    int err;///for error when entering the matrix A
    
    int error = 0;///for a degenerate matrix
    
    if (argc < 5 || argc > 6) {
        cout << "no arguments" << endl;
        return 1;
    }
    else {
        n = atoi(argv[1]);
        m = atoi(argv[2]);
        k = atoi(argv[3]);
        thread_count = atoi(argv[4]);
    }
    
    if (k > 4 || k < 0 || n < 1 || m < 1 || thread_count < 1) {
        cout << "arguments are wrong" << '\n';
        return 1;
    }
    
    ARGS * args;///argument for function SolveMatrix
    
    pthread_t * threads;
    
    double * A;///matrix

    double * B;///right part

    double * X;///solution

    double * x;///vector based on the reflection matrix

    double * A_i;///the column of the matrix A that is used to build the reflection matrix

    double * A1;///clone of matrix A
    
    double * B1;///clone of matrix B

    double * TMP;
    
    try{
        A = new double [n * n];
        B = new double [n];
        X = new double [n];
        x = new double [n];
        A_i = new double [n];
        A1 = new double[n * n];
        B1 = new double[n];
        TMP = new double [n];
        threads = new pthread_t [thread_count];
        args = new ARGS [thread_count];
    }
    catch (bad_alloc&) {
        cout << "memory allocation failure" << '\n';
        return - 1;
    }
    
    err = InputMatrix(A, n, k, argc, argv);
    if (err == 1)
        return 1;
    
    ///the begining of solving the matrix
    ///initialize vector B
    CreateB(A, B, n);
    
    cout << '\n' << "Matrix A:" << '\n';
    PrintMatrix(A, n, n, m);
    
    cout << "Vector В:" << '\n';
    PrintMatrix(B, n, 1, m);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A1[i * n + j] = A[i * n + j];
        }
    }
    for (int i = 0; i < n; i++) {
        B1[i] = B[i];
    }
    
    for (int i = 0; i < thread_count; i++) {
        args[i].A = A;
        args[i].B = B;
        args[i].X = X;
        args[i].n = n;
        args[i].thread_num = i;
        args[i].total_threads = thread_count;
        args[i].k = k;
        args[i].x = x;
        args[i].A_i = A_i;
        args[i].TMP = TMP;
        args[i].error_catch = & error;
    }
    
    for (int i = 0; i < thread_count; i++) {
        if (pthread_create (threads + i, 0, SolveMatrix, args + i)) {
            cout << "the thread with number " << i << " has not been created" << '\n';
            delete [] A;
            delete [] B;
            delete [] X;
            delete [] x;
            delete [] A_i;
            delete [] A1;
            delete [] B1;
            delete []args;
            delete []threads;
            return - 1;
        }
    }
    
    for (int i = 0; i < thread_count; i++) {
        if (pthread_join(threads[i], 0)) {
            cout << "the thread with number " << i << " has not been started" << '\n';
            delete [] A;
            delete [] B;
            delete [] X;
            delete [] x;
            delete [] A_i;
            delete [] A1;
            delete [] B1;
            delete []args;
            delete []threads;
            return - 1;
        }
    }
    
    time = args[0].time_thread;
        
    for (int i = 1; i < thread_count; i++) {
        if (time < args[i].time_thread)
            time = args[i].time_thread;
    }

    if(error) {
        cout << "no solution, the matrix is degenerate" << '\n';
        delete [] A;
        delete [] B;
        delete [] X;
        delete [] x;
        delete [] A_i;
        delete [] A1;
        delete [] B1;
        delete [] TMP;
        delete []args;
        delete []threads;
        return - 1;
    }
    
    time = args[0].time_thread;
        
    for (int i = 1; i < thread_count; i++) {
        if (time < args[i].time_thread)
            time = args[i].time_thread;
    }
    
    cout << "Solution:" << '\n';
    PrintMatrix(X, n, 1, m);
    
    cout << "Discrepancy ---> " << Discrepancy(A1, B1, X, n) << '\n' << '\n';
    
    cout << "Accuracy ------> " << Accuracy(X, n) << '\n';
    
    cerr << '\n' << "Milliseconds (" << thread_count << " threads) --> " << (int)(time * 1000) << '\n' << '\n';

    delete [] X;
    delete [] A;
    delete [] B;
    delete [] x;
    delete [] A_i;
    delete [] A1;
    delete [] B1;
    delete [] TMP;
    delete []args;
    delete []threads;
    
    return 0;
}
