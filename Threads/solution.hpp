//#pragma once

void UToA(double * A, double * x, int i, int tem, int n, int start, int finish);

double FindNorm(double * X, int n);

double MatrixNorm(double * A, int n);

void FindReflectionMatrix(int n, double * A_i, double * x);

void SolveMatrix(double * A, double * B, double * X, int n, int thread_num, int thread_count, int k,
                 double * x, double * A_i, double * TMP, int * error);

void synchronize(int total_threads);
