#pragma once

#include <stdlib.h>
#include <cmath>

//#define EPS 0.000000001 //погрешность сравнения double

void Solve(int n,double* A,double* b,double* x, double EPS, int thread_num,int total_threads, int* err);
//int TA(int i, int j, double* A,double* b, int n,double EPS, int thread_num,int total_threads);
