#pragma once

#include <stdlib.h>
#include <cmath>

//#define EPS 0.000001 //погрешность сравнения double

int Solve(int n,double* A,double* x, double EPS, double* Q_cos, double* Q_sin);
int TA(int i, int j, double* A,int n,int m, double EPS,double* Q_cos, double* Q_sin);
int T2(int i, int j, double* A,int n, double EPS,double* Q_cos, double* Q_sin);
