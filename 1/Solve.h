#pragma once

#include <stdlib.h>
#include <cmath>

#define EPS 0.000000001 //погрешность сравнения double

int Solve(int n,double* A,double* b,double* x);
int TA(int i, int j, double* A,double* b, int n);
