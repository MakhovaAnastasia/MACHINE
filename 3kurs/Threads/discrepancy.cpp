#include "discrepancy.hpp"
#include "solution.hpp"
#include "input.hpp"
#include <iostream>
#include <math.h>

double Discrepancy(double * A, double * B, double * X, int n) {
    double sum;
    for (int i = 0; i < n; i++) {
        sum = 0;
        for (int j = 0; j < n; j++) {
            sum += A[i * n + j] * X[j];
        }
        B[i] -= sum;
    }
    return FindNorm(B, n);
}

double Accuracy(double * X, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
        if (i % 2 == 0)
            sum += (X[i] - 1) * (X[i] - 1);
        else
            sum += X[i] * X[i];
    }
    return sqrt(sum);
}
