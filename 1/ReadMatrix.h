#pragma once

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <float.h>

using namespace std;

int ReadMatrix(double*A, int N, int K, string FileName);
void Read_from_file(double*A, int N, string FileName);
void Read_by_func(double*A, int N, int K);
double f(int k, int n, int i, int j);