#pragma once

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <float.h>

using namespace std;

int ReadMatrix(double*A, int N, int K, string FileName);
int Read_from_file(double*A, int N, string FileName);
double f(int K,int N,int i,int j);
bool isValid(string input);

int PrintMatrix(double* M, int l, int n, int m);
