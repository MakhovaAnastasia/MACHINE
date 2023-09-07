#include <iostream>
#include <fstream>

using namespace std;

double f(int k, int n, int i, int j) {
    i++;
    j++;
    switch (k) {
        case 1:
            //if (j == 1)
                //return 0;
            return n - max(i, j) + 1;
            break;
        case 2:
            return max(i, j);
            break;
        case 3:
            return abs(i - j);
            break;
        case 4:
            return 1 / (static_cast<double>(i + j) - 1);
            break;
        default:
            cerr << "argument k is wrong";
            exit(1);
    }
}

void PrintMatrix(double * A, int l, int n, int m) {
    bool first_row = true;
    cout << '\n';
    for (int i = 0; i < l && i < m; i++) {
        if (!first_row)
            cout << '\n';
        for (int j = 0; j < n && j < m; j++) {
            cout << scientific << A[i * n + j] << " ";///using scientific floating-point notation
        }
        first_row = false;
    }
    cout << '\n' << '\n';
}

int InputMatrix(double * A, int n, int k, int argc, const char * argv[]) {
    if(k == 0) {
        if (argc < 6) {
            cout << "please, enter a file" << '\n';
            exit(1);
        }
        string filename(argv[5]);
        ifstream fin(filename);
        if (!fin) {
            cout << "no such file in the directory" << '\n';
            exit(1);
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                fin >> A[i * n + j];
                if (!fin) {
                    cout << "bad file" << '\n';
                    exit(1);
                }
            }
        }
    }
    else {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                A[i * n + j] = f(k, n, i, j);
        }
    }
    return 0;
}

void CreateB(double * A, double * B, int n) {
    double sum_B = 0;
    for (int i = 0; i < n; i++){
        for (int t = 0; t < (n + 1)/2; t++){
            sum_B += A[i * n + 2 * t];
        }
        B[i] = sum_B;
        sum_B = 0;
    }
}
