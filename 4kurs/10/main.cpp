#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <float.h>
#include <math.h>
#include <iomanip>
#include <algorithm>

using namespace std;


//plot '2.txt' using 1:2 with linespoints, '2.txt' using 1:3 with lines

int main(void)
{
    int A[3] = {1,10,1000};
    int N[4] = {1, 2, 3, 6};

    cout<<setw(15)<<"# "<<setw(15)<<"E1 "<<setw(15)<<"E2 "<<setw(15)<<"E3 "<<setw(15)<<"E6 "<<setw(15)<<"m "<<setw(15)<<"A "<<endl;
    for(int i = 1; i<= 6; i++)//#
    {
        for(int k = 0; k < 3; k++)//Ak
        {
            cout<<setw(15)<<i;
            for(int j = 0; j < 4; j++)//Ej
            {
                cout<<setw(15)<<N[j];
            }
            cout<<setw(15)<<"m "<<setw(15)<<A[k]<<endl;
        }
    }

/*
    //проверим скалярное произведение
    double    h = (1/(double)N);
    cout<<setprecision(10)<<fixed;

    for(int i = 1; i < N; i++) //k = 1...N-1
    {
        for(int j = 1; j < N; j++)
        {
            double res = 0;
            for(int n = 0; n < N; n++ )
            {
                 res+=y_k_n(i,n,N)*y_k_n(j,n,N);
            }
            cout<<setw(15)<<res*h;
        }
        cout<<endl;

    }
*/
    return 0;
}




