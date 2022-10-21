
#include "ReadMatrix.h"

int main()
{
    string filename;
    int n = -1; //размер матрицы
    int m = -1; //размер вывода
    int k = -1; // формула
    cin>>n>>m>>k;

    if(k == 0)
    {
        if(!getline(cin, filename))
        {
            cout<<"-2"<<endl;
            return 0;
        }
    }
    
    ReadMatrix(n,k, filename);


    return 0;
}