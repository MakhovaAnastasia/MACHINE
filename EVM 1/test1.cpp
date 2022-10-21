#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <float.h>

using namespace std;

int Find_local_max(const string& filename);
int compare(double a, double b);

int main()
{
    string filename;
    if(!getline(cin, filename))
    {
        cout<<"-2"<<endl;
        return 0;
    }
    cout<<Find_local_max(filename)<<endl;
    return 0;
}

int Find_local_max(const string& filename)
{
    int num_of_local_max = 0; //число максимумов
    double prev = -DBL_MAX; //предыдущее число
    double new_number = -DBL_MAX; //новое число
    bool may_be_local_max = 0; //если =1 и new_number < local_max, то мы нашли локальный максимум
    
	ifstream in(filename);
	if (in.is_open())
    {
        while (!in.eof())
        {
            in>>new_number;
            //if(in.fail())
            //{
            //    return  -4; //проблема с чтением
            //}
            int comp = compare(new_number,prev);
            if(comp == 1) // new_number > prev
            {
                may_be_local_max = 1;
            }
            else if((comp == -1)&&(may_be_local_max == 1))
            {
                num_of_local_max++;
                may_be_local_max = 0;
            }
            prev = new_number;
        }
        if(may_be_local_max == 1)//последнее число больше предыдущего
        {
            num_of_local_max++;
        }
    }
    else{
        return -1;//не открылся файл
    }
    in.close();
    //if(in.fail())
    //{
     //   return -3; //файл не закрылся
    //}
    return num_of_local_max;
}

int compare(double a, double b)
{
    if(a > b)
        return 1; //>
    if(a < b)
        return -1;  //<
    return 0;  //=
}
