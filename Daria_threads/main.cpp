#include <stdio.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include "matrix.hpp"
#include "invert.hpp"
#include <sys/time.h>
#include <pthread.h>

using namespace std;

long double get_time();

typedef struct
{
    double *a;
    double *a_inv;
    double *x;
    int n;
    int thread_num;
    int threads_count;
    int *continue_flag;
    int *return_flag;
    long double time;
} Args;

void *invert(void *Arg)
{
    Args *arg = (Args*)Arg;
    long double t;
    
    synchronize(arg->threads_count);
    t = get_time();

    invert(arg->a, arg->a_inv, arg->x, arg->n, arg->thread_num, arg->threads_count, arg->continue_flag, arg->return_flag);
    
    synchronize(arg->threads_count);
    arg->time = get_time() - t;

    return NULL;
}

int main(int argc, char **argv)
{
    int n, m, k, i;
    double *a;
    double *a_inv;
    double *x;
    int continue_flag = 1, return_flag = 1;
    char filename[120];
    FILE* fin = nullptr;
    long double t;
    int threads_count;
    pthread_t *threads;
    Args *args;
    int flag;
    
    if (argc < 5)
    {
        printf("Некорректный запуск программы. Правильный формат:\n./a.out n m threads k *filename (если k != 0)");
        return -1;
    }
    
    if (sscanf(argv[1], "%d", &n) != 1 || sscanf(argv[2], "%d", &m) != 1 || sscanf(argv[3], "%d", &threads_count) != 1 || sscanf(argv[4], "%d", &k) != 1)
    {
        printf("Данные запуска некорректны.\n");
        return -1;
    }
    
    if ((k == 0 && argc != 6) || (k != 0 && argc != 5))
    {
        printf("Данные запуска некорректны.\n");
        return -1;
    }
    
    if (n < 0 || m < 0 || m > n || k < 0 || k > 4 || threads_count < 1)
    {
        printf("Данные некорректны.\n");
        return -1;
    }

    if (k == 0)
    {
        if(sscanf(argv[5], "%s", filename) != 1)
        {
            printf("Данные запуска некорректны.\n");
            return -1;
        }
            
        fin = fopen(filename, "r");
        
        if (!fin)
        {
            printf("Файла не существует.\n");
            fclose(fin);
            return -2;
        }
    }
    
    try
    {
        a = new double [n*n];
        a_inv = new double [n*n];
        x = new double [n];
        args = new Args [threads_count];
        threads = new pthread_t [threads_count];
    }
    catch (bad_alloc&)
    {
        printf("Недостаточно памяти.\n");

        if (k == 0)
            fclose(fin);
        
        return -2;
    }
    
    flag = enter_matrix(a, n, k, fin);
    
    if (flag < 0)
    {
        printf("Матрица некорректна.\n");
        
        if (k == 0)
            fclose(fin);
        
        delete []a;
        delete []a_inv;
        delete []x;
        delete []args;
        delete []threads;
        
        return -2;
    }
    
    printf("\nИзначальная матрица:\n");
    print_matrix(a, n, m);
    
    for (i = 0; i < threads_count; ++i)
    {
        args[i].a = a;
        args[i].a_inv = a_inv;
        args[i].x = x;
        args[i].n = n;
        args[i].thread_num = i;
        args[i].threads_count = threads_count;
        args[i].continue_flag = &continue_flag;
        args[i].return_flag = &return_flag;
    }
    
    for (i = 0; i < threads_count; ++i)
    {
        if (pthread_create(threads+i, 0, invert, args+i))
        {
            printf("Поток не создался.\n");
            
            if (k == 0)
                fclose(fin);
            
            delete []a;
            delete []a_inv;
            delete []x;
            delete []args;
            delete []threads;
            
            return -1;
        }
    }
    
    for (i = 0; i < threads_count; ++i)
    {
        if (pthread_join(threads[i], 0))
        {
            printf("Поток не запустился");
            
            if (k == 0)
                fclose(fin);
            
            delete []a;
            delete []a_inv;
            delete []x;
            delete []args;
            delete []threads;
            
            return -1;
        }
    }
    
    if(!return_flag)
    {
        printf("Матрица вырождена.\n");
        
        if (k == 0)
            fclose(fin);
        
        delete []a;
        delete []a_inv;
        delete []x;
        delete []args;
        delete []threads;
        
        return -1;
    }
    
//    printf("%d, %LF\n", 0, args[0].time);
    t = args[0].time;
    
    for (i = 1; i < threads_count; ++i)
    {
//        printf("%d, %LF\n", i, args[i].time);
        if (t < args[i].time)
            t = args[i].time;
    }
    
    printf("\nОбратная матрица:\n");
    print_matrix(a_inv, n, m);
    printf("\nВремя: %Lf с.\n", t);
    
    if (k == 0)
        fseek(fin, 0, SEEK_SET);
    
    flag = enter_matrix(a, n, k, fin);
    
    printf("\nПогрешность: %10.3e\n", error_norm(a, a_inv, n));
    
    if (k == 0)
        fclose(fin);
    
    delete []a;
    delete []a_inv;
    delete []x;
    delete []args;
    delete []threads;

    return 0;
}

long double get_time()
{
    struct timeval t;
    gettimeofday(&t, 0);
    return t.tv_sec + t.tv_usec/1000000.0;
}
