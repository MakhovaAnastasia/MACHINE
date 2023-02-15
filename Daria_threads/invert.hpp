#ifndef invert_hpp
#define invert_hpp

#include <stdio.h>
#include <pthread.h>

void invert(double *a, double *a_inv, double *x, int n, int thread_num, int threads_count, int *continue_flag, int *return_flag);
void synchronize(int total_threads);

#endif /* invert_hpp */
