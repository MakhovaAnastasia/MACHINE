#include <iostream>
#include <ctime>
#include "gaus.hpp"
#include "for_thread.hpp"
#include "read_print.hpp"
#include <chrono>

int main(int argc, char **argv)
{
  int n, k, m, total_threads;
  std::string file_name;
  double *matrix;

  if (argc < 5 || argc > 6) {
    std::cout << "wrong argument's number\n" << "usage: 'prog n m 0 filename' or 'prog n m k t'" << std::endl;
    return -1;
  }
  else {
    n = atoi(argv[1]);
    m = atoi(argv[2]);
    k = atoi(argv[3]);
    if (k == 0) {
      if (argc == 6) {
        file_name = argv[4];
        total_threads = atoi(argv[5]);
      }
      else {
        std::cout << "if k == 0, u need to give file with matrix: 'prog n m 0 filename t'" << std::endl;
        return -2;
      }
    }
    else
      total_threads = atoi(argv[4]);
    if (n <= 0 || m <= 0 || k < 0 || k > 4 || total_threads < 1) {
      std::cout << "usage: 'prog n m 0 filename' or 'prog n m k', where n > 0, m > 0, 4 >= k > 0, t >= 1" << std::endl;
      return -2;
    }
  }

  matrix = new double[n * n];
  if (matrix == NULL) {
    std::cout << "can not allocate " << n << "*" << n << " input matrix" << std::endl;
    return -3;
  }
  if (k == 0) {
    int read_err = read_matrix(file_name.c_str(), matrix, n);
    switch (read_err) {
    case (-1):
      std::cout << "can not open file " << file_name.c_str() << std::endl;
      delete[]matrix;
      return -4;
    case (-2):
      std::cout << "can not read matrix from file " << file_name.c_str() << std::endl;
      delete[]matrix;
      return -5;
    case (-3):
      std::cout << "error with input array allocation" << std::endl;
      delete[]matrix;
      return -6;
    default:
      break;
    }
  }
  else
    fill_matrix(k, n, matrix);
  double *matrix_copy = new double[n * n];
  fill_matrix(k, n, matrix_copy);
  double *reverse = new double[n * n];
  if (reverse == NULL) {
    std::cout << "can not allocate " << n << "*" << n << " output matrix" << std::endl;
    delete[] matrix, reverse, matrix_copy;
    return -3;
  }
  fill_matrix(5, n, reverse);

  int *idxs = new int[n];
  for (int i = 0; i < n; i++)
    idxs[i] = i;
  pthread_t *threads = new pthread_t[total_threads];
  ARGS *args = new ARGS[total_threads];
  if (args == nullptr || threads == nullptr) {
    std::cout << "can not allocate threads or args" << std::endl;
    if (args) delete[]args;
    if (threads) delete[]threads;
    delete[] matrix, reverse, matrix_copy, idxs;
    return -7;
  }
  //int gaus_fprop_err = gaus_full_algorithm(matrix, reverse, n, idxs);
  for (int i = 0; i < total_threads; i++)
  {
    args[i].n = n;
    args[i].matrix = matrix;
    args[i].reverse = reverse;
    args[i].idxs = idxs;
    args[i].thread_idx = i;
    args[i].total_threads = total_threads;
  }
  double total_time = get_full_time();
  for (int i = 0; i < total_threads; i++)
    if (pthread_create(threads + i, 0, gaus_full_algorithm, args + i)) {
      std::cout << "error when creating thread " << i << std::endl;
      delete[] args, threads, matrix, reverse, matrix_copy, idxs;
      return -8;
    }
  for (int i = 0; i < total_threads; i++)
    if (pthread_join(threads[i], 0))
    {
      std::cout << "error in joining threads" << std::endl;
      if (threads) delete[]threads;
      delete[] args, matrix, reverse, matrix_copy, idxs;
      return -9;
    }
  double gaus_id_error = error(matrix_copy, matrix, n);
  int gaus_fprop_err = args[0].algo_error;
  if(gaus_fprop_err == -1) {
    std::cout << "det(matrix) = 0" << std::endl;
    delete[] args, threads, matrix, reverse, matrix_copy, idxs;
    return -10;
  }

  std::cout << "nerror: " << gaus_id_error << std::endl;

  std::cout << "total time: " << get_full_time() - total_time << " s" << std::endl;
  double thread_time = 0;
  for (int i = 0; i < total_threads; i++)
    thread_time += args->thread_time / (double)CLOCKS_PER_SEC;
  std::cout << "single thread time: " << args[0].thread_time << " s" << std::endl;
//   std::cout << "matrix_copy" << std::endl;
//   print_matrix(m, m, n, matrix_copy);
//   std::cout << "reverse matrix" << std::endl;
//   print_matrix(m, m, n, reverse);
//   std::cout << "E matrix" << std::endl;
//   print_matrix(m, m, n, matrix);
  delete[] args, threads, matrix, reverse, matrix_copy, idxs;
  return 0;
}

// size = [100, 300, 500, 750, 1000, 1250, 1500, 1750, 2000]
// time = [167, 436, 1297, 4056, 9105, 17359, 29770, 46996, 70696]
/*
  27.00 2.61
  4.63 2.97
  3.38 3.13
  2.37 2.24
  1.95 1.91
  1.73 1.71
  1.59 1.58
  1.49 1.50
*/

// size = [1000, 2000, 4000]
// time = [3642, 24976, 197951]
