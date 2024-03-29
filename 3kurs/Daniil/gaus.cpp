#include "gaus.hpp"
#include "read_print.hpp"
#include "for_thread.hpp"
#include <algorithm>
#define swap double a double b {double t = a; a = b; b = t;}
#define EPS 2e-16

inline static int divUp(int a, int b) { return 1 + (a - 1) / b; }

inline static int main_element(double *matrix, int n, int i)
{
  double max_el = fabs(matrix[i * n + i]);
  int max_idx = i * n + i;
  for(int k = i; k < n; k ++)
    for(int l = i; l < n; l ++)
      if (fabs(matrix[k * n + l]) > max_el) {
        max_el = fabs(matrix[k * n + l]);
        max_idx = k * n + l;
      }
  return max_idx;
}

inline static void replace_with_main_element(double *matrix, double *reverse, int n, int i, int main_i, int main_j)
{
  double rabotyga;
  double *matrix_i = matrix + i * n, *reverse_i = reverse + i * n;
  double *matrix_main_i = matrix + main_i * n, *reverse_main_i = reverse + main_i * n;
  if (main_i == i) {
    // òà æå ñòðîêà => ìåíÿåì ïåðåìåííûå, ò.å. ñòîëáöû
    for(int k = 0; k < n; k ++) {
      rabotyga = matrix[k * n + i];
      matrix[k * n + i] = matrix[k * n + main_j];
      matrix[k * n + main_j] = rabotyga;
      rabotyga = reverse[k * n + i];
      reverse[k * n + i] = reverse[k * n + main_j];
      reverse[k * n + main_j] = rabotyga;
    }
  }
  // òîò æå ñòîëáåö => ìåíÿåì ïîäñòðîêè â îáû÷íîé, è ñòðîêè â ïðèïèñàííîé
  else if (main_j == i) {
    for (int k = i; k < n; k++) {
      rabotyga = matrix_i[k];
      matrix_i[k] = matrix_main_i[k];
      matrix_main_i[k] = rabotyga;
    }
    for (int k = 0; k < n; k++) {
      rabotyga = reverse_i[k];
      reverse_i[k] = reverse_main_i[k];
      reverse_main_i[k] = rabotyga;
    }
  }
  else {
    // êîìáèíàöèÿ ïåðåñòàíîâîê
    for (int k = 0; k < n; k++) {
      rabotyga = matrix[k * n + i];
      matrix[k * n + i] = matrix[k * n + main_j];
      matrix[k * n + main_j] = rabotyga;
      rabotyga = reverse[k * n + i];
      reverse[k * n + i] = reverse[k * n + main_j];
      reverse[k * n + main_j] = rabotyga;
    }
    for (int k = i; k < n; k++) {
      rabotyga = matrix_i[k];
      matrix_i[k] = matrix_main_i[k];
      matrix_main_i[k] = rabotyga;
    }
    for (int k = 0; k < n; k++) {
      rabotyga = reverse_i[k];
      reverse_i[k] = reverse_main_i[k];
      reverse_main_i[k] = rabotyga;
    }
  }
}

int gaus_fprop(double *matrix, double *reverse, int n,
               int *idxs, int thread_idx, int total_threads)
{
  int i, j, k;
  int tmp;
  int main_el_idx, main_i, main_j;
  int fst, lst, threads_num, thread_range;
  double *matrix_i, *reverse_i, *matrix_j, *reverse_j;
  double a_ii_rev, m_ji;
  double m_nn, m_nn_rev, *m_i, *r_i;
  for(i = 0; i < n; i++) {
    if (thread_idx == 0) {
      main_el_idx = main_element(matrix, n, i);
      main_i = main_el_idx / n;
      main_j = main_el_idx % n;
      tmp = idxs[main_j];
      idxs[main_j] = idxs[i];
      idxs[i] = tmp;
      replace_with_main_element(matrix, reverse, n, i, main_i, main_j);
      if (fabs(matrix[i * n + i]) < EPS) {
        std::cout << "zero max on: " << i << " " << fabs(matrix[i * n + i]) << std::endl;
        return -1;
      }
      a_ii_rev = 1 / matrix[i * n + i];
      matrix_i = matrix + i * n;
      reverse_i = reverse + i * n;
      // äîìíîæåàåì ïîäñòðîêè â îáû÷íîé è ñòðîêè â ïðèïèñàííîé
      for (k = i; k < n; k++) {
        matrix_i[k] *= a_ii_rev;
      }
      for (k = 0; k < n; k++) {
        reverse_i[k] *= a_ii_rev;
      }
    }
    synchronize(total_threads);
    if (i < n - 1) {
      threads_num = std::min(n - i - 1, total_threads);
      thread_range = divUp(n - i - 1, threads_num);
      fst = i + 1 + thread_range * thread_idx;
      lst = std::min(n, fst + thread_range);
      for(j = fst; j < lst; j++) {
        matrix_i = matrix + i * n;
        matrix_j = matrix + j * n;
        reverse_j = reverse + j * n;
        reverse_i = reverse + i * n;
        //printf("%f\n", matrix_j[i]);
        m_ji = matrix_j[i];
        // çà÷èùàåì ïîäñòðîêè â îáû÷íîé ìàòðèöå è ñòðîêè â ïðèïèñàííîé
        for (k = i; k < n; k++) {
          matrix_j[k] -= m_ji * matrix_i[k];
        }
        for (k = 0; k < n; k++) {
          reverse_j[k] -= m_ji * reverse_i[k];
        }
      }
    }
    synchronize(total_threads);
  }
  return 0;
}

int gaus_bprop(double *matrix, double *reverse, int n,
               int thread_idx, int total_threads)
{
  int i, j, k;
  int threads_num, thread_range, fst, lst;
  for (i = n - 1; i > 0; i--) {
    threads_num = std::min(i + 1, total_threads);
    thread_range = divUp(i + 1, threads_num);
    fst = thread_range * thread_idx;
    lst = std::min(i, fst + thread_range);
    for (j = fst; j < lst; j++) {
      //double m_ij = matrix[j * n + i];
      for (k = 0; k < n; k++) {
        //matrix[j * n + k] -= matrix[i * n + k] * m_ij;
        reverse[j * n + k] -= reverse[i * n + k] * matrix[j * n + i];
      }
    }
    synchronize(total_threads);
  }
  return 0;
}

double error(double *matrix, double *reverse, int n) {
  double error = 0;
  for (int i = 0; i < n; i++) {
    double stolb_error = 0;
    for (int j = 0; j < n; j++) {
      double el_error = 0;
      for (int k = 0; k < n; k++)
        el_error += matrix[i * n + k] * reverse[k * n + j];
      if (i == j)
        el_error -= 1;
      stolb_error += fabs(el_error);
    }
    error = std::max(stolb_error, error);
  }
  return error;
}

void replace(double * pre_reverse, double * or_matrix, int n,
             int *idxs, int total_threads, int thread_idx) {
  int thread_range = divUp(n, total_threads);
  int fst = thread_range * thread_idx;
  int lst = std::min(n, fst + thread_range);
  for (int i = fst; i < lst; i++)
    for (int j = 0; j < n; j++)
      or_matrix[idxs[i] * n + idxs[j]] = pre_reverse[i * n + j];
}

void *gaus_full_algorithm(void *void_args) {
  ARGS *args = (ARGS*)void_args;
  args->thread_time = get_full_time();
  //InvMatrix(arg->n, arg->a, arg->x, arg->index, arg->my_rank, arg->total_threads);
  int gaus_fprop_err = gaus_fprop(args->matrix, args->reverse, args->n, args->idxs,
                                  args->thread_idx, args->total_threads);
  args[0].algo_error = gaus_fprop_err;

  if (gaus_fprop_err == 0) {
    gaus_bprop(args->matrix, args->reverse, args->n,
               args->thread_idx, args->total_threads);
    replace(args->reverse, args->matrix, args->n,
            args->idxs, args->total_threads, args->thread_idx);
  }
  synchronize(args->total_threads);
  args->thread_time = get_full_time() - args->thread_time;
  return nullptr;
}
//std::cout << "==================================" << std::endl;
//print_matrix(7, n, n, matrix);
//std::cout << "--------------------------------------------------------------------" << std::endl;
//print_matrix(7, n, n, reverse);
//std::cout << "--------------------------------------------------------------------" << std::endl;
//std::cout << "||||||||||||||||||||||||||||||||||" << std::endl;
