#ifndef __HELPER_H__
#define __HELPER_H__

void print_vector(double * vec, size_t J);
double max_vector(double * vec, size_t J);
double min_vector(double * vec, size_t J);
double sum_vector(double * vec, size_t J);
double subsum_vector(double * vec, size_t i, size_t j);
int cmp_vector(double * vec, size_t J, double val, double eps);
int fill_vector(double * vec, size_t J, double val);
void swap_ptr(double ** one, double ** two);
void print_progress(unsigned n, size_t N);

#endif
