#ifndef __HELPER_H__
#define __HELPER_H__

void print_vector(double * vec, size_t J);
double max_vector(double * vec, size_t J);
double min_vector(double * vec, size_t J);
double sum_vector(double * vec, size_t J);
int cmp_vector(double * vec, size_t J, double val, double eps);
int fill_vector(double * vec, size_t J, double val);
void print_progress(unsigned n, size_t N);

#endif
