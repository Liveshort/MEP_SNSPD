#include <stdio.h>
#include <stdlib.h>

void print_vector(double * vec, size_t J) {
    printf("[");
    for (unsigned j=0; j<J; ++j) {
        printf("%4.2E", vec[j]);
        if (j == J-1) continue;
        printf(" ");
    }
    printf("]\n");
}

void print_matrix(double * mat, size_t m) {
    printf("####### matrix (%zu x %zu):\n", m, m);
    for (unsigned j=0; j<m; ++j) {
        for (unsigned k=0; k<m; ++k) {
            printf("%7.2g ", mat[m*j+k]);
        }
        puts("");
    }
    puts("####### end matrix");
}

double max_vector(double * vec, size_t J) {
    double max = -1E99;

    for (unsigned j=0; j<J; ++j)
        if (vec[j] > max) max = vec[j];

    return max;
}

double min_vector(double * vec, size_t J) {
    double min = 1E99;

    for (unsigned j=0; j<J; ++j)
        if (vec[j] < min) min = vec[j];

    return min;
}

double sum_vector(double * vec, size_t J) {
    double sum = 0;

    for (unsigned j=0; j<J; ++j)
        sum += vec[j];

    return sum;
}

double subsum_vector(double * vec, size_t i, size_t j) {
    double subsum = 0;

    for (unsigned k=i; k<j; ++k)
        subsum += vec[k];

    return subsum;
}

int cmp_vector(double * vec, size_t J, double val, double eps) {
    for (unsigned j=0; j<J; ++j)
        if (abs(vec[j] - val) > eps) return 1;

    return 0;
}

int fill_vector(double * vec, size_t J, double val) {
    for (unsigned j=0; j<J; ++j)
        vec[j] = val;

    return 0;
}

// print progress
void print_progress(unsigned n, size_t N) {
    if ((n+1) % (N/100) == 0) {
        if ((n+1) / (N/100) != 1) printf("\r");
        printf("    Progress: [");
        for (unsigned i=0; i<50; ++i) {
            if (i < (n+1) / (N/50)) printf("#");
            else printf(" ");
        }
        printf("] %7d / %7zu", n+1, N);
        fflush(stdout);
    }
}
