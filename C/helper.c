#include <stdio.h>

void print_vector(double * vec, size_t J) {
    printf("[");
    for (unsigned j=0; j<J; ++j) {
        printf("%4.2E", vec[j]);
        if (j == J-1) continue;
        printf(" ");
    }
    printf("]\n");
}

double max_vector(double * vec, size_t J) {
    double max = -1E99;

    for (unsigned j=0; j<J; ++j) {
        if (vec[j] > max) {
            max = vec[j];
        }
    }

    return max;
}

double min_vector(double * vec, size_t J) {
    double min = 1E99;

    for (unsigned j=0; j<J; ++j) {
        if (vec[j] < min) {
            min = vec[j];
        }
    }

    return min;
}

double sum_vector(double * vec, size_t J) {
    double sum = 0;

    for (unsigned j=0; j<J; ++j) {
        sum += vec[j];
    }

    return sum;
}

// print progress
void print_progress(unsigned n, size_t N) {
    if ((n+1) % (N/100) == 0) {
        if ((n+1) / (N/100) != 1) printf("\r");
        printf("Progress: [");
        for (unsigned i=0; i<50; ++i) {
            if (i < (n+1) / (N/50)) printf("#");
            else printf(" ");
        }
        printf("] %7d / %7zu", n+1, N);
        fflush(stdout);
    }
}
