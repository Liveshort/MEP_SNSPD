#include <stdio.h>

void print_vector(double * vec, size_t J) {
    printf("[");
    for (unsigned j=0; j<J; ++j) {
        printf("%1.2E", vec[j]);
        if (j == J-1) continue;
        printf(" ");
    }
    printf("]\n");
}
