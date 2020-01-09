#include <stdlib.h>

#include "types.h"

// function that frees SimRes struct
void free_simres(SimRes * res) {
    // first free all the matrix contents
    for (unsigned i=0; i<res->numberOfI; ++i) {
        free(res->I[i]);
    }
    for (unsigned r=0; r<res->numberOfR; ++r) {
        free(res->R[r]);
    }
    for (unsigned t=0; t<res->numberOfT; ++t) {
        for (unsigned n=0; n<res->N; ++n) {
            free(res->T[t][n]);
        }
        free(res->T[t]);
    }
    free(res->T);
    free(res->I);
    free(res->R);
    free(res->dX);
    free(res);

    return;
}

// function that frees SimData struct
void free_simdata(SimData * simData) {
    free(simData);

    return;
}
