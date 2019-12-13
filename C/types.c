#include <stdlib.h>

#include "types.h"

// function that frees SimRes struct
void free_simres(SimRes * simRes) {
    // first free all the matrix contents
    for (unsigned n=0; n<simRes->N; ++n) {
        free(simRes->T[n]);
        free(simRes->I[n]);
        free(simRes->R[n]);
    }
    free(simRes->T);
    free(simRes->I);
    free(simRes->R);
    free(simRes);

    return;
}

// function that frees SimRes struct
void free_simdata(SimData * simData) {
    free(simData);

    return;
}
