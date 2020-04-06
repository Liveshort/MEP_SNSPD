#include <stdlib.h>

#include "types.h"

// some global constants
const double Kb = 1.3806503E-23;      // boltzmann constant
const double Lorentz = 2.45E-8;       // Lorentz number

// function that frees SimRes struct
void free_simres(SimRes * res) {
    // first free all the matrix contents
    for (unsigned i=0; i<res->numberOfI; ++i)
        free(res->I[i]);
    for (unsigned r=0; r<res->numberOfR; ++r)
        free(res->R[r]);
    for (unsigned v=0; v<res->numberOfC; ++v)
        free(res->V_c[v]);
    for (unsigned t=0; t<res->numberOfT; ++t) {
        for (unsigned n=0; n<res->N/res->timeskip; ++n)
            free(res->T[t][n]);
        free(res->T[t]);
    }
    free(res->J);
    free(res->T);
    free(res->I);
    free(res->R);
    free(res->V_c);
    free(res->dX);
    free(res->I_b);
    free(res);

    return;
}

// function that frees SimData struct
void free_simdata(SimData * data) {
    free(data);

    return;
}
