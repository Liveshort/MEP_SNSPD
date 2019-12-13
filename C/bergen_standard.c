#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "types.h"
#include "helper.h"
#include "thermal.h"

int run_bergen_standard(SimRes * res, SimData * data, double dX, double dt) {
    // first locally save some important parameters that we will need all the time
    size_t J = data->J;
    size_t N = data->N;
    double ** T = res->T;
    double ** I = res->I;
    double ** R = res->R;

    // set up initial values at t = 0
    for (unsigned j=0; j<J; ++j) {
        T[0][j] = data->T_sub;
    }
    // determine halfway point and set up a beginning hotspot
    unsigned halfway = J/2;
    unsigned initHS_segs = (unsigned) (data->initHS_l_std/dX) + 1;
    // check if there is a nonzero number of segments
    if (initHS_segs < 2) {
        puts("Number of segments in initial hot-spot smaller than 2.\nReturning empty result with error code 2 (wrong initial hot-spot size)...");
        res->exitValue = 2;
        return 2;
    }
    for (unsigned j=halfway - initHS_segs/2; j<halfway + initHS_segs/2; ++j) {
        T[0][j] = data->initHS_T_std;
    }

    print_vector(T[0], J);
}
