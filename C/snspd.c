#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "types.h"
#include "snspd.h"
#include "helper.h"
#include "bergen_standard.h"

// function that coordinates the overall simulation. Data comes into this function from python
//     or whatever, is processed by the library, and is then returned to the user as a result
//     struct.
// runType is used as follows:
//   - 0: the snspd standard model is simulated. it assumes nothing but the snspd and a load
//            resistor connected after a capacitor
SimRes * run_snspd_simulation(SimData * data, int runType) {
    printf("Runtype %d\n", runType);
    // first locally save some important parameters that we will need all the time
    size_t J = data->J;
    size_t N = data->N;

    // create the simulation result struct and allocate sufficient memory
    //   number of time samples N
    //   number of spacial samples J
    SimRes * res = calloc(1, sizeof(SimRes));
    res->J = J;
    res->N = N;
    res->T = calloc(N, sizeof(double *));
    res->I = calloc(N, sizeof(double *));
    res->R = calloc(N, sizeof(double *));
    for (unsigned n=0; n<N; ++n) {
        res->T[n] = calloc(J, sizeof(double));
        res->I[n] = calloc(data->numberOfI, sizeof(double));
        res->R[n] = calloc(data->numberOfR, sizeof(double));
    }

    // calculate delta x and delta t
    double dX = data->wireLength / (J - 1);
    double dt = data->tMax / (N - 1);

    switch(runType) {
        case 0:
            run_bergen_standard(res, data, dX, dt);
            break;
        default:
            printf("Unknown runtype %d...\nReturning empty result with error 1 (wrong runType)...", runType);
            res->exitValue = 1;
            return res;
    }

    res->exitValue = 0;

    return res;
}
