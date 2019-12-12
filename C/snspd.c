#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "snspd.h"
#include "thermal.h"

//int run_bergen_standard(SimRes * simRes, SimData * data, double dX, double dt) {
//    // first locally save some important parameters that we will need all the time
//    size_t J = data->J;
//    size_t N = data->N;
//
//    // set up initial values at t = 0
//    for (unsigned j=0; j<J; ++j) {
//        T[0][j] = data->T_sub;
//    }
//    // determine halfway point and set up a beginning hotspot
//    unsigned halfway = J/2;
//    unsigned initHS_segs = (unsigned) data->initHS_l_std/dX
//    for (unsigned j=halfway - initHS_segs); j<halfway + initHS_segs; ++j) {
//        T[0][j] = data->T_hs;
//    }
//}

// function that coordinates the overall simulation. Data comes into this function from python
//     or whatever, is processed by the library, and is then returned to the user as a result
//     struct.
// runType is used as follows:
//   - 0: the snspd standard model is simulated. it assumes nothing but the snspd and a load
//            resistor connected after a capacitor
SimRes * run_snspd_simulation(SimData * data, int runType) {
    printf("Runtype %d", runType);
    // first locally save some important parameters that we will need all the time
    size_t J = data->J;
    size_t N = data->N;

    // create the simulation result struct and allocate sufficient memory
    //   number of time samples N
    //   number of spacial samples J
    SimRes * simRes = calloc(1, sizeof(SimRes));
    simRes->J = J;
    simRes->N = N;
    simRes->T = "mooie tekst";
    double ** T = calloc(N, sizeof(double *));
    double ** I = calloc(N, sizeof(double *));
    double ** R = calloc(N, sizeof(double *));
    for (unsigned n=0; n<N; ++n) {
        T[n] = calloc(J, sizeof(double));
        I[n] = calloc(data->numberOfI, sizeof(double));
        R[n] = calloc(data->numberOfR, sizeof(double));
    }

    // calculate delta x and delta t
    double dX = data->wireLength / (J - 1);
    double dt = data->tMax / (N - 1);

    switch(runType) {
        case 0:
            printf("Some input received: %zu, %zu, %zu, %zu, %f", J, N, data->numberOfI, data->numberOfR, data->wireLength);
            break;
        default:
            printf("Unknown runtype %d...\nReturning empty result with error 1...", runType);
            simRes->exitValue = 1;
            return simRes;
    }

    simRes->exitValue = 0;

    return simRes;
}

int main(int argc, char * argv[]) {
    puts("test");

    exit(0);
}
