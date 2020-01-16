#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "types.h"
#include "snspd.h"
#include "helper.h"
#include "yang.h"
#include "yang_parallel.h"
#include "two_stage_waterfall_res.h"

// function that coordinates the overall simulation. Data comes into this function from python
//     or whatever, is processed by the library, and is then returned to the user as a result
//     struct.
// runType is used as follows:
//   - 0: the snspd standard model is simulated. it assumes nothing but the snspd and a load
//            resistor connected after a capacitor
//   - 1: the snspd standard model with a parallel high pass filter is simulated. it assumes an
//            snspd with a parallel resistor and inductor, connected to a load resistor and
//            capacitor.
SimRes * run_snspd_simulation(SimData * data, int runType) {
    printf("Runtype %d\n", runType);
    size_t J0, J1;
    // first locally save some important parameters that we will need all the time
    J0 = data->J0;
    if (data->numberOfT > 1) J1 = data->J1;
    size_t N = data->N;
    size_t NT = data->N/data->timeskip;
    size_t NE = data->N*data->ETratio;

    // create the simulation result struct and allocate sufficient memory
    //   number of time samples N
    //   number of spacial samples J
    SimRes * res = calloc(1, sizeof(SimRes));
    res->runType = data->runType;
    res->J = calloc(data->numberOfT, sizeof(size_t));
    res->J[0] = J0;
    if (data->numberOfT > 1) res->J[1] = J1;
    res->N = N;

    res->numberOfT = data->numberOfT;
    res->numberOfI = data->numberOfI;
    res->numberOfR = data->numberOfR;
    res->numberOfC = data->numberOfC;

    res->timeskip = data->timeskip;
    res->ETratio = data->ETratio;

    res->T = calloc(data->numberOfT, sizeof(double **));
    for (unsigned t=0; t<data->numberOfT; ++t) {
        res->T[t] = calloc(NT, sizeof(double *));
        for (unsigned n=0; n<NT; ++n)
            res->T[t][n] = calloc(res->J[t], sizeof(double));
    }

    res->I = calloc(data->numberOfI, sizeof(double *));
    for (unsigned i=0; i<data->numberOfI; ++i)
        res->I[i] = calloc(NE, sizeof(double));

    res->R = calloc(data->numberOfR, sizeof(double *));
    for (unsigned r=0; r<data->numberOfR; ++r)
        res->R[r] = calloc(NE, sizeof(double));

    res->V_c = calloc(data->numberOfC, sizeof(double *));
    for (unsigned v=0; v<data->numberOfC; ++v)
        res->V_c[v] = calloc(NE, sizeof(double));

    // forward bias currents
    res->I_b = calloc(res->numberOfT, sizeof(double));
    if (runType == 0 || runType == 1) res->I_b[0] = data->I_b_std;
    if (runType == 2 || runType == 3) {
        res->I_b[0] = data->I_b0_wtf;
        res->I_b[1] = data->I_b1_wtf;
    }

    // calculate delta x and delta t
    res->dX = calloc(res->numberOfT, sizeof(double));
    res->dX[0] = data->wireLength / (J0 - 1);
    if (data->numberOfT > 1) res->dX[1] = data->wireLength_1 / (J1 - 1);
    res->dt = data->tMax / (N - 1);

    switch(runType) {
        case 0:
            run_yang(res, data, res->dX[0], res->dt, J0, N, NT, NE);
            break;
        case 1:
            run_yang_parallel(res, data, res->dX[0], res->dt, J0, N, NT, NE);
            break;
        case 2:
            run_two_stage_waterfall_res(res, data, res->dX[0], res->dX[1], res->dt, J0, J1, N, NT, NE);
            break;
        default:
            printf("Unknown runtype %d...\nReturning empty result with error 1 (wrong runtype)...", runType);
            res->exitValue = 1;
            return res;
    }

    res->exitValue = 0;

    return res;
}
