#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "types.h"
#include "snspd.h"
#include "thermal.h"

int main(int argc, char * argv[]) {
    puts("We're running");

    SimData * simData = calloc(1, sizeof(SimData));
    simData->J = 10;
    simData->N = 10;
    simData->numberOfI = 1;
    simData->numberOfR = 1;
    simData->wireLength = 1.5E-6;
    simData->wireThickness = 4E-9;
    simData->wireWidth = 100E-9;
    simData->tMax = 1E-9;
    simData->T_c = 10.5;
    simData->I_c0 = 20E-6;
    simData->c_p = 9800;
    simData->c_e = 2400;
    simData->alpha = 8E5;
    simData->T_sub = 2;
    simData->R_L_std = 50;
    simData->C_m_std = 100E-9;
    simData->I_b_std = 16.5E-6;
    simData->initHS_l_std = 0.5E-6;//15E-9;
    simData->initHS_T_std = 8;
    simData->rho_norm_std = simData->wireThickness * 600;
    simData->L_w_std = 808E-9;

    SimRes * simRes = run_snspd_simulation(simData, 0);

    free_simres(simRes);
    free_simdata(simData);

    exit(0);
}
