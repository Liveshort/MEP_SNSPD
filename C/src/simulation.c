#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "snspd.h"
#include "thermal.h"

int main(int argc, char * argv[]) {
    SimData * simData = calloc(1, sizeof(SimData));
    simData->J = 20;
    simData->N = 10;
    simData->numberOfI = 2;
    simData->numberOfR = 2;
    simData->wireLength = 3.56E-6;

    SimRes * simRes = run_snspd_simulation(simData, 0);

    exit(0);
}
