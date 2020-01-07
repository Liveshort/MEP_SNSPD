#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "types.h"
#include "snspd.h"
#include "helper.h"

#include "linalg.h"

int main(int argc, char * argv[]) {
    // Put in the experiment data and run the experiment
    SimData * data = calloc(1, sizeof(SimData));
    data->J = 1000;
    data->N = 1000;
    data->numberOfT = 1;
    data->numberOfI = 1;
    data->numberOfR = 1;
    data->wireLength = 1.5E-6;//1.5E-8;//1.5E-6;
    data->wireThickness = 4E-9;
    data->wireWidth = 100E-9;
    data->tMax = 1E-10;
    data->T_c = 10.5;
    data->I_c0 = 20E-6;
    data->c_p = 9800;
    data->c_e = 2400;
    data->alpha = 8E5;
    data->T_sub = 2;
    data->T_sub_eps = 0.001;
    data->R_L_std = 50;
    data->C_m_std = 100E-9;
    data->I_b_std = 16.5E-6;
    data->initHS_l_std = 15E-9;//0.25E-8;//15E-9;
    data->initHS_T_std = 8;
    data->rho_norm_std = data->wireThickness * 600;
    data->L_w_std = 808E-9;
    data->T_ref_std = 10;

    SimRes * res = run_snspd_simulation(data, 0);

    FILE * fp = fopen("../simres/T.bin", "wb");
    for (unsigned n=0; n<res->N; ++n) {
        fwrite(res->T[0][n], sizeof(double), res->J, fp);
    }
    fclose(fp);

    fp = fopen("../simres/I.bin", "wb");
    fwrite(res->I[0], sizeof(double), res->N, fp);
    fclose(fp);

    fp = fopen("../simres/R.bin", "wb");
    fwrite(res->R[0], sizeof(double), res->N, fp);
    fclose(fp);

    free_simres(res);
    free_simdata(data);

    //double * diag = calloc(5, sizeof(double));
    //double * offdiag = calloc(5, sizeof(double));
    //double * rhs = calloc(5, sizeof(double));
//
    //diag[0] = 3;
    //diag[1] = 3;
    //diag[2] = 6;
    //diag[3] = 5;
    //diag[4] = 3;
    ////diag[5] = 3;
//
    //offdiag[0] = 0.8;
    //offdiag[1] = 0.5;
    //offdiag[2] = 0.2;
    //offdiag[3] = 0.4;
    //offdiag[4] = 0.5;
    ////offdiag[5] = 1;
//
    //rhs[0] = 10;
    //rhs[1] = 4;
    //rhs[2] = 3;
    //rhs[3] = 5;
    //rhs[4] = 8;
    ////rhs[5] = 6;
//
    //double * res = calloc(5, sizeof(double));
//
    //TDM_solve(res, 5, diag, offdiag, rhs);
    //print_vector(res, 5);

    exit(0);
}
