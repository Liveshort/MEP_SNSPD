#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <lapacke.h>

#include "types.h"
#include "helper.h"
#include "thermal.h"
#include "transmission.h"

int advance_time_electric_yang_parallel(double * I0_np1, double * I1_np1, double * V_c_np1, double I0_n, double I1_n, double V_c_n, double X1, double X2, double Y, double R_w_n, double R_w_np1, double R_p, double R_L, double R_s, double I_b) {
    // set up matrix and vector for Ax = b calculation
    lapack_int n, nrhs, info;
    n = 3; nrhs = 1;

    double * A = calloc(n*n, sizeof(double));
    double * b = calloc(n*nrhs, sizeof(double));
    lapack_int * ipiv = calloc(n, sizeof(lapack_int));

    A[0] = X1 + R_L + R_w_np1 + R_s;
    A[1] = R_L;
    A[2] = -1;

    A[3] = R_L;
    A[4] = X2 + R_L + R_p;
    A[5] = -1;

    A[6] = Y;
    A[7] = Y;
    A[8] = 1;

    b[0] = V_c_n + 2*R_L*I_b + (X1 - R_L - R_w_n - R_s)*I0_n - R_L*I1_n;
    b[1] = V_c_n + 2*R_L*I_b + (X2 - R_L - R_p)*I1_n - R_L*I0_n;
    b[2] = V_c_n + Y*(2*I_b - I0_n - I1_n);

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, n, ipiv, b, nrhs);

    if (info != 0) {
        puts("Error encountered in matrix calculation of electrical model...\nExiting with code 3.");
        return 3;
    }

    *I0_np1 = b[0];
    *I1_np1 = b[1];
    *V_c_np1 = b[2];

    free(A);
    free(b);
    free(ipiv);

    return 0;
}

// main function that runs the simulation
int run_yang_parallel(SimRes * res, SimData * data, double dX, double dt, size_t J, size_t N, size_t NE, size_t NTL) {
    // first locally save some important parameters that we will need all the time
    double ** T = res->T[0];
    double * I0 = res->I[0];
    double * I1 = res->I[1];
    double * R = res->R[0];
    double * V_c = res->V_c[0];
    double * Iload = res->I[2];

    // set up vectors to temporarily save the current and next T
    // this saves a lot of memory for larger simulations
    double * T_stash_1 = calloc(J, sizeof(double));
    double * T_stash_2 = calloc(J, sizeof(double));
    double * T_prev = T[0];
    double * T_curr = T[1];

    // set up initial thermal values at t = 1, add a steady state time step at t = 0
    fill_vector(T_prev, J, data->T_sub);
    fill_vector(T_curr, J, data->T_sub);
    // determine halfway point and set up a beginning hotspot at t = 1
    unsigned halfway = J/2;
    unsigned initHS_segs = (unsigned) (data->initHS_l_std/dX) + 1;
    // check if there is a nonzero number of segments
    if (initHS_segs < 2) {
        puts("Number of segments in initial hot-spot smaller than 2.\nReturning empty result with error code 2 (wrong initial hot-spot size)...");
        res->exitValue = 2;
        return 2;
    }
    for (unsigned j=halfway - initHS_segs/2; j<halfway + initHS_segs/2; ++j) {
        T_curr[j] = data->initHS_T_std;
    }

    // perform some math voodoo to work out what bias current we need to get the
    //     target current that the user provided in the nanowire if target is provided
    // if bias is provided, just use that
    double r, v_I0, v_I1;
    r = data->R_p_parallel / (data->R_s_std + data->R_p_parallel);
    if (fabs(data->I_b_std) < 1e-12) {
        v_I0 = data->I_t_std;
        v_I1 = (1-r)*data->I_t_std/r;
    } else {
        v_I0 = r*data->I_b_std;
        v_I1 = (1-r)*data->I_b_std;
    }
    for (unsigned n=0; n<=data->timeskip; ++n) {
        // set up initial current through the snspd in steady state (t = 0)
        I0[n] = v_I0;
        I1[n] = v_I1;
        // set up initial voltage drop over the capacitor
        V_c[n] = data->R_p_parallel*data->R_s_std / (data->R_s_std + data->R_p_parallel)*(v_I0 + v_I1);
    }
    // put right value in the results
    res->I_b[0] = v_I0 + v_I1;

    // prepare model parameters for estimating alpha, kappa and c
    // these parameters are considered partially state and temperature dependent
    // formulae for these can be found in Yang
    double DeltaRef = 2.15*Kb*data->T_c*(1 - (data->T_ref_std/data->T_c)*(data->T_ref_std/data->T_c));
    double A = data->c_e*exp(DeltaRef/(data->T_ref_std*Kb));
    double gamma = A/(2.43*data->T_c);
    double B = data->alpha/(pow(data->T_ref_std, 3));

    printf("    DeltaRef: %e\n    A:     %e\n    gamma: %e\n    B:     %e\n\n", DeltaRef, A, gamma, B);

    // define the resistance of a segment of wire in the normal state
    double R_seg = data->rho_norm_std*dX/(data->wireWidth*data->wireThickness);
    // declare the nanowire resistance and current density
    for (unsigned n=0; n<=data->timeskip; ++n)
        R[0] = 0;
    double currentDensity_w = 0;

    // allocate space for the state and temperature dependent variables for each time step
    double * alpha_n = calloc(J, sizeof(double));
    double * kappa_n = calloc(J, sizeof(double));
    double * c_n = calloc(J, sizeof(double));
    double * rho_seg_n = calloc(J, sizeof(double));
    double * R_seg_n = calloc(J, sizeof(double));

    // set up two characteristic numbers for the electrical calculations
    double X1 = (2*data->L_w_std)/dt;
    double X2 = (2*data->L_p_parallel)/dt;
    double Y = dt/(2*data->C_m_std);

    // set a flag to check if done
    int flagDone = 0;

    // main time loop
    for (unsigned n=data->timeskip+1; n<NE; ++n) {
        // print progress
        print_progress(n, NE);

        // advance the thermal model to the next time step after the initial step
        if (n > data->timeskip+1 && n < N) {
            if (cmp_vector(T_prev, J, data->T_sub, data->T_sub_eps) || !data->allowOpt)
                advance_time_thermal(T_prev, T_curr, J, data->T_sub, alpha_n, c_n, rho_seg_n, kappa_n, data->wireThickness, currentDensity_w, dt, dX);
            else {
                fill_vector(T_curr, J, data->T_sub);
                flagDone = 1;
            }
        }

        if (!flagDone && n < N) {
            // first update the thermal values used in the differential equation,
            //     the targets are included as the first five parameters
            update_thermal_values(alpha_n, kappa_n, c_n, rho_seg_n, R_seg_n, T_curr, J, A, B, gamma, data->T_c, I0[n-1], data->I_c0, data->rho_norm_std, data->c_p, data->T_ref_std, R_seg);
            // update the current nanowire resistance
            R[n] = sum_vector(R_seg_n, J);
        } else {
            R[n] = 0;
        }

        // update the current density through the nanowire
        currentDensity_w = I0[n-1]/(data->wireWidth*data->wireThickness);
        // update the electric values
        advance_time_electric_yang_parallel(&I0[n], &I1[n], &V_c[n], I0[n-1], I1[n-1], V_c[n-1], X1, X2, Y, R[n-1], R[n], data->R_p_parallel, data->R_L_std, data->R_s_std, v_I0 + v_I1);
        Iload[n] = v_I0 + v_I1 - I0[n] - I1[n];

        // shuffle the T pointers around so the old and new timestep don't point to the same array
        T_prev = T_curr;
        if (n % data->timeskip == 0)
            T_curr = T[n/data->timeskip];
        else {
            if (T_prev == T_stash_1)
                T_curr = T_stash_2;
            else
                T_curr = T_stash_1;
        }
    }

    // free allocated space
    free(T_stash_1);
    free(T_stash_2);
    free(alpha_n);
    free(kappa_n);
    free(c_n);
    free(rho_seg_n);
    free(R_seg_n);

    // transmission line loop
    if (data->simTL)
        sim_transmission_line(data, res, NE, NTL, 1);

    // print result
    puts("\n    Simulation completed.");
    res->exitValue = 0;
    return 0;
}
