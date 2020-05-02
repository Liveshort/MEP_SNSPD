#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <lapacke.h>
#include <omp.h>

#include "types.h"
#include "helper.h"
#include "thermal.h"
#include "transmission.h"

// advances the electric model one timestep using BLAS matrix algebra.
// you need the LAPACK for this, or work out the matrix logic on your own (not recommended)
int advance_time_electric_wtf_3s(double * I0_np1, double * I1_np1, double * I2_np1, double * I3_np1, double * I4_np1, double * I5_np1, double * V_c_np1, double I0_n, double I1_n, double I2_n, double I3_n, double I4_n, double I5_n, double V_c_n, double XW0, double XP0, double XW1, double XP1, double XW2, double XP2, double X12, double XM, double Y, double R_w0_n, double R_w0_np1, double R_w1_n, double R_w1_np1, double R_w2_n, double R_w2_np1, double R_L, double R_01, double R_12, double R_small, double R_s0, double R_s1, double R_s2, double R_p0, double R_p1, double R_p2, double I_b0, double I_b1, double I_b2) {
    // set up matrix and vector for Ax = b calculation
    lapack_int n, nrhs, info;
    n = 7; nrhs = 1;

    double Q_ap = R_small + R_L + XM + X12 + R_12 + R_01;
    double Q_bp = R_small + R_L + XM + X12 + R_12;
    double Q_gp = R_small + R_L + XM;

    double I_alp = I0_n + I1_n;
    double I_bet = I0_n + I1_n + I2_n + I3_n;
    double I_gam = I0_n + I1_n + I2_n + I3_n + I4_n + I5_n;
    double I_b_bet = I_b0 + I_b1;
    double I_b_gam = I_b0 + I_b1 + I_b2;

    double * A = calloc(n*n, sizeof(double));
    double * b = calloc(n*nrhs, sizeof(double));
    lapack_int * ipiv = calloc(n, sizeof(lapack_int));

    A[0] = XW0 + R_w0_np1 + R_s0 + Q_ap;
    A[1] = Q_ap;
    A[2] = Q_bp;
    A[3] = Q_bp;
    A[4] = Q_gp;
    A[5] = Q_gp;
    A[6] = -1;

    A[7] = Q_ap;
    A[8] = XP0 + R_p0 + Q_ap;
    A[9] = Q_bp;
    A[10] = Q_bp;
    A[11] = Q_gp;
    A[12] = Q_gp;
    A[13] = -1;

    A[14] = Q_bp;
    A[15] = Q_bp;
    A[16] = XW1 + R_w1_np1 + R_s1 + Q_bp;
    A[17] = Q_bp;
    A[18] = Q_gp;
    A[19] = Q_gp;
    A[20] = -1;

    A[21] = Q_bp;
    A[22] = Q_bp;
    A[23] = Q_bp;
    A[24] = XP1 + R_p1 + Q_bp;
    A[25] = Q_gp;
    A[26] = Q_gp;
    A[27] = -1;

    A[28] = Q_gp;
    A[29] = Q_gp;
    A[30] = Q_gp;
    A[31] = Q_gp;
    A[32] = XW2 + R_w2_np1 + R_s2 + Q_gp;
    A[33] = Q_gp;
    A[34] = -1;

    A[35] = Q_gp;
    A[36] = Q_gp;
    A[37] = Q_gp;
    A[38] = Q_gp;
    A[39] = Q_gp;
    A[40] = XP2 + R_p2 + Q_gp;
    A[41] = -1;

    A[42] = Y;
    A[43] = Y;
    A[44] = Y;
    A[45] = Y;
    A[46] = Y;
    A[47] = Y;
    A[48] = 1;

    b[0] = V_c_n + 2*I_b0*R_01 + 2*I_b_bet*R_12 + 2*I_b_gam*(R_L + R_small) - I0_n*(R_w0_n + R_s0 - XW0) - I_alp*R_01 - I_bet*(R_12 - X12)  - I_gam*(R_small + R_L - XM);
    b[1] = V_c_n + 2*I_b0*R_01 + 2*I_b_bet*R_12 + 2*I_b_gam*(R_L + R_small) - I1_n*(R_p0 - XP0) - I_alp*R_01 - I_bet*(R_12 - X12) - I_gam*(R_small + R_L - XM);
    b[2] = V_c_n + 2*I_b_bet*R_12 + 2*I_b_gam*(R_L + R_small) - I2_n*(R_w1_n + R_s1 - XW1) - I_bet*(R_12 - X12) - I_gam*(R_small + R_L - XM);
    b[3] = V_c_n + 2*I_b_bet*R_12 + 2*I_b_gam*(R_L + R_small) - I3_n*(R_p1 - XP1) - I_bet*(R_12 - X12) - I_gam*(R_small + R_L - XM);
    b[4] = V_c_n + 2*I_b_gam*(R_L + R_small) - I4_n*(R_w2_n + R_s2 - XW2) - I_gam*(R_small + R_L - XM);
    b[5] = V_c_n + 2*I_b_gam*(R_L + R_small) - I5_n*(R_p2 - XP2) - I_gam*(R_small + R_L - XM);
    b[6] = V_c_n + Y*(2*I_b_gam - I_gam);

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, n, ipiv, b, nrhs);

    if (info != 0) {
        printf("Error encountered in matrix calculation of electrical model (info = %d)...\nExiting with code 3.\n", info);
        exit(3);
    }

    *I0_np1 = b[0];
    *I1_np1 = b[1];
    *I2_np1 = b[2];
    *I3_np1 = b[3];
    *I4_np1 = b[4];
    *I5_np1 = b[5];
    *V_c_np1 = b[6];

    free(A);
    free(b);
    free(ipiv);

    return 0;
}

// calculates the bias currents for certain target currents
int calculate_bias_from_target_currents_wtf_3s(double * v_I_b0, double * v_I_b1, double * v_I_b2, double I_t0, double I_t1, double I_t2, double R_s0, double R_s1, double R_s2, double R_01, double R_12, double R_p0, double R_p1, double R_p2) {
    // set up matrix and vector for Ax = b calculation
    lapack_int n, nrhs, info;
    n = 3; nrhs = 1;

    double * A = calloc(n*n, sizeof(double));
    double * b = calloc(n*nrhs, sizeof(double));
    lapack_int * ipiv = calloc(n, sizeof(lapack_int));

    double R_v0 = R_p0*R_s0/(R_p0 + R_s0);
    double R_v1 = R_p1*R_s1/(R_p1 + R_s1);
    double R_v2 = R_p2*R_s2/(R_p2 + R_s2);

    double R_z0 = R_01 + ((R_12 + R_v2)*R_v1) / (R_12 + R_v1 + R_v2);
    double R_z1 = (R_12 + R_v2)*(R_01 + R_v0)/(R_12 + R_v2 + R_01 + R_v0);
    double R_z2 = R_12 + ((R_01 + R_v0)*R_v1) / (R_01 + R_v0 + R_v1);

    double I_v0 = I_t0/R_p0*(R_p0 + R_s0);
    double I_v1 = I_t1/R_p1*(R_p1 + R_s1);
    double I_v2 = I_t2/R_p2*(R_p2 + R_s2);

    A[0] = R_z0/(R_v0 + R_z0);
    A[1] = R_v0/(R_v0 + R_01) * (R_z0 - R_01)/(R_z0 + R_v0);
    A[2] = R_v2/(R_z2 + R_v2) * (R_z2 - R_12)/R_z2 * R_v1/(R_v1 + R_v0 + R_01) * R_v0/(R_v0 + R_01);

    A[3] = R_v0/(R_v0 + R_z0) * (R_12 + R_v2)/(R_12 + R_v2 + R_v1);
    A[4] = R_z1/(R_z1 + R_v1);
    A[5] = R_v2/(R_v2 + R_z2) * (R_01 + R_v0)/(R_01 + R_v0 + R_v1);

    A[6] = R_v0/(R_z0 + R_v0) * (R_z0 - R_01)/R_z0 * R_v1/(R_v1 + R_v2 + R_12) * R_v2/(R_v2 + R_12);
    A[7] = R_v2/(R_v2 + R_12) * (R_z2 - R_12)/(R_z2 + R_v2);
    A[8] = R_z2/(R_z2 + R_v2);

    b[0] = I_v0;
    b[1] = I_v1;
    b[2] = I_v2;

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, n, ipiv, b, nrhs);

    if (info != 0) {
        puts("Error encountered in matrix calculation of bias currents...\nExiting with code 4.");
        exit(4);
    }

    *v_I_b0 = b[0];
    *v_I_b1 = b[1];
    *v_I_b2 = b[2];

    free(A);
    free(b);
    free(ipiv);

    return 0;
}

int calculate_initial_currents_wtf_3s(double v_I_b0, double v_I_b1, double v_I_b2, double * v_I_t0, double * v_I_t1, double * v_I_t2, double * v_I_t3, double * v_I_t4, double * v_I_t5, double R_s0, double R_s1, double R_s2, double R_p0, double R_p1, double R_p2, double R_01, double R_12) {
    double R_v0 = R_p0*R_s0/(R_p0 + R_s0);
    double R_v1 = R_p1*R_s1/(R_p1 + R_s1);
    double R_v2 = R_p2*R_s2/(R_p2 + R_s2);
    double R_z0 = R_01 + ((R_12 + R_v2)*R_v1) / (R_12 + R_v1 + R_v2);
    double R_z1 = (R_12 + R_v2)*(R_01 + R_v0)/(R_12 + R_v2 + R_01 + R_v0);
    double R_z2 = R_12 + ((R_01 + R_v0)*R_v1) / (R_01 + R_v0 + R_v1);
    // add an edge case for when the series resistor is zero
    if (fabs(R_s0) < 1e-6) {
        *v_I_t0 = R_v2/R_z2 * (R_z2 - R_12)/R_z2 * R_v1/(R_v1 + R_v0 + R_01) * R_v0/(R_v0 + R_01)*v_I_b2 + R_v0/(R_v0 + R_01) * (R_z0 - R_01)/(R_z0 + R_v1)*v_I_b1 + R_z0/(R_z0 + R_v0)*v_I_b0;
        *v_I_t1 = 0;
    } else {
        *v_I_t0 = R_p0/(R_s0 + R_p0)*(R_v2/R_z2 * (R_z2 - R_12)/R_z2 * R_v1/(R_v1 + R_v0 + R_01) * R_v0/(R_v0 + R_01)*v_I_b2 + R_v0/(R_v0 + R_01) * (R_z0 - R_01)/(R_z0 + R_v1)*v_I_b1 + R_z0/(R_z0 + R_v0)*v_I_b0);
        *v_I_t1 = R_s0/(R_s0 + R_p0)*(R_v2/R_z2 * (R_z2 - R_12)/R_z2 * R_v1/(R_v1 + R_v0 + R_01) * R_v0/(R_v0 + R_01)*v_I_b2 + R_v0/(R_v0 + R_01) * (R_z0 - R_01)/(R_z0 + R_v1)*v_I_b1 + R_z0/(R_z0 + R_v0)*v_I_b0);
    }
    if (fabs(R_s1) < 1e-6) {
        *v_I_t2 = R_v0/(R_v0 + R_z0) * (R_12 + R_v2)/(R_12 + R_v2 + R_v1)*v_I_b0 + R_z1/(R_z1 + R_v1)*v_I_b1 + R_v2/(R_v2 + R_z2) * (R_01 + R_v0)/(R_01 + R_v0 + R_v1)*v_I_b2;
        *v_I_t3 = 0;
    } else {
        *v_I_t2 = R_p1/(R_s1 + R_p1)*(R_v0/(R_v0 + R_z0) * (R_12 + R_v2)/(R_12 + R_v2 + R_v1)*v_I_b0 + R_z1/(R_z1 + R_v1)*v_I_b1 + R_v2/(R_v2 + R_z2) * (R_01 + R_v0)/(R_01 + R_v0 + R_v1)*v_I_b2);
        *v_I_t3 = R_s1/(R_s1 + R_p1)*(R_v0/(R_v0 + R_z0) * (R_12 + R_v2)/(R_12 + R_v2 + R_v1)*v_I_b0 + R_z1/(R_z1 + R_v1)*v_I_b1 + R_v2/(R_v2 + R_z2) * (R_01 + R_v0)/(R_01 + R_v0 + R_v1)*v_I_b2);
    }
    if (fabs(R_s2) < 1e-6) {
        *v_I_t4 = R_v0/(R_z0 + R_v0) * (R_z0 - R_01)/R_z0 * R_v1/(R_v1 + R_v2 + R_12) * R_v2/(R_v2 + R_12)*v_I_b0 + R_v2/(R_v2 + R_12) * (R_z2 - R_12)/(R_z2 + R_v2)*v_I_b1 + R_z2/(R_z2 + R_v2)*v_I_b2;
        *v_I_t5 = 0;
    } else {
        *v_I_t4 = R_p2/(R_s2 + R_p2)*(R_v0/(R_z0 + R_v0) * (R_z0 - R_01)/R_z0 * R_v1/(R_v1 + R_v2 + R_12) * R_v2/(R_v2 + R_12)*v_I_b0 + R_v2/(R_v2 + R_12) * (R_z2 - R_12)/(R_z2 + R_v2)*v_I_b1 + R_z2/(R_z2 + R_v2)*v_I_b2);
        *v_I_t5 = R_s2/(R_s2 + R_p2)*(R_v0/(R_z0 + R_v0) * (R_z0 - R_01)/R_z0 * R_v1/(R_v1 + R_v2 + R_12) * R_v2/(R_v2 + R_12)*v_I_b0 + R_v2/(R_v2 + R_12) * (R_z2 - R_12)/(R_z2 + R_v2)*v_I_b1 + R_z2/(R_z2 + R_v2)*v_I_b2);
    }

    return 0;
}

// main function that runs the simulation
int run_waterfall_3s_res(SimRes * res, SimData * data, double * dX, double dt, size_t * J, size_t N, size_t NE, size_t NTL) {
    // first locally save some important parameters that we will need all the time
    double *** T = calloc(data->numberOfT, sizeof(double **));
    for (unsigned j=0; j<data->numberOfT; ++j)
        T[j] = res->T[j];

    double ** I = calloc(data->numberOfI, sizeof(double *));
    for (unsigned j=0; j<data->numberOfI; ++j)
        I[j] = res->I[j];

    double ** R = calloc(data->numberOfR, sizeof(double *));
    for (unsigned j=0; j<data->numberOfR; ++j)
        R[j] = res->R[j];

    double * V_c = res->V_c[0];
    double * Iload = res->I[6];

    double * wireWidth = calloc(data->numberOfT, sizeof(double));
    wireWidth[0] = data->wireWidth;
    wireWidth[1] = data->wireWidth_1;
    wireWidth[2] = data->wireWidth_2;

    double * wireThickness = calloc(data->numberOfT, sizeof(double));
    wireThickness[0] = data->wireThickness;
    wireThickness[1] = data->wireThickness_1;
    wireThickness[2] = data->wireThickness_2;

    // set up vectors to temporarily save the current and next T
    // this saves a lot of memory for larger simulations
    double ** T_stash = calloc(2*data->numberOfT, sizeof(double *));
    for (unsigned j=0; j<data->numberOfT; ++j) {
        T_stash[2*j] = calloc(J[j], sizeof(double));
        T_stash[2*j+1] = calloc(J[j], sizeof(double));
    }
    double ** T_prev = calloc(data->numberOfT, sizeof(double *));
    double ** T_curr = calloc(data->numberOfT, sizeof(double *));
    for (unsigned j=0; j<data->numberOfT; ++j) {
        T_prev[j] = T[j][0];
        T_curr[j] = T[j][1];
        // set up initial thermal values at t = 1, add a steady state time step at t = 0
        fill_vector(T_prev[j], J[j], data->T_sub);
        fill_vector(T_curr[j], J[j], data->T_sub);
    }

    // determine halfway point and set up a beginning hotspot at t = 1
    unsigned halfway = J[0]/2;
    unsigned initHS_segs = (unsigned) (data->initHS_l_wtf/dX[0]) + 1;
    // check if there is a nonzero number of segments
    if (initHS_segs < 2) {
        puts("Number of segments in initial hot-spot smaller than 2.\nReturning empty result with error code 2 (wrong initial hot-spot size)...");
        res->exitValue = 2;
        exit(2);
    }
    for (unsigned j=halfway - initHS_segs/2; j<halfway + initHS_segs/2; ++j) {
        T_curr[0][j] = data->initHS_T_wtf;
    }

    // calculate the correct bias current
    double v_I_b0, v_I_b1, v_I_b2, v_I_t0, v_I_t1, v_I_t2, v_I_t3, v_I_t4, v_I_t5;
    if (fabs(data->I_b0_wtf) < 1e-12 && fabs(data->I_b1_wtf) < 1e-12 && fabs(data->I_b2_wtf) < 1e-12) {
        calculate_bias_from_target_currents_wtf_3s(&v_I_b0, &v_I_b1, &v_I_b2, data->I_t0_wtf, data->I_t1_wtf, data->I_t2_wtf, data->R_s0_wtf, data->R_s1_wtf, data->R_s2_wtf, data->R_01_wtf, data->R_12_wtf, data->R_p0_wtf, data->R_p1_wtf, data->R_p2_wtf);
    } else {
        v_I_b0 = data->I_b0_wtf;
        v_I_b1 = data->I_b1_wtf;
        v_I_b2 = data->I_b2_wtf;
    }
    // print bias currents
    printf("    bias currents:    %4.2e %4.2e %4.2e\n", v_I_b0, v_I_b1, v_I_b2);
    // put the right currents in the results
    res->I_b[0] = v_I_b0;
    res->I_b[1] = v_I_b1;
    res->I_b[2] = v_I_b2;

    // determine initial conditions
    calculate_initial_currents_wtf_3s(v_I_b0, v_I_b1, v_I_b2, &v_I_t0, &v_I_t1, &v_I_t2, &v_I_t3, &v_I_t4, &v_I_t5, data->R_s0_wtf, data->R_s1_wtf, data->R_s2_wtf, data->R_p0_wtf, data->R_p1_wtf, data->R_p2_wtf, data->R_01_wtf, data->R_12_wtf);
    for (unsigned n=0; n<=data->timeskip; ++n) {
        // set up initial current through the snspd in steady state (t = 0)
        I[0][n] = v_I_t0;
        I[1][n] = v_I_t1;
        I[2][n] = v_I_t2;
        I[3][n] = v_I_t3;
        I[4][n] = v_I_t4;
        I[5][n] = v_I_t5;
        // set up initial voltage drop over the capacitor
        V_c[n] = (v_I_t4+v_I_t5)*data->R_p2_wtf*data->R_s2_wtf/(data->R_p2_wtf + data->R_s2_wtf);
    }
    // print initial currents
    printf("    initial currents: %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e\n\n", I[0][0], I[1][0], I[2][0], I[3][0], I[4][0], I[5][0]);
    // set up vectors for the critical currents to simulate nanowire impurities
    double ** I_c = calloc(data->numberOfT, sizeof(double *));
    for (unsigned j=0; j<data->numberOfT; ++j) {
        I_c[j] = calloc(J[j], sizeof(double));

        double curr_I_c;
        if (j == 0) curr_I_c = data->I_c0_wtf;
        else if (j == 1) curr_I_c = data->I_c1_wtf;
        else if (j == 2) curr_I_c = data->I_c2_wtf;

        for (unsigned k=0; k<J[j]; ++k)
            I_c[j][k] = curr_I_c + data->impurityOffset + data->impuritySpread*(((double) rand() / (double) RAND_MAX));

        I_c[j][J[j]/2] = curr_I_c;
    }

    // prepare model parameters for estimating alpha, kappa and c
    // these parameters are considered partially state and temperature dependent
    // formulae for these can be found in Yang
    double DeltaRef = 2.15*Kb*data->T_c*(1 - (data->T_ref_wtf/data->T_c)*(data->T_ref_wtf/data->T_c));
    double A = data->c_e*exp(DeltaRef/(data->T_ref_wtf*Kb));
    double gamma = A/(2.43*data->T_c);
    double B = data->alpha/(pow(data->T_ref_wtf, 3));

    printf("    Delta: %e\n    A:     %e\n    gamma: %e\n    B:     %e\n\n", DeltaRef, A, gamma, B);

    // determine surface ratio between the cross sections of the wires
    double surfaceRatio10 = wireThickness[1]*wireWidth[1]/wireThickness[0]/wireWidth[0];
    double surfaceRatio21 = wireThickness[2]*wireWidth[2]/wireThickness[0]/wireWidth[0];
    // define the resistance of a segment of wire in the normal state
    double * R_seg = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; j++)
        R_seg[j] = data->rho_norm_wtf*dX[j]/(wireWidth[j]*wireThickness[j]);
    // declare the nanowire resistance and current density
    for (unsigned n=0; n<=data->timeskip; ++n)
        for (unsigned j=0; j<data->numberOfR; ++j)
            R[j][n] = 0;
    double * currentDensity_w = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT; ++j)
        currentDensity_w[j] = 0;

    // allocate space for the state and temperature dependent variables for each time step
    double ** alpha_n = calloc(data->numberOfT, sizeof(double *));
    double ** kappa_n = calloc(data->numberOfT, sizeof(double *));
    double ** c_n = calloc(data->numberOfT, sizeof(double *));
    double ** rho_seg_n = calloc(data->numberOfT, sizeof(double *));
    double ** R_seg_n = calloc(data->numberOfT, sizeof(double *));
    for (unsigned j=0; j<data->numberOfT; ++j) {
        alpha_n[j] = calloc(J[j], sizeof(double));
        kappa_n[j] = calloc(J[j], sizeof(double));
        c_n[j] = calloc(J[j], sizeof(double));
        rho_seg_n[j] = calloc(J[j], sizeof(double));
        R_seg_n[j] = calloc(J[j], sizeof(double));
    }

    // set up two characteristic numbers for the electrical calculations
    double XW0 = (2*data->L_w0_wtf)/dt;
    double XP0 = (2*data->L_p0_wtf)/dt;
    double XW1 = (2*data->L_w1_wtf)/dt;
    double XP1 = (2*data->L_p1_wtf)/dt;
    double XW2 = (2*data->L_w2_wtf)/dt;
    double XP2 = (2*data->L_p2_wtf)/dt;
    double X12 = (2*data->L_12_wtf)/dt;
    double XM = (2*data->L_m_wtf)/dt;
    double Y = dt/(2*data->C_m_wtf);

    // set a flag to check if done
    int flagDone = 0;

    // main time loop
    for (unsigned n=data->timeskip+1; n<NE; ++n) {
        // print progress
        print_progress(n, NE);

        // advance the thermal model to the next time step after the initial step
        if (n > data->timeskip+1 && n < N) {
            if (!data->allowOpt || cmp_vector(T_prev[0], J[0], data->T_sub, data->T_sub_eps) || cmp_vector(T_prev[1], J[1], data->T_sub, data->T_sub_eps) || cmp_vector(T_prev[2], J[2], data->T_sub, data->T_sub_eps))
                #pragma omp parallel for
                for (unsigned j=0; j<data->numberOfT; ++j)
                    advance_time_thermal(T_prev[j], T_curr[j], J[j], data->T_sub, alpha_n[j], c_n[j], rho_seg_n[j], kappa_n[j], wireThickness[j], currentDensity_w[j], dt, dX[j]);
            else {
                for (unsigned j=0; j<data->numberOfT; ++j)
                    fill_vector(T_curr[j], J[j], data->T_sub);
                flagDone = 1;
            }
        }

        if (!flagDone && n < N) {
            // first update the thermal values used in the differential equation,
            //     the targets are included as the first five parameters
            update_thermal_values(alpha_n[0], kappa_n[0], c_n[0], rho_seg_n[0], R_seg_n[0], T_curr[0], J[0], A, B, gamma, data->T_c, I[0][n-1], I_c[0], data->rho_norm_wtf, data->c_p, data->T_ref_wtf, R_seg[0]);
            update_thermal_values(alpha_n[1], kappa_n[1], c_n[1], rho_seg_n[1], R_seg_n[1], T_curr[1], J[1], A, B, gamma, data->T_c, I[2][n-1], I_c[1], data->rho_norm_wtf/surfaceRatio10, data->c_p, data->T_ref_wtf, R_seg[1]);
            update_thermal_values(alpha_n[2], kappa_n[2], c_n[2], rho_seg_n[2], R_seg_n[2], T_curr[2], J[2], A, B, gamma, data->T_c, I[4][n-1], I_c[2], data->rho_norm_wtf/surfaceRatio21, data->c_p, data->T_ref_wtf, R_seg[2]);
            // update the current nanowire resistance
            for (unsigned j=0; j<data->numberOfR; ++j)
                R[j][n] = sum_vector(R_seg_n[j], J[j]);
        } else {
            for (unsigned j=0; j<data->numberOfR; ++j)
                R[j][n] = 0;
        }

        // update the current density through the nanowire
        for (unsigned j=0; j<data->numberOfT; ++j)
            currentDensity_w[j] = I[2*j][n-1]/(wireWidth[j]*wireThickness[j]);
        // update the electric values
        advance_time_electric_wtf_3s(&I[0][n], &I[1][n], &I[2][n], &I[3][n], &I[4][n], &I[5][n], &V_c[n], I[0][n-1], I[1][n-1], I[2][n-1], I[3][n-1], I[4][n-1], I[5][n-1], V_c[n-1], XW0, XP0, XW1, XP1, XW2, XP2, X12, XM, Y, R[0][n-1], R[0][n], R[1][n-1], R[1][n], R[2][n-1], R[2][n], data->R_L_wtf, data->R_01_wtf, data->R_12_wtf, data->R_small_wtf, data->R_s0_wtf, data->R_s1_wtf, data->R_s2_wtf, data->R_p0_wtf, data->R_p1_wtf, data->R_p2_wtf, v_I_b0, v_I_b1, v_I_b2);
        Iload[n] = v_I_b0 + v_I_b1 + v_I_b2 - I[0][n] - I[1][n] - I[2][n] - I[3][n] - I[4][n] - I[5][n];

        // shuffle the T pointers around so the old and new timestep don't point to the same array
        for (unsigned j=0; j<data->numberOfT; ++j)
            T_prev[j] = T_curr[j];
        if (n % data->timeskip == 0 && n < N) {
            for (unsigned j=0; j<data->numberOfT; ++j)
                T_curr[j] = T[j][n/data->timeskip];
        } else {
            for (unsigned j=0; j<data->numberOfT; ++j) {
                if (T_prev[j] == T_stash[2*j])
                    T_curr[j] = T_stash[2*j+1];
                else
                    T_curr[j] = T_stash[2*j];
            }
        }
    }

    // free allocated space
    for (unsigned j=0; j<2*data->numberOfT; ++j)
        free(T_stash[j]);
    free(T_stash);

    free(T_prev);
    free(T_curr);

    for (unsigned j=0; j<data->numberOfT; ++j) {
        free(alpha_n[j]);
        free(kappa_n[j]);
        free(c_n[j]);
        free(rho_seg_n[j]);
        free(R_seg_n[j]);
        free(I_c[j]);
    }
    free(alpha_n);
    free(kappa_n);
    free(c_n);
    free(rho_seg_n);
    free(R_seg_n);
    free(I_c);

    free(T);
    free(I);
    free(R);
    free(wireWidth);
    free(wireThickness);

    free(R_seg);
    free(currentDensity_w);

    // transmission line loop
    if (data->simTL)
        sim_transmission_line(data, res, NE, NTL, 3);

    // print result
    puts("\n    Simulation completed.");
    res->exitValue = 0;
    return 0;
}
