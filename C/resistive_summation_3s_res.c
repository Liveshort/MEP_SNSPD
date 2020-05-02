#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <lapacke.h>

#include "types.h"
#include "helper.h"
#include "thermal.h"
#include "transmission.h"
#include "omp.h"

// advances the electric model one timestep using BLAS matrix algebra.
// you need the LAPACK for this, or work out the matrix logic on your own (not recommended)
int advance_time_electric_res_sum_3s(unsigned n, double ** I, double ** V_c, double XW0, double XP0, double XW1, double XP1, double XW2, double XP2, double X12, double XM, double Y, double ** R_w, double R_L, double R_01, double R_12, double R_small, double R_s0, double R_s1, double R_s2, double R_p0, double R_p1, double R_p2, double I_b0, double I_b1, double I_b2, double R_k0, double R_k1, double R_k2, double R_k3, double R_Lp) {
    // set up matrix and vector for Ax = b calculation
    lapack_int m, mrhs, info;
    m = 28; mrhs = 1;

    double * Q_ap = calloc(4, sizeof(double));
    double * Q_bp = calloc(4, sizeof(double));
    double * Q_gp = calloc(4, sizeof(double));

    Q_ap[0] = R_k0 + R_small + R_L + XM + X12 + R_12 + R_01;
    Q_bp[0] = R_k0 + R_small + R_L + XM + X12 + R_12;
    Q_gp[0] = R_k0 + R_small + R_L + XM;
    Q_ap[1] = R_k1 + R_small + R_L + XM + X12 + R_12 + R_01;
    Q_bp[1] = R_k1 + R_small + R_L + XM + X12 + R_12;
    Q_gp[1] = R_k1 + R_small + R_L + XM;
    Q_ap[2] = R_k2 + R_small + R_L + XM + X12 + R_12 + R_01;
    Q_bp[2] = R_k2 + R_small + R_L + XM + X12 + R_12;
    Q_gp[2] = R_k2 + R_small + R_L + XM;
    Q_ap[3] = R_k3 + R_small + R_L + XM + X12 + R_12 + R_01;
    Q_bp[3] = R_k3 + R_small + R_L + XM + X12 + R_12;
    Q_gp[3] = R_k3 + R_small + R_L + XM;

    double * I_alp = calloc(4, sizeof(double));
    double * I_bet = calloc(4, sizeof(double));
    double * I_gam = calloc(4, sizeof(double));

    for (unsigned k=0; k<4; ++k) {
        I_alp[k] = I[6*k][n-1] + I[6*k+1][n-1];
        I_bet[k] = I[6*k][n-1] + I[6*k+1][n-1] + I[6*k+2][n-1] + I[6*k+3][n-1];
        I_gam[k] = I[6*k][n-1] + I[6*k+1][n-1] + I[6*k+2][n-1] + I[6*k+3][n-1] + I[6*k+4][n-1] + I[6*k+5][n-1];
    }

    double I_sum = I_gam[0] + I_gam[1] + I_gam[2] + I_gam[3];

    double * Rend = calloc(4, sizeof(double));
    Rend[0] = R_k0 + R_small + R_L;
    Rend[1] = R_k1 + R_small + R_L;
    Rend[2] = R_k2 + R_small + R_L;
    Rend[3] = R_k3 + R_small + R_L;

    double I_b_bet = I_b0 + I_b1;
    double I_b_gam = I_b0 + I_b1 + I_b2;

    double * A = calloc(m*m, sizeof(double));
    double * b = calloc(m*mrhs, sizeof(double));
    lapack_int * ipiv = calloc(m, sizeof(lapack_int));

    for (unsigned k=0; k<4; ++k) {
        for (unsigned j=0; j<6*k; ++j)
            A[m*6*k+j] = R_L;
        A[m*6*k+6*k+0] = XW0 + R_w[3*k][n] + R_s0 + Q_ap[k];
        A[m*6*k+6*k+1] = Q_ap[k];
        A[m*6*k+6*k+2] = Q_bp[k];
        A[m*6*k+6*k+3] = Q_bp[k];
        A[m*6*k+6*k+4] = Q_gp[k];
        A[m*6*k+6*k+5] = Q_gp[k];
        for (unsigned j=6*(k+1); j<24; ++j)
            A[m*6*k+j] = R_L;
        A[m*6*k+24+k] = -1;

        for (unsigned j=0; j<6*k; ++j)
            A[m*(6*k+1)+j] = R_L;
        A[m*(6*k+1)+6*k+0] = Q_ap[k];
        A[m*(6*k+1)+6*k+1] = XP0 + R_p0 + Q_ap[k];
        A[m*(6*k+1)+6*k+2] = Q_bp[k];
        A[m*(6*k+1)+6*k+3] = Q_bp[k];
        A[m*(6*k+1)+6*k+4] = Q_gp[k];
        A[m*(6*k+1)+6*k+5] = Q_gp[k];
        for (unsigned j=6*(k+1); j<24; ++j)
            A[m*(6*k+1)+j] = R_L;
        A[m*(6*k+1)+24+k] = -1;

        for (unsigned j=0; j<6*k; ++j)
            A[m*(6*k+2)+j] = R_L;
        A[m*(6*k+2)+6*k+0] = Q_bp[k];
        A[m*(6*k+2)+6*k+1] = Q_bp[k];
        A[m*(6*k+2)+6*k+2] = XW1 + R_w[3*k+1][n] + R_s1 + Q_bp[k];
        A[m*(6*k+2)+6*k+3] = Q_bp[k];
        A[m*(6*k+2)+6*k+4] = Q_gp[k];
        A[m*(6*k+2)+6*k+5] = Q_gp[k];
        for (unsigned j=6*(k+1); j<24; ++j)
            A[m*(6*k+2)+j] = R_L;
        A[m*(6*k+2)+24+k] = -1;

        for (unsigned j=0; j<6*k; ++j)
            A[m*(6*k+3)+j] = R_L;
        A[m*(6*k+3)+6*k+0] = Q_bp[k];
        A[m*(6*k+3)+6*k+1] = Q_bp[k];
        A[m*(6*k+3)+6*k+2] = Q_bp[k];
        A[m*(6*k+3)+6*k+3] = XP1 + R_p1 + Q_bp[k];
        A[m*(6*k+3)+6*k+4] = Q_gp[k];
        A[m*(6*k+3)+6*k+5] = Q_gp[k];
        for (unsigned j=6*(k+1); j<24; ++j)
            A[m*(6*k+3)+j] = R_L;
        A[m*(6*k+3)+24+k] = -1;

        for (unsigned j=0; j<6*k; ++j)
            A[m*(6*k+4)+j] = R_L;
        A[m*(6*k+4)+6*k+0] = Q_gp[k];
        A[m*(6*k+4)+6*k+1] = Q_gp[k];
        A[m*(6*k+4)+6*k+2] = Q_gp[k];
        A[m*(6*k+4)+6*k+3] = Q_gp[k];
        A[m*(6*k+4)+6*k+4] = XW2 + R_w[3*k+2][n] + R_s2 + Q_gp[k];
        A[m*(6*k+4)+6*k+5] = Q_gp[k];
        for (unsigned j=6*(k+1); j<24; ++j)
            A[m*(6*k+4)+j] = R_L;
        A[m*(6*k+4)+24+k] = -1;

        for (unsigned j=0; j<6*k; ++j)
            A[m*(6*k+5)+j] = R_L;
        A[m*(6*k+5)+6*k+0] = Q_gp[k];
        A[m*(6*k+5)+6*k+1] = Q_gp[k];
        A[m*(6*k+5)+6*k+2] = Q_gp[k];
        A[m*(6*k+5)+6*k+3] = Q_gp[k];
        A[m*(6*k+5)+6*k+4] = Q_gp[k];
        A[m*(6*k+5)+6*k+5] = XP2 + R_p2 + Q_gp[k];
        for (unsigned j=6*(k+1); j<24; ++j)
            A[m*(6*k+5)+j] = R_L;
        A[m*(6*k+5)+24+k] = -1;

        A[(24+k)*m+6*k+0] = Y;
        A[(24+k)*m+6*k+1] = Y;
        A[(24+k)*m+6*k+2] = Y;
        A[(24+k)*m+6*k+3] = Y;
        A[(24+k)*m+6*k+4] = Y;
        A[(24+k)*m+6*k+5] = Y;
        A[(24+k)*m+24+k] = 1;
    }

    for (unsigned k=0; k<4; ++k) {
        b[6*k] = V_c[k][n-1] + 2*I_b0*R_01 + 2*I_b_bet*R_12 + 2*I_b_gam*Rend[k] + 2*3*I_b_gam*R_L - I[6*k][n-1]*(R_w[3*k][n-1] + R_s0 - XW0) - I_alp[k]*R_01 - I_bet[k]*(R_12 - X12)  - I_gam[k]*(Rend[k] - R_L - XM) - I_sum*R_L;
        b[6*k+1] = V_c[k][n-1] + 2*I_b0*R_01 + 2*I_b_bet*R_12 + 2*I_b_gam*Rend[k] + 2*3*I_b_gam*R_L - I[6*k+1][n-1]*(R_p0 - XP0) - I_alp[k]*R_01 - I_bet[k]*(R_12 - X12) - I_gam[k]*(Rend[k] - R_L - XM) - I_sum*R_L;
        b[6*k+2] = V_c[k][n-1] + 2*I_b_bet*R_12 + 2*I_b_gam*Rend[k] + 2*3*I_b_gam*R_L - I[6*k+2][n-1]*(R_w[3*k+1][n-1] + R_s1 - XW1) - I_bet[k]*(R_12 - X12) - I_gam[k]*(Rend[k] - R_L - XM) - I_sum*R_L;
        b[6*k+3] = V_c[k][n-1] + 2*I_b_bet*R_12 + 2*I_b_gam*Rend[k] + 2*3*I_b_gam*R_L - I[6*k+3][n-1]*(R_p1 - XP1) - I_bet[k]*(R_12 - X12) - I_gam[k]*(Rend[k] - R_L - XM) - I_sum*R_L;
        b[6*k+4] = V_c[k][n-1] + 2*I_b_gam*Rend[k] + 2*3*I_b_gam*R_L - I[6*k+4][n-1]*(R_w[3*k+2][n-1] + R_s2 - XW2) - I_gam[k]*(Rend[k] - R_L - XM) - I_sum*R_L;
        b[6*k+5] = V_c[k][n-1] + 2*I_b_gam*Rend[k] + 2*3*I_b_gam*R_L - I[6*k+5][n-1]*(R_p2 - XP2) - I_gam[k]*(Rend[k] - R_L - XM) - I_sum*R_L;

        b[24+k] = V_c[k][n-1] + Y*(2*I_b_gam - I_gam[k]);
    }


    //print_matrix(A, m);
    //print_vector(b, m);

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, m, mrhs, A, m, ipiv, b, mrhs);

    if (info != 0) {
        printf("Error encountered in matrix calculation of electrical model (info = %d)...\nExiting with code 3.\n", info);
        exit(3);
    }

    for (unsigned j=0; j<24; ++j)
        I[j][n] = b[j];
    for (unsigned k=0; k<4; ++k)
        V_c[k][n] = b[24+k];

    free(A);
    free(b);
    free(ipiv);

    free(Q_ap);
    free(Q_bp);
    free(Q_gp);
    free(I_alp);
    free(I_bet);
    free(I_gam);

    free(Rend);

    return 0;
}

// calculates the bias currents for certain target currents
int calculate_bias_from_target_currents_res_sum_3s(double * v_I_b0, double * v_I_b1, double * v_I_b2, double I_t0, double I_t1, double I_t2, double R_s0, double R_s1, double R_s2, double R_01, double R_12, double R_p0, double R_p1, double R_p2) {
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

int calculate_initial_currents_res_sum_3s(double v_I_b0, double v_I_b1, double v_I_b2, double * v_I_t0, double * v_I_t1, double * v_I_t2, double * v_I_t3, double * v_I_t4, double * v_I_t5, double R_s0, double R_s1, double R_s2, double R_p0, double R_p1, double R_p2, double R_01, double R_12) {
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
int run_resistive_summation_3s_res(SimRes * res, SimData * data, double * dX, double dt, size_t * J, size_t N, size_t NE, size_t NTL) {
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

    double ** V_c = calloc(data->numberOfC, sizeof(double *));
    for (unsigned j=0; j<data->numberOfC; ++j)
        V_c[j] = res->V_c[j];
    double * Iload = res->I[24];

    double * wireWidth = calloc(data->numberOfT, sizeof(double));
    for (unsigned k=0; k<4; ++k) {
        wireWidth[3*k] = data->wireWidth;
        wireWidth[3*k+1] = data->wireWidth_1;
        wireWidth[3*k+2] = data->wireWidth_2;
    }

    double * wireThickness = calloc(data->numberOfT, sizeof(double));
    for (unsigned k=0; k<4; ++k) {
        wireThickness[3*k] = data->wireThickness;
        wireThickness[3*k+1] = data->wireThickness_1;
        wireThickness[3*k+2] = data->wireThickness_2;
    }

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
        calculate_bias_from_target_currents_res_sum_3s(&v_I_b0, &v_I_b1, &v_I_b2, data->I_t0_wtf, data->I_t1_wtf, data->I_t2_wtf, data->R_s0_wtf, data->R_s1_wtf, data->R_s2_wtf, data->R_01_wtf, data->R_12_wtf, data->R_p0_wtf, data->R_p1_wtf, data->R_p2_wtf);
    } else {
        v_I_b0 = data->I_b0_wtf;
        v_I_b1 = data->I_b1_wtf;
        v_I_b2 = data->I_b2_wtf;
    }
    // print bias currents
    printf("    bias currents:    %4.2e %4.2e %4.2e\n", v_I_b0, v_I_b1, v_I_b2);
    // put the right currents in the results
    for (unsigned k=0; k<4; ++k) {
        res->I_b[3*k] = v_I_b0;
        res->I_b[3*k+1] = v_I_b1;
        res->I_b[3*k+2] = v_I_b2;
    }
    // set up vectors for the critical currents to simulate nanowire impurities
    double ** I_c = calloc(data->numberOfT, sizeof(double *));
    for (unsigned j=0; j<data->numberOfT; ++j) {
        I_c[j] = calloc(J[j], sizeof(double));

        double curr_I_c;
        if (j % 3 == 0) curr_I_c = data->I_c0_wtf;
        else if (j % 3 == 1) curr_I_c = data->I_c1_wtf;
        else if (j % 3 == 2) curr_I_c = data->I_c2_wtf;

        for (unsigned k=0; k<J[j]; ++k)
            I_c[j][k] = curr_I_c + data->impurityOffset + data->impuritySpread*(((double) rand() / (double) RAND_MAX));

        I_c[j][J[j]/2] = curr_I_c;
    }

    // determine initial conditions
    calculate_initial_currents_res_sum_3s(v_I_b0, v_I_b1, v_I_b2, &v_I_t0, &v_I_t1, &v_I_t2, &v_I_t3, &v_I_t4, &v_I_t5, data->R_s0_wtf, data->R_s1_wtf, data->R_s2_wtf, data->R_p0_wtf, data->R_p1_wtf, data->R_p2_wtf, data->R_01_wtf, data->R_12_wtf);
    for (unsigned n=0; n<=data->timeskip; ++n) {
        // set up initial current through the snspd in steady state (t = 0)
        for (unsigned k=0; k<4; ++k) {
            I[6*k][n] = v_I_t0;
            I[6*k+1][n] = v_I_t1;
            I[6*k+2][n] = v_I_t2;
            I[6*k+3][n] = v_I_t3;
            I[6*k+4][n] = v_I_t4;
            I[6*k+5][n] = v_I_t5;
            // set up initial voltage drop over the capacitor
            V_c[k][n] = (v_I_t4+v_I_t5)*data->R_p2_wtf*data->R_s2_wtf/(data->R_p2_wtf + data->R_s2_wtf);
        }
    }
    // print initial currents
    printf("    initial currents: %4.2e %4.2e %4.2e %4.2e %4.2e %4.2e\n\n", I[0][0], I[1][0], I[2][0], I[3][0], I[4][0], I[5][0]);

    // prepare model parameters for estimating alpha, kappa and c
    // these parameters are considered partially state and temperature dependent
    // formulae for these can be found in Yang
    double DeltaRef = 2.15*Kb*data->T_c*(1 - (data->T_ref_wtf/data->T_c)*(data->T_ref_wtf/data->T_c));
    double A = data->c_e*exp(DeltaRef/(data->T_ref_wtf*Kb));
    double gamma = A/(2.43*data->T_c);
    double B = data->alpha/(pow(data->T_ref_wtf, 3));

    printf("    Delta: %e\n    A:     %e\n    gamma: %e\n    B:     %e\n\n", DeltaRef, A, gamma, B);

    // determine surface ratio between the cross sections of the wires
    double * surfaceRatio = calloc(data->numberOfT, sizeof(double));
    for (unsigned j=0; j<data->numberOfT/4; j++)
        for (unsigned k=0; k<4; k++)
            surfaceRatio[3*k+j] = wireThickness[3*k+j]*wireWidth[3*k+j]/wireThickness[3*k]/wireWidth[3*k];
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

    // main time loop
    for (unsigned n=data->timeskip+1; n<NE; ++n) {
        // trigger a second pixel somewhere down the line
        if (n == (unsigned) (0.5E-9/dt)) {
            halfway = J[4]/2;
            initHS_segs = (unsigned) (data->initHS_l_wtf/dX[4]) + 1;
            // check if there is a nonzero number of segments
            if (initHS_segs < 2) {
                puts("Number of segments in second hot-spot smaller than 2.\nReturning empty result with error code 2 (wrong initial hot-spot size)...");
                res->exitValue = 2;
                exit(2);
            }
            for (unsigned j=halfway - initHS_segs/2; j<halfway + initHS_segs/2; ++j) {
                T_prev[4][j] = data->initHS_T_wtf;
            }
        }

        // print progress
        print_progress(n, NE);

        // advance the thermal model to the next time step after the initial step
        if (n > data->timeskip+1 && n < N)
            #pragma omp parallel for
            for (unsigned j=0; j<data->numberOfT; ++j)
                advance_time_thermal(T_prev[j], T_curr[j], J[j], data->T_sub, alpha_n[j], c_n[j], rho_seg_n[j], kappa_n[j], wireThickness[j], currentDensity_w[j], dt, dX[j]);

        if (n < N) {
            // first update the thermal values used in the differential equation,
            //     the targets are included as the first five parameters
            #pragma omp parallel for
            for (unsigned j=0; j<data->numberOfT; ++j) {
                update_thermal_values(alpha_n[j], kappa_n[j], c_n[j], rho_seg_n[j], R_seg_n[j], T_curr[j], J[j], A, B, gamma, data->T_c, I[2*j][n-1], I_c[j], data->rho_norm_wtf/surfaceRatio[j], data->c_p, data->T_ref_wtf, R_seg[j]);
            }
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
        advance_time_electric_res_sum_3s(n, I, V_c, XW0, XP0, XW1, XP1, XW2, XP2, X12, XM, Y, R, data->R_L_wtf, data->R_01_wtf, data->R_12_wtf, data->R_small_wtf, data->R_s0_wtf, data->R_s1_wtf, data->R_s2_wtf, data->R_p0_wtf, data->R_p1_wtf, data->R_p2_wtf, v_I_b0, v_I_b1, v_I_b2, data->R_k0_rsm, data->R_k1_rsm, data->R_k2_rsm, data->R_k3_rsm, data->R_Lp_rsm);

        double sumI = 0;
        for (unsigned j=0; j<data->numberOfI-3; ++j)
            sumI += I[j][n];
        Iload[n] = sum_vector(res->I_b, data->numberOfT) - sumI;

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
    free(surfaceRatio);

    // transmission line loop
    if (data->simTL)
        sim_transmission_line(data, res, NE, NTL, 3);

    // print result
    puts("\n    Simulation completed.");
    res->exitValue = 0;
    return 0;
}
