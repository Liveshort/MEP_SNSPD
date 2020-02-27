#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <lapacke.h>

#include "types.h"
#include "helper.h"
#include "thermal.h"

// returns the critical current for a segment of a given temperature T
static inline double I_cT_2_wtf_res_par(double I_c0, double T, double T_c) {
    return (I_c0 * (1 - T/T_c*T/T_c)*(1 - T/T_c*T/T_c));
}

// updates the alpha, kappa, c and conductivity values for all segments
// takes into account state dependence etc, specific to the model used in [Yang]
int update_thermal_values_2_wtf_res_par(double * alpha_n, double * kappa_n, double * c_n, double * rho_seg_n, double * R_seg_n, double * T_n, size_t J, double A, double B, double gamma, double T_c, double I_n, double I_c0, double rho_norm, double c_p, double T_ref, double R_seg) {
    // loop over all segments at the current timestep
    for (unsigned j=0; j<J; ++j) {
        // alpha is taken to be state independent, not strictly true, but for more, see [Yang]
        alpha_n[j] = B * pow(T_n[j], 3);
        // model values for kappa, c and rho for normal state
        if (T_n[j] > T_c || I_n > I_cT_2_wtf_res_par(I_c0, T_n[j], T_c)) {
            kappa_n[j] = Lorentz*T_n[j]/rho_norm;
            c_n[j] = gamma*T_n[j] + c_p*pow(T_n[j]/T_ref, 3);
            rho_seg_n[j] = rho_norm;
            R_seg_n[j] = R_seg;
        }
        // model values for kappa c and rho for superconducting state
        else {
            kappa_n[j] = Lorentz*T_n[j]/rho_norm * T_n[j]/T_c;
            double Delta = 2.15*T_c*Kb*(1 - (T_n[j]/T_c)*(T_n[j]/T_c));
            c_n[j] = A*exp(-Delta/(Kb*T_n[j])) + c_p*pow(T_n[j]/T_ref, 3);
            rho_seg_n[j] = 0;
            R_seg_n[j] = 0;
        }
    }

    return 0;
}

// advances the electric model one timestep using BLAS matrix algebra.
// you need the LAPACK for this, or work out the matrix logic on your own (not recommended)
int advance_time_electric_2_wtf_res_par(double * I0_np1, double * I1_np1, double * I2_np1, double * I3_np1, double * V_c_np1, double I0_n, double I1_n, double I2_n, double I3_n, double V_c_n, double XW0, double XP0, double XW1, double XP1, double XM, double Y, double R_w0_n, double R_w0_np1, double R_w1_n, double R_w1_np1, double R_L, double R_01, double R_small, double R_s0, double R_s1, double R_p0, double R_p1, double I_b0, double I_b1) {
    // set up matrix and vector for Ax = b calculation
    lapack_int n, nrhs, info;
    n = 5; nrhs = 1;

    double Q = R_small + R_L + XM;

    double I_alp = I0_n + I1_n;
    double I_bet = I0_n + I1_n + I2_n + I3_n;
    double I_b = I_b0 + I_b1;

    double * A = calloc(n*n, sizeof(double));
    double * b = calloc(n*nrhs, sizeof(double));
    lapack_int * ipiv = calloc(n, sizeof(lapack_int));

    A[0] = XW0 + Q + R_01 + R_w0_np1 + R_s0;
    A[1] = Q + R_01;
    A[2] = Q;
    A[3] = Q;
    A[4] = -1;

    A[5] = Q + R_01;
    A[6] = XP0 + Q + R_01 + R_p0;
    A[7] = Q;
    A[8] = Q;
    A[9] = -1;

    A[10] = Q;
    A[11] = Q;
    A[12] = XW1 + Q + R_w1_np1 + R_s1;
    A[13] = Q;
    A[14] = -1;

    A[15] = Q;
    A[16] = Q;
    A[17] = Q;
    A[18] = XP1 + Q + R_p1;
    A[19] = -1;

    A[20] = Y;
    A[21] = Y;
    A[22] = Y;
    A[23] = Y;
    A[24] = 1;

    b[0] = V_c_n + 2*I_b0*R_01 + 2*I_b*(R_L + R_small) - I0_n*(R_w0_n + R_s0 - XW0) - I_alp*R_01 - I_bet*(R_small + R_L - XM);
    b[1] = V_c_n + 2*I_b0*R_01 + 2*I_b*(R_L + R_small) - I1_n*(R_p0 - XP0) - I_alp*R_01 - I_bet*(R_small + R_L - XM);
    b[2] = V_c_n + 2*I_b*(R_L + R_small) - I2_n*(R_w1_n + R_s1 - XW1) - I_bet*(R_small + R_L - XM);
    b[3] = V_c_n + 2*I_b*(R_L + R_small) - I3_n*(R_p1 - XP1) - I_bet*(R_small + R_L - XM);
    b[4] = V_c_n + Y*(2*I_b - I_bet);

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, n, ipiv, b, nrhs);

    if (info != 0) {
        puts("Error encountered in matrix calculation of electrical model...\nExiting with code 3.");
        exit(3);
    }

    *I0_np1 = b[0];
    *I1_np1 = b[1];
    *I2_np1 = b[2];
    *I3_np1 = b[3];
    *V_c_np1 = b[4];

    free(A);
    free(b);
    free(ipiv);

    return 0;
}

// calculates the bias currents for certain target currents
int calculate_bias_from_target_currents_wtf_res_par(double * v_I_b0, double * v_I_b1, double I_t0, double I_t1, double R_s0, double R_s1, double R_01, double R_p0, double R_p1) {
    // set up matrix and vector for Ax = b calculation
    lapack_int n, nrhs, info;
    n = 2; nrhs = 1;

    double * A = calloc(n*n, sizeof(double));
    double * b = calloc(n*nrhs, sizeof(double));
    lapack_int * ipiv = calloc(n, sizeof(lapack_int));

    double R_v0 = R_p0*R_s0/(R_p0 + R_s0);
    double R_v1 = R_p1*R_s1/(R_p1 + R_s1);
    double R_tot = R_01 + R_v0 + R_v1;

    double I_v0 = I_t0/R_p0*(R_p0 + R_s0);
    double I_v1 = I_t1/R_p1*(R_p1 + R_s1);

    A[0] = (R_01 + R_v1)/R_tot;
    A[1] = R_v1/R_tot;
    A[2] = R_v0/R_tot;
    A[3] = (R_01 + R_v0)/R_tot;

    b[0] = I_v0;
    b[1] = I_v1;

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, n, ipiv, b, nrhs);

    if (info != 0) {
        puts("Error encountered in matrix calculation of bias currents...\nExiting with code 4.");
        return 4;
    }

    *v_I_b0 = b[0];
    *v_I_b1 = b[1];

    free(A);
    free(b);
    free(ipiv);

    return 0;
}

// main function that runs the simulation
int run_two_stage_waterfall_res_par(SimRes * res, SimData * data, double dX0, double dX1, double dt, size_t J0, size_t J1, size_t N, size_t NT, size_t NE) {
    // first locally save some important parameters that we will need all the time
    double ** T0 = res->T[0];
    double ** T1 = res->T[1];
    double * I0 = res->I[0];
    double * I1 = res->I[1];
    double * I2 = res->I[2];
    double * I3 = res->I[3];
    double * R0 = res->R[0];
    double * R1 = res->R[1];
    double * V_c = res->V_c[0];

    // set up vectors to temporarily save the current and next T
    // this saves a lot of memory for larger simulations
    double * T_stash_1 = calloc(J0, sizeof(double));
    double * T_stash_2 = calloc(J0, sizeof(double));
    double * T_stash_3 = calloc(J1, sizeof(double));
    double * T_stash_4 = calloc(J1, sizeof(double));
    double * T0_prev = T0[0];
    double * T0_curr = T0[1];
    double * T1_prev = T1[0];
    double * T1_curr = T1[1];

    // set up initial thermal values at t = 1, add a steady state time step at t = 0
    fill_vector(T0_prev, J0, data->T_sub);
    fill_vector(T0_curr, J0, data->T_sub);
    fill_vector(T1_prev, J1, data->T_sub);
    fill_vector(T1_curr, J1, data->T_sub);
    // determine halfway point and set up a beginning hotspot at t = 1
    unsigned halfway = J0/2;
    unsigned initHS_segs = (unsigned) (data->initHS_l_wtf/dX0) + 1;
    // check if there is a nonzero number of segments
    if (initHS_segs < 2) {
        puts("Number of segments in initial hot-spot smaller than 2.\nReturning empty result with error code 2 (wrong initial hot-spot size)...");
        res->exitValue = 2;
        return 2;
    }
    for (unsigned j=halfway - initHS_segs/2; j<halfway + initHS_segs/2; ++j) {
        T0_curr[j] = data->initHS_T_wtf;
    }

    // calculate the correct bias current
    double v_I_b0, v_I_b1, v_I_t0, v_I_t1, v_I_t2, v_I_t3;
    if (fabs(data->I_b0_wtf) < 1e-12 && fabs(data->I_b1_wtf) < 1e-12) {
        calculate_bias_from_target_currents_wtf_res_par(&v_I_b0, &v_I_b1, data->I_t0_wtf, data->I_t1_wtf, data->R_s0_wtf, data->R_s1_wtf, data->R_01_wtf, data->R_p0_wtf, data->R_p1_wtf);
    } else {
        v_I_b0 = data->I_b0_wtf;
        v_I_b1 = data->I_b1_wtf;
    }
    printf("bias currents %4.2e %4.2e\n", v_I_b0, v_I_b1);
    // determine initial conditions
    double R_v0 = data->R_p0_wtf*data->R_s0_wtf/(data->R_p0_wtf + data->R_s0_wtf);
    double R_v1 = data->R_p1_wtf*data->R_s1_wtf/(data->R_p1_wtf + data->R_s1_wtf);
    double R_tot = data->R_01_wtf + R_v0 + R_v1;
    // add an edge case for when the series resistor is zero
    if (fabs(data->R_s0_wtf) < 1e-6) {
        v_I_t0 = R_v1/R_tot*v_I_b1 + (data->R_01_wtf + R_v1)/R_tot*v_I_b0;
        v_I_t1 = 0;
    } else {
        v_I_t0 = data->R_p0_wtf/(data->R_s0_wtf + data->R_p0_wtf)*(R_v1/R_tot*v_I_b1 + (data->R_01_wtf + R_v1)/R_tot*v_I_b0);
        v_I_t1 = data->R_s0_wtf/(data->R_s0_wtf + data->R_p0_wtf)*(R_v1/R_tot*v_I_b1 + (data->R_01_wtf + R_v1)/R_tot*v_I_b0);
    }
    if (fabs(data->R_s1_wtf) < 1e-6) {
        v_I_t2 = R_v0/R_tot*v_I_b0 + (data->R_01_wtf + R_v0)/R_tot*v_I_b1;
        v_I_t3 = 0;
    } else {
        v_I_t2 = data->R_p1_wtf/(data->R_s1_wtf + data->R_p1_wtf)*(R_v0/R_tot*v_I_b0 + (data->R_01_wtf + R_v0)/R_tot*v_I_b1);
        v_I_t3 = data->R_s1_wtf/(data->R_s1_wtf + data->R_p1_wtf)*(R_v0/R_tot*v_I_b0 + (data->R_01_wtf + R_v0)/R_tot*v_I_b1);
    }
    for (unsigned n=0; n<=data->timeskip; ++n) {
        // set up initial current through the snspd in steady state (t = 0)
        I0[n] = v_I_t0;
        I1[n] = v_I_t1;
        I2[n] = v_I_t2;
        I3[n] = v_I_t3;
        // set up initial voltage drop over the capacitor
        V_c[n] = (v_I_t2+v_I_t3)*R_v1;
    }
    // put the right currents in the results
    res->I_b[0] = v_I_b0;
    res->I_b[1] = v_I_b1;

    // prepare model parameters for estimating alpha, kappa and c
    // these parameters are considered partially state and temperature dependent
    // formulae for these can be found in Yang
    double DeltaRef = 2.15*Kb*data->T_c*(1 - (data->T_ref_wtf/data->T_c)*(data->T_ref_wtf/data->T_c));
    double A = data->c_e*exp(DeltaRef/(data->T_ref_wtf*Kb));
    double gamma = A/(2.43*data->T_c);
    double B = data->alpha/(pow(data->T_ref_wtf, 3));

    printf("DeltaRef: %e\nA:     %e\ngamma: %e\nB:     %e\n", DeltaRef, A, gamma, B);

    // determine surface ratio between the cross sections of the wires
    double surfaceRatio10 = data->wireThickness_1*data->wireWidth_1/data->wireThickness/data->wireWidth;
    // define the resistance of a segment of wire in the normal state
    double R_seg0 = data->rho_norm_wtf*dX0/(data->wireWidth*data->wireThickness);
    double R_seg1 = data->rho_norm_wtf/surfaceRatio10*dX1/(data->wireWidth*data->wireThickness);
    // declare the nanowire resistance and current density
    for (unsigned n=0; n<=data->timeskip; ++n) {
        R0[n] = 0;
        R1[n] = 0;
    }
    double currentDensity_w0 = 0;
    double currentDensity_w1 = 0;

    // allocate space for the state and temperature dependent variables for each time step
    double * alpha0_n = calloc(J0, sizeof(double));
    double * kappa0_n = calloc(J0, sizeof(double));
    double * c0_n = calloc(J0, sizeof(double));
    double * rho_seg0_n = calloc(J0, sizeof(double));
    double * R_seg0_n = calloc(J0, sizeof(double));

    double * alpha1_n = calloc(J1, sizeof(double));
    double * kappa1_n = calloc(J1, sizeof(double));
    double * c1_n = calloc(J1, sizeof(double));
    double * rho_seg1_n = calloc(J1, sizeof(double));
    double * R_seg1_n = calloc(J1, sizeof(double));

    // set up two characteristic numbers for the electrical calculations
    double XW0 = (2*data->L_w0_wtf)/dt;
    double XP0 = (2*data->L_p0_wtf)/dt;
    double XW1 = (2*data->L_w1_wtf)/dt;
    double XP1 = (2*data->L_p1_wtf)/dt;
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
            if (!data->allowOpt || cmp_vector(T0_prev, J0, data->T_sub, data->T_sub_eps) || cmp_vector(T1_prev, J1, data->T_sub, data->T_sub_eps)) {
                advance_time_thermal(T0_prev, T0_curr, J0, data->T_sub, alpha0_n, c0_n, rho_seg0_n, kappa0_n, data->wireThickness, currentDensity_w0, dt, dX0);
                advance_time_thermal(T1_prev, T1_curr, J1, data->T_sub, alpha1_n, c1_n, rho_seg1_n, kappa1_n, data->wireThickness_1, currentDensity_w1, dt, dX1);
            } else {
                fill_vector(T0_curr, J0, data->T_sub);
                fill_vector(T1_curr, J1, data->T_sub);
                flagDone = 1;
            }
        }

        if (!flagDone && n < N) {
            // first update the thermal values used in the differential equation,
            //     the targets are included as the first five parameters
            update_thermal_values_2_wtf_res_par(alpha0_n, kappa0_n, c0_n, rho_seg0_n, R_seg0_n, T0_curr, J0, A, B, gamma, data->T_c, I0[n-1], data->I_c0_wtf, data->rho_norm_wtf, data->c_p, data->T_ref_wtf, R_seg0);
            update_thermal_values_2_wtf_res_par(alpha1_n, kappa1_n, c1_n, rho_seg1_n, R_seg1_n, T1_curr, J1, A, B, gamma, data->T_c, I2[n-1], data->I_c1_wtf, data->rho_norm_wtf/surfaceRatio10, data->c_p, data->T_ref_wtf, R_seg1);
            // update the current nanowire resistance
            R0[n] = sum_vector(R_seg0_n, J0);
            R1[n] = sum_vector(R_seg1_n, J1);
        } else {
            R0[n] = 0;
            R1[n] = 0;
        }

        // update the current density through the nanowire
        currentDensity_w0 = I0[n-1]/(data->wireWidth*data->wireThickness);
        currentDensity_w1 = I2[n-1]/(data->wireWidth_1*data->wireThickness_1);
        // update the electric values
        advance_time_electric_2_wtf_res_par(&I0[n], &I1[n], &I2[n], &I3[n], &V_c[n], I0[n-1], I1[n-1], I2[n-1], I3[n-1], V_c[n-1], XW0, XP0, XW1, XP1, XM, Y, R0[n-1], R0[n], R1[n-1], R1[n], data->R_L_wtf, data->R_01_wtf, data->R_small_wtf, data->R_s0_wtf, data->R_s1_wtf, data->R_p0_wtf, data->R_p1_wtf, v_I_b0, v_I_b1);

        // shuffle the T pointers around so the old and new timestep don't point to the same array
        T0_prev = T0_curr;
        T1_prev = T1_curr;
        if (n % data->timeskip == 0 && n < N) {
            T0_curr = T0[n/data->timeskip];
            T1_curr = T1[n/data->timeskip];
        } else {
            if (T0_prev == T_stash_1)
                T0_curr = T_stash_2;
            else
                T0_curr = T_stash_1;

            if (T1_prev == T_stash_3)
                T1_curr = T_stash_4;
            else
                T1_curr = T_stash_3;
        }
    }

    // free allocated space
    free(T_stash_1);
    free(T_stash_2);
    free(T_stash_3);
    free(T_stash_4);

    free(alpha0_n);
    free(kappa0_n);
    free(c0_n);
    free(rho_seg0_n);
    free(R_seg0_n);

    free(alpha1_n);
    free(kappa1_n);
    free(c1_n);
    free(rho_seg1_n);
    free(R_seg1_n);

    // print result
    puts("\nSimulation completed.");
    res->exitValue = 0;
    return 0;
}
