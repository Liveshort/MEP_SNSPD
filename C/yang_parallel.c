#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <lapacke.h>

#include "types.h"
#include "helper.h"
#include "thermal.h"

// returns the critical current for a segment of a given temperature T
static inline double I_cT_yang_parallel(double I_c0, double T, double T_c) {
    return (I_c0 * (1 - T/T_c*T/T_c)*(1 - T/T_c*T/T_c));
}

// updates the alpha, kappa, c and conductivity values for all segments
// takes into account state dependence etc, specific to the model used in [Yang]
int update_thermal_values_yang_parallel(double * alpha_n, double * kappa_n, double * c_n, double * rho_seg_n, double * R_seg_n, double * T_n, size_t J, double A, double B, double gamma, double T_c, double I_n, double I_c0, double rho_norm, double c_p, double T_ref, double R_seg) {
    // loop over all segments at the current timestep
    for (unsigned j=0; j<J; ++j) {
        // alpha is taken to be state independent, not strictly true, but for more, see [Yang]
        alpha_n[j] = B * pow(T_n[j], 3);
        // model values for kappa, c and rho for normal state
        if (T_n[j] > T_c || I_n > I_cT_yang_parallel(I_c0, T_n[j], T_c)) {
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

int advance_time_electric_yang_parallel(double * I1_np1, double * I2_np1, double * V_c_np1, double I1_n, double I2_n, double V_c_n, double X1, double X2, double Y, double R_w_n, double R_w_np1, double R_p, double R_L, double R_s, double I_b) {
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
    A[4] = X2 + R_L;
    A[5] = -1;

    A[6] = Y;
    A[7] = Y;
    A[8] = 1;

    b[0] = V_c_n + 2*R_L*I_b + (X1 - R_L - R_w_n - R_s)*I1_n - R_L*I2_n;
    b[1] = V_c_n + 2*R_L*I_b + (X2 - R_L - 2*R_p)*I2_n - R_L*I1_n;
    b[2] = V_c_n + Y*(2*I_b - I1_n - I2_n);

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, n, ipiv, b, nrhs);

    if (info != 0) {
        puts("Error encountered in matrix calculation of electrical model...\nExiting with code 3.");
        return 3;
    }

    *I1_np1 = b[0];
    *I2_np1 = b[1];
    *V_c_np1 = b[2];

    free(A);
    free(b);
    free(ipiv);

    return 0;
}

// main function that runs the simulation
int run_yang_parallel(SimRes * res, SimData * data, double dX, double dt) {
    // first locally save some important parameters that we will need all the time
    size_t J = data->J;
    size_t N = data->N;
    size_t NE = data->N*data->ETratio;
    double ** T = res->T[0];
    double * I1 = res->I[0];
    double * I2 = res->I[1];
    double * R = res->R[0];
    double * V_c = res->V_c[0];

    // set up initial thermal values at t = 1, add a steady state time step at t = 0
    for (unsigned j=0; j<J; ++j) {
        T[0][j] = data->T_sub;
        T[1][j] = data->T_sub;
    }
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
        T[1][j] = data->initHS_T_std;
    }
    // set up initial current through the snspd in steady state (t = 0)
    I1[0] = data->R_p_parallel*data->I_b_std / (data->R_s_std + data->R_p_parallel);
    I2[0] = data->R_s_std*data->I_b_std / (data->R_s_std + data->R_p_parallel);
    // set up initial voltage drop over the capacitor
    V_c[0] = (data->R_s_std*data->R_p_parallel)*data->I_b_std / (data->R_s_std + data->R_p_parallel);

    // prepare model parameters for estimating alpha, kappa and c
    // these parameters are considered partially state and temperature dependent
    // formulae for these can be found in Yang
    double Delta = 2.15*Kb*data->T_c*(1 - (data->T_ref_std/data->T_c)*(data->T_ref_std/data->T_c));
    double A = data->c_e*exp(Delta/(data->T_ref_std*Kb));
    double gamma = A/(2.43*data->T_c);
    double B = data->alpha/(pow(data->T_ref_std, 3));

    printf("Delta: %e\nA:     %e\ngamma: %e\nB:     %e\n", Delta, A, gamma, B);

    // define the resistance of a segment of wire in the normal state
    double R_seg = data->rho_norm_std*dX/(data->wireWidth*data->wireThickness);
    // declare the nanowire resistance and current density
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
    int flag = 0;

    // main time loop
    for (unsigned n=1; n<NE; ++n) {
        // print progress
        print_progress(n, NE);
        // advance the thermal model to the next time step after the initial step
        if (n > 1 && n < N) {
            if (cmp_vector(T[n-1], J, data->T_sub, data->T_sub_eps) || !data->allowOpt)
                advance_time_thermal(T[n-1], T[n], J, data->T_sub, alpha_n, c_n, rho_seg_n, kappa_n, data->wireThickness, currentDensity_w, dt, dX);
            else {
                fill_vector(T[n], J, data->T_sub);
                flag = 1;
            }
        }

        if (!flag && n < N) {
            // first update the thermal values used in the differential equation,
            //     the targets are included as the first five parameters
            update_thermal_values_yang_parallel(alpha_n, kappa_n, c_n, rho_seg_n, R_seg_n, T[n], J, A, B, gamma, data->T_c, I1[n-1], data->I_c0, data->rho_norm_std, data->c_p, data->T_ref_std, R_seg);
            // update the current nanowire resistance
            R[n] = sum_vector(R_seg_n, J);
        } else {
            R[n] = 0;
        }

        // update the current density through the nanowire
        currentDensity_w = I1[n-1]/(data->wireWidth*data->wireThickness);
        // update the electric values
        advance_time_electric_yang_parallel(&I1[n], &I2[n], &V_c[n], I1[n-1], I2[n-1], V_c[n-1], X1, X2, Y, R[n-1], R[n], data->R_p_parallel, data->R_L_std, data->R_s_std, data->I_b_std);
    }

    // free allocated space
    free(alpha_n);
    free(kappa_n);
    free(c_n);
    free(rho_seg_n);
    free(R_seg_n);

    // print result
    puts("\nDone.");
    res->exitValue = 0;
    return 0;
}
