#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <lapacke.h>

#include "types.h"
#include "helper.h"
#include "thermal.h"
#include "lapacke_example_aux.h"

// returns the critical current for a segment of a given temperature T
double I_cT(double I_c0, double T, double T_c) {
    //puts("critical current:");
    //printf("%4.2e %4.2e %4.2e %4.2e\n", I_c0, T, T_c, I_c0 * (1 - T/T_c*T/T_c)*(1 - T/T_c*T/T_c));
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
        if (T_n[j] > T_c || I_n > I_cT(I_c0, T_n[j], T_c)) {
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

    //print_vector(alpha_n, J);
    //print_vector(kappa_n, J);
    //print_vector(c_n, J);
    //print_vector(rho_seg_n, J);
    //puts("");

    return 0;
}

int advance_time_electric_yang_parallel(double * I_np1, double * V_c_np1, double I_n, double V_c_n, double X, double Y, double R_w_n, double R_w_np1, double R_L, double I_b) {
    // set up matrix and vector for Ax = b calculation
    lapack_int n, nrhs, info;
    n = 2; nrhs = 1;

    double * A = calloc(n*n, sizeof(double));
    double * b = calloc(n*nrhs, sizeof(double));
    lapack_int * ipiv = calloc(n, sizeof(lapack_int));

    A[0] = X + R_L + R_w_np1;
    A[1] = -1;
    A[2] = Y;
    A[3] = 1;

    b[0] = V_c_n + 2*R_L*I_b + (X - R_L - R_w_n)*I_n;
    b[1] = V_c_n + Y*(2*I_b - I_n);

    //print_matrix_rowmajor( "Entry Matrix A", n, n, A, n );
    //print_matrix_rowmajor( "Right Rand Side b", n, nrhs, b, nrhs );

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, n, ipiv, b, nrhs);

    if (info != 0) {
        puts("Error encountered in matrix calculation of electrical model...\nExiting with code 3.");
        return 3;
    }

    //print_matrix_rowmajor( "Solution", n, nrhs, b, nrhs );

    *I_np1 = b[0];
    *V_c_np1 = b[1];

    free(A);
    free(b);
    free(ipiv);

    return 0;
}

int advance_time_electric_basic_yang_parallel(double * I_np1, double * V_c_np1, double I_n, double V_c_n, double X, double Y, double R_w_n, double R_w_np1, double R_L, double I_b) {
    double k1 = X*I_n + 2*R_L*I_b - R_L*I_n - R_w_n*I_n + V_c_n;
    double k2 = Y*(2*I_b - I_n) + V_c_n;

    *I_np1 = (k1 + k2) / (X + R_L + R_w_np1 + Y);
    *V_c_np1 = k2 - Y*I_n;

    return 0;
}

// main function that runs the simulation
int run_yang_parallel(SimRes * res, SimData * data, double dX, double dt) {
    // first locally save some important parameters that we will need all the time
    size_t J = data->J;
    size_t N = data->N;
    double ** T = res->T[0];
    double * I = res->I[0];
    double * R = res->R[0];

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
    I[0] = data->I_b_std;

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
    double V_c_nm1 = 0;
    double V_c_n;

    // allocate space for the state and temperature dependent variables for each time step
    double * alpha_n = calloc(J, sizeof(double));
    double * kappa_n = calloc(J, sizeof(double));
    double * c_n = calloc(J, sizeof(double));
    double * rho_seg_n = calloc(J, sizeof(double));
    double * R_seg_n = calloc(J, sizeof(double));

    // set up two characteristic numbers for the electrical calculations
    double X = (2*data->L_w_std)/dt;
    double Y = dt/(2*data->C_m_std);

    // main time loop
    for (unsigned n=1; n<N; ++n) {
        // print progress
        print_progress(n, N);
        // advance the thermal model to the next time step after the initial step
        if (n > 1) {
            advance_time_thermal(T[n-1], T[n], J, data->T_sub, alpha_n, c_n, rho_seg_n, kappa_n, data->wireThickness, currentDensity_w, dt, dX);
        }

        // TODO: optimization with T_sub_eps

        // first update the thermal values used in the differential equation,
        //     the targets are included as the first five parameters
        update_thermal_values_yang_parallel(alpha_n, kappa_n, c_n, rho_seg_n, R_seg_n, T[n], J, A, B, gamma, data->T_c, I[n-1], data->I_c0, data->rho_norm_std, data->c_p, data->T_ref_std, R_seg);

        //puts("R_seg");
        //print_vector(R_seg_n, J);

        // update the current nanowire resistance
        R[n] = sum_vector(R_seg_n, J);
        // update the current density through the nanowire
        currentDensity_w = I[n-1]/(data->wireWidth*data->wireThickness);
        // update the electric values
        advance_time_electric_yang_parallel(&I[n], &V_c_n, I[n-1], V_c_nm1, X, Y, R[n-1], R[n], data->R_L_std, data->I_b_std);
        // shift the data into the right position for the next loop
        //printf("V_c: %4.2e\n", V_c_n);
        V_c_nm1 = V_c_n;

        //puts("Current I, and resistance R:");
        //print_vector(I, 10);
        //print_vector(R, 10);
        //puts("");
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
