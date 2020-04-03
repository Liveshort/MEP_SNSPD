#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <cblas.h>
#include <lapacke.h>

#include "helper.h"
#include "types.h"

const double c0 = 299792458;         // speed of light
const double R_r_const = 50;

// function that fills up the static part of the transmission line matrix, to be supplemented by A[0] with R_r
int fill_transmission_matrix_varRr(double * A, size_t NTL, double XT, double YT, double R_L) {
    unsigned rl = 2*NTL+1;

    A[NTL+1] = -1;

    A[rl] = XT;
    A[rl + 1] = XT;
    A[rl + NTL+1] = 1;
    A[rl + NTL+2] = -1;

    for (unsigned i=2; i<NTL; i++) {
        A[i*rl + i] = XT;
        A[i*rl + NTL+i-1] = -1;
        A[i*rl + NTL+i] = 2;
        A[i*rl + NTL+i+1] = -1;
    }

    for (unsigned i=0; i<=NTL; i++)
        A[NTL*rl + i] = R_L;
    A[(NTL+1)*rl - 1] = 1;

    for (unsigned i=1; i<=NTL; i++) {
        A[(NTL+i)*rl + i] = -YT;
        A[(NTL+i)*rl + NTL+i] = 1;
    }

    //for (unsigned i=0; i<rl; i++) {
    //    for (unsigned j=0; j<rl; j++) {
    //        printf("%6.1f ", A[rl*i + j]);
    //    }
    //    puts("");
    //}

    return 0;
}

// function that fills up the transmission line matrix with static R_r and inverts to obtain a 'fast' matrix A
int fill_transmission_matrix_constRr(double * A, size_t NTL, double XT, double YT, double R_L) {
    unsigned rl = 2*NTL+1;

    A[0] = XT + R_r_const;
    A[NTL+1] = -1;

    A[rl] = XT;
    A[rl + 1] = XT;
    A[rl + NTL+1] = 1;
    A[rl + NTL+2] = -1;

    for (unsigned i=2; i<NTL; i++) {
        A[i*rl + i] = XT;
        A[i*rl + NTL+i-1] = -1;
        A[i*rl + NTL+i] = 2;
        A[i*rl + NTL+i+1] = -1;
    }

    for (unsigned i=0; i<=NTL; i++)
        A[NTL*rl + i] = R_L;
    A[(NTL+1)*rl - 1] = 1;

    for (unsigned i=1; i<=NTL; i++) {
        A[(NTL+i)*rl + i] = -YT;
        A[(NTL+i)*rl + NTL+i] = 1;
    }

    //for (unsigned i=0; i<rl; i++) {
    //    for (unsigned j=0; j<rl; j++) {
    //        printf("%6.1f ", A[rl*i + j]);
    //    }
    //    puts("");
    //}

    // set up matrix and vector for matrix inversion (matrix is static, yay!)
    lapack_int info;
    lapack_int * ipiv = calloc(rl, sizeof(lapack_int));

    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, rl, rl, A, rl, ipiv);

    if (info != 0) {
        puts("Error encountered in matrix calculation of transmission line model...\nExiting with code 7.");
        exit(7);
    }

    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, rl, A, rl, ipiv);

    if (info != 0) {
        puts("Error encountered in matrix calculation of transmission line model...\nExiting with code 7.");
        exit(7);
    }

    //for (unsigned i=0; i<rl; i++) {
    //    for (unsigned j=0; j<rl; j++) {
    //        printf("%8.1g ", A[rl*i + j]);
    //    }
    //    puts("");
    //}

    return 0;
}

// function that fills up the transmission line matrix with static R_r and inverts to obtain a 'fast' matrix A
int fill_transmission_matrix_noRr(double * A, size_t NTL, double XT, double YT, double R_L) {
    unsigned rl = 2*NTL;

    A[0] = XT;
    A[NTL] = 1;
    A[NTL+1] = -1;

    for (unsigned i=1; i<NTL-1; i++) {
        A[i*rl + i] = XT;
        A[i*rl + NTL+i-1] = -1;
        A[i*rl + NTL+i] = 2;
        A[i*rl + NTL+i+1] = -1;
    }

    for (unsigned i=0; i<NTL; i++)
        A[(NTL-1)*rl + i] = R_L;
    A[NTL*rl - 1] = 1;

    for (unsigned i=0; i<NTL; i++) {
        A[(NTL+i)*rl + i] = -YT;
        A[(NTL+i)*rl + NTL+i] = 1;
    }

    //for (unsigned i=0; i<rl; i++) {
    //    for (unsigned j=0; j<rl; j++) {
    //        printf("%6.1f ", A[rl*i + j]);
    //    }
    //    puts("");
    //}

    // set up matrix and vector for matrix inversion (matrix is static, yay!)
    lapack_int info;
    lapack_int * ipiv = calloc(rl, sizeof(lapack_int));

    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, rl, rl, A, rl, ipiv);

    if (info != 0) {
        puts("Error encountered in matrix calculation of transmission line model...\nExiting with code 7.");
        exit(7);
    }

    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, rl, A, rl, ipiv);

    if (info != 0) {
        puts("Error encountered in matrix calculation of transmission line model...\nExiting with code 7.");
        exit(7);
    }

    //for (unsigned i=0; i<rl; i++) {
    //    for (unsigned j=0; j<rl; j++) {
    //        printf("%8.1g ", A[rl*i + j]);
    //    }
    //    puts("");
    //}

    return 0;
}

// advance time for variable R_r by solving the Ax=b problem every time step (time consuming)
int advance_time_transmission_varRr(double * I_tl_np1, double * V_c_tl_np1, double * I_tl_n, double * V_c_tl_n, double * A, double * A_blueprint, double * b, lapack_int * ipiv, size_t NTL, double XT, double YT, double R_L, double R_r, double I_b_tl) {
    // set up matrix and vector for Ax = b calculation
    lapack_int info;
    //lapack_int * ipiv = calloc(2*NTL+1, sizeof(lapack_int));

    // copy the blueprint matrix
    memcpy(A, A_blueprint, (2*NTL+1)*(2*NTL+1)*sizeof(double));

    // update first matrix element
    A[0] = R_r + XT;

    //for (unsigned i=0; i<(2*NTL+1); i++) {
    //    for (unsigned j=0; j<(2*NTL+1); j++) {
    //        printf("%6.1f ", A[(2*NTL+1)*i + j]);
    //    }
    //    puts("");
    //}

    // update rhs of equation
    b[0] = (XT-R_r)*I_tl_n[0] + V_c_tl_n[0];
    b[1] = XT*(I_tl_n[0] + I_tl_n[1]) - V_c_tl_n[0] + V_c_tl_n[1];
    for (unsigned i=2; i<NTL; i++)
        b[i] = XT*I_tl_n[i] + V_c_tl_n[i-2] - 2*V_c_tl_n[i-1] + V_c_tl_n[i];

    b[NTL] = -V_c_tl_n[NTL-1] + 2*I_b_tl*R_L - sum_vector(I_tl_n, NTL+1)*R_L;

    for (unsigned i=0; i<NTL; i++)
        b[NTL+1 + i] = V_c_tl_n[i] + YT*I_tl_n[i+1];

    //for (unsigned j=0; j<(2*NTL+1); j++) {
    //    printf("%4.2e\n", b[j]);
    //}

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, 2*NTL+1, 1, A, 2*NTL+1, ipiv, b, 1);

    if (info != 0) {
        puts("Error encountered in matrix calculation of transmission line model...\nExiting with code 7.");
        exit(7);
    }

    // put the currents and voltages in the right spot
    for (unsigned i=0; i<=NTL; i++)
        I_tl_np1[i] = b[i];
    for (unsigned i=0; i<NTL; i++)
        V_c_tl_np1[i] = b[NTL+1 + i];

    //printf("%4.2e\n", I_tl_np1[0] - I_tl_n[0]);
    //for (unsigned j=0; j<NTL; j++) {
    //    printf("%4.2e %4.2e\n", I_tl_np1[j+1], V_c_tl_np1[j]);
    //}

    //free(ipiv);

    return 0;
}

// advance time for static R_r by multiplying inv(A)b to obtain x (faster, but still time consuming)
int advance_time_transmission_constRr(double * I_tl_np1, double * V_c_tl_np1, double * I_tl_n, double * V_c_tl_n, double * A, double * b, size_t NTL, double XT, double YT, double R_L, double I_b_tl_np1, double I_b_tl_n) {
    // set up matrix and vector for Ax = b calculation
    lapack_int info;
    double * tmp = calloc(2*NTL+1, sizeof(double));

    // update rhs of equation
    b[0] = (XT - R_r_const)*I_tl_n[0] + V_c_tl_n[0] + XT*(I_b_tl_np1 - I_b_tl_n);
    b[1] = XT*(I_tl_n[0] + I_tl_n[1]) - V_c_tl_n[0] + V_c_tl_n[1] + XT*(I_b_tl_np1 - I_b_tl_n);
    for (unsigned i=2; i<NTL; i++)
        b[i] = XT*I_tl_n[i] + V_c_tl_n[i-2] - 2*V_c_tl_n[i-1] + V_c_tl_n[i];

    b[NTL] = -V_c_tl_n[NTL-1] + 2*I_b_tl_np1*R_L - sum_vector(I_tl_n, NTL+1)*R_L;

    for (unsigned i=0; i<NTL; i++)
        b[NTL+1 + i] = V_c_tl_n[i] + YT*I_tl_n[i+1];

    //printf("%4.2e\n", I_b_tl);

    //puts("before mv:");
    //printf("%4.2e\n", b[0]);
    //for (unsigned j=0; j<NTL; j++) {
    //    printf("%4.2e %4.2e\n", b[j+1], b[j+1+NTL]);
    //}

    cblas_dgemv(CblasRowMajor, CblasNoTrans, 2*NTL+1, 2*NTL+1, 1, A, 2*NTL+1, b, 1, 0, tmp, 1);

    // put the currents and voltages in the right spot
    for (unsigned i=0; i<NTL+1; i++)
        I_tl_np1[i] = tmp[i];
    for (unsigned i=0; i<NTL; i++)
        V_c_tl_np1[i] = tmp[NTL+1 + i];

    //puts("after mv:");
    //printf("%4.2e\n", I_tl_np1[0]);
    //for (unsigned j=0; j<NTL; j++) {
    //    printf("%4.2e %4.2e\n", I_tl_np1[j+1], V_c_tl_np1[j]);
    //}

    free(tmp);

    return 0;
}

// advance time for static R_r by multiplying inv(A)b to obtain x (faster, but still time consuming)
int advance_time_transmission_noRr(double * I_tl_np1, double * V_c_tl_np1, double * I_tl_n, double * V_c_tl_n, double * A, double * b, size_t NTL, double XT, double YT, double R_L, double I_b_tl_np1, double I_b_tl_n) {
    // set up matrix and vector for Ax = b calculation
    lapack_int info;
    double * tmp = calloc(2*NTL, sizeof(double));

    // update rhs of equation
    //printf("vals %4.2e %4.2e %4.2e %4.2e %4.2e\n", I_tl_n[0], V_c_tl_n[0], V_c_tl_n[1], I_b_tl_np1, I_b_tl_n);
    b[0] = XT*I_tl_n[0] - V_c_tl_n[0] + V_c_tl_n[1] + XT*(I_b_tl_np1 - I_b_tl_n);
    for (unsigned i=1; i<NTL-1; i++)
        b[i] = XT*I_tl_n[i] + V_c_tl_n[i-1] - 2*V_c_tl_n[i] + V_c_tl_n[i+1];
    b[NTL-1] = -V_c_tl_n[NTL-1] + 2*I_b_tl_np1*R_L - sum_vector(I_tl_n, NTL)*R_L;

    for (unsigned i=0; i<NTL; i++)
        b[NTL + i] = V_c_tl_n[i] + YT*I_tl_n[i];

    //printf("%4.2e\n", I_b_tl);

    //puts("before mv:");
    //for (unsigned j=0; j<NTL; j++) {
    //    printf("%4.2e %4.2e\n", b[j], b[j+NTL]);
    //}

    cblas_dgemv(CblasRowMajor, CblasNoTrans, 2*NTL, 2*NTL, 1, A, 2*NTL, b, 1, 0, tmp, 1);

    // put the currents and voltages in the right spot
    for (unsigned i=0; i<NTL; i++)
        I_tl_np1[i] = tmp[i];
    for (unsigned i=0; i<NTL; i++)
        V_c_tl_np1[i] = tmp[NTL + i];

    //puts("after mv:");
    //for (unsigned j=0; j<NTL; j++) {
    //    printf("%4.2e %4.2e\n", I_tl_np1[j], V_c_tl_np1[j]);
    //}

    free(tmp);

    return 0;
}

// function that simulates transmission line without feedback into the system
int sim_transmission_line(SimData * data, SimRes * res, size_t NE, size_t NTL, size_t numberOfIb) {
    double delay = data->LTL/(data->VF*c0);
    unsigned steps = (unsigned) (delay/res->dt);

    if (steps >= NE) {
        printf("\nError: estimated delay (%4.2e [s] ~~ %u steps) is larger than simulation space (%ld)...\nTreating transmission line type as 0 instead", delay, steps, NE);

        return -1;
    }

    double C_T = delay/(50*NTL);
    double L_T = (50*delay)/NTL;

    if (data->simTL == DELAYTRNS) {
        double * Iload = res->I[data->numberOfI-3];

        printf("\nTransmission line delay: %4.2e s ==> %u steps\nComputing...", delay, steps);
        for (unsigned n=NE-1; n>=steps; --n)
            Iload[n] = Iload[n-steps];
        for (unsigned n=0; n<steps; ++n)
            Iload[n] = 0;
    } else if (data->simTL == VARRPLCTRNS) {
        printf("\nTransmission line properties: C = %4.2e F, L = %4.2e H", C_T, L_T);
        // calculate transmission line constants
        double XT = (2*L_T)/res->dt;
        double YT = res->dt/(2*C_T);

        // set up the currents for the transmission line
        double * I_tl_prev = calloc(NTL+1, sizeof(double));
        double * I_tl_curr = calloc(NTL+1, sizeof(double));
        double * V_c_tl_prev = calloc(NTL, sizeof(double));
        double * V_c_tl_curr = calloc(NTL, sizeof(double));
        double * Iload = res->I[data->numberOfI-3];
        double * Itlsum = res->I[data->numberOfI-2];
        double * Itl = res->I[data->numberOfI-1];

        // get the total bias current
        double v_I = sum_vector(res->I_b, numberOfIb);
        // set up initial condition for the transmission line model
        I_tl_prev[0] = v_I;
        // set up the transmission matrix
        double * AT = calloc((2*NTL+1)*(2*NTL+1), sizeof(double));
        double * AT_tmp = calloc((2*NTL+1)*(2*NTL+1), sizeof(double));
        double * bT = calloc(2*NTL+1, sizeof(double));
        fill_transmission_matrix_varRr(AT, NTL, XT, YT, data->R_L_std);
        lapack_int * ipiv = calloc(2*NTL+1, sizeof(lapack_int));

        puts("\nTransmission line loop:");
        for (unsigned n=data->timeskip+1; n<NE; ++n) {
            // print progress
            print_progress(n, NE);

            // update the transmission line
            double Isum_n = 0;
            for (unsigned i=0; i<res->numberOfI-3; ++i)
                Isum_n += res->I[i][n];

            double R_r = data->R_L_std*(v_I/Isum_n - 1);
            advance_time_transmission_varRr(I_tl_curr, V_c_tl_curr, I_tl_prev, V_c_tl_prev, AT_tmp, AT, bT, ipiv, NTL, XT, YT, data->R_L_std, R_r, v_I);
            Itlsum[n] = subsum_vector(I_tl_curr, 1, NTL+1);
            Itl[n] = Isum_n - I_tl_curr[0];
            Iload[n] = v_I - I_tl_curr[0] - Itlsum[n];

            // shuffle the transmission currents around for the next timestep
            swap_ptr((void **) &I_tl_prev, (void **) &I_tl_curr);
            swap_ptr((void **) &V_c_tl_prev, (void **) &V_c_tl_curr);
        }

        free(I_tl_prev);
        free(I_tl_curr);
        free(V_c_tl_prev);
        free(V_c_tl_curr);
        free(AT);
        free(AT_tmp);
        free(bT);
        free(ipiv);
    } else if (data->simTL == CONSTRPLCTRNS) {
        printf("\nTransmission line properties: C = %4.2e F, L = %4.2e H", C_T, L_T);
        // calculate transmission line constants
        double XT = (2*L_T)/res->dt;
        double YT = res->dt/(2*C_T);

        // set up the currents for the transmission line
        double * I_tl_prev = calloc(NTL+1, sizeof(double));
        double * I_tl_curr = calloc(NTL+1, sizeof(double));
        double * V_c_tl_prev = calloc(NTL, sizeof(double));
        double * V_c_tl_curr = calloc(NTL, sizeof(double));
        double I_b_tl_prev, I_b_tl_curr;
        double * Iload = res->I[data->numberOfI-3];
        double * Itlsum = res->I[data->numberOfI-2];
        double * Itl = res->I[data->numberOfI-1];

        // get the total bias current
        double v_I = sum_vector(res->I_b, numberOfIb);
        // set up initial condition for the transmission line model
        I_b_tl_prev = 0;
        // set up the transmission matrix
        double * AT = calloc((2*NTL+1)*(2*NTL+1), sizeof(double));
        double * bT = calloc(2*NTL+1, sizeof(double));
        fill_transmission_matrix_constRr(AT, NTL, XT, YT, data->R_L_std);

        puts("\nTransmission line loop:");
        for (unsigned n=data->timeskip+1; n<NE; ++n) {
            // print progress
            print_progress(n, NE);

            // update the transmission line
            double Isum_n = 0;
            for (unsigned i=0; i<res->numberOfI-3; ++i)
                Isum_n += res->I[i][n];

            I_b_tl_curr = 2*(v_I - Isum_n);
            advance_time_transmission_constRr(I_tl_curr, V_c_tl_curr, I_tl_prev, V_c_tl_prev, AT, bT, NTL, XT, YT, data->R_L_std, I_b_tl_curr, I_b_tl_prev);
            Itlsum[n] = subsum_vector(I_tl_curr, 1, NTL+1);
            Itl[n] = 0;
            Iload[n] = I_b_tl_curr - I_tl_curr[0] - Itlsum[n];

            // shuffle the transmission currents around for the next timestep
            swap_ptr((void **) &I_tl_prev, (void **) &I_tl_curr);
            swap_ptr((void **) &V_c_tl_prev, (void **) &V_c_tl_curr);
            swap_dbl(&I_b_tl_curr, &I_b_tl_prev);
        }

        free(I_tl_prev);
        free(I_tl_curr);
        free(V_c_tl_prev);
        free(V_c_tl_curr);
        free(AT);
        free(bT);
    } else if (data->simTL == NORPLCTRNS) {
        printf("\nTransmission line properties: C = %4.2e F, L = %4.2e H", C_T, L_T);
        // calculate transmission line constants
        double XT = (2*L_T)/res->dt;
        double YT = res->dt/(2*C_T);

        // set up the currents for the transmission line
        double * I_tl_prev = calloc(NTL, sizeof(double));
        double * I_tl_curr = calloc(NTL, sizeof(double));
        double * V_c_tl_prev = calloc(NTL, sizeof(double));
        double * V_c_tl_curr = calloc(NTL, sizeof(double));
        double I_b_tl_prev, I_b_tl_curr;
        double * Iload = res->I[data->numberOfI-3];
        double * Itlsum = res->I[data->numberOfI-2];
        double * Itl = res->I[data->numberOfI-1];

        // get the total bias current
        double v_I = sum_vector(res->I_b, numberOfIb);
        // set up initial condition for the transmission line model
        I_b_tl_prev = 0;
        // set up the transmission matrix
        double * AT = calloc(2*NTL*2*NTL, sizeof(double));
        double * bT = calloc(2*NTL, sizeof(double));
        fill_transmission_matrix_noRr(AT, NTL, XT, YT, data->R_L_std);

        puts("\nTransmission line loop:");
        for (unsigned n=data->timeskip+1; n<NE; ++n) {
            // print progress
            print_progress(n, NE);

            // update the transmission line
            double Isum_n = 0;
            for (unsigned i=0; i<res->numberOfI-3; ++i)
                Isum_n += res->I[i][n];

            I_b_tl_curr = v_I - Isum_n;
            advance_time_transmission_noRr(I_tl_curr, V_c_tl_curr, I_tl_prev, V_c_tl_prev, AT, bT, NTL, XT, YT, data->R_L_std, I_b_tl_curr, I_b_tl_prev);
            Itlsum[n] = sum_vector(I_tl_curr, NTL);
            Itl[n] = 0;
            Iload[n] = I_b_tl_curr - Itlsum[n];

            // shuffle the transmission currents around for the next timestep
            swap_ptr((void **) &I_tl_prev, (void **) &I_tl_curr);
            swap_ptr((void **) &V_c_tl_prev, (void **) &V_c_tl_curr);
            swap_dbl(&I_b_tl_curr, &I_b_tl_prev);
        }

        free(I_tl_prev);
        free(I_tl_curr);
        free(V_c_tl_prev);
        free(V_c_tl_curr);
        free(AT);
        free(bT);
    }

    return 0;
}
