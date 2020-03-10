#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>

#include "helper.h"
#include "types.h"

// function that fills up the static part of the transmission line matrix when reflections are accounted for
int fill_transmission_matrix_rfl(double * A, size_t NTL, double XT, double YT, double R_L) {
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

// function that fills up the static part of the transmission line matrix without reflections
int fill_transmission_matrix_norfl(double * A, size_t NTL, double XT, double YT, double R_L) {
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

    return 0;
}

int advance_time_transmission_rfl(double * I_tl_np1, double * V_c_tl_np1, double * I_tl_n, double * V_c_tl_n, double * A, double * A_blueprint, double * b, size_t NTL, double XT, double YT, double R_L, double R_r, double I_b_tl) {
    // set up matrix and vector for Ax = b calculation
    lapack_int info;
    lapack_int * ipiv = calloc(2*NTL+1, sizeof(lapack_int));

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

    printf("%4.2e\n", I_tl_np1[0] - I_tl_n[0]);
    //for (unsigned j=0; j<NTL; j++) {
    //    printf("%4.2e %4.2e\n", I_tl_np1[j+1], V_c_tl_np1[j]);
    //}

    free(ipiv);

    return 0;
}
