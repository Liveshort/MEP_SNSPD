#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <lapacke.h>

#include "types.h"
#include "helper.h"
#include "thermal.h"
#include "lapacke_example_aux.h"

double Kb = 1.3806503E-23;      // boltzmann constant
double Lorentz = 2.45E-8;       // Lorentz number

// returns the critical current for a segment of a given temperature T
double I_cT(double I_c0, double T, double T_c) {
    return (I_c0 * (1 - T/T_c*T/T_c)*(1 - T/T_c*T/T_c));
}

void test_lapack() {
    lapack_int n, nrhs, lda, ldb, info;
    int i, j;
    double *A, *b;
    lapack_int *ipiv;
    n = 5; nrhs = 1;

    lda=n, ldb=nrhs;
    A = (double *)malloc(n*n*sizeof(double)) ;
    if (A==NULL){ printf("error of memory allocation\n"); exit(0); }
    b = (double *)malloc(n*nrhs*sizeof(double)) ;
    if (b==NULL){ printf("error of memory allocation\n"); exit(0); }
    ipiv = (lapack_int *)malloc(n*sizeof(lapack_int)) ;
    if (ipiv==NULL){ printf("error of memory allocation\n"); exit(0); }

    for( i = 0; i < n; i++ ) {
            for( j = 0; j < n; j++ ) A[i*lda+j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
    }

    for(i=0;i<n*nrhs;i++)
        b[i] = ((double) rand()) / ((double) RAND_MAX) - 0.5;

    /* Print Entry Matrix */
    print_matrix_rowmajor( "Entry Matrix A", n, n, A, lda );
    /* Print Right Rand Side */
    print_matrix_rowmajor( "Right Rand Side b", n, nrhs, b, ldb );
    printf( "\n" );

    info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n, nrhs, A, lda, ipiv,
                    b, ldb );
    /* Check for the exact singularity */
    if( info > 0 ) {
            printf( "The diagonal element of the triangular factor of A,\n" );
            printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
            printf( "the solution could not be computed.\n" );
            exit( 1 );
    }
    if (info <0) exit( 1 );

    /* Print solution */
    print_matrix_rowmajor( "Solution", n, nrhs, b, ldb );
    /* Print details of LU factorization */
    print_matrix_rowmajor( "Details of LU factorization", n, n, A, lda );
    /* Print pivot indices */
    print_vector_lpk( "Pivot indices", n, ipiv );

    free(A);
    free(b);
    free(ipiv);
}

// updates the alpha, kappa, c and conductivity values for all segments
// takes into account state dependence etc, specific to the model used in [Yang]
int update_thermal_values(double * alpha_n, double * kappa_n, double * c_n, double * rho_seg_n, double * R_seg_n, double * T_n, size_t J, double A, double B, double gamma, double T_c, double I_n, double I_c0, double rho_norm, double c_p, double T_ref, double R_seg) {
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

    print_vector(kappa_n, J);
    print_vector(alpha_n, J);
    print_vector(c_n, J);
    print_vector(rho_seg_n, J);
    puts("");

    return 0;
}

// main function that runs the simulation
int run_yang(SimRes * res, SimData * data, double dX, double dt) {
    // first locally save some important parameters that we will need all the time
    size_t J = data->J;
    size_t N = data->N;
    double ** T = res->T[0];
    double * I = res->I[0];
    double * R = res->R[0];

    // set up initial thermal values at t = 0
    for (unsigned j=0; j<J; ++j) {
        T[0][j] = data->T_sub;
    }
    // determine halfway point and set up a beginning hotspot
    unsigned halfway = J/2;
    unsigned initHS_segs = (unsigned) (data->initHS_l_std/dX) + 1;
    // check if there is a nonzero number of segments
    if (initHS_segs < 2) {
        puts("Number of segments in initial hot-spot smaller than 2.\nReturning empty result with error code 2 (wrong initial hot-spot size)...");
        res->exitValue = 2;
        return 2;
    }
    for (unsigned j=halfway - initHS_segs/2; j<halfway + initHS_segs/2; ++j) {
        T[0][j] = data->initHS_T_std;
    }
    // set up initial current through the snspd
    I[0] = data->I_b_std;

    // prepare model parameters for estimating alpha, kappa and c
    // these parameters are considered partially state and temperature dependent
    // formulae for these can be found in Yang
    double Delta = 2.15*Kb*data->T_c*(1 - (data->T_ref_std/data->T_c)*(data->T_ref_std/data->T_c));
    double A = data->c_e*exp(Delta/(data->T_ref_std*Kb));
    double gamma = A/(2.43*data->T_c);
    double B = data->alpha/(pow(data->T_ref_std, 3));

    // define the resistance of a segment of wire in the normal state
    double R_seg = data->rho_norm_std*dX/(data->wireWidth*data->wireThickness);
    // declare the nanowire resistance and current density
    double R_w;
    double currentDensity_w;

    // allocate space for the state and temperature dependent variables for each time step
    double * alpha_n = calloc(J, sizeof(double));
    double * kappa_n = calloc(J, sizeof(double));
    double * c_n = calloc(J, sizeof(double));
    double * rho_seg_n = calloc(J, sizeof(double));
    double * R_seg_n = calloc(J, sizeof(double));

    // main time loop
    for (unsigned n=0; n<2; ++n) {

        // TODO: optimization with T_sub_eps

        // first update the thermal values used in the differential equation,
        //     the targets are included as the first five parameters
        update_thermal_values(alpha_n, kappa_n, c_n, rho_seg_n, R_seg_n, T[n], J, A, B, gamma, data->T_c, I[n], data->I_c0, data->rho_norm_std, data->c_p, data->T_ref_std, R_seg);
        // update the current nanowire resistance
        R_w = sum_vector(R_seg_n, J);
        // update the current density through the nanowire
        currentDensity_w = I[n]/(data->wireWidth*data->wireThickness);
        // advance the thermal model to the next time step
        advance_time(T[n], T[n+1], J, data->T_sub, alpha_n, c_n, rho_seg_n, kappa_n, data->wireThickness, currentDensity_w, dt, dX);
    }

    print_vector(T[0], J);

    test_lapack();

    // free allocated space
    free(alpha_n);
    free(kappa_n);
    free(c_n);
    free(rho_seg_n);
    free(R_seg_n);
}
