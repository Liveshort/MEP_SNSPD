#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "linalg.h"
#include "helper.h"

// function that advances time a single time step. It sets up the tridiagonal matrix and hands
//     it to the solver in the linalg code.
int advance_time_thermal(double * T_n, double * T_np1, size_t J, double T_sub, double * alpha_n, double * c_n, double * rho_seg_n, double * kappa_n, double wireThickness, double currentDensity_w, double dt, double dX) {
    // allocate some space for the matrix defining diagonals and off diagonals and the right hand side
    double * lhsDiag_n = calloc(J, sizeof(double));
    double * lhsOffDiag_n = calloc(J, sizeof(double));
    double * rhs_n = calloc(J, sizeof(double));
    double r, w, q;

    // fill in the boundaries
    lhsDiag_n[0] = lhsDiag_n[J-1] = 1;
    rhs_n[0] = rhs_n[J-1] = T_sub;
    // fill in the rest
    for (unsigned j=1; j<J-1; ++j) {
        // calculate some convenient constants per segment
        r = kappa_n[j]*dt / (2*dX*dX*c_n[j]);
        w = alpha_n[j]*dt / (2*c_n[j]*wireThickness);
        q = currentDensity_w*currentDensity_w*rho_seg_n[j]*dt/c_n[j] + 2*w*T_sub;
        // fill up the matrix diagonal, off diagonal and rhs
        lhsDiag_n[j] = 1 + w + 2*r;
        lhsOffDiag_n[j] = -r;
        rhs_n[j] = r*(T_n[j-1] + T_n[j+1]) + (1 - w - 2*r)*T_n[j] + q;
    }
    //puts("Matrix values");
    //print_vector(lhsOffDiag_n, J);
    //print_vector(lhsDiag_n, J);
    //print_vector(rhs_n, J);

    // solve the tridiagonal matrix
    TDM_solve(T_np1, J, lhsDiag_n, lhsOffDiag_n, rhs_n);

    //puts("T_n and T_np1");
    //print_vector(T_n, J);
    //print_vector(T_np1, J);
    //puts("kappa, alpha, c, rho_seg");

    free(lhsDiag_n);
    free(lhsOffDiag_n);
    free(rhs_n);

    return 0;
}
