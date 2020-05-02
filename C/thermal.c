#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "linalg.h"
#include "helper.h"
#include "types.h"

// returns the critical current for a segment of a given temperature T
static inline double I_cT(double I_c0, double T, double T_c) {
    return (I_c0 * (1 - T/T_c*T/T_c)*(1 - T/T_c*T/T_c));
}

// updates the alpha, kappa, c and conductivity values for all segments
// takes into account state dependence etc, specific to the model used in [Yang]
int update_thermal_values(double * alpha_n, double * kappa_n, double * c_n, double * rho_seg_n, double * R_seg_n, double * T_n, size_t J, double A, double B, double gamma, double T_c, double I_n, double * I_c0, double rho_norm, double c_p, double T_ref, double R_seg) {
    // loop over all segments at the current timestep
    for (unsigned j=0; j<J; ++j) {
        // alpha is taken to be state independent, not strictly true, but for more, see [Yang]
        alpha_n[j] = B * pow(T_n[j], 3);
        // model values for kappa, c and rho for normal state
        if (T_n[j] > T_c || I_n > I_cT(I_c0[j], T_n[j], T_c)) {
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
