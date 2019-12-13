#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "types.h"
#include "helper.h"
#include "thermal.h"

double Kb = 1.3806503E-23;      // boltzmann constant
double Lorentz = 2.45E-8;       // Lorentz number

// returns the critical current for a segment of a given temperature T
double I_cT(double I_c0, double T, double T_c) {
    return (I_c0 * (1 - T/T_c*T/T_c)*(1 - T/T_c*T/T_c));
}

// updates the alpha, kappa, c and conductivity values for all segments
// takes into account state dependence etc, specific to the model used in [Yang]
int update_thermal_values(double * alpha_n, double * kappa_n, double * c_n, double * rho_seg_n, double * T_n, size_t J, double A, double B, double gamma, double T_c, double I, double I_c0, double rho_norm, double c_p, double T_ref) {
    // loop over all segments at the current timestep
    for (unsigned j=0; j<J; ++j) {
        // alpha is taken to be state independent, not strictly true, but for more, see [Yang]
        alpha_n[j] = B * pow(T_n[j], 3);
        // model values for kappa, c and rho for normal state
        if (T_n[j] > T_c || I > I_cT(I_c0, T_n[j], T_c)) {
            kappa_n[j] = Lorentz*T_n[j]/rho_norm;
            c_n[j] = gamma*T_n[j] + c_p*pow(T_n[j]/T_ref, 3);
            rho_seg_n[j] = rho_norm;
        }
        // model values for kappa c and rho for superconducting state
        else {
            kappa_n[j] = Lorentz*T_n[j]/rho_norm * T_n[j]/T_c;
            double Delta = 2.15*T_c*(1 - (T_n[j]/T_c)*(T_n[j]/T_c));
            c_n[j] = A*exp(-Delta/(Kb*T_n[j])) + c_p*pow(T_n[j]/T_ref, 3);
            rho_seg_n[j] = 0;
        }
    }

    return 0;
}

// main function that runs the simulation
int run_yang(SimRes * res, SimData * data, double dX, double dt) {
    // first locally save some important parameters that we will need all the time
    size_t J = data->J;
    size_t N = data->N;
    double ** T = res->T;
    double ** I = res->I;
    double ** R = res->R;

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
    I[0][0] = data->I_b_std;

    // prepare model parameters for estimating alpha, kappa and c
    // these parameters are considered partially state and temperature dependent
    // formulae for these can be found in Yang
    double Delta = 2.15*Kb*data->T_c*(1 - (data->T_ref_std/data->T_c)*(data->T_ref_std/data->T_c));
    double A = data->c_e*exp(Delta/(data->T_ref_std*Kb));
    double gamma = A/(2.43*data->T_c);
    double B = data->alpha/(pow(data->T_ref_std, 3));

    // allocate space for the state and temperature dependent variables for each time step
    double * alpha_n = calloc(J, sizeof(double));
    double * kappa_n = calloc(J, sizeof(double));
    double * c_n = calloc(J, sizeof(double));
    double * rho_seg_n = calloc(J, sizeof(double));

    print_vector(T[0], J);
}
