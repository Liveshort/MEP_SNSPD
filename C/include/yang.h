#ifndef __YANG_H__
#define __YANG_H__

double I_cT(double I_c0, double T, double T_c);
int run_yang(SimRes * res, SimData * data, double dX, double dt);
int update_thermal_values(double * alpha_n, double * kappa_n, double * c_n, double * rho_seg_n, double * R_seg_n, double * T_n, size_t J, double A, double B, double gamma, double T_c, double * I_n, double I_c0, double rho_norm, double c_p, double T_ref, double R_seg);

#endif
