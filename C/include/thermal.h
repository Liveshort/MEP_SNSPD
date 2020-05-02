#ifndef __THERMAL_H__
#define __THERMAL_H__

int advance_time_thermal(double * T_n, double * T_np1, size_t J, double T_sub, double * alpha_n, double * c_n, double * rho_seg_n, double * kappa_n, double wireThickness, double currentDensity_w, double dt, double dX);

int update_thermal_values(double * alpha_n, double * kappa_n, double * c_n, double * rho_seg_n, double * R_seg_n, double * T_n, size_t J, double A, double B, double gamma, double T_c, double I_n, double * I_c0, double rho_norm, double c_p, double T_ref, double R_seg);

#endif
