#ifndef __THERMAL_H__
#define __THERMAL_H__

int advance_time_thermal(double * T_n, double * T_np1, size_t J, double T_sub, double * alpha_n, double * c_n, double * rho_seg_n, double * kappa_n, double wireThickness, double currentDensity_w, double dt, double dX);

#endif
