#ifndef __TRANSMISSION_H__
#define __TRANSMISSION_H__

// some global constants
const double c0;            // speed of light

int sim_transmission_line(SimData * data, SimRes * res, size_t NE, size_t NTL, size_t numberOfIb);
//int fill_transmission_matrix_rfl(double * A, size_t NTL, double XT, double YT, double R_L);
//int fill_transmission_matrix_norfl(double * A, size_t NTL, double XT, double YT, double R_L);
//int fill_transmission_matrix_norfl_eq(double * A, size_t NTL, double XT, double YT, double R_L);
//int advance_time_transmission_rfl(double * I_tl_np1, double * V_c_tl_np1, double * I_tl_n, double * V_c_tl_n, double * A, double * A_blueprint, double * b, lapack_int * ipiv, size_t NTL, double XT, double YT, double R_L, double R_r, double I_b_tl);
//int advance_time_transmission_norfl(double * I_tl_np1, double * V_c_tl_np1, double * I_tl_n, double * V_c_tl_n, double * A, double * b, size_t NTL, double XT, double YT, double R_L, double I_b_tl);
//int advance_time_transmission_norfl_eq(double * I_tl_np1, double * V_c_tl_np1, double * I_tl_n, double * V_c_tl_n, double * A, double * b, lapack_int * ipiv, size_t NTL, double XT, double YT, double R_L, double I_b_tl);

#endif
