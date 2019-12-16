#ifndef __THERMAL_H__
#define __THERMAL_H__

int TDM_solve(double * T_np1, size_t J, double * lhsDiag_n, double * lhsOffDiag_n, double * rhs);

#endif
