#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <lapacke.h>

#include "types.h"
#include "helper.h"
#include "thermal.h"

// returns the critical current for a segment of a given temperature T
static inline double I_cT_3_wtf_res_par(double I_c0, double T, double T_c) {
    return (I_c0 * (1 - T/T_c*T/T_c)*(1 - T/T_c*T/T_c));
}

// updates the alpha, kappa, c and conductivity values for all segments
// takes into account state dependence etc, specific to the model used in [Yang]
int update_thermal_values_3_wtf_res_par(double * alpha_n, double * kappa_n, double * c_n, double * rho_seg_n, double * R_seg_n, double * T_n, size_t J, double A, double B, double gamma, double T_c, double I_n, double I_c0, double rho_norm, double c_p, double T_ref, double R_seg) {
    // loop over all segments at the current timestep
    for (unsigned j=0; j<J; ++j) {
        // alpha is taken to be state independent, not strictly true, but for more, see [Yang]
        alpha_n[j] = B * pow(T_n[j], 3);
        // model values for kappa, c and rho for normal state
        if (T_n[j] > T_c || I_n > I_cT_3_wtf_res_par(I_c0, T_n[j], T_c)) {
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

// advances the electric model one timestep using BLAS matrix algebra.
// you need the LAPACK for this, or work out the matrix logic on your own (not recommended)
int advance_time_electric_3_wtf_res_par(double * I0_np1, double * I1_np1, double * I2_np1, double * I3_np1, double * I4_np1, double * I5_np1, double * V_c_np1, double * I_trns_np1, double * V_c_trns_np1, double I0_n, double I1_n, double I2_n, double I3_n, double I4_n, double I5_n, double V_c_n, double * I_trns_n, double * V_c_trns_n, double XW0, double XP0, double XW1, double XP1, double XW2, double XP2, double X01, double XM, double XT, double YM, double YT, double R_w0_n, double R_w0_np1, double R_w1_n, double R_w1_np1, double R_w2_n, double R_w2_np1, double R_L, double R_01, double R_12, double R_small, double R_s0, double R_s1, double R_s2, double R_p0, double R_p1, double R_p2, double I_b0, double I_b1, double I_b2) {
    // set up matrix and vector for Ax = b calculation
    lapack_int n, nrhs, info;
    n = 27; nrhs = 1;
    unsigned NT = 10;

    double Q_ap = R_small + R_L + XM + X01 + R_12 + R_01;
    double Q_bp = R_small + R_L + XM + X01 + R_12;
    double Q_gp = R_small + R_L + XM;

    double I_alp = I0_n + I1_n;
    double I_bet = I0_n + I1_n + I2_n + I3_n;
    double I_gam = I0_n + I1_n + I2_n + I3_n + I4_n + I5_n;
    double I_b_bet = I_b0 + I_b1;
    double I_b_gam = I_b0 + I_b1 + I_b2;

    double * A = calloc(n*n, sizeof(double));
    double * b = calloc(n*nrhs, sizeof(double));
    lapack_int * ipiv = calloc(n, sizeof(lapack_int));

    A[0] = XW0 + R_w0_np1 + R_s0 + Q_ap + NT*XT;
    A[1] = Q_ap + NT*XT;
    A[2] = Q_bp + NT*XT;
    A[3] = Q_bp + NT*XT;
    A[4] = Q_gp + NT*XT;
    A[5] = Q_gp + NT*XT;
    for (unsigned i=1; i<=NT; i++)
        A[5+i] = R_L + (NT-i)*XT;
    A[16] = -1;

    A[27] = Q_ap + NT*XT;
    A[28] = XP0 + R_p0 + Q_ap + NT*XT;
    A[29] = Q_bp + NT*XT;
    A[30] = Q_bp + NT*XT;
    A[31] = Q_gp + NT*XT;
    A[32] = Q_gp + NT*XT;
    for (unsigned i=1; i<=NT; i++)
        A[32+i] = R_L + (NT-i)*XT;
    A[43] = -1;

    A[54] = Q_bp + NT*XT;
    A[55] = Q_bp + NT*XT;
    A[56] = XW1 + R_w1_np1 + R_s1 + Q_bp + NT*XT;
    A[57] = Q_bp + NT*XT;
    A[58] = Q_gp + NT*XT;
    A[59] = Q_gp + NT*XT;
    for (unsigned i=1; i<=NT; i++)
        A[59+i] = R_L + (NT-i)*XT;
    A[70] = -1;

    A[81] = Q_bp + NT*XT;
    A[82] = Q_bp + NT*XT;
    A[83] = Q_bp + NT*XT;
    A[84] = XP1 + R_p1 + Q_bp + NT*XT;
    A[85] = Q_gp + NT*XT;
    A[86] = Q_gp + NT*XT;
    for (unsigned i=1; i<=NT; i++)
        A[86+i] = R_L + (NT-i)*XT;
    A[97] = -1;

    A[108] = Q_gp + NT*XT;
    A[109] = Q_gp + NT*XT;
    A[110] = Q_gp + NT*XT;
    A[111] = Q_gp + NT*XT;
    A[112] = XW2 + R_w2_np1 + R_s2 + Q_gp + NT*XT;
    A[113] = Q_gp + NT*XT;
    for (unsigned i=1; i<=NT; i++)
        A[113+i] = R_L + (NT-i)*XT;
    A[124] = -1;

    A[135] = Q_gp + NT*XT;
    A[136] = Q_gp + NT*XT;
    A[137] = Q_gp + NT*XT;
    A[138] = Q_gp + NT*XT;
    A[139] = Q_gp + NT*XT;
    A[140] = XP2 + R_p2 + Q_gp + NT*XT;
    for (unsigned i=1; i<=NT; i++)
        A[140+i] = R_L + (NT-i)*XT;
    A[151] = -1;

    for (unsigned i=1; i<=6; i++)
        A[161+i] = R_L + (NT-1)*XT;
    for (unsigned i=1; i<=10; i++)
        A[167+i] = R_L + (NT-i)*XT;
    A[179] = 1;

    for (unsigned i=1; i<=7; i++)
        A[188+i] = R_L + (NT-2)*XT;
    for (unsigned i=2; i<=10; i++)
        A[194+i] = R_L + (NT-i)*XT;
    A[207] = 1;

    for (unsigned i=1; i<=8; i++)
        A[215+i] = R_L + (NT-3)*XT;
    for (unsigned i=3; i<=10; i++)
        A[221+i] = R_L + (NT-i)*XT;
    A[235] = 1;

    for (unsigned i=1; i<=9; i++)
        A[242+i] = R_L + (NT-4)*XT;
    for (unsigned i=4; i<=10; i++)
        A[248+i] = R_L + (NT-i)*XT;
    A[263] = 1;

    for (unsigned i=1; i<=10; i++)
        A[269+i] = R_L + (NT-5)*XT;
    for (unsigned i=5; i<=10; i++)
        A[275+i] = R_L + (NT-i)*XT;
    A[291] = 1;

    for (unsigned i=1; i<=11; i++)
        A[296+i] = R_L + (NT-6)*XT;
    for (unsigned i=6; i<=10; i++)
        A[302+i] = R_L + (NT-i)*XT;
    A[319] = 1;

    for (unsigned i=1; i<=12; i++)
        A[323+i] = R_L + (NT-7)*XT;
    for (unsigned i=7; i<=10; i++)
        A[329+i] = R_L + (NT-i)*XT;
    A[347] = 1;

    for (unsigned i=1; i<=13; i++)
        A[350+i] = R_L + (NT-8)*XT;
    for (unsigned i=8; i<=10; i++)
        A[356+i] = R_L + (NT-i)*XT;
    A[375] = 1;

    for (unsigned i=1; i<=14; i++)
        A[377+i] = R_L + (NT-9)*XT;
    for (unsigned i=9; i<=10; i++)
        A[383+i] = R_L + (NT-i)*XT;
    A[403] = 1;

    for (unsigned i=1; i<=16; i++)
        A[404+i] = R_L;
    A[431] = 1;

    A[432] = YM;
    A[433] = YM;
    A[434] = YM;
    A[435] = YM;
    A[436] = YM;
    A[437] = YM;
    A[448] = 1;

    A[465] = -YT;
    A[476] = 1;

    A[493] = -YT;
    A[504] = 1;

    A[521] = -YT;
    A[532] = 1;

    A[549] = -YT;
    A[560] = 1;

    A[577] = -YT;
    A[588] = 1;

    A[605] = -YT;
    A[616] = 1;

    A[633] = -YT;
    A[644] = 1;

    A[661] = -YT;
    A[672] = 1;

    A[689] = -YT;
    A[700] = 1;

    A[717] = -YT;
    A[728] = 1;

    //for (int i=0; i<27; i++) {
    //    for (int j=0; j<27; j++) {
    //        printf("%+8.2g ", A[27*i+j]);
    //    }
    //    puts("");
    //}

    b[0] = V_c_n + 2*I_b0*R_01 + 2*I_b_bet*R_12 + 2*I_b_gam*(R_L + R_small) - I0_n*(R_w0_n + R_s0 - XW0) - I_alp*R_01 - I_bet*(R_12 - X01) - I_gam*(R_small + R_L - XM) + I_gam*NT*XT + I_trns_n[0]*(9*XT - R_L) + I_trns_n[1]*(8*XT - R_L) + I_trns_n[2]*(7*XT - R_L) + I_trns_n[3]*(6*XT - R_L) + I_trns_n[4]*(5*XT - R_L) + I_trns_n[5]*(4*XT - R_L) + I_trns_n[6]*(3*XT - R_L) + I_trns_n[7]*(2*XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(R_L);
    b[1] = V_c_n + 2*I_b0*R_01 + 2*I_b_bet*R_12 + 2*I_b_gam*(R_L + R_small) - I1_n*(R_p0 - XP0) - I_alp*R_01 - I_bet*(R_12 - X01) - I_gam*(R_small + R_L - XM) + I_gam*NT*XT + I_trns_n[0]*(9*XT - R_L) + I_trns_n[1]*(8*XT - R_L) + I_trns_n[2]*(7*XT - R_L) + I_trns_n[3]*(6*XT - R_L) + I_trns_n[4]*(5*XT - R_L) + I_trns_n[5]*(4*XT - R_L) + I_trns_n[6]*(3*XT - R_L) + I_trns_n[7]*(2*XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(R_L);
    b[2] = V_c_n + 2*I_b_bet*R_12 + 2*I_b_gam*(R_L + R_small) - I2_n*(R_w1_n + R_s1 - XW1) - I_bet*(R_12 - X01) - I_gam*(R_small + R_L - XM) + I_gam*NT*XT + I_trns_n[0]*(9*XT - R_L) + I_trns_n[1]*(8*XT - R_L) + I_trns_n[2]*(7*XT - R_L) + I_trns_n[3]*(6*XT - R_L) + I_trns_n[4]*(5*XT - R_L) + I_trns_n[5]*(4*XT - R_L) + I_trns_n[6]*(3*XT - R_L) + I_trns_n[7]*(2*XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(R_L);
    b[3] = V_c_n + 2*I_b_bet*R_12 + 2*I_b_gam*(R_L + R_small) - I3_n*(R_p1 - XP1) - I_bet*(R_12 - X01) - I_gam*(R_small + R_L - XM) + I_gam*NT*XT + I_trns_n[0]*(9*XT - R_L) + I_trns_n[1]*(8*XT - R_L) + I_trns_n[2]*(7*XT - R_L) + I_trns_n[3]*(6*XT - R_L) + I_trns_n[4]*(5*XT - R_L) + I_trns_n[5]*(4*XT - R_L) + I_trns_n[6]*(3*XT - R_L) + I_trns_n[7]*(2*XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(R_L);
    b[4] = V_c_n + 2*I_b_gam*(R_L + R_small) - I4_n*(R_w2_n + R_s2 - XW2) - I_gam*(R_small + R_L - XM) + I_gam*NT*XT + I_trns_n[0]*(9*XT - R_L) + I_trns_n[1]*(8*XT - R_L) + I_trns_n[2]*(7*XT - R_L) + I_trns_n[3]*(6*XT - R_L) + I_trns_n[4]*(5*XT - R_L) + I_trns_n[5]*(4*XT - R_L) + I_trns_n[6]*(3*XT - R_L) + I_trns_n[7]*(2*XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(R_L);
    b[5] = V_c_n + 2*I_b_gam*(R_L + R_small) - I5_n*(R_p2 - XP2) - I_gam*(R_small + R_L - XM) + I_gam*NT*XT + I_trns_n[0]*(9*XT - R_L) + I_trns_n[1]*(8*XT - R_L) + I_trns_n[2]*(7*XT - R_L) + I_trns_n[3]*(6*XT - R_L) + I_trns_n[4]*(5*XT - R_L) + I_trns_n[5]*(4*XT - R_L) + I_trns_n[6]*(3*XT - R_L) + I_trns_n[7]*(2*XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(R_L);

    b[6] = -V_c_trns_n[0] + 2*I_b_gam*R_L + I_gam*(9*XT - R_L) + I_trns_n[0]*(9*XT - R_L) + I_trns_n[1]*(8*XT - R_L) + I_trns_n[2]*(7*XT - R_L) + I_trns_n[3]*(6*XT - R_L) + I_trns_n[4]*(5*XT - R_L) + I_trns_n[5]*(4*XT - R_L) + I_trns_n[6]*(3*XT - R_L) + I_trns_n[7]*(2*XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(-R_L);
    b[7] = -V_c_trns_n[1] + 2*I_b_gam*R_L + I_gam*(8*XT - R_L) + I_trns_n[0]*(8*XT - R_L) + I_trns_n[1]*(8*XT - R_L) + I_trns_n[2]*(7*XT - R_L) + I_trns_n[3]*(6*XT - R_L) + I_trns_n[4]*(5*XT - R_L) + I_trns_n[5]*(4*XT - R_L) + I_trns_n[6]*(3*XT - R_L) + I_trns_n[7]*(2*XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(-R_L);
    b[8] = -V_c_trns_n[2] + 2*I_b_gam*R_L + I_gam*(7*XT - R_L) + I_trns_n[0]*(7*XT - R_L) + I_trns_n[1]*(7*XT - R_L) + I_trns_n[2]*(7*XT - R_L) + I_trns_n[3]*(6*XT - R_L) + I_trns_n[4]*(5*XT - R_L) + I_trns_n[5]*(4*XT - R_L) + I_trns_n[6]*(3*XT - R_L) + I_trns_n[7]*(2*XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(-R_L);
    b[9] = -V_c_trns_n[3] + 2*I_b_gam*R_L + I_gam*(6*XT - R_L) + I_trns_n[0]*(6*XT - R_L) + I_trns_n[1]*(6*XT - R_L) + I_trns_n[2]*(6*XT - R_L) + I_trns_n[3]*(6*XT - R_L) + I_trns_n[4]*(5*XT - R_L) + I_trns_n[5]*(4*XT - R_L) + I_trns_n[6]*(3*XT - R_L) + I_trns_n[7]*(2*XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(-R_L);
    b[10] = -V_c_trns_n[4] + 2*I_b_gam*R_L + I_gam*(5*XT - R_L) + I_trns_n[0]*(5*XT - R_L) + I_trns_n[1]*(5*XT - R_L) + I_trns_n[2]*(5*XT - R_L) + I_trns_n[3]*(5*XT - R_L) + I_trns_n[4]*(5*XT - R_L) + I_trns_n[5]*(4*XT - R_L) + I_trns_n[6]*(3*XT - R_L) + I_trns_n[7]*(2*XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(-R_L);
    b[11] = -V_c_trns_n[5] + 2*I_b_gam*R_L + I_gam*(4*XT - R_L) + I_trns_n[0]*(4*XT - R_L) + I_trns_n[1]*(4*XT - R_L) + I_trns_n[2]*(4*XT - R_L) + I_trns_n[3]*(4*XT - R_L) + I_trns_n[4]*(4*XT - R_L) + I_trns_n[5]*(4*XT - R_L) + I_trns_n[6]*(3*XT - R_L) + I_trns_n[7]*(2*XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(-R_L);
    b[12] = -V_c_trns_n[6] + 2*I_b_gam*R_L + I_gam*(3*XT - R_L) + I_trns_n[0]*(3*XT - R_L) + I_trns_n[1]*(3*XT - R_L) + I_trns_n[2]*(3*XT - R_L) + I_trns_n[3]*(3*XT - R_L) + I_trns_n[4]*(3*XT - R_L) + I_trns_n[5]*(3*XT - R_L) + I_trns_n[6]*(3*XT - R_L) + I_trns_n[7]*(2*XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(-R_L);
    b[13] = -V_c_trns_n[7] + 2*I_b_gam*R_L + I_gam*(2*XT - R_L) + I_trns_n[0]*(2*XT - R_L) + I_trns_n[1]*(2*XT - R_L) + I_trns_n[2]*(2*XT - R_L) + I_trns_n[3]*(2*XT - R_L) + I_trns_n[4]*(2*XT - R_L) + I_trns_n[5]*(2*XT - R_L) + I_trns_n[6]*(2*XT - R_L) + I_trns_n[7]*(2*XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(-R_L);
    b[14] = -V_c_trns_n[8] + 2*I_b_gam*R_L + I_gam*(XT - R_L) + I_trns_n[0]*(XT - R_L) + I_trns_n[1]*(XT - R_L) + I_trns_n[2]*(XT - R_L) + I_trns_n[3]*(XT - R_L) + I_trns_n[4]*(XT - R_L) + I_trns_n[5]*(XT - R_L) + I_trns_n[6]*(XT - R_L) + I_trns_n[7]*(XT - R_L) + I_trns_n[8]*(XT - R_L) + I_trns_n[9]*(-R_L);
    b[15] = -V_c_trns_n[9] + 2*I_b_gam*R_L + I_gam*(- R_L) + I_trns_n[0]*(- R_L) + I_trns_n[1]*(- R_L) + I_trns_n[2]*(- R_L) + I_trns_n[3]*(- R_L) + I_trns_n[4]*(- R_L) + I_trns_n[5]*(- R_L) + I_trns_n[6]*(- R_L) + I_trns_n[7]*(- R_L) + I_trns_n[8]*(- R_L) + I_trns_n[9]*(-R_L);

    b[16] = V_c_n + YM*(2*I_b_gam - I_gam);

    b[17] = V_c_trns_n[0] + YT*I_trns_n[0];
    b[18] = V_c_trns_n[1] + YT*I_trns_n[1];
    b[19] = V_c_trns_n[2] + YT*I_trns_n[2];
    b[20] = V_c_trns_n[3] + YT*I_trns_n[3];
    b[21] = V_c_trns_n[4] + YT*I_trns_n[4];
    b[22] = V_c_trns_n[5] + YT*I_trns_n[5];
    b[23] = V_c_trns_n[6] + YT*I_trns_n[6];
    b[24] = V_c_trns_n[7] + YT*I_trns_n[7];
    b[25] = V_c_trns_n[8] + YT*I_trns_n[8];
    b[26] = V_c_trns_n[9] + YT*I_trns_n[9];

    //for (int j=0; j<27; j++) {
    //    printf("%4.2e\n", b[j]);
    //}
    //puts("");

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, n, ipiv, b, nrhs);

    if (info != 0) {
        printf("Error encountered in matrix calculation of electrical model (info = %d)...\nExiting with code 3.\n", info);
        exit(3);
    }

    *I0_np1 = b[0];
    *I1_np1 = b[1];
    *I2_np1 = b[2];
    *I3_np1 = b[3];
    *I4_np1 = b[4];
    *I5_np1 = b[5];
    I_trns_np1[0] = b[6];
    I_trns_np1[1] = b[7];
    I_trns_np1[2] = b[8];
    I_trns_np1[3] = b[9];
    I_trns_np1[4] = b[10];
    I_trns_np1[5] = b[11];
    I_trns_np1[6] = b[12];
    I_trns_np1[7] = b[13];
    I_trns_np1[8] = b[14];
    I_trns_np1[9] = b[15];
    *V_c_np1 = b[16];
    V_c_trns_np1[0] = b[17];
    V_c_trns_np1[1] = b[18];
    V_c_trns_np1[2] = b[19];
    V_c_trns_np1[3] = b[20];
    V_c_trns_np1[4] = b[21];
    V_c_trns_np1[5] = b[22];
    V_c_trns_np1[6] = b[23];
    V_c_trns_np1[7] = b[24];
    V_c_trns_np1[8] = b[25];
    V_c_trns_np1[9] = b[26];

    //for (int j=0; j<27; j++) {
    //    printf("%4.2e\n", b[j]);
    //}
    //puts("");

    free(A);
    free(b);
    free(ipiv);

    return 0;
}

// calculates the bias currents for certain target currents
int calculate_bias_from_target_currents_3_wtf_res_par(double * v_I_b0, double * v_I_b1, double * v_I_b2, double I_t0, double I_t1, double I_t2, double R_s0, double R_s1, double R_s2, double R_01, double R_12, double R_p0, double R_p1, double R_p2) {
    // set up matrix and vector for Ax = b calculation
    lapack_int n, nrhs, info;
    n = 3; nrhs = 1;

    double * A = calloc(n*n, sizeof(double));
    double * b = calloc(n*nrhs, sizeof(double));
    lapack_int * ipiv = calloc(n, sizeof(lapack_int));

    double R_v0 = R_p0*R_s0/(R_p0 + R_s0);
    double R_v1 = R_p1*R_s1/(R_p1 + R_s1);
    double R_v2 = R_p2*R_s2/(R_p2 + R_s2);

    double R_z0 = R_01 + ((R_12 + R_v2)*R_v1) / (R_12 + R_v1 + R_v2);
    double R_z1 = (R_12 + R_v2)*(R_01 + R_v0)/(R_12 + R_v2 + R_01 + R_v0);
    double R_z2 = R_12 + ((R_01 + R_v0)*R_v1) / (R_01 + R_v0 + R_v1);

    double I_v0 = I_t0/R_p0*(R_p0 + R_s0);
    double I_v1 = I_t1/R_p1*(R_p1 + R_s1);
    double I_v2 = I_t2/R_p2*(R_p2 + R_s2);

    A[0] = R_z0/(R_v0 + R_z0);
    A[1] = R_v0/(R_v0 + R_01) * (R_z0 - R_01)/(R_z0 + R_v0);
    A[2] = R_v2/(R_z2 + R_v2) * (R_z2 - R_12)/R_z2 * R_v1/(R_v1 + R_v0 + R_01) * R_v0/(R_v0 + R_01);

    A[3] = R_v0/(R_v0 + R_z0) * (R_12 + R_v2)/(R_12 + R_v2 + R_v1);
    A[4] = R_z1/(R_z1 + R_v1);
    A[5] = R_v2/(R_v2 + R_z2) * (R_01 + R_v0)/(R_01 + R_v0 + R_v1);

    A[6] = R_v0/(R_z0 + R_v0) * (R_z0 - R_01)/R_z0 * R_v1/(R_v1 + R_v2 + R_12) * R_v2/(R_v2 + R_12);
    A[7] = R_v2/(R_v2 + R_12) * (R_z2 - R_12)/(R_z2 + R_v2);
    A[8] = R_z2/(R_z2 + R_v2);

    b[0] = I_v0;
    b[1] = I_v1;
    b[2] = I_v2;

    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, A, n, ipiv, b, nrhs);

    if (info != 0) {
        puts("Error encountered in matrix calculation of bias currents...\nExiting with code 4.");
        exit(4);
    }

    *v_I_b0 = b[0];
    *v_I_b1 = b[1];
    *v_I_b2 = b[2];

    free(A);
    free(b);
    free(ipiv);

    return 0;
}

// main function that runs the simulation
int run_three_stage_waterfall_res_par(SimRes * res, SimData * data, double dX0, double dX1, double dX2, double dt, size_t J0, size_t J1, size_t J2, size_t N, size_t NT, size_t NE) {
    // first locally save some important parameters that we will need all the time
    double ** T0 = res->T[0];
    double ** T1 = res->T[1];
    double ** T2 = res->T[2];
    double * I0 = res->I[0];
    double * I1 = res->I[1];
    double * I2 = res->I[2];
    double * I3 = res->I[3];
    double * I4 = res->I[4];
    double * I5 = res->I[5];
    double * R0 = res->R[0];
    double * R1 = res->R[1];
    double * R2 = res->R[2];
    double * V_c = res->V_c[0];

    double * I_trns_prev = calloc(10, sizeof(double));
    double * I_trns_curr = calloc(10, sizeof(double));
    double * V_c_trns_prev = calloc(10, sizeof(double));
    double * V_c_trns_curr = calloc(10, sizeof(double));

    double * Itrns = res->I[6];

    // set up vectors to temporarily save the current and next T
    // this saves a lot of memory for larger simulations
    double * T_stash_1 = calloc(J0, sizeof(double));
    double * T_stash_2 = calloc(J0, sizeof(double));
    double * T_stash_3 = calloc(J1, sizeof(double));
    double * T_stash_4 = calloc(J1, sizeof(double));
    double * T_stash_5 = calloc(J2, sizeof(double));
    double * T_stash_6 = calloc(J2, sizeof(double));
    double * T0_prev = T0[0];
    double * T0_curr = T0[1];
    double * T1_prev = T1[0];
    double * T1_curr = T1[1];
    double * T2_prev = T2[0];
    double * T2_curr = T2[1];

    // set up initial thermal values at t = 1, add a steady state time step at t = 0
    fill_vector(T0_prev, J0, data->T_sub);
    fill_vector(T0_curr, J0, data->T_sub);
    fill_vector(T1_prev, J1, data->T_sub);
    fill_vector(T1_curr, J1, data->T_sub);
    fill_vector(T2_prev, J2, data->T_sub);
    fill_vector(T2_curr, J2, data->T_sub);
    // determine halfway point and set up a beginning hotspot at t = 1
    unsigned halfway = J0/2;
    unsigned initHS_segs = (unsigned) (data->initHS_l_wtf/dX0) + 1;
    // check if there is a nonzero number of segments
    if (initHS_segs < 2) {
        puts("Number of segments in initial hot-spot smaller than 2.\nReturning empty result with error code 2 (wrong initial hot-spot size)...");
        res->exitValue = 2;
        return 2;
    }
    for (unsigned j=halfway - initHS_segs/2; j<halfway + initHS_segs/2; ++j) {
        T0_curr[j] = data->initHS_T_wtf;
    }

    // calculate the correct bias current
    double v_I_b0, v_I_b1, v_I_b2, v_I_t0, v_I_t1, v_I_t2, v_I_t3, v_I_t4, v_I_t5;
    if (fabs(data->I_b0_wtf) < 1e-12 && fabs(data->I_b1_wtf) < 1e-12 && fabs(data->I_b2_wtf) < 1e-12) {
        calculate_bias_from_target_currents_3_wtf_res_par(&v_I_b0, &v_I_b1, &v_I_b2, data->I_t0_wtf, data->I_t1_wtf, data->I_t2_wtf, data->R_s0_wtf, data->R_s1_wtf, data->R_s2_wtf, data->R_01_wtf, data->R_12_wtf, data->R_p0_wtf, data->R_p1_wtf, data->R_p2_wtf);
    } else {
        v_I_b0 = data->I_b0_wtf;
        v_I_b1 = data->I_b1_wtf;
        v_I_b2 = data->I_b2_wtf;
    }
    printf("bias currents %4.2e %4.2e %4.2e\n", v_I_b0, v_I_b1, v_I_b2);
    // determine initial conditions
    double R_v0 = data->R_p0_wtf*data->R_s0_wtf/(data->R_p0_wtf + data->R_s0_wtf);
    double R_v1 = data->R_p1_wtf*data->R_s1_wtf/(data->R_p1_wtf + data->R_s1_wtf);
    double R_v2 = data->R_p2_wtf*data->R_s2_wtf/(data->R_p2_wtf + data->R_s2_wtf);
    double R_01 = data->R_01_wtf;
    double R_12 = data->R_12_wtf;
    double R_z0 = R_01 + ((R_12 + R_v2)*R_v1) / (R_12 + R_v1 + R_v2);
    double R_z1 = (R_12 + R_v2)*(R_01 + R_v0)/(R_12 + R_v2 + R_01 + R_v0);
    double R_z2 = R_12 + ((R_01 + R_v0)*R_v1) / (R_01 + R_v0 + R_v1);
    // add an edge case for when the series resistor is zero
    if (fabs(data->R_s0_wtf) < 1e-6) {
        v_I_t0 = R_v2/R_z2 * (R_z2 - R_12)/R_z2 * R_v1/(R_v1 + R_v0 + R_01) * R_v0/(R_v0 + R_01)*v_I_b2 + R_v0/(R_v0 + R_01) * (R_z0 - R_01)/(R_z0 + R_v1)*v_I_b1 + R_z0/(R_z0 + R_v0)*v_I_b0;
        v_I_t1 = 0;
    } else {
        v_I_t0 = data->R_p0_wtf/(data->R_s0_wtf + data->R_p0_wtf)*(R_v2/R_z2 * (R_z2 - R_12)/R_z2 * R_v1/(R_v1 + R_v0 + R_01) * R_v0/(R_v0 + R_01)*v_I_b2 + R_v0/(R_v0 + R_01) * (R_z0 - R_01)/(R_z0 + R_v1)*v_I_b1 + R_z0/(R_z0 + R_v0)*v_I_b0);
        v_I_t1 = data->R_s0_wtf/(data->R_s0_wtf + data->R_p0_wtf)*(R_v2/R_z2 * (R_z2 - R_12)/R_z2 * R_v1/(R_v1 + R_v0 + R_01) * R_v0/(R_v0 + R_01)*v_I_b2 + R_v0/(R_v0 + R_01) * (R_z0 - R_01)/(R_z0 + R_v1)*v_I_b1 + R_z0/(R_z0 + R_v0)*v_I_b0);
    }
    if (fabs(data->R_s1_wtf) < 1e-6) {
        v_I_t2 = R_v0/(R_v0 + R_z0) * (R_12 + R_v2)/(R_12 + R_v2 + R_v1)*v_I_b0 + R_z1/(R_z1 + R_v1)*v_I_b1 + R_v2/(R_v2 + R_z2) * (R_01 + R_v0)/(R_01 + R_v0 + R_v1)*v_I_b2;
        v_I_t3 = 0;
    } else {
        v_I_t2 = data->R_p1_wtf/(data->R_s1_wtf + data->R_p1_wtf)*(R_v0/(R_v0 + R_z0) * (R_12 + R_v2)/(R_12 + R_v2 + R_v1)*v_I_b0 + R_z1/(R_z1 + R_v1)*v_I_b1 + R_v2/(R_v2 + R_z2) * (R_01 + R_v0)/(R_01 + R_v0 + R_v1)*v_I_b2);
        v_I_t3 = data->R_s1_wtf/(data->R_s1_wtf + data->R_p1_wtf)*(R_v0/(R_v0 + R_z0) * (R_12 + R_v2)/(R_12 + R_v2 + R_v1)*v_I_b0 + R_z1/(R_z1 + R_v1)*v_I_b1 + R_v2/(R_v2 + R_z2) * (R_01 + R_v0)/(R_01 + R_v0 + R_v1)*v_I_b2);
    }
    if (fabs(data->R_s2_wtf) < 1e-6) {
        v_I_t4 = R_v0/(R_z0 + R_v0) * (R_z0 - R_01)/R_z0 * R_v1/(R_v1 + R_v2 + R_12) * R_v2/(R_v2 + R_12)*v_I_b0 + R_v2/(R_v2 + R_12) * (R_z2 - R_12)/(R_z2 + R_v2)*v_I_b1 + R_z2/(R_z2 + R_v2)*v_I_b2;
        v_I_t5 = 0;
    } else {
        v_I_t4 = data->R_p2_wtf/(data->R_s2_wtf + data->R_p2_wtf)*(R_v0/(R_z0 + R_v0) * (R_z0 - R_01)/R_z0 * R_v1/(R_v1 + R_v2 + R_12) * R_v2/(R_v2 + R_12)*v_I_b0 + R_v2/(R_v2 + R_12) * (R_z2 - R_12)/(R_z2 + R_v2)*v_I_b1 + R_z2/(R_z2 + R_v2)*v_I_b2);
        v_I_t5 = data->R_s2_wtf/(data->R_s2_wtf + data->R_p2_wtf)*(R_v0/(R_z0 + R_v0) * (R_z0 - R_01)/R_z0 * R_v1/(R_v1 + R_v2 + R_12) * R_v2/(R_v2 + R_12)*v_I_b0 + R_v2/(R_v2 + R_12) * (R_z2 - R_12)/(R_z2 + R_v2)*v_I_b1 + R_z2/(R_z2 + R_v2)*v_I_b2);
    }
    for (unsigned n=0; n<=data->timeskip; ++n) {
        // set up initial current through the snspd in steady state (t = 0)
        I0[n] = v_I_t0;
        I1[n] = v_I_t1;
        I2[n] = v_I_t2;
        I3[n] = v_I_t3;
        I4[n] = v_I_t4;
        I5[n] = v_I_t5;
        // set up initial voltage drop over the capacitor
        V_c[n] = (v_I_t4+v_I_t5)*R_v2;
    }
    printf("%4.2e ", I0[0]);
    printf("%4.2e ", I1[0]);
    printf("%4.2e ", I2[0]);
    printf("%4.2e ", I3[0]);
    printf("%4.2e ", I4[0]);
    printf("%4.2e\n", I5[0]);
    // put the right currents in the results
    res->I_b[0] = v_I_b0;
    res->I_b[1] = v_I_b1;
    res->I_b[2] = v_I_b2;

    // prepare model parameters for estimating alpha, kappa and c
    // these parameters are considered partially state and temperature dependent
    // formulae for these can be found in Yang
    double DeltaRef = 2.15*Kb*data->T_c*(1 - (data->T_ref_wtf/data->T_c)*(data->T_ref_wtf/data->T_c));
    double A = data->c_e*exp(DeltaRef/(data->T_ref_wtf*Kb));
    double gamma = A/(2.43*data->T_c);
    double B = data->alpha/(pow(data->T_ref_wtf, 3));

    printf("DeltaRef: %e\nA:     %e\ngamma: %e\nB:     %e\n", DeltaRef, A, gamma, B);

    // determine surface ratio between the cross sections of the wires
    double surfaceRatio10 = data->wireThickness_1*data->wireWidth_1/data->wireThickness/data->wireWidth;
    double surfaceRatio21 = data->wireThickness_2*data->wireWidth_2/data->wireThickness/data->wireWidth;
    // define the resistance of a segment of wire in the normal state
    double R_seg0 = data->rho_norm_wtf*dX0/(data->wireWidth*data->wireThickness);
    double R_seg1 = data->rho_norm_wtf/surfaceRatio10*dX1/(data->wireWidth*data->wireThickness);
    double R_seg2 = data->rho_norm_wtf/surfaceRatio21*dX2/(data->wireWidth*data->wireThickness);
    // declare the nanowire resistance and current density
    for (unsigned n=0; n<=data->timeskip; ++n) {
        R0[n] = 0;
        R1[n] = 0;
        R2[n] = 0;
    }
    double currentDensity_w0 = 0;
    double currentDensity_w1 = 0;
    double currentDensity_w2 = 0;

    // allocate space for the state and temperature dependent variables for each time step
    double * alpha0_n = calloc(J0, sizeof(double));
    double * kappa0_n = calloc(J0, sizeof(double));
    double * c0_n = calloc(J0, sizeof(double));
    double * rho_seg0_n = calloc(J0, sizeof(double));
    double * R_seg0_n = calloc(J0, sizeof(double));

    double * alpha1_n = calloc(J1, sizeof(double));
    double * kappa1_n = calloc(J1, sizeof(double));
    double * c1_n = calloc(J1, sizeof(double));
    double * rho_seg1_n = calloc(J1, sizeof(double));
    double * R_seg1_n = calloc(J1, sizeof(double));

    double * alpha2_n = calloc(J2, sizeof(double));
    double * kappa2_n = calloc(J2, sizeof(double));
    double * c2_n = calloc(J2, sizeof(double));
    double * rho_seg2_n = calloc(J2, sizeof(double));
    double * R_seg2_n = calloc(J2, sizeof(double));

    // set up two characteristic numbers for the electrical calculations
    double XW0 = (2*data->L_w0_wtf)/dt;
    double XP0 = (2*data->L_p0_wtf)/dt;
    double XW1 = (2*data->L_w1_wtf)/dt;
    double XP1 = (2*data->L_p1_wtf)/dt;
    double XW2 = (2*data->L_w2_wtf)/dt;
    double XP2 = (2*data->L_p2_wtf)/dt;
    double X01 = (2*data->L_01_wtf)/dt;
    double XM = (2*data->L_m_wtf)/dt;
    double YM = dt/(2*data->C_m_wtf);
    double C_T = 2E-12;
    double L_T = 5E-9;
    double XT = (2*L_T)/dt;
    double YT = dt/(2*C_T);

    // set a flag to check if done
    int flagDone = 0;

    // main time loop
    for (unsigned n=data->timeskip+1; n<NE; ++n) {
        // print progress
        print_progress(n, NE);

        // advance the thermal model to the next time step after the initial step
        if (n > data->timeskip+1 && n < N) {
            if (!data->allowOpt || cmp_vector(T0_prev, J0, data->T_sub, data->T_sub_eps) || cmp_vector(T1_prev, J1, data->T_sub, data->T_sub_eps) || cmp_vector(T2_prev, J2, data->T_sub, data->T_sub_eps)) {
                advance_time_thermal(T0_prev, T0_curr, J0, data->T_sub, alpha0_n, c0_n, rho_seg0_n, kappa0_n, data->wireThickness, currentDensity_w0, dt, dX0);
                advance_time_thermal(T1_prev, T1_curr, J1, data->T_sub, alpha1_n, c1_n, rho_seg1_n, kappa1_n, data->wireThickness_1, currentDensity_w1, dt, dX1);
                advance_time_thermal(T2_prev, T2_curr, J2, data->T_sub, alpha2_n, c2_n, rho_seg2_n, kappa2_n, data->wireThickness_2, currentDensity_w2, dt, dX2);
            } else {
                fill_vector(T0_curr, J0, data->T_sub);
                fill_vector(T1_curr, J1, data->T_sub);
                fill_vector(T2_curr, J2, data->T_sub);
                flagDone = 1;
            }
        }

        if (!flagDone && n < N) {
            // first update the thermal values used in the differential equation,
            //     the targets are included as the first five parameters
            update_thermal_values_3_wtf_res_par(alpha0_n, kappa0_n, c0_n, rho_seg0_n, R_seg0_n, T0_curr, J0, A, B, gamma, data->T_c, I0[n-1], data->I_c0_wtf, data->rho_norm_wtf, data->c_p, data->T_ref_wtf, R_seg0);
            update_thermal_values_3_wtf_res_par(alpha1_n, kappa1_n, c1_n, rho_seg1_n, R_seg1_n, T1_curr, J1, A, B, gamma, data->T_c, I2[n-1], data->I_c1_wtf, data->rho_norm_wtf/surfaceRatio10, data->c_p, data->T_ref_wtf, R_seg1);
            update_thermal_values_3_wtf_res_par(alpha2_n, kappa2_n, c2_n, rho_seg2_n, R_seg2_n, T2_curr, J2, A, B, gamma, data->T_c, I4[n-1], data->I_c2_wtf, data->rho_norm_wtf/surfaceRatio21, data->c_p, data->T_ref_wtf, R_seg2);
            // update the current nanowire resistance
            R0[n] = sum_vector(R_seg0_n, J0);
            R1[n] = sum_vector(R_seg1_n, J1);
            R2[n] = sum_vector(R_seg2_n, J2);
        } else {
            R0[n] = 0;
            R1[n] = 0;
            R2[n] = 0;
        }

        // update the current density through the nanowire
        currentDensity_w0 = I0[n-1]/(data->wireWidth*data->wireThickness);
        currentDensity_w1 = I2[n-1]/(data->wireWidth_1*data->wireThickness_1);
        currentDensity_w2 = I4[n-1]/(data->wireWidth_2*data->wireThickness_2);
        // update the electric values
        advance_time_electric_3_wtf_res_par(&I0[n], &I1[n], &I2[n], &I3[n], &I4[n], &I5[n], &V_c[n], I_trns_curr, V_c_trns_curr, I0[n-1], I1[n-1], I2[n-1], I3[n-1], I4[n-1], I5[n-1], V_c[n-1], I_trns_prev, V_c_trns_prev, XW0, XP0, XW1, XP1, XW2, XP2, X01, XM, XT, YM, YT, R0[n-1], R0[n], R1[n-1], R1[n], R2[n-1], R2[n], data->R_L_wtf, data->R_01_wtf, data->R_12_wtf, data->R_small_wtf, data->R_s0_wtf, data->R_s1_wtf, data->R_s2_wtf, data->R_p0_wtf, data->R_p1_wtf, data->R_p2_wtf, v_I_b0, v_I_b1, v_I_b2);

        // shuffle the T pointers around so the old and new timestep don't point to the same array
        T0_prev = T0_curr;
        T1_prev = T1_curr;
        T2_prev = T2_curr;
        if (n % data->timeskip == 0 && n < N) {
            T0_curr = T0[n/data->timeskip];
            T1_curr = T1[n/data->timeskip];
            T2_curr = T2[n/data->timeskip];
        } else {
            if (T0_prev == T_stash_1)
                T0_curr = T_stash_2;
            else
                T0_curr = T_stash_1;

            if (T1_prev == T_stash_3)
                T1_curr = T_stash_4;
            else
                T1_curr = T_stash_3;

            if (T2_prev == T_stash_5)
                T2_curr = T_stash_6;
            else
                T2_curr = T_stash_5;
        }

        Itrns[n] = sum_vector(I_trns_curr, 10);

        // update transmission currents and voltages
        double * tmp = I_trns_prev;
        I_trns_prev = I_trns_curr;
        I_trns_curr = tmp;
        tmp = V_c_trns_prev;
        V_c_trns_prev = V_c_trns_curr;
        V_c_trns_curr = tmp;
    }

    // free allocated space
    free(T_stash_1);
    free(T_stash_2);
    free(T_stash_3);
    free(T_stash_4);
    free(T_stash_5);
    free(T_stash_6);

    free(alpha0_n);
    free(kappa0_n);
    free(c0_n);
    free(rho_seg0_n);
    free(R_seg0_n);

    free(alpha1_n);
    free(kappa1_n);
    free(c1_n);
    free(rho_seg1_n);
    free(R_seg1_n);

    free(alpha2_n);
    free(kappa2_n);
    free(c2_n);
    free(rho_seg2_n);
    free(R_seg2_n);

    free(I_trns_prev);
    free(I_trns_curr);
    free(V_c_trns_prev);
    free(V_c_trns_curr);

    // print result
    puts("\nSimulation completed.");
    res->exitValue = 0;
    return 0;
}
