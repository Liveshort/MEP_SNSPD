#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#include <sys/sysinfo.h>

#include "types.h"
#include "snspd.h"
#include "helper.h"
#include "linalg.h"

int collect_data(FILE * fp, SimData * data) {
    // create throwaway variable for text from the setup file we don't use
    char dump[500];

    // scan runtype
    if (fscanf(fp, "%d;%2000[^\n]\n", &data->runType, dump) < 1) exit(6);

    // set up number of vectors needed based on the runtype
    if (data->runType == 0) {
        data->numberOfT = 1;
        data->numberOfI = 4;
        data->numberOfR = 1;
        data->numberOfC = 1;
    } else if (data->runType == 1) {
        data->numberOfT = 1;
        data->numberOfI = 5;
        data->numberOfR = 1;
        data->numberOfC = 1;
    } else if (data->runType == 2) {
        data->numberOfT = 2;
        data->numberOfI = 5;
        data->numberOfR = 2;
        data->numberOfC = 1;
    } else if (data->runType == 3) {
        data->numberOfT = 2;
        data->numberOfI = 5;
        data->numberOfR = 2;
        data->numberOfC = 2;
    } else if (data->runType == 4) {
        data->numberOfT = 2;
        data->numberOfI = 7;
        data->numberOfR = 2;
        data->numberOfC = 1;
    } else if (data->runType == 5) {
        data->numberOfT = 2;
        data->numberOfI = 7;
        data->numberOfR = 2;
        data->numberOfC = 2;
    } else if (data->runType == 6) {
        data->numberOfT = 3;
        data->numberOfI = 9;
        data->numberOfR = 3;
        data->numberOfC = 1;
    } else if (data->runType == 7) {
        data->numberOfT = 3;
        data->numberOfI = 9;
        data->numberOfR = 3;
        data->numberOfC = 3;
    } else if (data->runType == 10) {
        data->numberOfT = 12;
        data->numberOfI = 27;
        data->numberOfR = 12;
        data->numberOfC = 4;
    }

    // scan data runtype 0 and 1
    if (data->runType == 0 || data->runType == 1) {
        if (fscanf(fp, "%d;%2000[^\n]\n", &data->allowOpt, dump) < 1) exit(6);
        if (fscanf(fp, "%d;%2000[^\n]\n", &data->simTL, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->J0, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->N, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->tMax, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->timeskip, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->ETratio, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->impurityOffset, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->impuritySpread, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->NTL, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->VF, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->LTL, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireLength, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireThickness, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireWidth, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_c, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_c0, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->c_p, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->c_e, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->alpha, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_sub, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_sub_eps, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_L_std, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_s_std, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->C_m_std, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_b_std, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_t_std, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->initHS_l_std, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->initHS_T_std, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->rho_norm_std, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_w_std, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_ref_std, dump) < 1) exit(6);
    }

    if (data->runType == 1) {
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_p_parallel, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_p_parallel, dump) < 1) exit(6);
    }

    // scan data runtype 2-5
    if (data->runType >= 2 && data->runType <= 5) {
        if (fscanf(fp, "%d;%2000[^\n]\n", &data->allowOpt, dump) < 1) exit(6);
        if (fscanf(fp, "%d;%2000[^\n]\n", &data->simTL, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->J0, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->J1, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->N, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->tMax, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->timeskip, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->ETratio, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->impurityOffset, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->impuritySpread, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->NTL, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->VF, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->LTL, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireLength, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireThickness, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireWidth, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireLength_1, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireThickness_1, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireWidth_1, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_c, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->c_p, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->c_e, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->alpha, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_sub, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_sub_eps, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_L_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_s0_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_s1_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_small_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_01_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->C_01_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->C_m_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_b0_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_b1_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_t0_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_t1_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->initHS_l_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->initHS_T_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->rho_norm_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_w0_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_w1_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_ref_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_c0_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_c1_wtf, dump) < 1) exit(6);
    }

    if (data->runType == 4 || data->runType == 5) {
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_p0_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_p1_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_p0_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_p1_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_m_wtf, dump) < 1) exit(6);
    }

    // scan data runtype 6-7
    if (data->runType >= 6 && data->runType <= 7 || data->runType == 10) {
        if (fscanf(fp, "%d;%2000[^\n]\n", &data->allowOpt, dump) < 1) exit(6);
        if (fscanf(fp, "%d;%2000[^\n]\n", &data->simTL, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->J0, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->J1, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->J2, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->N, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->tMax, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->timeskip, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->ETratio, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->impurityOffset, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->impuritySpread, dump) < 1) exit(6);
        if (fscanf(fp, "%zu;%2000[^\n]\n", &data->NTL, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->VF, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->LTL, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireLength, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireThickness, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireWidth, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireLength_1, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireThickness_1, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireWidth_1, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireLength_2, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireThickness_2, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireWidth_2, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_c, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->c_p, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->c_e, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->alpha, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_sub, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_sub_eps, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_L_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_s0_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_s1_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_s2_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_small_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_01_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_12_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->C_01_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->C_12_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->C_m_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_b0_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_b1_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_b2_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_t0_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_t1_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_t2_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->initHS_l_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->initHS_T_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->rho_norm_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_w0_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_w1_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_w2_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_ref_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_c0_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_c1_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_c2_wtf, dump) < 1) exit(6);

        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_p0_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_p1_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_p2_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_p0_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_p1_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_p2_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_m_wtf, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_12_wtf, dump) < 1) exit(6);
    }

    if (data->runType == 10) {
        if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_k0_rsm, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_k1_rsm, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_k2_rsm, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_k3_rsm, dump) < 1) exit(6);
        if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_Lp_rsm, dump) < 1) exit(6);
    }

    fclose(fp);

    return 0;
}

int write_results(char * outputPath, FILE * fp, SimRes * res) {
    // generate correct filenames
    const char Tbin[] = "T.bin";
    const char Ibin[] = "I.bin";
    const char Rbin[] = "R.bin";
    const char Cbin[] = "V_c.bin";
    const char paramInfo[] = "param.info";
    char *TFilename = calloc(strlen(outputPath) + strlen(Tbin) + 1, sizeof(char));
    char *IFilename = calloc(strlen(outputPath) + strlen(Ibin) + 1, sizeof(char));
    char *RFilename = calloc(strlen(outputPath) + strlen(Rbin) + 1, sizeof(char));
    char *CFilename = calloc(strlen(outputPath) + strlen(Cbin) + 1, sizeof(char));
    char *paramInfoFilename = calloc(strlen(outputPath) + strlen(paramInfo) + 1, sizeof(char));
    snprintf(TFilename, strlen(outputPath) + strlen(Tbin) + 1, "%s%s", outputPath, Tbin);
    snprintf(IFilename, strlen(outputPath) + strlen(Ibin) + 1, "%s%s", outputPath, Ibin);
    snprintf(RFilename, strlen(outputPath) + strlen(Rbin) + 1, "%s%s", outputPath, Rbin);
    snprintf(CFilename, strlen(outputPath) + strlen(Cbin) + 1, "%s%s", outputPath, Cbin);
    snprintf(paramInfoFilename, strlen(outputPath) + strlen(paramInfo) + 1, "%s%s", outputPath, paramInfo);

    // write data to binary files
    fp = fopen(TFilename, "wb");
    for (unsigned q=0; q<res->numberOfT; ++q) {
        for (unsigned n=0; n<res->N/res->timeskip; ++n)
        fwrite(res->T[q][n], sizeof(double), res->J[q], fp);
    }
    fclose(fp);

    fp = fopen(IFilename, "wb");
    for (unsigned q=0; q<res->numberOfI; ++q) {
        for (unsigned n=0; n<res->N*res->ETratio; n += res->timeskip/10)
            fwrite(&res->I[q][n], sizeof(double), 1, fp);
    }
    fclose(fp);

    fp = fopen(RFilename, "wb");
    for (unsigned q=0; q<res->numberOfR; ++q) {
        for (unsigned n=0; n<res->N*res->ETratio; n += res->timeskip/10)
            fwrite(&res->R[q][n], sizeof(double), 1, fp);
    }
    fclose(fp);

    fp = fopen(CFilename, "wb");
    for (unsigned q=0; q<res->numberOfC; ++q) {
        for (unsigned n=0; n<res->N*res->ETratio; n += res->timeskip/10)
            fwrite(&res->V_c[q][n], sizeof(double), 1, fp);
    }
    fclose(fp);

    // write simulation parameters to file for readout
    fp = fopen(paramInfoFilename, "w");
    fprintf(fp, "%40s; %d\n", "runtype of the simulation", res->runType);
    fprintf(fp, "%40s; %d\n", "transmission simulation type", res->tlType);
    fprintf(fp, "%40s; ", "J (# of spatial elements)");
    for (unsigned i=0; i<res->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%zu", res->J[i]);
    }
    fprintf(fp, "\n%40s; %zu\n", "N (# of temporal elements)", res->N);
    fprintf(fp, "%40s; %zu\n", "timeskip factor", res->timeskip);
    fprintf(fp, "%40s; %zu\n", "electrical / thermal time ratio", res->ETratio);
    fprintf(fp, "%40s; %zu\n", "# of nanowires", res->numberOfT);
    fprintf(fp, "%40s; %zu\n", "# of currents", res->numberOfI);
    fprintf(fp, "%40s; %zu\n", "# of resistances", res->numberOfR);
    fprintf(fp, "%40s; %zu\n", "# of capacitor voltages", res->numberOfC);
    fprintf(fp, "%40s; ", "I_b [A]");
    for (unsigned i=0; i<res->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%8.6e", res->I_b[i]);
    }
    fprintf(fp, "\n%40s; ", "dX [m]");
    for (unsigned i=0; i<res->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%8.6e", res->dX[i]);
    }
    fprintf(fp, "\n%40s; %8.6e\n", "dt [s]", res->dt);
    fclose(fp);

    free(TFilename);
    free(IFilename);
    free(RFilename);
    free(CFilename);
    free(paramInfoFilename);
}

int main(int argc, char * argv[]) {
    puts("########################################################################################\n################################### SNSPD SIMULATION ###################################\n########################################################################################\n");

    // check and open setup file
    if (argc != 3) {
        puts("Wrong number of arguments...\nUsage: ./simulation inputFile outputFolder\nExiting...");
        exit(4);
    }

    printf("    %d processors detected. %d processors available.\n", get_nprocs_conf(), get_nprocs());
    printf("    OMP_NUM_THREADS set to %s, using that amount of threads.\n\n", getenv("OMP_NUM_THREADS"));

    printf("    Input file:       %s\n", argv[1]);
    printf("    Output folder:    %s\n\n", argv[2]);

    FILE * fp;

    if ((fp = fopen(argv[1], "r")) == NULL) {
        printf("Something went wrong trying to open the file \"%s\"...\nExiting...\n", argv[1]);
        exit(5);
    }

    // put in the experiment data
    SimData * data = calloc(1, sizeof(SimData));
    collect_data(fp, data);

    // check if some of the parameters make sense
    if (data->N % data->timeskip != 0 || data->N <= 0) {
        printf("Error: invalid N (%ld), not larger than 0 or N is not an integer multiple of timeskip (%ld)...\nExiting with code 8.\n", data->N, data->timeskip);
        exit(8);
    }

    if (data->timeskip % 10 != 0 || data->timeskip > 100000) {
        printf("Error: timeskip is not a multiple of 10...\nExiting with code 9.\n");
        exit(9);
    }

    // run simulation
    SimRes * res = run_snspd_simulation(data, data->runType);

    puts("    Writing files...");
    write_results(argv[2], fp, res);
    puts("    Done.\n");

    free_simres(res);
    free_simdata(data);

    puts("########################################################################################\n#################################### END SIMULATION ####################################\n########################################################################################");

    exit(0);
}
