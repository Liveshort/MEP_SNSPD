#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

#include "types.h"
#include "snspd.h"
#include "helper.h"

#include "linalg.h"

int main(int argc, char * argv[]) {
    // check and open setup file
    if (argc != 3) {
        puts("Wrong number of arguments...\nUsage: ./simulation inputFile outputFolder\nExiting...");
        exit(4);
    }

    printf("Input file:    %s\n", argv[1]);
    printf("Output folder: %s\n", argv[2]);

    FILE * fp;

    if ((fp = fopen(argv[1], "r")) == NULL) {
        printf("Something went wrong trying to open the file \"%s\"...\nExiting...\n", argv[1]);
        exit(5);
    }

    // create throwaway variable for text from the setup file we don't use
    char dump[500];

    // put in the experiment data
    SimData * data = calloc(1, sizeof(SimData));
    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);
    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->J, dump) < 1) exit(6);
    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->N, dump) < 1) exit(6);
    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->numberOfT, dump) < 1) exit(6);
    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->numberOfI, dump) < 1) exit(6);
    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->numberOfR, dump) < 1) exit(6);
    if (fscanf(fp, "%zu;%2000[^\n]\n", &data->timeskip, dump) < 1) exit(6);

    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireLength, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireThickness, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->wireWidth, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->tMax, dump) < 1) exit(6);

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
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->C_m_std, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->I_b_std, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->initHS_l_std, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->initHS_T_std, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->rho_norm_std, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_w_std, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->T_ref_std, dump) < 1) exit(6);

    if (fscanf(fp, "%2000[^\n]\n", dump) < 1) exit(6);

    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->L_p_parallel, dump) < 1) exit(6);
    if (fscanf(fp, "%lf;%2000[^\n]\n", &data->R_p_parallel, dump) < 1) exit(6);

    // run simulation
    SimRes * res = run_snspd_simulation(data, 0);

    // generate correct filenames
    const char Tbin[] = "T.bin";
    const char Ibin[] = "I.bin";
    const char Rbin[] = "R.bin";
    const char paramInfo[] = "param.info";
    char *TFilename = malloc(strlen(argv[2]) + strlen(Tbin) + 1);
    char *IFilename = malloc(strlen(argv[2]) + strlen(Ibin) + 1);
    char *RFilename = malloc(strlen(argv[2]) + strlen(Rbin) + 1);
    char *paramInfoFilename = malloc(strlen(argv[2]) + strlen(paramInfo) + 1);
    snprintf(TFilename, strlen(argv[2]) + strlen(Tbin) + 1, "%s%s", argv[2], Tbin);
    snprintf(IFilename, strlen(argv[2]) + strlen(Tbin) + 1, "%s%s", argv[2], Ibin);
    snprintf(RFilename, strlen(argv[2]) + strlen(Tbin) + 1, "%s%s", argv[2], Rbin);
    snprintf(paramInfoFilename, strlen(argv[2]) + strlen(paramInfo) + 1, "%s%s", argv[2], paramInfo);

    // write data to binary files
    fp = fopen(TFilename, "wb");
    for (unsigned n=0; n<res->N; n += res->timeskip)
        fwrite(res->T[0][n], sizeof(double), res->J, fp);
    fclose(fp);

    fp = fopen(IFilename, "wb");
    for (unsigned n=0; n<res->N; n += res->timeskip)
        fwrite(&res->I[0][n], sizeof(double), 1, fp);
    fclose(fp);

    fp = fopen(RFilename, "wb");
    for (unsigned n=0; n<res->N; n += res->timeskip)
        fwrite(&res->R[0][n], sizeof(double), 1, fp);
    fclose(fp);

    // write simulation parameters to file for readout
    fp = fopen(paramInfoFilename, "w");
    fprintf(fp, "%40s; %zu\n", "J (# of spatial elements)", res->J);
    fprintf(fp, "%40s; %zu\n", "N (# of temporal elements)", res->N);
    fprintf(fp, "%40s; %zu\n", "timeskip factor", res->timeskip);
    fprintf(fp, "%40s; %zu\n", "# of nanowires", res->numberOfT);
    fprintf(fp, "%40s; %zu\n", "# of currents", res->numberOfI);
    fprintf(fp, "%40s; %zu\n", "# of resistances", res->numberOfR);
    fprintf(fp, "%40s; ", "dX [m]");
    for (unsigned i=0; i<res->numberOfT; ++i) {
        if (i > 0) fprintf(fp, "; ");
        fprintf(fp, "%8.6e\n", res->dX[i]);
    }
    fprintf(fp, "%40s; %8.6e\n", "dt [s]", res->dt);
    fclose(fp);

    free_simres(res);
    free_simdata(data);

    //double * diag = calloc(5, sizeof(double));
    //double * offdiag = calloc(5, sizeof(double));
    //double * rhs = calloc(5, sizeof(double));
//
    //diag[0] = 3;
    //diag[1] = 3;
    //diag[2] = 6;
    //diag[3] = 5;
    //diag[4] = 3;
    ////diag[5] = 3;
//
    //offdiag[0] = 0.8;
    //offdiag[1] = 0.5;
    //offdiag[2] = 0.2;
    //offdiag[3] = 0.4;
    //offdiag[4] = 0.5;
    ////offdiag[5] = 1;
//
    //rhs[0] = 10;
    //rhs[1] = 4;
    //rhs[2] = 3;
    //rhs[3] = 5;
    //rhs[4] = 8;
    ////rhs[5] = 6;
//
    //double * res = calloc(5, sizeof(double));
//
    //TDM_solve(res, 5, diag, offdiag, rhs);
    //print_vector(res, 5);

    exit(0);
}
