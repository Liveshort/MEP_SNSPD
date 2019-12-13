#ifndef __TYPES_H__
#define __TYPES_H__

// structure that contains the result of the simulation of an snspd
typedef struct _simres {
    size_t J;               // number of spatial samples
    size_t N;               // number of time samples
    double ** T;            // temperature matrix T[time][space]
    double ** I;            // current vectors I[time][number]
    double ** R;            // resistance vectors R[time][number]
    int exitValue;          // 0 for succes, anything else for error
} SimRes;

// structure that contains the (input) data of the simulation of an snspd
typedef struct _simdata {
    // general info
    size_t J;               // number of spatial samples
    size_t N;               // number of time samples
    size_t numberOfI;       // number of currents that need to be measured
    size_t numberOfR;       // number of resistances that need to be measured
    // physical dimensions
    double wireLength;      // length of the nanowire (divided into J segments)
    double wireThickness;   // thickness of the nanowire
    double wireWidth;       // width of the nanowire
    double tMax;            // maximum time, sim will run from 0 to tMax
    // experiment specific data
    double T_c;             // critical temperature
    double I_c0;            // critical current at 0K
    double c_p;             // phonon specific heat
    double c_e;             // electron specific heat
    double alpha;           // thermal boundary conductivity
    double T_sub;           // substrate temperature
    // data specific to the standard model (runtype 0)
    double R_L_std;         // load resistor
    double C_m_std;         // dc port bias tee
    double I_b_std;         // bias current
    double initHS_l_std;    // initial hotspot length (to simulate a photon hit)
    double initHS_T_std;    // initial hotspot temperature
    double rho_norm_std;    // conductivity of the nanowire in normal state
    double L_w_std;         // kinetic inductance of the nanowire
    double T_ref_std;       // reference temperature for model parameters [Yang]
} SimData;

void free_simres(SimRes * simRes);
void free_simdata(SimData * simData);

#endif
