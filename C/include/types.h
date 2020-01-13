#ifndef __TYPES_H__
#define __TYPES_H__

// some global constants
const double Kb;            // boltzmann constant
const double Lorentz;       // Lorentz number

// structure that contains the result of the simulation of an snspd
typedef struct _simres {
    int runType;            // runtype of the simulation (0: yang, 1: yang_parallel)
    size_t J;               // number of spatial samples
    size_t N;               // number of time samples
    double * dX;            // delta X
    double dt;              // delta t
    size_t timeskip;        // factor to reduce timesteps that are used in the calculation, but not shown in the result
    size_t ETratio;         // ratio between time calculated for electrical model and thermal model
    size_t numberOfT;       // number of temperature vectors that need to be measured
    size_t numberOfI;       // number of currents that need to be measured
    size_t numberOfR;       // number of resistances that need to be measured
    size_t numberOfC;       // number of capacitors that need to be measured
    double *** T;           // temperature matrices T[number][time][space]
    double ** I;            // current vectors I[number][time]
    double ** R;            // resistance vectors R[number][time]
    double ** V_c;          // voltage vectors for capacitors
    int exitValue;          // 0 for succes, anything else for error
} SimRes;

// structure that contains the (input) data of the simulation of an snspd
typedef struct _simdata {
    // runtype
    int runType;            // runtype of the simulation (0: yang, 1: yang parallel)
    int allowOpt;           // if 1, optimizes thermal when wire has cooled to within T_sub_eps
    // general info
    size_t J;               // number of spatial samples
    size_t N;               // number of time samples
    double tMax;            // maximum time, sim will run from 0 to tMax
    size_t timeskip;        // factor to reduce timesteps that are used in the calculation, but not shown in the result
    size_t ETratio;         // ratio between time calculated for electrical model and thermal model
    size_t numberOfT;       // number of temperature vectors that need to be measured
    size_t numberOfI;       // number of currents that need to be measured
    size_t numberOfR;       // number of resistances that need to be measured
    size_t numberOfC;       // number of capacitors that need to be measured
    // physical dimensions detector wire
    double wireLength;      // length of the nanowire (divided into J segments)
    double wireThickness;   // thickness of the nanowire
    double wireWidth;       // width of the nanowire
    // experiment specific data
    double T_c;             // critical temperature
    double I_c0;            // critical current at 0K
    double c_p;             // phonon specific heat
    double c_e;             // electron specific heat
    double alpha;           // thermal boundary conductivity
    double T_sub;           // substrate temperature
    double T_sub_eps;       // sub temp epsilon, optimization strategy (detect steady state)
    // data specific to the standard model (runtype 0)
    double R_L_std;         // load resistor
    double R_s_std;         // series resistor
    double C_m_std;         // dc port bias tee
    double I_b_std;         // bias current
    double initHS_l_std;    // initial hotspot length (to simulate a photon hit)
    double initHS_T_std;    // initial hotspot temperature
    double rho_norm_std;    // conductivity of the nanowire in normal state
    double L_w_std;         // kinetic inductance of the nanowire
    double T_ref_std;       // reference temperature for model parameters [Yang]
    // data specific to the standard model with parallel filter (runtype 1)
    double L_p_parallel;    // kinetic inductance of the parallel inductor
    double R_p_parallel;    // resistance of the parallel resistor
} SimData;

void free_simres(SimRes * simRes);
void free_simdata(SimData * simData);

#endif
