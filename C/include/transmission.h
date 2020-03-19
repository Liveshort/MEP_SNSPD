#ifndef __TRANSMISSION_H__
#define __TRANSMISSION_H__

// some global constants
const double c0;            // speed of light

int sim_transmission_line(SimData * data, SimRes * res, size_t NE, size_t NTL, size_t numberOfIb);

#endif
