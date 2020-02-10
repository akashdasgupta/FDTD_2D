#ifndef SOFTSOURCE
#define SOFTSOURCE
#include <math.h>  
#include <cmath>
//Creates a gausian pulse pertubation in the electric field
void inject_soft_source(double Ez[], int iterations, double dt, double sigma, double t0);
// Creates a contimous sin wave
void inject_soft_source2(double Ez[], int iterations, double dt, double sigma, double t0);
#endif
