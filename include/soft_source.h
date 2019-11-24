#ifndef SOFTSOURCE
#define SOFTSOURCE

#include <math.h>  
# include <cmath>

void inject_soft_source(double Ez[], double Hx[], double Hy[], int iterations, double dx, double dy, double dt, double nsrc, double ersrc, double muxsrc, double muysrc, double sigma, double t0);

#endif
