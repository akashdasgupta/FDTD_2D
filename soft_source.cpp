#include <soft_source.h>
#include <math.h>  
# include <cmath>
#include <iostream>

/*
 * Creates a gausian pulse pertubation in the electric field
 * @param [out] Ez array holding Electric field pertubations
 * @param iterations total simulation iterations
 * @param dt timestep of simulation
 * @param sigma the width of the gausian pulse
 * @param t0 time delay before source injected
 */
void inject_soft_source(double Ez[], int iterations, double dt, double sigma, double t0)
{
    double t [iterations];
    for(int i=0; i<iterations; ++i)
    {
        t[i] = i * dt; // array of times 
    }
    
    for (int i=0; i<iterations; ++i)
    {
        Ez[i] = std::exp(-1 * ((t[i] - t0)/ sigma) * ((t[i] - t0)/ sigma)); 
    }
}

/*
 * Creates a contimous sin wave
 * @param [out] Ez array holding Electric field pertubations
 * @param iterations total simulation iterations
 * @param dt timestep of simulation
 * @param sigma strength of graded injection
 * @param t0 time delay before source injected
 */
void inject_soft_source2(double Ez[], int iterations, double dt, double sigma, double t0)
{
    double c{299792458};

    double t [iterations];
    for(int i=0; i<iterations; ++i)
    {
        t[i] = i * dt; // array of times
    }

    for (int i=0; i<iterations; ++i)
    {
        // tanh term ensures the source injection is gradual
        double wavelength {200e-9} // can change if needed
        Ez[i] = (0.5 * std::tanh(((t[i] - t0)/ sigma)) + 0.5) * std::sin(2 * 3.14159265359* c *t[i] /wavelength);
    }
}






