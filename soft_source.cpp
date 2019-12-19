#include <soft_source.h>
#include <math.h>  
# include <cmath>
#include <iostream>

void inject_soft_source(double Ez[], double Hx[], double Hy[], int iterations, double dx, double dy, double dt, double nsrc, double ersrc, double muxsrc, double muysrc, double sigma, double t0)
{
    double c{299792458};

    double t [iterations];
    for(int i=0; i<iterations; ++i)
    {
        t[i] = i * dt;
        //std::cout << t[i] << std::endl;
    }
    
    double deltahx {(nsrc * dy / (2*c)) + dt/2};
    double Hx_over_E {-1 * sqrt(ersrc / muxsrc)};
    
    double deltahy {(nsrc * dx / (2*c)) + dt/2};
    double Hy_over_E {-1 * sqrt(ersrc / muysrc)};
    
    for (int i=0; i<iterations; ++i)
    {
        //std::cout<< std::exp(-1 * ((t[i] - t0)/ sigma) * ((t[i] - t0)/ sigma)) << std::endl;
        Ez[i] = std::exp(-1 * ((t[i] - t0)/ sigma) * ((t[i] - t0)/ sigma));
//         Hx[i] = Hx_over_E * std::exp(-1 * ((t[i] - t0+ deltahx) / sigma) * ((t[i] - t0+ deltahx)/ sigma));
//         Hy[i] = Hy_over_E * std::exp(-1 * ((t[i] - t0+ deltahy) / sigma) * ((t[i] - t0+ deltahy)/ sigma));
        
    }
}



void inject_soft_source2(double Ez[], double Hx[], double Hy[], int iterations, double dx, double dy, double dt, double nsrc, double ersrc, double muxsrc, double muysrc, double sigma, double t0)
{
    double c{299792458};

    double t [iterations];
    for(int i=0; i<iterations; ++i)
    {
        t[i] = i * dt;
        //std::cout << t[i] << std::endl;
    }
    
    double deltahx {(nsrc * dy / (2*c)) + dt/2};
    double Hx_over_E {-1 * sqrt(ersrc / muxsrc)};
    
    double deltahy {(nsrc * dx / (2*c)) + dt/2};
    double Hy_over_E {-1 * sqrt(ersrc / muysrc)};
    
    for (int i=0; i<iterations; ++i)
    {
        //std::cout<< std::exp(-1 * ((t[i] - t0)/ sigma) * ((t[i] - t0)/ sigma)) << std::endl;
        Ez[i] = (0.5 * std::tanh(((t[i] - t0)/ sigma)) + 0.5) * std::sin(2 * 3.14159265359* c *t[i] / 520e-9);
//         Hx[i] = Hx_over_E * std::exp(-1 * ((t[i] - t0+ deltahx) / sigma) * ((t[i] - t0+ deltahx)/ sigma));
//         Hy[i] = Hy_over_E * std::exp(-1 * ((t[i] - t0+ deltahy) / sigma) * ((t[i] - t0+ deltahy)/ sigma));
        
    }
}






