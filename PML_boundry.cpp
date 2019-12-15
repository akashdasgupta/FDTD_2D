#include <PML_boundry.h>
#include <iostream>
#include <cmath>


double epsilon_naught{8.8541878128e-12};
double c{299792458};


void sigma_setter(double sigma_x[], double sigma_y[], int pml_size_x, int pml_size_y, double dt)
{
    for(int i=0; i<pml_size_x; ++i)
    {
        double factor{static_cast<double>(i+1) / static_cast<double>(pml_size_x)};
        //std::cout << factor << std::endl;

        sigma_x[i] = epsilon_naught * (factor*factor* factor) / (2 * dt);
        //std::cout << sigma_x[i] << std::endl;
    }
    
    for(int i=0; i<pml_size_y; ++i)
    {
        double factor{static_cast<double>(i / pml_size_y)};
        sigma_y[i] = epsilon_naught * (factor*factor* factor) / (2 * dt);
    }
}

void param_setter_x(PML_coefs Hx[], PML_coefs Hy[], PML_coefs Dz[], double sigma_x[], double sigma_y[], double mu_x, double mu_y, double dt, int size)
{
    for (int i=0; i<size; ++i)
    {
        Hx[i].init_coefs_Hx(sigma_x[i], 0, 1, dt);
        Hy[i].init_coefs_Hy(sigma_x[i], 0, 1, dt);
        Dz[i].init_coefs_Dz(sigma_x[i], 0, dt);
    }
}


void PML_coefs::init_coefs_Hx(double sigma_x, double sigma_y, double mu_x, double dt)
{
    double m0 {(1/dt) + (sigma_y / (2 * epsilon_naught))};

    m_m1 = ((1/dt) - (sigma_y / (2 * epsilon_naught))) / m0;
    m_m2 = -c / (mu_x * m0);
    m_m3 = - (c * dt * sigma_x) / (m0 * epsilon_naught* mu_x);
    m_m4 = 0;
}



void PML_coefs::init_coefs_Hy(double sigma_x, double sigma_y, double mu_y, double dt)
{
    double m0 {(1/dt) + (sigma_x / (2 * epsilon_naught)) };
    
    m_m1 = ((1/dt) - (sigma_x / (2 * epsilon_naught))) / m0;
    m_m2 = -c / (mu_y * m0);
    m_m3 = - (c * dt * sigma_y) / (m0 * epsilon_naught *mu_y);
    m_m4 = 0;
}



void PML_coefs::init_coefs_Dz(double sigma_x, double sigma_y, double dt)
{
    double m0 {(1/dt) + ((sigma_x + sigma_y) / (2 * epsilon_naught)) + ((sigma_x *sigma_y * dt)/(4 * epsilon_naught * epsilon_naught))};
    
    m_m1 = ((1/dt) - ((sigma_x + sigma_y) / (2 * epsilon_naught)) - ((sigma_x *sigma_y * dt)/(4 * epsilon_naught * epsilon_naught))) / m0;
    m_m2 =  c / m0;
    m_m3 =  0;
    m_m4 = - (dt * sigma_y * sigma_x) / (m0 * epsilon_naught * epsilon_naught);
}

std::ostream& operator<< (std::ostream &out, const PML_coefs &coef)
{
    out << "m1: " << coef.m_m1 << '\n'
        << "m2: " << coef.m_m2 << '\n'
        << "m3: " << coef.m_m3 << '\n'
        << "m4: " << coef.m_m4 << '\n';
    return out;
}
