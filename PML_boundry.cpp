#include <PML_boundry.h>
#include <iostream>
#include <cmath>


double epsilon_naught{8.8541878128e-12};
double c{299792458};

int index(int i,int j, int Nx){return ((i*Nx) + j);}

void sigma_setter(double sigma_x[], double sigma_y[], int pml_size, double dt)
{
    for(int i=0; i<pml_size; ++i)
    {
        double factor{static_cast<double>(i+1) / static_cast<double>(pml_size)};
        //std::cout << factor << std::endl;

        sigma_x[i] = epsilon_naught * (factor*factor* factor) / (2 * dt);
        //std::cout << sigma_x[i] << std::endl;
    }
    
    for(int i=0; i<pml_size; ++i)
    {
        double factor{static_cast<double>(i / pml_size)};
        sigma_y[i] = epsilon_naught * (factor*factor* factor) / (2 * dt);
    }
}

void param_setter_x(PML_coefs Hx[], PML_coefs Hy[], PML_coefs Dz[], double sigma_x[], double mu_x, double mu_y, double dt, int size)
{
    for (int i=0; i<size; ++i)
    {
        Hx[i].init_coefs_Hx(sigma_x[i], 0, mu_x, dt);
        Hy[i].init_coefs_Hy(sigma_x[i], 0, mu_y, dt);
        Dz[i].init_coefs_Dz(sigma_x[i], 0, dt);
    }
}

void param_setter_ytop(PML_coefs Hx[], PML_coefs Hy[], PML_coefs Dz[], double sigma_x[], double sigma_y[],  double mu_x, double mu_y, double dt, int pml_size, int Nx)
{
    for(int i=0; i<pml_size; ++i)
    {
        for(int j=0; j<pml_size; ++j)
        {
            Hx[index(i,j,Nx)].init_coefs_Hx(sigma_x[pml_size-(j+1)], sigma_y[pml_size-(i+1)], mu_x, dt);
            Hy[index(i,j,Nx)].init_coefs_Hy(sigma_x[pml_size-(j+1)], sigma_y[pml_size-(i+1)], mu_y, dt);
            Dz[index(i,j,Nx)].init_coefs_Dz(sigma_x[pml_size-(j+1)], sigma_y[pml_size-(i+1)], dt);
        }
        
        for(int j=pml_size; j<Nx-pml_size; ++j)
        {
            Hx[index(i,j,Nx)].init_coefs_Hx(0, sigma_y[pml_size-(i+1)], mu_x, dt);
            Hy[index(i,j,Nx)].init_coefs_Hy(0, sigma_y[pml_size-(i+1)], mu_y, dt);
            Dz[index(i,j,Nx)].init_coefs_Dz(0, sigma_y[pml_size-(i+1)], dt);
        }
        
        for(int j=Nx-pml_size; j<Nx; ++j)
        {
            Hx[index(i,j,Nx)].init_coefs_Hx(sigma_x[j-Nx+pml_size], sigma_y[pml_size-(i+1)], mu_x, dt);
            Hy[index(i,j,Nx)].init_coefs_Hy(sigma_x[j-Nx+pml_size], sigma_y[pml_size-(i+1)], mu_y, dt);
            Dz[index(i,j,Nx)].init_coefs_Dz(sigma_x[j-Nx+pml_size], sigma_y[pml_size-(i+1)], dt);
        }
    }
}

void param_setter_ybottom(PML_coefs Hx[], PML_coefs Hy[], PML_coefs Dz[], double sigma_x[], double sigma_y[],  double mu_x, double mu_y, double dt, int pml_size, int Nx)
{
    for(int i=0; i<pml_size; ++i)
    {
        for(int j=0; j<pml_size; ++j)
        {
            Hx[index(i,j,Nx)].init_coefs_Hx(sigma_x[pml_size-(j+1)], sigma_y[i], mu_x, dt);
            Hy[index(i,j,Nx)].init_coefs_Hy(sigma_x[pml_size-(j+1)], sigma_y[i], mu_y, dt);
            Dz[index(i,j,Nx)].init_coefs_Dz(sigma_x[pml_size-(j+1)], sigma_y[i], dt);
        }
        
        for(int j=pml_size; j<Nx-pml_size; ++j)
        {
            Hx[index(i,j,Nx)].init_coefs_Hx(0, sigma_y[i], mu_x, dt);
            Hy[index(i,j,Nx)].init_coefs_Hy(0, sigma_y[i], mu_y, dt);
            Dz[index(i,j,Nx)].init_coefs_Dz(0, sigma_y[i], dt);
        }
        
        for(int j=Nx-pml_size; j<Nx; ++j)
        {
            Hx[index(i,j,Nx)].init_coefs_Hx(sigma_x[j-Nx+pml_size], sigma_y[i], mu_x, dt);
            Hy[index(i,j,Nx)].init_coefs_Hy(sigma_x[j-Nx+pml_size], sigma_y[i], mu_y, dt);
            Dz[index(i,j,Nx)].init_coefs_Dz(sigma_x[j-Nx+pml_size], sigma_y[i], dt);
        }
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
