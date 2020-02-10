#include <PML_boundry.h>
#include <core_funcs.h>
#include <cmath>

// note these are not global so dont freak out!
double epsilon_naught{8.8541878128e-12};
double c{299792458};

/*
 *Creates an arrsy of conductivity, length of PML  size, which gradually increases. Gradual modeled here by a cubic function. Final conductivity given as epsilon_0 / 2dt 
 * @param [out] sigma_x conductivity in the right/left PML
 * @param [out] sigma_y conductivity in the top/bottom PML
 * @param pml_size size of PML
 * @param dt timestep
 */
void sigma_setter(double sigma_x[], double sigma_y[], int pml_size, double dt)
{
    // at the moment these are the same, but gives some flexibility in case you want one to be stronger than the other
    
    //sigma x loop:
    for(int i=0; i<pml_size; ++i)
    {
        double factor{static_cast<double>(i+1) / static_cast<double>(pml_size)};
        sigma_x[i] = epsilon_naught * (factor*factor* factor) / (2 * dt);
    }
    // sigma y loop:
    for(int i=0; i<pml_size; ++i)
    {
        double factor{static_cast<double>(i+1) / static_cast<double>(pml_size)};
        sigma_y[i] = epsilon_naught * (factor*factor* factor) / (2 * dt);
    }
}

/* Sets parameter array used by the bulk update functions. The sizejust neds to be 1 PML size,it can be re-ued for each row
 * @param [out] Hx array of coficient class for magnetic field in x
 * @param [out] Hy array of coficient class for magnetic field in y
 * @param [out] Dz array of coficient class for electric field in z
 * @param sigma x array of sigmas in the PML region
 * @param mux permibility in x
 * @param muy permibility in y
 * @param dt timestep
 * @paramsize PML size
 */
void param_setter_x(PML_coefs Hx[], PML_coefs Hy[], PML_coefs Dz[], double sigma_x[], double mu_x, double mu_y, double dt, int size)
{
    for (int i=0; i<size; ++i)
    {
        Hx[i].init_coefs_Hx(sigma_x[i], 0, mu_x, dt);
        Hy[i].init_coefs_Hy(sigma_x[i], 0, mu_y, dt);
        Dz[i].init_coefs_Dz(sigma_x[i], 0, dt);
    }
}

/* Sets parameter array used by the PML row update functions.In these ranks the entire region is PML, so unique coficient needed at each point
 * @param [out] Hx array of coficient class for magnetic field in x
 * @param [out] Hy array of coficient class for magnetic field in y
 * @param [out] Dz array of coficient class for electric field in z
 * @param sigma x array of sigmas in the PML region
 * @param sigma y array of sigmas in the PML region
 * @param mux permibility in x
 * @param muy permibility in y
 * @param dt timestep
 * @param pml_size PML size
 * @param ysize number of rows
 * @param Nx rowlength
 */
void param_setter(PML_coefs Hx[], PML_coefs Hy[], PML_coefs Dz[], double sigma_x[], double sigma_y[],  double mu_x, double mu_y, double dt, int pml_size,int y_size, int Nx)
{
    for(int i=0; i<y_size; ++i)
    {
        // x left
        for(int j=0; j<pml_size; ++j)
        {
            Hx[index(i,j,Nx)].init_coefs_Hx(sigma_x[pml_size-(j+1)], sigma_y[i], mu_x, dt);
            Hy[index(i,j,Nx)].init_coefs_Hy(sigma_x[pml_size-(j+1)], sigma_y[i], mu_y, dt);
            Dz[index(i,j,Nx)].init_coefs_Dz(sigma_x[pml_size-(j+1)], sigma_y[i], dt);
        }
        // y PML only
        for(int j=pml_size; j<Nx-pml_size; ++j)
        {
            Hx[index(i,j,Nx)].init_coefs_Hx(0, sigma_y[i], mu_x, dt);
            Hy[index(i,j,Nx)].init_coefs_Hy(0, sigma_y[i], mu_y, dt);
            Dz[index(i,j,Nx)].init_coefs_Dz(0, sigma_y[i], dt);
        }
        // x right
        for(int j=Nx-pml_size; j<Nx; ++j)
        {
            Hx[index(i,j,Nx)].init_coefs_Hx(sigma_x[j-Nx+pml_size], sigma_y[i], mu_x, dt);
            Hy[index(i,j,Nx)].init_coefs_Hy(sigma_x[j-Nx+pml_size], sigma_y[i], mu_y, dt);
            Dz[index(i,j,Nx)].init_coefs_Dz(sigma_x[j-Nx+pml_size], sigma_y[i], dt);
        }
    }
}

/*
 * Initilises a PML_coefs class to be Hx type (holds coficients for magnetic field)
 * @param sigma_x conductivity value at point (x pml type)
 * @param sigma_y conductivity value at point (y pml type)
 * @param mu_x mu for x component
 * @param dt timestep
 */
void PML_coefs::init_coefs_Hx(double sigma_x, double sigma_y, double mu_x, double dt)
{
    double m0 {(1/dt) + (sigma_y / (2 * epsilon_naught))};

    m_m1 = ((1/dt) - (sigma_y / (2 * epsilon_naught))) / m0;
    m_m2 = -c / (mu_x * m0);
    m_m3 = - (c * dt * sigma_x) / (m0 * epsilon_naught* mu_x);
    m_m4 = 0;
}

/*
 * Initilises a PML_coefs class to be Hy type (holds coficients for magnetic field)
 * @param sigma_x conductivity value at point (x pml type)
 * @param sigma_y conductivity value at point (y pml type)
 * @param mu_y mu for y component
 * @param dt timestep
 */
void PML_coefs::init_coefs_Hy(double sigma_x, double sigma_y, double mu_y, double dt)
{
    double m0 {(1/dt) + (sigma_x / (2 * epsilon_naught)) };
    
    m_m1 = ((1/dt) - (sigma_x / (2 * epsilon_naught))) / m0;
    m_m2 = -c / (mu_y * m0);
    m_m3 = - (c * dt * sigma_y) / (m0 * epsilon_naught *mu_y);
    m_m4 = 0;
}

/*
 * Initilises a PML_coefs class to be Dz type (holds coficients for electric field)
 * @param sigma_x conductivity value at point (x pml type)
 * @param sigma_y conductivity value at point (y pml type)
 * @param dt timestep
 */
void PML_coefs::init_coefs_Dz(double sigma_x, double sigma_y, double dt)
{
    double m0 {(1/dt) + ((sigma_x + sigma_y) / (2 * epsilon_naught)) + ((sigma_x *sigma_y * dt)/(4 * epsilon_naught * epsilon_naught))};
    
    m_m1 = ((1/dt) - ((sigma_x + sigma_y) / (2 * epsilon_naught)) - ((sigma_x *sigma_y * dt)/(4 * epsilon_naught * epsilon_naught))) / m0;
    m_m2 =  c / m0;
    m_m3 =  0;
    m_m4 = - (dt * sigma_y * sigma_x) / (m0 * epsilon_naught * epsilon_naught);
}

