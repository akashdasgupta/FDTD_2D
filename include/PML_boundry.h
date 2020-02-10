#ifndef PMLBOUNDRY
#define PMLBOUNDRY
/*Holds the 4 possible coefficients for the PML update funcs. Can initilise but afterwards immutable.*/
class PML_coefs
{
private:
    // 4 cofiveints as member func:
    double m_m1{0};
    double m_m2{0};
    double m_m3{0};
    double m_m4{0};
public:
    // innits:
    void init_coefs_Hx(double sigma_x, double sigma_y, double mu_x, double dt);
    void init_coefs_Hy(double sigma_x, double sigma_y, double mu_y, double dt);
    void init_coefs_Dz(double sigma_x, double sigma_y, double dt);
    
    // some getters:
    double Getm1(){return m_m1;}
    double Getm2(){return m_m2;}
    double Getm3(){return m_m3;}
    double Getm4(){return m_m4;}

};
// Creates an arrsy of conductivity, length of PML  size, which gradually increases. Gradual modeled here by a cubic function. Final conductivity given as epsilon_0 / 2dt 
void sigma_setter(double sigma_x[], double sigma_y[], int pml_size, double dt);
// Sets parameter array used by the bulk update functions. The sizejust neds to be 1 PML size,it can be re-ued for each row
void param_setter_x(PML_coefs Hx[], PML_coefs Hy[], PML_coefs Dz[], double sigma_x[], double mu_x, double mu_y, double dt, int size);

void param_setter(PML_coefs Hx[], PML_coefs Hy[], PML_coefs Dz[], double sigma_x[], double sigma_y[],  double mu_x, double mu_y, double dt, int pml_size,int y_size, int Nx);
#endif
