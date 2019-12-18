#ifndef PMLBOUNDRY
#define PMLBOUNDRY
#include<iostream>

class PML_coefs
{
private:
    double m_m1{0};
    double m_m2{0};
    double m_m3{0};
    double m_m4{0};


public:
    void init_coefs_Hx(double sigma_x, double sigma_y, double mu_x, double dt);
    void init_coefs_Hy(double sigma_x, double sigma_y, double mu_y, double dt);
    void init_coefs_Dz(double sigma_x, double sigma_y, double dt);
    
    friend std::ostream& operator<< (std::ostream &out, const PML_coefs &coef);
    
//     void Setm1(double value);
//     void Setm2(double value);
//     void Setm3(double value);
//     void Setm4(double value);
        
    double Getm1(){return m_m1;}
    double Getm2(){return m_m2;}
    double Getm3(){return m_m3;}
    double Getm4(){return m_m4;}



};
int index(int i,int j, int Nx);
void sigma_setter(double sigma_x[], double sigma_y[], int pml_size, double dt);
void param_setter_x(PML_coefs Hx[], PML_coefs Hy[], PML_coefs Dz[], double sigma_x[], double mu_x, double mu_y, double dt, int size);
void param_setter_ytop(PML_coefs Hx[], PML_coefs Hy[], PML_coefs Dz[], double sigma_x[], double sigma_y[],  double mu_x, double mu_y, double dt,  int pml_size, int Nx);
void param_setter_ybottom(PML_coefs Hx[], PML_coefs Hy[], PML_coefs Dz[], double sigma_x[], double sigma_y[],  double mu_x, double mu_y, double dt, int pml_size, int Nx);
#endif
