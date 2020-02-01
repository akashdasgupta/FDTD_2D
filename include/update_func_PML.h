#ifndef UPDATE
#define UPDATE

// Updates magnetic field on all points of simulation space, for region within the 'top' PML
void update_H_pml(double Ez[], double Hx[], double Hy[], int Nx, int Ny, double dx, double dy, PML_coefs coefs_Hx[], PML_coefs coefs_Hy[], double ICHx[], double ICHy[]);

// Updates electric field on all points of simulation space, for region within the 'top' PML 
void update_E_pml(double Ez[], double Dz[], double Hx[], double Hy[], int Nx, int Ny, double ep[], double dx, double dy, double dt, PML_coefs coefs_Dz[], double IDz[]);

// Updates magnetic field on all points of simulation space, outside region of the 'top' and 'bottom' PML
void update_H_bulk(double Ez[], double Hx[], double Hy[], int Nx, int Ny, double mux[], double muy[], double dx, double dy, int pml_size, PML_coefs coefs_Hx[], PML_coefs coefs_Hy[], double ICHx[], double ICHy[]);

// Updates electric field on all points of simulation space, outside region of the 'top' and 'bottom' PML
void update_E_bulk(double Ez[], double Dz[], double Hx[], double Hy[], int Nx, int Ny, double ep[], double dx, double dy, double dt,int pml_size, PML_coefs coefs_Dz[], double IDz[]);


#endif
