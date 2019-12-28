#ifndef UPDATE
#define UPDATE

void update_E_master_top(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt, PML_coefs coefs_Dz[], double IDz[]);

void update_H_master_top(Point space[], int Nx, int Ny, double mux[], double muy[], double dx, double dy, PML_coefs coefs_Hx[], PML_coefs coefs_Hy[], double ICHx[], double ICHy[]);

void update_H_worker(Point space[], int Nx, int Ny, double mux[], double muy[], double dx, double dy, int pml_size, PML_coefs coefs_Hx[], PML_coefs coefs_Hy[], double IHx[], double IHy[], double ICHx[], double ICHy[]);

void update_E_worker(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt,int pml_size, PML_coefs coefs_Dz[], double IDz[], double ICDz[]);

void update_H_master_bottom(Point space[], int Nx, int Ny, double mux[], double muy[], double dx, double dy, PML_coefs coefs_Hx[], PML_coefs coefs_Hy[], double ICHx[], double ICHy[]);

void update_E_master_bottom(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt, PML_coefs coefs_Dz[], double IDz[]);

int index(int i,int j, int Nx);

#endif
