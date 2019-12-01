#ifndef UPDATE
#define UPDATE

void update_H_master(Point space[], int Nx, int Ny, double mux[], double muy[], double dx, double dy);
void update_E_master(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt);
void update_H_worker(Point space[], int Nx, int Ny, double mux[], double muy[], double dx, double dy, Point row_bellow[]);
void update_E_worker(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt, Point row_above[]);

int index(int i,int j, int Nx);

#endif
