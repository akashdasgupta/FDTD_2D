#ifndef UPDATE
#define UPDATE

void update_H(Point space[], int &Nx, int &Ny, double mux[],double muy[], double dx, double dy);
void update_E(Point space[], int &Nx, int &Ny, double ep[], double dx, double dy, double dt);
int index(int i,int j, int Nx);

#endif
