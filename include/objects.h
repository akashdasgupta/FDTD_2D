#ifndef OBJECTS
#define OBJECTS
//  Modifies permitivity grid to create planar convex lense
void lense_col(double ep[], int Nx, int Ny, int center_x, int center_y, int radius, int width);
// Modifies permitivity grid to create biconvex lense
void lense_foc(double ep[], int Nx, int Ny, int center_x, int center_y, int radius, int width);
// Modifies permitivity grid to create a double slit
void double_slit(double ep[], int Nx, int Ny, int height, int center_x, int slit_seperation, int width);
// Updates magnetic field on all points of simulation space, outside region of the 'top' and 'bottom' PML 
void tunnelling(double ep[], int Nx, int Ny, int width);
#endif
