#include <math.h>
#include <objects.h>
#include <core_funcs.h>

/*
 * Modifies permitivity grid to create planar convex lense
 * @param [out] ep the permitivity grid to act on
 * @param Nx Gridsixe (x)
 * @param Ny Gridsize (y)
 * @param center_x position in x of center of lens axis
 * @ param center_x position in y of center of lens axis
 * @param radius radius of curvature of the lens
 * @param width the lens diameter
 */
void lense_col(double ep[], int Nx, int Ny, int center_x, int center_y, int radius, int width)
{
    double rad_center_y{(center_y) + sqrt((radius * radius ) -(width * width / 4))};
    double distance{}; 
    
    for(int i=0; i<Ny; ++i)
    {
        for(int j=0; j<Nx; ++j)
        {
            distance = sqrt(((j - center_x)*(j - center_x)) +
                            ((i - rad_center_y)*(i - rad_center_y)));
            if (i <= center_y && distance <= radius)
            {
                ep[index(i,j,Nx)] = 4.7; //rel. of pyrex
            }
            else
            {
                ep[index(i,j,Nx)] = 1;
            }
        }
    }  
}

/*
 * Modifies permitivity grid to create biconvex lense
 * @param [out] ep the permitivity grid to act on
 * @param Nx Gridsixe (x)
 * @param Ny Gridsize (y)
 * @param center_x position in x of center of lens axis
 * @ param center_x position in y of center of lens axis
 * @param radius radius of curvature of the lens
 * @param width the lens diameter
 */
void lense_foc(double ep[], int Nx, int Ny, int center_x, int center_y, int radius, int width)
{
    double rad_center_under_y{(center_y) + sqrt((radius * radius ) -(width * width / 4))};
    double rad_center_above_y{(center_y) - sqrt((radius * radius ) -(width * width / 4))};
    double distance_from_under{}; 
    double distance_from_above{}; 
    
    for(int i=0; i<Ny; ++i)
    {
        for(int j=0; j<Nx; ++j)
        {
            distance_from_under = sqrt(((j - center_x)*(j - center_x)) +
                            ((i - rad_center_under_y)*(i - rad_center_under_y)));
            distance_from_above = sqrt(((j - center_x)*(j - center_x)) +
                            ((i - rad_center_above_y)*(i - rad_center_above_y)));
            if (distance_from_above <= radius && distance_from_under <= radius)
            {
                ep[index(i,j,Nx)] = 4.7; //rel. of pyrex
            }
            else
            {
                ep[index(i,j,Nx)] = 1;
            }
        }
    }  
}

/*
 * Modifies permitivity grid to create a double slit
 * @param [out] ep the permitivity grid to act on
 * @param Nx Gridsixe (x)
 * @param Ny Gridsize (y)
 * @param height y position of wall
 * @param center_x position in x of center of double slit
 * @param slit_seperation seperation of the double slits
 * @param width width of the slits
 */
void double_slit(double ep[], int Nx, int Ny, int height, int center_x, int slit_seperation, int width)
{
    // creates a wall of high permitivity, with  gaps between boundry 1 to 2, 3 to 4
    int boundry_1{center_x - (slit_seperation + width) /2};
    int boundry_2{center_x - (slit_seperation - width) /2};
    int boundry_3{center_x + (slit_seperation - width) /2};
    int boundry_4{center_x + (slit_seperation + width) /2};
    
    // starts off with 1s grid:
    for(int i=0; i<Ny; ++i)
    {
        for(int j=0; j<Nx; ++j)
        {
            ep[index(i,j,Nx)] = 1.0;
        }
    }
    
    // creats wall:
    for (int j = 0; j<boundry_1; ++j)
    {
        ep[index(height,j,Nx)] = 1e20;
    }
    
    for (int j = boundry_2; j<boundry_3; ++j)
    {
        ep[index(height,j,Nx)] = 1e20;
    }
    
    for (int j = boundry_4; j<Nx; ++j)
    {
        ep[index(height,j,Nx)] = 1e20;
    }
}

/*
 * Modifies permitivity grid to create tunneling setup: whole grid @ refractive index of root 2.5, with 45 deg gap of index 1 diagonally cutting
 * @param [out] ep the permitivity grid to act on
 * @param Nx Gridsixe (x)
 * @param Ny Gridsize (y)
 * @param width width of gap
 */
void tunnelling(double ep[], int Nx, int Ny, int width)
{
    for (int i=0;i<Ny;++i)
    {
        for(int j=0;j<Nx;++j)
        {
            if (i< j || i>(width+j))
            {
                ep[index(i,j,Nx)] = 2.5;
            }
            else
            {
                ep[index(i,j,Nx)] = 1;
            }
        }
    }
}
