#include <math.h>
#include <objects.h>
#include <core_funcs.h>

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

void double_slit(double ep[], int Nx, int Ny, int height, int center_x, int slit_seperation, int width)
{
    int boundry_1{center_x - (slit_seperation + width) /2};
    int boundry_2{center_x - (slit_seperation - width) /2};
    int boundry_3{center_x + (slit_seperation - width) /2};
    int boundry_4{center_x + (slit_seperation + width) /2};
    
    for(int i=0; i<Ny; ++i)
    {
        for(int j=0; j<Nx; ++j)
        {
            ep[index(i,j,Nx)] = 1.0;
        }
    }
    
    std::cout <<height << " , " << boundry_1 << " , " <<boundry_2 << " , " <<boundry_3 << " , " <<boundry_4 << '\n';
    
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
