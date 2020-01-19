#include <math.h>
#include <objects.h>
#include <PML_boundry.h>

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

// void lense_foc(double ep[], int Nx, int Ny, int center_x, int center_y, int radius, int width)
// {
//     double rad_center_under_y{(center_y) + sqrt((radius * radius ) -(width * width / 4))};
//     double rad_center_above_y{(center_y) - sqrt((radius * radius ) -(width * width / 4))};
//     double distance_from_under{}; 
//     double distance_from_above{}; 
//     
//     for(int i=0; i<Ny; ++i)
//     {
//         for(int j=0; j<Nx; ++j)
//         {
//             distance_from_under = sqrt(((j - center_x)*(j - center_x)) +
//                             ((i - rad_center_under_y)*(i - rad_center_under_y)));
//             distance_from_above = sqrt(((j - center_x)*(j - center_x)) +
//                             ((i - rad_center_above_y)*(i - rad_center_above_y)));
//             if (distance_from_above <= radius && distance_from_under <= radius)
//             {
//                 ep[index(i,j,Nx)] = 4.7; //rel. of pyrex
//             }
//             else
//             {
//                 ep[index(i,j,Nx)] = 1;
//             }
//         }
//     }  
// }
