#include "Ez_point.h"
#include <update_func_metalic.h>
#include <iostream>
#include <fstream>
#include <string>

void SaveToFile(int size, Point row[], std::string name)
{
    std::fstream fs;
    fs.open(name, std::fstream::in | std::fstream::out | std::fstream::app);
    for(int i; i<size; ++i)
    {
        fs << row[i];
    }

}

int main()
{
    double c{299792458};

    std::string savefile_path{"data"};
    int Nx{100};
    int Ny{100};
    int iterations{100};
    double dx {10e-9};
    double dy {10e-9};
    double dt {10e-9 / c};
    
    Point *sim_space = new Point[Ny*Nx];
    
    double *ep = new double[Ny*Nx];
    double *mux = new double[Ny*Nx];
    double *muy = new double[Ny*Nx];
    
    
    for (int i=0; i<Nx; ++i)
    {
        for (int j=0; j<Ny; ++j)
        {
            ep[i*j] = 1.0;
            mux[i*j] = c * dt;
            muy[i*j] = c * dt;
        }
    }
    
    for (int t=1; t<=iterations; ++t)
    {
        update_H(sim_space, Nx, Ny, mux, muy, dx, dy)
        update_E(sim_space, Nx, Ny, ep,dx, dy)
    }
    
    SaveToFile(Nx*Ny, sim_space, "data/test_array.txt");

    delete sim_space;
    delete mux;
    delete muy;
    delete ep;
    
    return 0;
}
