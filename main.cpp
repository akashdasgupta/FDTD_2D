#include "Ez_point.h"
#include "soft_source.h"
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
    int iterations{1000};
    double dx {1e-5};
    double dy {1e-5};
    double dt {dx / (2 * c)};
    
    double source_sigma{10*dt};
    double t0 {6 * source_sigma};

    
    Point *sim_space = new Point[Ny*Nx];
    
    double *ep = new double[Ny*Nx];
    double *mux = new double[Ny*Nx];
    double *muy = new double[Ny*Nx];
    
    double *sourceE = new double[iterations];
    double *sourceHx = new double[iterations];
    double *sourceHy = new double[iterations];
    
    
    for (int i=0; i<Nx; ++i)
    {
        for (int j=0; j<Ny; ++j)
        {
            ep[index(i,j,Nx)] = 1.0;
            mux[index(i,j,Nx)] = c * dt;
            muy[index(i,j,Nx)] = c * dt;
        }
    }
    
    //inject_soft_source(sourceE,sourceHx,sourceHy,iterations,dx,dy,dt, 1, 1, 1,1, source_sigma, t0);
    
//    for (int i=0; i<iterations; ++i)
 //   {
 //       std::cout << sourceE[i] << std::endl;
 //   }
    int source_center{50};
    //sim_space[index(source_center,source_center, Nx)].SetEz(1);
    
    for (int t=1; t<=iterations; ++t)
    {
     sim_space[index(source_center,source_center, Nx)].InjectEz(sourceE[t]);
     sim_space[index(source_center,source_center, Nx)].InjectDz(sourceE[t]);
     sim_space[index(source_center,source_center, Nx)].InjectHx(sourceHx[t]);
     sim_space[index(source_center,source_center, Nx)].InjectHy(sourceHy[t]);
//         
        
    update_H(sim_space, Nx, Ny, mux, muy, dx, dy);
    update_E(sim_space, Nx, Ny, ep,dx, dy, dt);
    SaveToFile(Nx*Ny, sim_space, "data/"+std::to_string(t)+".txt");
    }
    


    delete sim_space;
    delete mux;
    delete muy;
    delete ep;
    
    delete sourceE;
    delete sourceHx;
    delete sourceHy;
    
    return 0;
}
