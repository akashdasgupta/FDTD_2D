#include "Ez_point.h"

int index(int i,int j, int Nx){return ((i*Nx) + j);}

void update_H_master(Point space[], int Nx, int Ny, double mux[], double muy[], double dx, double dy)
{
    double CEx{}; //curl of Ez in x
    double CEy{}; //curl of Ez in y
    
    double new_Hx{};
    double new_Hy{};
    
    
    for (int i=0; i<Ny-1; ++i)
    {
        for (int j=0; j<Nx-1; ++j)
        {
            CEx = (space[index(i+1,j,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dy;
            CEy = (space[index(i,j+1,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dx;
            
            new_Hx = space[index(i,j,Nx)].GetHx() - mux[index(i,j,Nx)] * CEx;
            new_Hy = space[index(i,j,Nx)].GetHy() + muy[index(i,j,Nx)] * CEy;
            
            space[index(i,j,Nx)].SetHx(new_Hx);
            space[index(i,j,Nx)].SetHy(new_Hy);
            
        }
        int j{Ny-1};
        CEx = (space[index(i+1,j,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dy;
        CEy = (0 - space[index(i,j,Nx)].GetEz()) / dx;
        
        new_Hx = space[index(i,j,Nx)].GetHx() - mux[index(i,j,Nx)] * CEx;
        new_Hy = space[index(i,j,Nx)].GetHy() + muy[index(i,j,Nx)] * CEy;
        
        space[index(i,j,Nx)].SetHx(new_Hx);
        space[index(i,j,Nx)].SetHy(new_Hy);
    }
    
    int i{Ny-1};
    
    for (int j; j<Nx-1; ++j)
        {
            CEx = (0 - space[index(i,j,Nx)].GetEz()) / dy;
            CEy = (space[index(i,j+1,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dx;
            
            new_Hx = space[index(i,j,Nx)].GetHx() - mux[index(i,j,Nx)] * CEx;
            new_Hy = space[index(i,j,Nx)].GetHy() + muy[index(i,j,Nx)] * CEy;
            
            space[index(i,j,Nx)].SetHx(new_Hx);
            space[index(i,j,Nx)].SetHy(new_Hy);
            
        }
    int j{Nx-1};
    CEx = (0 - space[index(i,j,Nx)].GetEz()) / dy;
    CEy = (0 - space[index(i,j,Nx)].GetEz()) / dx;
    
    new_Hx = space[index(i,j,Nx)].GetHx() - mux[index(i,j,Nx)] * CEx;
    new_Hy = space[index(i,j,Nx)].GetHy() + muy[index(i,j,Nx)] * CEy;
    
    space[index(i,j,Nx)].SetHx(new_Hx);
    space[index(i,j,Nx)].SetHy(new_Hy);
}

void update_E_master(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt)
{
    double CHz{};
    double new_Dz{};
    double c{299792458};
    
    // i=0, j=0, xlow and y low coinside
    CHz = ((space[0].GetHy()) / dx ) - ((space[0].GetHx()) / dy );
    new_Dz = space[0].GetDz() + (c * dt * CHz);
    space[0].SetDz(new_Dz);
    space[0].SetEz(new_Dz / ep[index(0,0,Nx)]);

    // Xlow only, i=0, j != 0
    int i{0};
    for (int j=1; j<=Nx-1; ++j)
    {
        CHz = ((space[index(i,j,Nx)].GetHy() - space[index(i,j-1,Nx)].GetHy()) / dx ) - ((space[index(i,j,Nx)].GetHx()) / dy );
        new_Dz = space[index(i,j,Nx)].GetDz() + (c * dt * CHz);
        space[index(i,j,Nx)].SetDz(new_Dz);
        space[index(i,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
    }
       
    
    for (int i=1; i<=Ny-1; ++i)
    {
        //xlow only, i!=0, j=0
        int j{0};
        CHz = ((space[index(i,j,Nx)].GetHy()) / dx ) - ((space[index(i,j,Nx)].GetHx() - space[index(i-1,j,Nx)].GetHx()) / dy );
        new_Dz = space[index(i,j,Nx)].GetDz() + (c * dt * CHz);
        space[index(i,j,Nx)].SetDz(new_Dz);
        space[index(i,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);

        
        for (int j=1; j<=Nx-1; ++j)
        {
            CHz = ((space[index(i,j,Nx)].GetHy() - space[index(i,j-1,Nx)].GetHy()) / dx ) - ((space[index(i,j,Nx)].GetHx() - space[index(i-1,j,Nx)].GetHx()) / dy );
            new_Dz = space[index(i,j,Nx)].GetDz() + (c * dt * CHz);
            space[index(i,j,Nx)].SetDz(new_Dz);
            space[index(i,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
    }
        
}



void update_H_worker(Point space[], int Nx, int Ny, double mux[], double muy[], double dx, double dy, Point row_bellow[])
{
    double CEx{}; //curl of Ez in x
    double CEy{}; //curl of Ez in y
    
    double new_Hx{};
    double new_Hy{};
    
    
    for (int i=0; i<Ny-1; ++i)
    {
        for (int j=0; j<Nx-1; ++j)
        {
            CEx = (space[index(i+1,j,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dy;
            CEy = (space[index(i,j+1,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dx;
            
            new_Hx = space[index(i,j,Nx)].GetHx() - mux[index(i,j,Nx)] * CEx;
            new_Hy = space[index(i,j,Nx)].GetHy() + muy[index(i,j,Nx)] * CEy;
            
            space[index(i,j,Nx)].SetHx(new_Hx);
            space[index(i,j,Nx)].SetHy(new_Hy);
            
        }
        int j{Ny-1};
        CEx = (space[index(i+1,j,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dy;
        CEy = (0 - space[index(i,j,Nx)].GetEz()) / dx;
        
        new_Hx = space[index(i,j,Nx)].GetHx() - mux[index(i,j,Nx)] * CEx;
        new_Hy = space[index(i,j,Nx)].GetHy() + muy[index(i,j,Nx)] * CEy;
        
        space[index(i,j,Nx)].SetHx(new_Hx);
        space[index(i,j,Nx)].SetHy(new_Hy);
    }
    
    int i{Ny-1};
    
    for (int j; j<Nx-1; ++j)
        {
            CEx = (row_bellow[j].GetEz() - space[index(i,j,Nx)].GetEz()) / dy;
            CEy = (space[index(i,j+1,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dx;
            
            new_Hx = space[index(i,j,Nx)].GetHx() - mux[index(i,j,Nx)] * CEx;
            new_Hy = space[index(i,j,Nx)].GetHy() + muy[index(i,j,Nx)] * CEy;
            
            space[index(i,j,Nx)].SetHx(new_Hx);
            space[index(i,j,Nx)].SetHy(new_Hy);
            
        }
    int j{Nx-1};
    CEx = (row_bellow[j].GetEz() - space[index(i,j,Nx)].GetEz()) / dy;
    CEy = (0 - space[index(i,j,Nx)].GetEz()) / dx;
    
    new_Hx = space[index(i,j,Nx)].GetHx() - mux[index(i,j,Nx)] * CEx;
    new_Hy = space[index(i,j,Nx)].GetHy() + muy[index(i,j,Nx)] * CEy;
    
    space[index(i,j,Nx)].SetHx(new_Hx);
    space[index(i,j,Nx)].SetHy(new_Hy);
    
}

void update_E_worker(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt, Point row_above[])
{
    double CHz{};
    double new_Dz{};
    double c{299792458};
    
    int j{0};
    // i=0, j=0, xlow and y low coinside
    CHz = ((space[0].GetHy()) / dx ) - ((space[0].GetHx() - row_above[j].GetHx()) / dy );
    new_Dz = space[0].GetDz() + (c * dt * CHz);
    space[0].SetDz(new_Dz);
    space[0].SetEz(new_Dz / ep[index(0,0,Nx)]);
    
    int i{0};

    // Xlow only, i=0, j != 0
    for (int j=1; j<=Nx-1; ++j)
    {
        CHz = ((space[index(i,j,Nx)].GetHy() - space[index(i,j-1,Nx)].GetHy()) / dx ) - ((space[index(i,j,Nx)].GetHx() - row_above[j].GetHx()) / dy );
        new_Dz = space[index(i,j,Nx)].GetDz() + (c * dt * CHz);
        space[index(i,j,Nx)].SetDz(new_Dz);
        space[index(i,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
    }
       
    
    for (int i=1; i<=Ny-1; ++i)
    {
        //ylow only, i!=0, j=0
        int j{0};
        CHz = ((space[index(i,j,Nx)].GetHy()) / dx ) - ((space[index(i,j,Nx)].GetHx() - space[index(i-1,j,Nx)].GetHx()) / dy );
        new_Dz = space[index(i,j,Nx)].GetDz() + (c * dt * CHz);
        space[index(i,j,Nx)].SetDz(new_Dz);
        space[index(i,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);

        
        for (int j=1; j<=Nx-1; ++j)
        {
            CHz = ((space[index(i,j,Nx)].GetHy() - space[index(i,j-1,Nx)].GetHy()) / dx ) - ((space[index(i,j,Nx)].GetHx() - space[index(i-1,j,Nx)].GetHx()) / dy );
            new_Dz = space[index(i,j,Nx)].GetDz() + (c * dt * CHz);
            space[index(i,j,Nx)].SetDz(new_Dz);
            space[index(i,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
    }
        
}












