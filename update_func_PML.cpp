#include "Ez_point.h"
#include <PML_boundry.h>
#include <iostream>



void update_H_master_top(Point space[], int Nx, int Ny, double mux[], double muy[], double dx, double dy, PML_coefs coefs_Hx[], PML_coefs coefs_Hy[], double ICHx[], double ICHy[])
{
    double CEx{}; //curl of Ez in x
    double CEy{}; //curl of Ez in y
    
    double new_Hx{};
    double new_Hy{};
    
    for (int i=0; i<Ny; ++i)
    {
        for (int j=0; j<Nx-1; ++j)
        {
            CEx = (space[index(i+2,j,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dy;
            CEy = -(space[index(i+1,j+1,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dx;
            

            ICHx[index(i,j,Nx)] +=CEx;
            ICHy[index(i,j,Nx)] +=CEy;
                        
            new_Hx = coefs_Hx[index(i,j,Nx)].Getm1() * space[index(i+1,j,Nx)].GetHx() 
            + coefs_Hx[index(i,j,Nx)].Getm2() * CEx
             + coefs_Hx[index(i,j,Nx)].Getm3() * ICHx[index(i,j,Nx)];
            
            
            new_Hy = coefs_Hy[index(i,j,Nx)].Getm1()* space[index(i+1,j,Nx)].GetHy() 
            + coefs_Hy[index(i,j,Nx)].Getm2() *CEy
            + coefs_Hy[index(i,j,Nx)].Getm3()*ICHy[index(i,j,Nx)];

            
            space[index(i+1,j,Nx)].SetHx(new_Hx);
            space[index(i+1,j,Nx)].SetHy(new_Hy);
        }

        {
            int j{Nx-1};
            CEx = (space[index(i+2,j,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dy;
            CEy = -(0 - space[index(i+1,j,Nx)].GetEz()) / dx;
            

            ICHx[index(i,j,Nx)] +=CEx;
            ICHy[index(i,j,Nx)] +=CEy;
                        
            new_Hx = coefs_Hx[index(i,j,Nx)].Getm1() * space[index(i+1,j,Nx)].GetHx() 
            + coefs_Hx[index(i,j,Nx)].Getm2() * CEx
             + coefs_Hx[index(i,j,Nx)].Getm3() * ICHx[index(i,j,Nx)];
            
            
            new_Hy = coefs_Hy[index(i,j,Nx)].Getm1()* space[index(i+1,j,Nx)].GetHy() 
            + coefs_Hy[index(i,j,Nx)].Getm2() *CEy
            + coefs_Hy[index(i,j,Nx)].Getm3()*ICHy[index(i,j,Nx)];

            
            space[index(i+1,j,Nx)].SetHx(new_Hx);
            space[index(i+1,j,Nx)].SetHy(new_Hy);

            
        }
    }
}

void update_E_master_top(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt, PML_coefs coefs_Dz[], double IDz[])
{
    double CHz{};
    double new_Dz{};
    double c{299792458};

    {
    int i{0};
    {
    int j{0};
    CHz = ((space[index(i+1,j,Nx)].GetHy()) / dx ) - ((space[index(i+1,j,Nx)].GetHx() - 0) / dy );
    
    IDz[index(i,j,Nx)] +=space[index(i+1,j,Nx)].GetDz();


    new_Dz = coefs_Dz[index(i,j,Nx)].Getm1() * space[index(i+1,j,Nx)].GetDz() 
    + coefs_Dz[index(i,j,Nx)].Getm2()* CHz 
    + coefs_Dz[index(i,j,Nx)].Getm4() * IDz[index(i,j,Nx)] ;

    
    space[index(i+1,j,Nx)].SetDz(new_Dz);
    space[index(i+1,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
    }
    
    for (int j=1; j<=Nx-1; ++j)
    {
        CHz = ((space[index(i+1,j,Nx)].GetHy() - space[index(i+1,j,Nx)].GetHy()) / dx ) - ((space[index(i+1,j,Nx)].GetHx() - 0) / dy );
        
        IDz[index(i,j,Nx)] +=space[index(i+1,j,Nx)].GetDz();

    
        new_Dz = coefs_Dz[index(i,j,Nx)].Getm1() * space[index(i+1,j,Nx)].GetDz() 
        + coefs_Dz[index(i,j,Nx)].Getm2()* CHz 
         + coefs_Dz[index(i,j,Nx)].Getm4() * IDz[index(i,j,Nx)] ;

        space[index(i+1,j,Nx)].SetDz(new_Dz);
        space[index(i+1,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
    }
    }
    
    for (int i=1; i<Ny; ++i)
    {
        {
        //ylow only, i!=0, j=0
        int j{0};
        CHz = ((space[index(i+1,j,Nx)].GetHy()) / dx ) - ((space[index(i+1,j,Nx)].GetHx() - space[index(i+1,j,Nx)].GetHx()) / dy );
        
        IDz[index(i,j,Nx)] +=space[index(i+1,j,Nx)].GetDz();
        
         new_Dz = coefs_Dz[index(i,j,Nx)].Getm1() * space[index(i+1,j,Nx)].GetDz() 
         + coefs_Dz[index(i,j,Nx)].Getm2()* CHz 
         + coefs_Dz[index(i,j,Nx)].Getm4() * IDz[index(i,j,Nx)] ;
        
        
        space[index(i+1,j,Nx)].SetDz(new_Dz);
        space[index(i+1,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
        
        for (int j=1; j<=Nx-1; ++j)
        {
            CHz = ((space[index(i+1,j,Nx)].GetHy() - space[index(i+1,j-1,Nx)].GetHy()) / dx ) - ((space[index(i+1,j,Nx)].GetHx() - space[index(i,j,Nx)].GetHx()) / dy );
            
            IDz[index(i,j,Nx)] +=space[index(i+1,j,Nx)].GetDz();

        
            new_Dz = coefs_Dz[index(i,j,Nx)].Getm1() * space[index(i+1,j,Nx)].GetDz() 
            + coefs_Dz[index(i,j,Nx)].Getm2()* CHz
             + coefs_Dz[index(i,j,Nx)].Getm4() * IDz[index(i,j,Nx)] ;


            space[index(i+1,j,Nx)].SetDz(new_Dz);
            space[index(i+1,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
    }        
}



void update_H_worker(Point space[], int Nx, int Ny, double mux[], double muy[], double dx, double dy, int pml_size, PML_coefs coefs_Hx[], PML_coefs coefs_Hy[], double IHx[], double IHy[], double ICHx[], double ICHy[])
{
    double CEx{}; //curl of Ez in x
    double CEy{}; //curl of Ez in y
    
    double new_Hx{};
    double new_Hy{};
    double cdt{5e-10};
    
    
    for (int i=0; i<Ny; ++i)
    {
        for (int j=0; j<pml_size; ++j)
        {
            CEx = (space[index(i+2,j,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dy;
            CEy = -(space[index(i+1,j+1,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dx;
            

            ICHx[index(i,j,pml_size)] +=CEx;
            ICHy[index(i,j,pml_size)] +=CEy;
                        
            new_Hx = coefs_Hx[pml_size-j-1].Getm1() * space[index(i+1,j,Nx)].GetHx() 
            + coefs_Hx[pml_size-j-1].Getm2() * CEx
             + coefs_Hx[pml_size-j-1].Getm3() * ICHx[index(i,j,pml_size)];
            
            
            new_Hy = coefs_Hy[pml_size-j-1].Getm1()* space[index(i+1,j,Nx)].GetHy() 
            + coefs_Hy[pml_size-j-1].Getm2() *CEy
            + coefs_Hy[pml_size-j-1].Getm3()*ICHy[index(i,j,pml_size)];

            
            space[index(i+1,j,Nx)].SetHx(new_Hx);
            space[index(i+1,j,Nx)].SetHy(new_Hy);
            
        }

        for (int j=pml_size; j<Nx-pml_size; ++j)
        {
            CEx = (space[index(i+2,j,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dy;
            CEy = (space[index(i+1,j+1,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dx;
            
            new_Hx = space[index(i+1,j,Nx)].GetHx() - mux[index(i,j,Nx)] * CEx;
            new_Hy = space[index(i+1,j,Nx)].GetHy() + muy[index(i,j,Nx)] * CEy;
            //std::cout << new_Hy << '\n';
            
            space[index(i+1,j,Nx)].SetHx(new_Hx);
            space[index(i+1,j,Nx)].SetHy(new_Hy);
            
        }
        
        for (int j=Nx-pml_size; j<Nx-1; ++j)
        {
            CEx = (space[index(i+2,j,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dy;
            CEy = -(space[index(i+1,j+1,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dx;
            
            ICHx[index(i+Ny,j-Nx+pml_size,pml_size)] +=CEx;
            ICHy[index(i+Ny,j-Nx+pml_size,pml_size)] +=CEy;
            
            
            new_Hx = coefs_Hx[j-(Nx-pml_size)].Getm1() * space[index(i+1,j,Nx)].GetHx() 
            + coefs_Hx[j-(Nx-pml_size)].Getm2() * CEx
            + coefs_Hx[j-(Nx-pml_size)].Getm3() * ICHx[index(i+Ny,j-Nx+pml_size,pml_size)];
            
            
            new_Hy = coefs_Hy[j-(Nx-pml_size)].Getm1()* space[index(i+1,j,Nx)].GetHy() 
            + coefs_Hy[j-(Nx-pml_size)].Getm2() * CEy 
             + coefs_Hy[j-(Nx-pml_size)].Getm3()*ICHy[index(i+Ny,j-Nx+pml_size,pml_size)];

            space[index(i+1,j,Nx)].SetHx(new_Hx);
            space[index(i+1,j,Nx)].SetHy(new_Hy);
            
        }
        
        {
        int j {Nx-1};
        
        CEx = (space[index(i+2,j,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dy;
        CEy = -(0-space[index(i+1,j,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dx;
        
        ICHx[index(i+Ny,j-Nx+pml_size,pml_size)] +=CEx;
        ICHy[index(i+Ny,j-Nx+pml_size,pml_size)] +=CEy;
    
        
        new_Hx = coefs_Hx[j-(Nx-pml_size)].Getm1() * space[index(i+1,j,Nx)].GetHx() 
        + coefs_Hx[j-(Nx-pml_size)].Getm2() * CEx 
         + coefs_Hx[j-(Nx-pml_size)].Getm3() * ICHx[index(i+Ny,j-Nx+pml_size,pml_size)];
        
        
        new_Hy = coefs_Hy[j-(Nx-pml_size)].Getm1()* space[index(i+1,j,Nx)].GetHy() 
        + coefs_Hy[j-(Nx-pml_size)].Getm2()* CEy
        + coefs_Hy[j-(Nx-pml_size)].Getm3()*ICHy[index(i+Ny,j-Nx+pml_size,pml_size)];
        
        space[index(i+1,j,Nx)].SetHx(new_Hx);
        space[index(i+1,j,Nx)].SetHy(new_Hy);
        }
    }
}

void update_E_worker(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt,int pml_size, PML_coefs coefs_Dz[], double IDz[], double ICDz[])
{
    double CHz{};
    double new_Dz{};
    double c{299792458};

    for (int i=0; i<Ny; ++i)
    {
        {
        //ylow only, i!=0, j=0
        int j{0};
        CHz = ((space[index(i+1,j,Nx)].GetHy()) / dx ) - ((space[index(i+1,j,Nx)].GetHx() - space[index(i,j,Nx)].GetHx()) / dy );
        
        IDz[index(i,j,pml_size)] +=space[index(i+1,j,Nx)].GetDz();
        
         new_Dz = coefs_Dz[pml_size-(j+1)].Getm1() * space[index(i+1,j,Nx)].GetDz() 
         + coefs_Dz[pml_size-(j+1)].Getm2()* CHz 
         + coefs_Dz[pml_size-(j+1)].Getm4() * IDz[index(i,j,pml_size)] ;
        
        
        space[index(i+1,j,Nx)].SetDz(new_Dz);
        space[index(i+1,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
    
        for (int j=1; j<pml_size; ++j)
        {
            CHz = ((space[index(i+1,j,Nx)].GetHy() - space[index(i+1,j-1,Nx)].GetHy()) / dx ) - ((space[index(i+1,j,Nx)].GetHx() - space[index(i,j,Nx)].GetHx()) / dy );
            IDz[index(i,j,pml_size)] +=space[index(i+1,j,Nx)].GetDz();


            new_Dz = coefs_Dz[pml_size-(j+1)].Getm1() * space[index(i+1,j,Nx)].GetDz() 
            + coefs_Dz[pml_size-(j+1)].Getm2()* CHz
             + coefs_Dz[pml_size-(j+1)].Getm4() * IDz[index(i,j,pml_size)] ;
            

            space[index(i+1,j,Nx)].SetDz(new_Dz);
            space[index(i+1,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
        
        for (int j=pml_size; j<Nx-pml_size; ++j)
        {
            CHz = ((space[index(i+1,j,Nx)].GetHy() - space[index(i+1,j-1,Nx)].GetHy()) / dx ) - ((space[index(i+1,j,Nx)].GetHx() - space[index(i,j,Nx)].GetHx()) / dy );
            new_Dz = space[index(i+1,j,Nx)].GetDz() + (c * dt * CHz);
            
            space[index(i+1,j,Nx)].SetDz(new_Dz);
            space[index(i+1,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
        
        for (int j=Nx-pml_size; j<=Nx-1; ++j)
        {
            CHz = ((space[index(i+1,j,Nx)].GetHy() - space[index(i+1,j-1,Nx)].GetHy()) / dx ) - ((space[index(i+1,j,Nx)].GetHx() - space[index(i,j,Nx)].GetHx()) / dy );
            
            IDz[index(i+Ny,j-Nx+pml_size,pml_size)] +=space[index(i+1,j,Nx)].GetDz();

        
            new_Dz = coefs_Dz[j-(Nx-pml_size)].Getm1() * space[index(i+1,j,Nx)].GetDz() 
            + coefs_Dz[j-(Nx-pml_size)].Getm2()* CHz
             + coefs_Dz[j-(Nx-pml_size)].Getm4() * IDz[index(i+Ny,j-Nx+pml_size,pml_size)] ;


            space[index(i+1,j,Nx)].SetDz(new_Dz);
            space[index(i+1,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
    }        
}


void update_H_master_bottom(Point space[], int Nx, int Ny, double mux[], double muy[], double dx, double dy, PML_coefs coefs_Hx[], PML_coefs coefs_Hy[], double ICHx[], double ICHy[])
{
    double CEx{}; //curl of Ez in x
    double CEy{}; //curl of Ez in y
    
    double new_Hx{};
    double new_Hy{};
    
    for (int i=0; i<Ny; ++i)
    {
        for (int j=0; j<Nx-1; ++j)
        {
            CEx = (space[index(i+2,j,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dy;
            CEy = -(space[index(i+1,j+1,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dx;
            

            ICHx[index(i,j,Nx)] +=CEx;
            ICHy[index(i,j,Nx)] +=CEy;
                        
            new_Hx = coefs_Hx[index(i,j,Nx)].Getm1() * space[index(i+1,j,Nx)].GetHx() 
            + coefs_Hx[index(i,j,Nx)].Getm2() * CEx
             + coefs_Hx[index(i,j,Nx)].Getm3() * ICHx[index(i,j,Nx)];
            
            
            new_Hy = coefs_Hy[index(i,j,Nx)].Getm1()* space[index(i+1,j,Nx)].GetHy() 
            + coefs_Hy[index(i,j,Nx)].Getm2() *CEy
            + coefs_Hy[index(i,j,Nx)].Getm3()*ICHy[index(i,j,Nx)];

            
            space[index(i+1,j,Nx)].SetHx(new_Hx);
            space[index(i+1,j,Nx)].SetHy(new_Hy);
        }
        {
            int j{Nx-1};
            CEx = (space[index(i+2,j,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dy;
            CEy = -(0 - space[index(i+1,j,Nx)].GetEz()) / dx;
            

            ICHx[index(i,j,Nx)] +=CEx;
            ICHy[index(i,j,Nx)] +=CEy;
                        
            new_Hx = coefs_Hx[index(i,j,Nx)].Getm1() * space[index(i+1,j,Nx)].GetHx() 
            + coefs_Hx[index(i,j,Nx)].Getm2() * CEx
             + coefs_Hx[index(i,j,Nx)].Getm3() * ICHx[index(i,j,Nx)];
            
            
            new_Hy = coefs_Hy[index(i,j,Nx)].Getm1()* space[index(i+1,j,Nx)].GetHy() 
            + coefs_Hy[index(i,j,Nx)].Getm2() *CEy
            + coefs_Hy[index(i,j,Nx)].Getm3()*ICHy[index(i,j,Nx)];

            
            space[index(i+1,j,Nx)].SetHx(new_Hx);
            space[index(i+1,j,Nx)].SetHy(new_Hy);

            
        }
    }
    {
        int i{Ny-1};
        for (int j=0; j<Nx-1; ++j)
        {
            CEx = (0 - space[index(i+1,j,Nx)].GetEz()) / dy;
            CEy = -(space[index(i+1,j+1,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dx;
            

            ICHx[index(i,j,Nx)] +=CEx;
            ICHy[index(i,j,Nx)] +=CEy;
                        
            new_Hx = coefs_Hx[index(i,j,Nx)].Getm1() * space[index(i+1,j,Nx)].GetHx() 
            + coefs_Hx[index(i,j,Nx)].Getm2() * CEx
                + coefs_Hx[index(i,j,Nx)].Getm3() * ICHx[index(i,j,Nx)];
            
            
            new_Hy = coefs_Hy[index(i,j,Nx)].Getm1()* space[index(i+1,j,Nx)].GetHy() 
            + coefs_Hy[index(i,j,Nx)].Getm2() *CEy
            + coefs_Hy[index(i,j,Nx)].Getm3()*ICHy[index(i,j,Nx)];

            
            space[index(i+1,j,Nx)].SetHx(new_Hx);
            space[index(i+1,j,Nx)].SetHy(new_Hy);
        }
    
    {
        int j{Nx-1};
        CEx = (0- space[index(i+1,j,Nx)].GetEz()) / dy;
        CEy = -(0 - space[index(i+1,j,Nx)].GetEz()) / dx;
        

        ICHx[index(i,j,Nx)] +=CEx;
        ICHy[index(i,j,Nx)] +=CEy;
                    
        new_Hx = coefs_Hx[index(i,j,Nx)].Getm1() * space[index(i+1,j,Nx)].GetHx() 
        + coefs_Hx[index(i,j,Nx)].Getm2() * CEx
            + coefs_Hx[index(i,j,Nx)].Getm3() * ICHx[index(i,j,Nx)];
        
        
        new_Hy = coefs_Hy[index(i,j,Nx)].Getm1()* space[index(i+1,j,Nx)].GetHy() 
        + coefs_Hy[index(i,j,Nx)].Getm2() *CEy
        + coefs_Hy[index(i,j,Nx)].Getm3()*ICHy[index(i,j,Nx)];

        
        space[index(i+1,j,Nx)].SetHx(new_Hx);
        space[index(i+1,j,Nx)].SetHy(new_Hy);

            
        }
    }
}

void update_E_master_bottom(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt, PML_coefs coefs_Dz[], double IDz[])
{
    double CHz{};
    double new_Dz{};
    double c{299792458};

 
    for (int i=0; i<Ny; ++i)
    {
        {
        //ylow only, i!=0, j=0
        int j{0};
        CHz = ((space[index(i+1,j,Nx)].GetHy()) / dx ) - ((space[index(i+1,j,Nx)].GetHx() - space[index(i+1,j,Nx)].GetHx()) / dy );
        
        IDz[index(i,j,Nx)] +=space[index(i+1,j,Nx)].GetDz();
        
         new_Dz = coefs_Dz[index(i,j,Nx)].Getm1() * space[index(i+1,j,Nx)].GetDz() 
         + coefs_Dz[index(i,j,Nx)].Getm2()* CHz 
         + coefs_Dz[index(i,j,Nx)].Getm4() * IDz[index(i,j,Nx)] ;
        
        
        space[index(i+1,j,Nx)].SetDz(new_Dz);
        space[index(i+1,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
        
        for (int j=1; j<=Nx-1; ++j)
        {
            CHz = ((space[index(i+1,j,Nx)].GetHy() - space[index(i+1,j-1,Nx)].GetHy()) / dx ) - ((space[index(i+1,j,Nx)].GetHx() - space[index(i,j,Nx)].GetHx()) / dy );
            
            IDz[index(i,j,Nx)] +=space[index(i,j,Nx)].GetDz();

        
            new_Dz = coefs_Dz[index(i,j,Nx)].Getm1() * space[index(i+1,j,Nx)].GetDz() 
            + coefs_Dz[index(i,j,Nx)].Getm2()* CHz
             + coefs_Dz[index(i,j,Nx)].Getm4() * IDz[index(i,j,Nx)] ;


            space[index(i+1,j,Nx)].SetDz(new_Dz);
            space[index(i+1,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
    }        
}









