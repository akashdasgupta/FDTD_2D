#include "Ez_point.h"
#include <PML_boundry.h>
#include <iostream>

int index(int i,int j, int Nx){return ((i*Nx) + j);}

void update_H_master_top(Point space[], int Nx, int Ny, double mux[], double muy[], double dx, double dy,  Point row_bellow[], int pml_size, PML_coefs coefs_Hx[], PML_coefs coefs_Hy[],  PML_coefs coefs_Hx_top[],  PML_coefs coefs_Hy_top[], double ICHx[], double ICHy[])
{
    
}

void update_E_master_top(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt, int pml_size, PML_coefs coefs_Dz[], double IDz[], double ICDz[])
{

}



void update_H_worker(Point space[], int Nx, int Ny, double mux[], double muy[], double dx, double dy, Point row_bellow[], int pml_size, PML_coefs coefs_Hx[], PML_coefs coefs_Hy[], double IHx[], double IHy[], double ICHx[], double ICHy[])
{
    double CEx{}; //curl of Ez in x
    double CEy{}; //curl of Ez in y
    
    double new_Hx{};
    double new_Hy{};
    double cdt{5e-10};
    
    
    for (int i=0; i<Ny-1; ++i)
    {
        for (int j=0; j<pml_size; ++j)
        {
            CEx = (space[index(i+1,j,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dy;
            CEy = -(space[index(i,j+1,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dx;
            

            ICHx[index(i,j,pml_size)] +=CEx;
            ICHy[index(i,j,pml_size)] +=CEy;
                        
            new_Hx = coefs_Hx[pml_size-j-1].Getm1() * space[index(i,j,Nx)].GetHx() 
            + coefs_Hx[pml_size-j-1].Getm2() * CEx
             + coefs_Hx[pml_size-j-1].Getm3() * ICHx[index(i,j,pml_size)];
            
            
            new_Hy = coefs_Hy[pml_size-j-1].Getm1()* space[index(i,j,Nx)].GetHy() 
            + coefs_Hy[pml_size-j-1].Getm2() *CEy
            + coefs_Hy[pml_size-j-1].Getm3()*ICHy[index(i,j,pml_size)];

            
            space[index(i,j,Nx)].SetHx(new_Hx);
            space[index(i,j,Nx)].SetHy(new_Hy);
            
        }

        for (int j=pml_size; j<Nx-pml_size; ++j)
        {
            CEx = (space[index(i+1,j,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dy;
            CEy = (space[index(i,j+1,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dx;
            
            new_Hx = space[index(i,j,Nx)].GetHx() - mux[index(i,j,Nx)] * CEx;
            new_Hy = space[index(i,j,Nx)].GetHy() + muy[index(i,j,Nx)] * CEy;
            
            space[index(i,j,Nx)].SetHx(new_Hx);
            space[index(i,j,Nx)].SetHy(new_Hy);
            
        }
        
        for (int j=Nx-pml_size; j<Nx-1; ++j)
        {
            CEx = (space[index(i+1,j,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dy;
            CEy = -(space[index(i,j+1,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dx;
            
            ICHx[index(i+Ny,j-Nx+pml_size,pml_size)] +=CEx;
            ICHy[index(i+Ny,j-Nx+pml_size,pml_size)] +=CEy;
            
            
            new_Hx = coefs_Hx[j-(Nx-pml_size)].Getm1() * space[index(i,j,Nx)].GetHx() 
            + coefs_Hx[j-(Nx-pml_size)].Getm2() * CEx
            + coefs_Hx[j-(Nx-pml_size)].Getm3() * ICHx[index(i+Ny,j-Nx+pml_size,pml_size)];
            
            
            new_Hy = coefs_Hy[j-(Nx-pml_size)].Getm1()* space[index(i,j,Nx)].GetHy() 
            + coefs_Hy[j-(Nx-pml_size)].Getm2() * CEy 
             + coefs_Hy[j-(Nx-pml_size)].Getm3()*ICHy[index(i+Ny,j-Nx+pml_size,pml_size)];

            space[index(i,j,Nx)].SetHx(new_Hx);
            space[index(i,j,Nx)].SetHy(new_Hy);
            
        }
        
        {
        int j {Nx-1};
        
        CEx = (space[index(i+1,j,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dy;
        CEy = -(0-space[index(i,j,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dx;
        
        ICHx[index(i+Ny,j-Nx+pml_size,pml_size)] +=CEx;
        ICHy[index(i+Ny,j-Nx+pml_size,pml_size)] +=CEy;
    
        
        new_Hx = coefs_Hx[j-(Nx-pml_size)].Getm1() * space[index(i,j,Nx)].GetHx() 
        + coefs_Hx[j-(Nx-pml_size)].Getm2() * CEx 
         + coefs_Hx[j-(Nx-pml_size)].Getm3() * ICHx[index(i+Ny,j-Nx+pml_size,pml_size)];
        
        
        new_Hy = coefs_Hy[j-(Nx-pml_size)].Getm1()* space[index(i,j,Nx)].GetHy() 
        + coefs_Hy[j-(Nx-pml_size)].Getm2()* CEy
        + coefs_Hy[j-(Nx-pml_size)].Getm3()*ICHy[index(i+Ny,j-Nx+pml_size,pml_size)];
        
        space[index(i,j,Nx)].SetHx(new_Hx);
        space[index(i,j,Nx)].SetHy(new_Hy);
        }
    }
    
    {
    int i{Ny-1};
    
    for (int j=0; j<pml_size; ++j)
    {
        CEx = (row_bellow[j].GetEz() - space[index(i,j,Nx)].GetEz()) / dy;
        CEy = -(space[index(i,j+1,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dx;
        
        ICHx[index(i,j,pml_size)] +=CEx;
        ICHy[index(i,j,pml_size)] +=CEy;
                
        new_Hx = coefs_Hx[pml_size-(j+1)].Getm1() * space[index(i,j,Nx)].GetHx() 
        + coefs_Hx[pml_size-(j+1)].Getm2() * CEx 
         + coefs_Hx[pml_size-(j+1)].Getm3() * ICHx[index(i,j,pml_size)];
        
        
        new_Hy = coefs_Hy[pml_size-(j+1)].Getm1()* space[index(i,j,Nx)].GetHy() 
        + coefs_Hy[pml_size-(j+1)].Getm2() * CEy 
         + coefs_Hy[pml_size-(j+1)].Getm3()*ICHy[index(i,j,pml_size)];
        
        space[index(i,j,Nx)].SetHx(new_Hx);
        space[index(i,j,Nx)].SetHy(new_Hy);
        
    }

    for (int j=pml_size; j<Nx-pml_size; ++j)
    {
        CEx = (row_bellow[j].GetEz() - space[index(i,j,Nx)].GetEz()) / dy;
        CEy = (space[index(i,j+1,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dx;
        
        new_Hx = space[index(i,j,Nx)].GetHx() - mux[index(i,j,Nx)] * CEx;
        new_Hy = space[index(i,j,Nx)].GetHy() + muy[index(i,j,Nx)] * CEy;
        
        space[index(i,j,Nx)].SetHx(new_Hx);
        space[index(i,j,Nx)].SetHy(new_Hy);
        
    }
    
    for (int j=Nx-pml_size; j<Nx-1; ++j)
    {
        CEx = (row_bellow[j].GetEz() - space[index(i,j,Nx)].GetEz()) / dy;
        CEy = -(space[index(i,j+1,Nx)].GetEz() - space[index(i,j,Nx)].GetEz()) / dx;
        
        ICHx[index(i+Ny,j-Nx+pml_size,pml_size)] +=CEx;
        ICHy[index(i+Ny,j-Nx+pml_size,pml_size)] +=CEy;
        
        
        new_Hx = coefs_Hx[j-(Nx-pml_size)].Getm1() * space[index(i,j,Nx)].GetHx() 
        + coefs_Hx[j-(Nx-pml_size)].Getm2() * CEx 
         + coefs_Hx[j-(Nx-pml_size)].Getm3() * ICHx[index(i+Ny,j-Nx+pml_size,pml_size)];
        
        
        new_Hy = coefs_Hy[j-(Nx-pml_size)].Getm1()* space[index(i,j,Nx)].GetHy() 
        + coefs_Hy[j-(Nx-pml_size)].Getm2()* CEy 
         + coefs_Hy[j-(Nx-pml_size)].Getm3()*ICHy[index(i+Ny,j-Nx+pml_size,pml_size)];
        
        space[index(i,j,Nx)].SetHx(new_Hx);
        space[index(i,j,Nx)].SetHy(new_Hy);

    }
    
    
    {
    int j {Nx-1};
    CEx = (row_bellow[j].GetEz() - space[index(i,j,Nx)].GetEz()) / dy;
    CEy = -(0 - space[index(i,j,Nx)].GetEz()) / dx;
    
    ICHx[index(i+Ny,j-Nx+pml_size,pml_size)] +=CEx;
    ICHy[index(i+Ny,j-Nx+pml_size,pml_size)] +=CEy;
    
    new_Hx = coefs_Hx[j-(Nx-pml_size)].Getm1() * space[index(i,j,Nx)].GetHx() 
    + coefs_Hx[j-(Nx-pml_size)].Getm2() * CEx 
     + coefs_Hx[j-(Nx-pml_size)].Getm3() * ICHx[index(i+Ny,j-Nx+pml_size,pml_size)];
    
    
    new_Hy = coefs_Hy[j-(Nx-pml_size)].Getm1()* space[index(i,j,Nx)].GetHy() 
    + coefs_Hy[j-(Nx-pml_size)].Getm2() * CEy 
    + coefs_Hy[j-(Nx-pml_size)].Getm3()*ICHy[index(i+Ny,j-Nx+pml_size,pml_size)];
    
    
    space[index(i,j,Nx)].SetHx(new_Hx);
    space[index(i,j,Nx)].SetHy(new_Hy);
    }
    }

}

void update_E_worker(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt, Point row_above[],int pml_size, PML_coefs coefs_Dz[], double IDz[], double ICDz[])
{
    double CHz{};
    double new_Dz{};
    double c{299792458};

    {
    int i{0};
    {
    int j{0};
    CHz = ((space[index(i,j,Nx)].GetHy()) / dx ) - ((space[index(i,j,Nx)].GetHx() - row_above[j].GetHx()) / dy );
    
    IDz[index(i,j,pml_size)] +=space[index(i,j,Nx)].GetDz();


    new_Dz = coefs_Dz[pml_size-(j+1)].Getm1() * space[index(i,j,Nx)].GetDz() 
    + coefs_Dz[pml_size-(j+1)].Getm2()* CHz 
    + coefs_Dz[pml_size-(j+1)].Getm4() * IDz[index(i,j,pml_size)] ;

    
    space[index(i,j,Nx)].SetDz(new_Dz);
    space[index(i,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
    }
    
    for (int j=1; j<pml_size; ++j)
    {
        CHz = ((space[index(i,j,Nx)].GetHy() - space[index(i,j-1,Nx)].GetHy()) / dx ) - ((space[index(i,j,Nx)].GetHx() - row_above[j].GetHx()) / dy );
        IDz[index(i,j,pml_size)] +=space[index(i,j,Nx)].GetDz();


        new_Dz = coefs_Dz[pml_size-(j+1)].Getm1() * space[index(i,j,Nx)].GetDz() 
        + coefs_Dz[pml_size-(j+1)].Getm2()* CHz 
         + coefs_Dz[pml_size-(j+1)].Getm4() * IDz[index(i,j,pml_size)] ;

        
        
        space[index(i,j,Nx)].SetDz(new_Dz);
        space[index(i,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
    }
    
    for (int j=pml_size; j<Nx-pml_size; ++j)
    {
        CHz = ((space[index(i,j,Nx)].GetHy() - space[index(i,j-1,Nx)].GetHy()) / dx ) - ((space[index(i,j,Nx)].GetHx() - row_above[j].GetHx()) / dy );
        new_Dz = space[index(i,j,Nx)].GetDz() + (c * dt * CHz);
        
        
        space[index(i,j,Nx)].SetDz(new_Dz);
        space[index(i,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
    }
    
    for (int j=Nx-pml_size; j<=Nx-1; ++j)
    {
        CHz = ((space[index(i,j,Nx)].GetHy() - space[index(i,j-1,Nx)].GetHy()) / dx ) - ((space[index(i,j,Nx)].GetHx() - row_above[j].GetHx()) / dy );
        
        IDz[index(i+Ny,j-Nx+pml_size,pml_size)] +=space[index(i,j,Nx)].GetDz();

    
        new_Dz = coefs_Dz[j-(Nx-pml_size)].Getm1() * space[index(i,j,Nx)].GetDz() 
        + coefs_Dz[j-(Nx-pml_size)].Getm2()* CHz 
         + coefs_Dz[j-(Nx-pml_size)].Getm4() * IDz[index(i+Ny,j-Nx+pml_size,pml_size)] ;

        space[index(i,j,Nx)].SetDz(new_Dz);
        space[index(i,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
    }
    }
    
    for (int i=1; i<Ny; ++i)
    {
        {
        //ylow only, i!=0, j=0
        int j{0};
        CHz = ((space[index(i,j,Nx)].GetHy()) / dx ) - ((space[index(i,j,Nx)].GetHx() - space[index(i-1,j,Nx)].GetHx()) / dy );
        
        IDz[index(i,j,pml_size)] +=space[index(i,j,Nx)].GetDz();
        
         new_Dz = coefs_Dz[pml_size-(j+1)].Getm1() * space[index(i,j,Nx)].GetDz() 
         + coefs_Dz[pml_size-(j+1)].Getm2()* CHz 
         + coefs_Dz[pml_size-(j+1)].Getm4() * IDz[index(i,j,pml_size)] ;
        
        
        space[index(i,j,Nx)].SetDz(new_Dz);
        space[index(i,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
    
        for (int j=1; j<pml_size; ++j)
        {
            CHz = ((space[index(i,j,Nx)].GetHy() - space[index(i,j-1,Nx)].GetHy()) / dx ) - ((space[index(i,j,Nx)].GetHx() - space[index(i-1,j,Nx)].GetHx()) / dy );
            IDz[index(i,j,pml_size)] +=space[index(i,j,Nx)].GetDz();


            new_Dz = coefs_Dz[pml_size-(j+1)].Getm1() * space[index(i,j,Nx)].GetDz() 
            + coefs_Dz[pml_size-(j+1)].Getm2()* CHz
             + coefs_Dz[pml_size-(j+1)].Getm4() * IDz[index(i,j,pml_size)] ;
            

            space[index(i,j,Nx)].SetDz(new_Dz);
            space[index(i,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
        
        for (int j=pml_size; j<Nx-pml_size; ++j)
        {
            CHz = ((space[index(i,j,Nx)].GetHy() - space[index(i,j-1,Nx)].GetHy()) / dx ) - ((space[index(i,j,Nx)].GetHx() - space[index(i-1,j,Nx)].GetHx()) / dy );
            new_Dz = space[index(i,j,Nx)].GetDz() + (c * dt * CHz);
            
            space[index(i,j,Nx)].SetDz(new_Dz);
            space[index(i,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
        
        for (int j=Nx-pml_size; j<=Nx-1; ++j)
        {
            CHz = ((space[index(i,j,Nx)].GetHy() - space[index(i,j-1,Nx)].GetHy()) / dx ) - ((space[index(i,j,Nx)].GetHx() - space[index(i-1,j,Nx)].GetHx()) / dy );
            
            IDz[index(i+Ny,j-Nx+pml_size,pml_size)] +=space[index(i,j,Nx)].GetDz();

        
            new_Dz = coefs_Dz[j-(Nx-pml_size)].Getm1() * space[index(i,j,Nx)].GetDz() 
            + coefs_Dz[j-(Nx-pml_size)].Getm2()* CHz
             + coefs_Dz[j-(Nx-pml_size)].Getm4() * IDz[index(i+Ny,j-Nx+pml_size,pml_size)] ;


            space[index(i,j,Nx)].SetDz(new_Dz);
            space[index(i,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
    }        
}












