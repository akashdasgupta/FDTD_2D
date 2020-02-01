#include <Ez_point.h>
#include <PML_boundry.h>
#include <core_funcs.h>
#include <iostream>
// #pragma omp parallel default(none) num_threads(4)

/*
* Updates magnetic field on all points of simulation space, for region within the 'top' PML 
* @param[out] space simulation space over which to work 
* @param Nx number of cols
* @param Ny number of rows
* @param dx finitestep in x
* @param dy finite step in y
* @param coefs_Hx Coefficients for Hx at each point (same dimentions as sim space)
* @param coefs_Hy Coefficients for Hy at each point (same dimentions as sim space)
* @param[out] ICHx Curl integrals for Hx (same dimentions as sim space)
* @param[out] ICHy Curl integrals for Hy (same dimentions as sim space) 
*/
void update_H_pml(Point space[], int Nx, int Ny, double dx, double dy, PML_coefs coefs_Hx[], PML_coefs coefs_Hy[], double ICHx[], double ICHy[])
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
        // Zero at the right edge (mirror like)
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

/*
* Updates electric field on all points of simulation space, for region within the 'top' PML 
* @param[out] space simulation space over which to work 
* @param Nx number of cols
* @param Ny number of rows
* @param ep permitivities on grid
* @param dx finitestep in x
* @param dy finite step in y
* @param coefs_Dz Coefficients for Dz at each point (same dimentions as sim space)
* @param[out] IDz integrals for Dz (same dimentions as sim space) 
*/

void update_E_pml(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt, PML_coefs coefs_Dz[], double IDz[])
{
    double CHz{};
    double new_Dz{};
    double c{299792458};

    for (int i=0; i<Ny; ++i)
    {
        {
        int j{0};
        CHz = ((space[index(i+1,j,Nx)].GetHy() - 0) / dx ) - ((space[index(i+1,j,Nx)].GetHx() - space[index(i,j,Nx)].GetHx()) / dy );
        
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

/*
* Updates magnetic field on all points of simulation space, outside region of the 'top' and 'bottom' PML 
* @param[out] space simulation space over which to work 
* @param Nx number of cols
* @param Ny number of rows
* @param mux permibility in x
* @param muy permibility in y
* @param ep permitivities on grid
* @param dx finitestep in x
* @param dy finite step in y
* @param dt timestep
* @param coefs_Hx Coefficients for Hx at pml region point (size of pml size)
* @param coefs_Hy Coefficients for Hy at pml region points (size of pml size)
* @param[out] ICHx Curl integrals for Hx (size of pml size)
* @param[out] ICHy Curl integrals for Hy (size of pml size)
*/
void update_H_bulk(Point space[], int Nx, int Ny, double mux[], double muy[], double dx, double dy, int pml_size, PML_coefs coefs_Hx[], PML_coefs coefs_Hy[], double ICHx[], double ICHy[])
{
    double CEx{}; //curl of Ez in x
    double CEy{}; //curl of Ez in y
    
    double new_Hx{};
    double new_Hy{};   
    
    for (int i=0; i<Ny; ++i)
    {
        // left PML
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

        // Bulk
        for (int j=pml_size; j<Nx-pml_size; ++j)
        {
            CEx = (space[index(i+2,j,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dy;
            CEy = (space[index(i+1,j+1,Nx)].GetEz() - space[index(i+1,j,Nx)].GetEz()) / dx;
            
            new_Hx = space[index(i+1,j,Nx)].GetHx() - mux[index(i,j,Nx)] * CEx;
            new_Hy = space[index(i+1,j,Nx)].GetHy() + muy[index(i,j,Nx)] * CEy;
            
            space[index(i+1,j,Nx)].SetHx(new_Hx);
            space[index(i+1,j,Nx)].SetHy(new_Hy);
            
        }
   
        //Right PML
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
        
        //Right edge (mirror like)
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

/*
* Updates electric field on all points of simulation space, outside region of the 'top' and 'bottom' PML 
* @param[out] space simulation space over which to work 
* @param Nx number of cols
* @param Ny number of rows
* @param ep permitivities on grid
* @param dx finitestep in x
* @param dy finite step in y
* @param dt timestep
* @param coefs_Dz Coefficients for Dz at each point (size of pml size)
* @param[out] IDz integrals for Dz (size of pml size)
*/

void update_E_bulk(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt,int pml_size, PML_coefs coefs_Dz[], double IDz[])
{
    double CHz{};
    double new_Dz{};
    double c{299792458};

    for (int i=0; i<Ny; ++i)
    {
        // Left most edge, mirrior like:
        {
        int j{0};
        CHz = ((space[index(i+1,j,Nx)].GetHy()) / dx ) - ((space[index(i+1,j,Nx)].GetHx() - space[index(i,j,Nx)].GetHx()) / dy );
        
        IDz[index(i,j,pml_size)] +=space[index(i+1,j,Nx)].GetDz();
        
         new_Dz = coefs_Dz[pml_size-(j+1)].Getm1() * space[index(i+1,j,Nx)].GetDz() 
         + coefs_Dz[pml_size-(j+1)].Getm2()* CHz 
         + coefs_Dz[pml_size-(j+1)].Getm4() * IDz[index(i,j,pml_size)] ;
        
        
        space[index(i+1,j,Nx)].SetDz(new_Dz);
        space[index(i+1,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
        
        // Left PML region :
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
        
        // Bulk :
        for (int j=pml_size; j<Nx-pml_size; ++j)
        {
            CHz = ((space[index(i+1,j,Nx)].GetHy() - space[index(i+1,j-1,Nx)].GetHy()) / dx ) - ((space[index(i+1,j,Nx)].GetHx() - space[index(i,j,Nx)].GetHx()) / dy );
            new_Dz = space[index(i+1,j,Nx)].GetDz() + (c * dt * CHz);
            
            space[index(i+1,j,Nx)].SetDz(new_Dz);
            space[index(i+1,j,Nx)].SetEz(new_Dz / ep[index(i,j,Nx)]);
        }
        
        // Right PML: 
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







