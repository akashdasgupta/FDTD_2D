#include <Ez_point.h>
#include <soft_source.h>
#include <PML_boundry.h>
#include <update_func_PML.h>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include<stdio.h>

void SaveToFile(int size, Point row[], std::string name)
{
    std::fstream fs;
    fs.open(name, std::fstream::in | std::fstream::out | std::fstream::app);
    for(int i; i<size; ++i)
    {
        fs << row[i];
    }

}

int main(int argc, char *argv[])
{
    // Initilise and find 
    //how many processes we can play with:
    MPI_Init (&argc, &argv);
    int size{};
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rank{};
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Status status{};
    

    //Makes data sendable via MPI:
    MPI_Datatype MPI_POINT;
    MPI_Type_contiguous(4, MPI_DOUBLE, &MPI_POINT); //initilise
    MPI_Type_commit(&MPI_POINT);//activate
    double c{299792458};

    std::string savefile_path{"data"};
    int Nx{100};
    int Ny{100};
    int iterations{3000};
    double dx {0.1e-8};
    double dy {0.1e-8};
    double dt {dx / (4*c)};
    
    double source_sigma{20*dt};
    double t0 {6 * source_sigma};
    
    int pml_size{Nx / 10};
   
    
    int worker_Ny_base = (Ny-(2*pml_size)) / (size- 1);
    int worker_Ny_remainder = (Ny-(2*pml_size)) % (size-1);
    
    int source_center{50};
    int rank_occuring {50 / worker_Ny_base};
    std::cout << rank_occuring << '\n';
    
    double *sigma_x = new double[pml_size];
    double *sigma_y = new double[pml_size];
    sigma_setter(sigma_x,sigma_y, pml_size,dt);
    
    PML_coefs *HxX_coefs = new PML_coefs[pml_size];
    PML_coefs *HyX_coefs = new PML_coefs[pml_size];
    PML_coefs *DzX_coefs = new PML_coefs[pml_size];
    
    param_setter_x(HxX_coefs, HyX_coefs, DzX_coefs,sigma_x, 1, 1, dt, pml_size);

    
    if (rank == 0)
    {
        
//          for (int i=0; i<pml_size; ++i)
//          {
//               std::cout << HyX_coefs[i] << '\n';
//          }
        
        
        Point *sim_space_top = new Point[pml_size*Nx];
        Point *sim_space_bottom = new Point[pml_size*Nx];

        double *ep = new double[pml_size*Nx];
        double *mux = new double[pml_size*Nx];
        double *muy = new double[pml_size*Nx];


        for (int i=0; i<pml_size; ++i)
        {
            for (int j=0; j<Nx; ++j)
            {
                ep[index(i,j,Nx)] = 1.0;
                mux[index(i,j,Nx)] = 1;
                muy[index(i,j,Nx)] = 1;
            }
        }
        
        PML_coefs *HxY_coefs_top = new PML_coefs[pml_size*Nx];
        PML_coefs *HyY_coefs_top = new PML_coefs[pml_size*Nx];
        PML_coefs *DzY_coefs_top = new PML_coefs[pml_size*Nx];
        param_setter_ytop(HxY_coefs_top,HyY_coefs_top,DzY_coefs_top,sigma_x, sigma_x,1,1,dt,pml_size, Nx);

        
        double *IDz_top = new double[Nx * pml_size];
        double *ICHx_top = new double[Nx * pml_size];
        double *ICHy_top = new double[Nx * pml_size];
        
        PML_coefs *HxY_coefs_bottom = new PML_coefs[pml_size*Nx];
        PML_coefs *HyY_coefs_bottom = new PML_coefs[pml_size*Nx];
        PML_coefs *DzY_coefs_bottom = new PML_coefs[pml_size*Nx];
        param_setter_ybottom(HxY_coefs_bottom,HyY_coefs_bottom,DzY_coefs_bottom,sigma_x, sigma_x,1,1,dt,pml_size, Nx);

        
        double *IDz_bottom = new double[Nx * pml_size];
        double *ICHx_bottom = new double[Nx * pml_size];
        double *ICHy_bottom = new double[Nx * pml_size];



        Point *bottom = new Point[Nx];
        Point *top = new Point[Nx];
        Point *row_from_above = new Point[Nx];
        Point *row_from_bellow = new Point[Nx];

        int above_me{size-1};
        int under_me{1};

        std::cout << "RANK: "<< rank << " ABOVE ME: "<< above_me << " UNDER ME: " << under_me << std::endl ;

        for (int t=1; t<iterations; ++t)
        {   
            for (int k=0; k<Nx; ++k)
            {
                top[k] = sim_space_bottom[k];
            }
            
            for (int k=index(pml_size-1, 0, Nx); k<index(pml_size-1, Nx, Nx); ++k)
            {
                bottom[k-index(pml_size-1, 0, Nx)] = sim_space_top[k];

            } //2e5 kai
            

            
            MPI_Sendrecv(top, Nx, MPI_POINT, above_me ,
                    t, row_from_bellow , Nx,
                    MPI_POINT, under_me,t,
                    MPI_COMM_WORLD, &status);

            MPI_Sendrecv(bottom, Nx, MPI_POINT, under_me ,
                        t+iterations, row_from_above, Nx,
                        MPI_POINT, above_me,t+iterations,
                        MPI_COMM_WORLD, &status);                

            update_H_master_top(sim_space_top, Nx, pml_size, mux, muy, dx, dy, row_from_bellow, HxY_coefs_top, HyY_coefs_top,ICHx_top, ICHy_top);
            update_E_master_top(sim_space_top, Nx, pml_size, ep,dx, dy, dt, DzY_coefs_top,IDz_top);
            
            update_H_master_bottom(sim_space_bottom, Nx, pml_size, mux, muy, dx, dy, HxY_coefs_bottom, HyY_coefs_bottom,ICHx_bottom, ICHy_bottom);
            update_E_master_bottom(sim_space_bottom, Nx, pml_size, ep,dx, dy, dt, row_from_above, DzY_coefs_bottom,IDz_bottom);

            
            
            
            SaveToFile(Nx*pml_size, sim_space_top, "data/"+std::to_string(rank)+"TOP.txt");
            SaveToFile(Nx*(pml_size), sim_space_bottom, "data/"+std::to_string(rank)+"BOTTOM.txt");
        }

            
            delete ep;
            delete mux ;
            delete muy;
            delete top;
            delete bottom;
            delete sim_space_top;
            delete sim_space_bottom;
            delete row_from_above;
            delete row_from_bellow;
    }

    else if (rank == rank_occuring)
    {
        int local_Ny{};
        if ((rank-1) < worker_Ny_remainder)
        {
            local_Ny = worker_Ny_base + 1;
        }
        else 
        {
            local_Ny = worker_Ny_base;
        }


        double *ep = new double[local_Ny*Nx];
        double *mux = new double[local_Ny*Nx];
        double *muy = new double[local_Ny*Nx];


        for (int i=0; i<local_Ny; ++i)
        {
            for (int j=0; j<Nx; ++j)
            {
                ep[index(i,j,Nx)] = 1.0;
                mux[index(i,j,Nx)] = c * dt;
                muy[index(i,j,Nx)] = c * dt;
            }
        }

        double *sourceE = new double[iterations];
        double *sourceHx = new double[iterations];
        double *sourceHy = new double[iterations];


        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};

        inject_soft_source2(sourceE, sourceHx, sourceHy,iterations,dx,dy,dt, 1, 1, 1,1, source_sigma, t0);

        Point *bottom = new Point[Nx];
        Point *top = new Point[Nx];
        Point *sim_space = new Point[local_Ny*Nx];
        Point *row_from_above = new Point[Nx];
        Point *row_from_bellow = new Point[Nx];

        double *IHx = new double[2 * local_Ny * pml_size];
        double *IHy = new double[2 * local_Ny * pml_size];
        double *IDz = new double[2 * local_Ny * pml_size];

        double *ICHx = new double[2 * local_Ny * pml_size];
        double *ICHy = new double[2 * local_Ny * pml_size];
        double *ICDz = new double[2 * local_Ny * pml_size];


        
        int source_center_y {local_Ny/2};
        std::cout << "RANK: "<< rank << " ABOVE ME: "<< above_me << " UNDER ME: " << under_me <<"I think my y size is :  " <<local_Ny << '\n';
       

        for (int t=1; t<iterations; ++t)
        {
            for (int k=0; k<Nx; ++k)
            {
                top[k] = sim_space[k];
            }

            for (int k=index(local_Ny-1, 0, Nx); k<index(local_Ny-1, Nx, Nx); ++k)
            {
                bottom[k-index(local_Ny-1, 0, Nx)] = sim_space[k];
            }
            sim_space[index(source_center_y,50, Nx)].InjectEz(sourceE[t]);
            sim_space[index(source_center_y,50, Nx)].InjectDz(sourceE[t]);

            MPI_Sendrecv(top, Nx, MPI_POINT, above_me ,
                        t, row_from_bellow , Nx,
                        MPI_POINT, under_me,t,
                        MPI_COMM_WORLD, &status);

            MPI_Sendrecv(bottom, Nx, MPI_POINT, under_me ,
                        t+iterations, row_from_above, Nx,
                        MPI_POINT, above_me,t+iterations,
                        MPI_COMM_WORLD, &status);    

            update_H_worker(sim_space, Nx, local_Ny, mux, muy, dx, dy, row_from_bellow, pml_size, HxX_coefs, HyX_coefs, IHx, IHy,ICHx, ICHy);

            update_E_worker(sim_space, Nx, local_Ny, ep,dx, dy, dt, row_from_above, pml_size, DzX_coefs, IDz, ICDz);


 
//             std::cout << "t= " << t << 
//             " Hx= "<<sim_space[index(local_Ny-1, 0, Nx)].GetHx()<<
//             " Hy= "<<sim_space[index(local_Ny-1, 0, Nx)].GetHy() <<
//             " Dz= "<<sim_space[index(local_Ny-1, 0, Nx)].GetDz()  << std::endl;
                
            SaveToFile(Nx*local_Ny, sim_space, "data/"+std::to_string(rank)+".txt");   
        }

        delete ep;
        delete mux ;
        delete muy;
        delete top;
        delete bottom;
        delete sim_space;
        delete row_from_above;
        delete row_from_bellow;

    }
    
    else 
    {
        int local_Ny {};
        if ((rank-1) < worker_Ny_remainder)
        {
            local_Ny = worker_Ny_base + 1;
        }
        else 
        {
            local_Ny = worker_Ny_base;
        }

        double *ep = new double[local_Ny*Nx];
        double *mux = new double[local_Ny*Nx];
        double *muy = new double[local_Ny*Nx];

        for (int j=0; j<Nx; ++j)
        {
            for (int i=0; i<local_Ny; ++i)
            {
                ep[index(i,j,Nx)] = 1.0;
                mux[index(i,j,Nx)] = c * dt;
                muy[index(i,j,Nx)] = c * dt;
            }
        }

        Point *bottom = new Point[Nx];
        Point *top = new Point[Nx];
        Point *sim_space = new Point[local_Ny*Nx];
        Point *row_from_above = new Point[Nx];
        Point *row_from_bellow = new Point[Nx];
        
        double *IHx = new double[2 * local_Ny * pml_size];
        double *IHy = new double[2 * local_Ny * pml_size];
        double *IDz = new double[2 * local_Ny * pml_size];

        double *ICHx = new double[2 * local_Ny * pml_size];
        double *ICHy = new double[2 * local_Ny * pml_size];
        double *ICDz = new double[2 * local_Ny * pml_size];

                

        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};

        std::cout << "RANK: "<< rank << " ABOVE ME: "<< above_me << " UNDER ME: " << under_me << "I think my y size is :  " <<local_Ny << '\n';


        for (int t=1; t<iterations; ++t)
        {            
            for (int k=0; k<Nx; ++k)
            {
                top[k] = sim_space[k];
            }

            for (int k=index(local_Ny-1, 0, Nx); k<index(local_Ny-1, Nx, Nx); ++k)
            {
                bottom[k-index(local_Ny-1, 0, Nx)] = sim_space[k];

            }

            MPI_Sendrecv(top, Nx, MPI_POINT, above_me ,
                    t, row_from_bellow , Nx,
                    MPI_POINT, under_me,t,
                    MPI_COMM_WORLD, &status);

            MPI_Sendrecv(bottom, Nx, MPI_POINT, under_me ,
                        t+iterations, row_from_above, Nx,
                        MPI_POINT, above_me,t+iterations,
                        MPI_COMM_WORLD, &status);

            update_H_worker(sim_space, Nx, local_Ny, mux, muy, dx, dy, row_from_bellow, pml_size, HxX_coefs, HyX_coefs, IHx, IHy,ICHx, ICHy);
            
            update_E_worker(sim_space, Nx, local_Ny, ep,dx, dy, dt, row_from_above, pml_size, DzX_coefs, IDz, ICDz); 
            
            SaveToFile(Nx*local_Ny, sim_space, "data/"+std::to_string(rank)+".txt");

        }

        delete ep;
        delete mux ;
        delete muy;
        delete top;
        delete bottom;
        delete sim_space;
        delete row_from_above;
        delete row_from_bellow;
    }
    
    delete sigma_x;
    delete sigma_y;
    delete HxX_coefs;
    delete HyX_coefs;
    delete DzX_coefs;
    
    MPI_Finalize();    
    return 0;
}
