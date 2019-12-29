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
    int Nx{200};
    int Ny{200};
    int iterations{10000};
    double dx {0.01e-8};
    double dy {0.01e-8};
    double dt {dx / (2*c)};
    
    double source_sigma{10*dt};
    double t0 {6 * source_sigma};
    
    int pml_size{Nx / 20};
    
    int which_frame{0};
    int num_frame[3000];
    for(int i=0;i<3000;++i)
    {
        num_frame[i] = i * iterations/3000;
    }
    num_frame[0] = 1;
   
    
    int worker_Ny_base = (Ny-(2*pml_size)) / (size- 1);
    int worker_Ny_remainder = (Ny-(2*pml_size)) % (size-1);
    
    int source_center{Nx/2};
    int rank_occuring {(Nx/2) / worker_Ny_base};
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
        
//           for (int i=0; i<pml_size; ++i)
//           {
//                std::cout << DzX_coefs[i] << '\n';
//           }
        
        
        
        Point *sim_space_top = new Point[(pml_size+2)*Nx];
        Point *sim_space_bottom = new Point[(pml_size+2)*Nx];

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

        int above_me{size-1};
        int under_me{1};

        std::cout << "RANK: "<< rank << " ABOVE ME: "<< above_me << " UNDER ME: " << under_me << std::endl ;

        for (int t=1; t<iterations; ++t)
        {   


            MPI_Sendrecv(&sim_space_bottom[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t, &sim_space_top[index(pml_size+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t,
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space_top[index(pml_size,0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+iterations, &sim_space_bottom[0], Nx,
            MPI_POINT, above_me,rank+t+iterations,
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);
            
            update_H_master_top(sim_space_top, Nx, pml_size, mux, muy, dx, dy, HxY_coefs_top, HyY_coefs_top,ICHx_top, ICHy_top);
            
            update_H_master_bottom(sim_space_bottom, Nx, pml_size, mux, muy, dx, dy, HxY_coefs_bottom, HyY_coefs_bottom,ICHx_bottom, ICHy_bottom);

            MPI_Sendrecv(&sim_space_bottom[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t+(2*iterations), &sim_space_top[index(pml_size+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t+(2*iterations),
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space_top[index(pml_size,0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+(3*iterations), &sim_space_bottom[0], Nx,
            MPI_POINT, above_me,rank+t+(3*iterations),
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);

            update_E_master_top(sim_space_top, Nx, pml_size, ep,dx, dy, dt, DzY_coefs_top,IDz_top);
            update_E_master_bottom(sim_space_bottom, Nx, pml_size, ep,dx, dy, dt, DzY_coefs_bottom,IDz_bottom);


            
            if(t == num_frame[which_frame])
            {
            SaveToFile(Nx*pml_size, &sim_space_top[Nx], "data/"+std::to_string(rank)+"TOP.txt");
            SaveToFile(Nx*(pml_size), &sim_space_bottom[Nx], "data/"+std::to_string(rank)+"BOTTOM.txt");
            which_frame += 1;
            }


        }

            
            delete ep;
            delete mux ;
            delete muy;
            delete sim_space_top;
            delete sim_space_bottom;
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
        
        Point *sim_space = new Point[(local_Ny+2)*Nx];


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
        
        double *IHx = new double[2 * local_Ny * pml_size];
        double *IHy = new double[2 * local_Ny * pml_size];
        double *IDz = new double[2 * local_Ny * pml_size];

        double *ICHx = new double[2 * local_Ny * pml_size];
        double *ICHy = new double[2 * local_Ny * pml_size];
        double *ICDz = new double[2 * local_Ny * pml_size];


        
        int source_center_y {local_Ny/2};
        std::cout << "RANK: "<< rank << " ABOVE ME: "<< above_me << " UNDER ME: " << under_me <<" I think my y size is :  " <<local_Ny << "AND..." <<'\n';
       

        for (int t=1; t<iterations; ++t)
        {
            //std::cout <<  sim_space << IDz[0] <<'\n';        
            sim_space[index(source_center_y,source_center, Nx)].InjectEz(sourceE[t]);
            sim_space[index(source_center_y,source_center, Nx)].InjectDz(sourceE[t]);

            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
                    above_me+t, &sim_space[index(local_Ny+1,0,Nx)] , Nx,
                    MPI_POINT, under_me,rank +t,
                    MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(local_Ny,0,Nx)], Nx, MPI_POINT, under_me ,
                        under_me+t+iterations, &sim_space[0], Nx,
                        MPI_POINT, above_me,rank+t+iterations,
                        MPI_COMM_WORLD, &status);
            MPI_Barrier(MPI_COMM_WORLD);

            
            update_H_worker(sim_space, Nx, local_Ny, mux, muy, dx, dy, pml_size, HxX_coefs, HyX_coefs, IHx, IHy,ICHx, ICHy);
            
            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t+(2*iterations), &sim_space[index(local_Ny+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t+(2*iterations),
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(local_Ny,0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+(3*iterations), &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+(3*iterations),
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);

            update_E_worker(sim_space, Nx, local_Ny, ep,dx, dy, dt, pml_size, DzX_coefs, IDz, ICDz);

            if(t == num_frame[which_frame])
            {                
            SaveToFile(Nx*local_Ny, &sim_space[Nx], "data/"+std::to_string(rank)+".txt"); 
            which_frame+=1;
            }

        }

        delete ep;
        delete mux ;
        delete muy;
        delete sim_space;

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

        Point *sim_space = new Point[(local_Ny+2)*Nx];
        
        double *IHx = new double[2 * local_Ny * pml_size];
        double *IHy = new double[2 * local_Ny * pml_size];
        double *IDz = new double[2 * local_Ny * pml_size];

        double *ICHx = new double[2 * local_Ny * pml_size];
        double *ICHy = new double[2 * local_Ny * pml_size];
        double *ICDz = new double[2 * local_Ny * pml_size];

        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};

        std::cout << "RANK: "<< rank << " ABOVE ME: "<< above_me << " UNDER ME: " << under_me << " I think my y size is :  " <<local_Ny << '\n';


        for (int t=1; t<iterations; ++t)
        {    
//             if (t==90)
//             {
//             std::cout << "RANK: " << rank << " I sent over: " << '\n' << sim_space[index(1,100,Nx)];
//             }

            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t, &sim_space[index(local_Ny+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t,
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(local_Ny,0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+iterations, &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+iterations,
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);

            update_H_worker(sim_space, Nx, local_Ny, mux, muy, dx, dy, pml_size, HxX_coefs, HyX_coefs, IHx, IHy,ICHx, ICHy);
            
            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t+(2*iterations), &sim_space[index(local_Ny+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t+(2*iterations),
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(local_Ny,0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+(3*iterations), &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+(3*iterations),
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);
            
            
            update_E_worker(sim_space, Nx, local_Ny, ep,dx, dy, dt, pml_size, DzX_coefs, IDz, ICDz); 
            
//             if (t==90)
//             {
//             std::cout << "RANK: " << rank << " I recieved: " << '\n' << sim_space[index(local_Ny+1,100,Nx)];
//             }


            
            if(t == num_frame[which_frame])
            {            
            SaveToFile(Nx*local_Ny, &sim_space[Nx], "data/"+std::to_string(rank)+".txt");
            which_frame+=1;
            }

        }

        delete ep;
        delete mux ;
        delete muy;
        delete sim_space;

    }
    
    delete sigma_x;
    delete sigma_y;
    delete HxX_coefs;
    delete HyX_coefs;
    delete DzX_coefs;
    
    MPI_Finalize();    
    return 0;
}
