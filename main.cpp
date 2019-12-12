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
    int iterations{600};
    double dx {0.1e-8};
    double dy {0.1e-8};
    double dt {dx / (2*c)};
    
    double source_sigma{10*dt};
    double t0 {6 * source_sigma};
    
    int pml_size_x{Nx / 10};
    int pml_size_y{Ny / 10};    
    
    int worker_Ny_base = Ny / size ;
    int worker_Ny_remainder = Ny % size;
    
    int source_center{50};
    int rank_occuring {50 / worker_Ny_base};
    
    double *sigma_x = new double[pml_size_x];
    double *sigma_y = new double[pml_size_y];
    sigma_setter(sigma_x,sigma_y, pml_size_x, pml_size_y,dt);
    
    PML_coefs *HxX_coefs = new PML_coefs[pml_size_x];
    PML_coefs *HyX_coefs = new PML_coefs[pml_size_x];
    PML_coefs *DzX_coefs = new PML_coefs[pml_size_x];
    
    param_setter(HxX_coefs, HyX_coefs, DzX_coefs,sigma_x,sigma_y, 1, 1, dt, pml_size_x);

    PML_coefs *HxY_coefs = new PML_coefs[pml_size_y];
    PML_coefs *HyY_coefs = new PML_coefs[pml_size_y];
    PML_coefs *DzY_coefs = new PML_coefs[pml_size_y];
    
    param_setter(HxY_coefs, HyY_coefs, DzY_coefs,sigma_x,sigma_y, 1, 1, dt, pml_size_x);
    
    if (rank == 0)
    {
        
        for (int i=0; i<pml_size_x; ++i)
        {
             std::cout << DzX_coefs[i] << '\n' << "Sigma_x:  " <<sigma_x[i] << '\n';
        }
        
        
        int local_Ny = worker_Ny_base;
        int split_size{local_Ny / 2};
        Point *sim_space_top = new Point[split_size*Nx];
        Point *sim_space_bottom = new Point[(local_Ny-split_size)*Nx];

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
        
        double *IHx = new double[2 * local_Ny * pml_size_x];
        double *IHy = new double[2 * local_Ny * pml_size_y];
        double *ICHx = new double[2 * local_Ny * pml_size_x];
        double *ICHy = new double[2 * local_Ny * pml_size_y];


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
            
            for (int k=index(split_size-1, 0, Nx); k<index(split_size-1, Nx, Nx); ++k)
            {
                bottom[k-index(split_size-1, 0, Nx)] = sim_space_top[k];

            } //2e5 kai
            

            
            MPI_Sendrecv(top, Nx, MPI_POINT, above_me ,
                    t, row_from_bellow , Nx,
                    MPI_POINT, under_me,t,
                    MPI_COMM_WORLD, &status);

            MPI_Sendrecv(bottom, Nx, MPI_POINT, under_me ,
                        t+iterations, row_from_above, Nx,
                        MPI_POINT, above_me,t+iterations,
                        MPI_COMM_WORLD, &status);                

            //update_H_worker(sim_space_top, Nx, split_size, mux, muy, dx, dy,row_from_bellow);
            //update_E_master(sim_space_top, Nx, split_size, ep,dx, dy, dt);
            
            //void update_H_worker(Point space[], int Nx, int Ny, double mux[], double muy[], double dx, double dy, Point row_bellow[], int pml_size, PML_coefs coefs_Hx[], PML_coefs coefs_Hy[], double IHx[], double IHy[], double ICHx[], double ICHy[]);

            
            //update_H_master(sim_space_bottom, Nx, local_Ny-split_size, mux, muy, dx, dy);
            //update_E_worker(sim_space_bottom, Nx, local_Ny-split_size, ep,dx, dy, dt,row_from_above);
            
            SaveToFile(Nx*split_size, sim_space_top, "data/"+std::to_string(rank)+"TOP.txt");
            SaveToFile(Nx*(local_Ny-split_size), sim_space_bottom, "data/"+std::to_string(rank)+"BOTTOM.txt");
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
        if (rank <= worker_Ny_remainder)
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

        double *IHx = new double[2 * local_Ny * pml_size_x];
        double *IHy = new double[2 * local_Ny * pml_size_y];
        double *IDz = new double[2 * local_Ny * pml_size_y];

        double *ICHx = new double[2 * local_Ny * pml_size_x];
        double *ICHy = new double[2 * local_Ny * pml_size_y];
        double *ICDz = new double[2 * local_Ny * pml_size_y];


        
        int source_center_y {12};
        std::cout << "RANK: "<< rank << " ABOVE ME: "<< above_me << " UNDER ME: " << under_me << std::endl ;


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
            sim_space[index(source_center_y,source_center, Nx)].InjectEz(sourceE[t]);
            sim_space[index(source_center_y,source_center, Nx)].InjectDz(sourceE[t]);

            MPI_Sendrecv(top, Nx, MPI_POINT, above_me ,
                        t, row_from_bellow , Nx,
                        MPI_POINT, under_me,t,
                        MPI_COMM_WORLD, &status);

            MPI_Sendrecv(bottom, Nx, MPI_POINT, under_me ,
                        t+iterations, row_from_above, Nx,
                        MPI_POINT, above_me,t+iterations,
                        MPI_COMM_WORLD, &status);        

            update_H_worker(sim_space, Nx, local_Ny, mux, muy, dx, dy, row_from_bellow, pml_size_x, HxX_coefs, HyX_coefs, IHx, IHy,ICHx, ICHy);
            //void update_E_worker(Point space[], int Nx, int Ny, double ep[], double dx, double dy, double dt, Point row_above[],int pml_size, PML_coefs coefs_Dz[], double IDz[], double ICDz[]);
            update_E_worker(sim_space, Nx, local_Ny, ep,dx, dy, dt, row_from_above, pml_size_x, DzX_coefs, IDz, ICDz); 
                
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
        if (rank <= worker_Ny_remainder)
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
        
        double *IHx = new double[2 * local_Ny * pml_size_x];
        double *IHy = new double[2 * local_Ny * pml_size_y];
        double *IDz = new double[2 * local_Ny * pml_size_y];

        double *ICHx = new double[2 * local_Ny * pml_size_x];
        double *ICHy = new double[2 * local_Ny * pml_size_y];
        double *ICDz = new double[2 * local_Ny * pml_size_y];

                

        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};

        std::cout << "RANK: "<< rank << " ABOVE ME: "<< above_me << " UNDER ME: " << under_me << std::endl ;

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

            update_H_worker(sim_space, Nx, local_Ny, mux, muy, dx, dy, row_from_bellow, pml_size_x, HxX_coefs, HyX_coefs, IHx, IHy,ICHx, ICHy);
            update_E_worker(sim_space, Nx, local_Ny, ep,dx, dy, dt, row_from_above, pml_size_x, DzX_coefs, IDz, ICDz); 
            
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
    delete HxY_coefs;
    delete HyY_coefs;
    delete DzY_coefs;
    
    MPI_Finalize();    
    return 0;
}
