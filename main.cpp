#include <Ez_point.h>
#include <soft_source.h>
#include <PML_boundry.h>
#include <update_func_PML.h>
#include <objects.h>
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

int rank_index_y(int i, int Ny, int *sizes, int rank_size)
{
    int rank_occuring {}; // [rank, Ny]
    for(int k=0; k<rank_size; ++k)
    {
        if (i-(sizes[k]-1) > 0)
        {
            i -= sizes[k]-1; 
        }
        else
        {
            rank_occuring = k;
            break;
        }
    }
    return rank_occuring;
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
    
    //Data saving folder:
    std::string savefile_path{"data"};
    
    //INITIAL PARAMS:   
    int Nx{200};
    int Ny{200};
    int iterations{1000};
    double lambda_min{200e-9};
    
    double dx {lambda_min / 20}; // 10 points per half wavelength
    double dy {lambda_min / 20};
    double dt {dx / (2*c)}; // best dispersion to resolution
    int pml_size{Nx / 5};
    int pml_ratio{3};
    
    // source parameters (TEMP!!!!!!)
    double source_sigma{10*dt};
    double t0 {6 * source_sigma};
    
    int frames_to_save{500};
    int which_frame{0};
    int num_frame[frames_to_save];
    for(int i=0;i<frames_to_save;++i)
    {
        num_frame[i] = i * iterations/frames_to_save;
    }
    num_frame[0] = 1;
   
    
//     int Ny_size_list[size]{};
    
    double density_pml = static_cast<double>(pml_ratio) 
                            * static_cast<double>(size) 
                            / (static_cast<double>(Ny) - ((1-static_cast<double>(pml_ratio)) 
                            * static_cast<double>(pml_size)));
    int n_pml = density_pml * pml_size;
    if (n_pml == 0)
    {
        n_pml = 1; // Cant have 0 processes in pml!!!
    }
    int n_bulk = size - (2* n_pml);
    
    int pml_base_size = pml_size / n_pml;
    int pml_remainder = pml_size % n_pml;
    
    int worker_Ny_base = (Ny-(2*pml_size)) / n_bulk;
    int worker_Ny_remainder = (Ny-(2*pml_size)) % n_bulk;
    

    int sizes[size]{};
    int c_sizes[size]{};
    
    for(int i=0; i<n_pml; i++)
    {
        if (i < pml_remainder)
        {
            sizes[i] = pml_base_size + 1;
        }
        else
        {
            sizes[i] = pml_base_size;
        }
    }
    for(int i=n_pml; i<size-n_pml; i++)
    {
        if (i-(n_pml) < worker_Ny_remainder)
        {
            sizes[i] = worker_Ny_base + 1;
        }
        else
        {
            sizes[i] = worker_Ny_base;
        }
    }
    for(int i=size-n_pml; i<size; i++)
    {
        if (i-(size-n_pml) < pml_remainder)
        {
            sizes[i] = pml_base_size + 1;
        }
        else
        {
            sizes[i] = pml_base_size;
        }
    }
    
    c_sizes[0] = sizes[0];
    for(int i=1; i<size; ++i)
    {
        c_sizes[i] = sizes[i] + c_sizes[i-1];
    }
    
    
    int source_position_x {Nx/2};
    int source_position_y {Ny/5 + 10};
    int rank_source =rank_index_y(source_position_y, Ny, sizes, size);
    
    int center_x {Nx / 2};
    int center_y {(Ny / 2)};
    int radius_fo_curvature{Ny /4 }; // why not...
    int width{Ny /3};
    
    double *global_ep = new double[Nx*Ny];
    lense_foc(global_ep,Nx,Ny,center_x, center_y, radius_fo_curvature, width);
    
    std::fstream fs1;
    fs1.open("data/epmap.txt", std::fstream::in | std::fstream::out | std::fstream::trunc);

    for (int i=0; i<Nx*Ny; ++i)
    {
        fs1 << global_ep[i] << '\n';
    }
    fs1.close();

    
    
    double * sigma_x = new double[pml_size];
    double *sigma_y = new double[pml_size];
    sigma_setter(sigma_x,sigma_y, pml_size,dt);
    
    PML_coefs *HxX_coefs = new PML_coefs[pml_size];
    PML_coefs *HyX_coefs = new PML_coefs[pml_size];
    PML_coefs *DzX_coefs = new PML_coefs[pml_size];
    
    param_setter_x(HxX_coefs, HyX_coefs, DzX_coefs,sigma_x, 1, 1, dt, pml_size);
    

    if (rank == 0)
    {
        
        double *local_sigma_y = new double[sizes[rank]];
        for (int i=0; i<sizes[rank]; i++)
        {
            local_sigma_y[i] = sigma_x[pml_size-(i+1)];
            std::cout << "RANK: "<< rank  << " I tried to access element " << pml_size-(i+1) << '\n';
        }
        
        Point *sim_space = new Point[(sizes[rank]+2)*Nx];

        double *ep = new double[sizes[rank]*Nx];
        double *mux = new double[sizes[rank]*Nx];
        double *muy = new double[sizes[rank]*Nx];


        for (int i=0; i<sizes[rank]; ++i)
        {
            for (int j=0; j<Nx; ++j)
            {
                ep[index(i,j,Nx)] = global_ep[index(i,j,Nx)];
                mux[index(i,j,Nx)] = 1;
                muy[index(i,j,Nx)] = 1;
            }
        }
        
        PML_coefs *HxY_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *HyY_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *DzY_coefs = new PML_coefs[sizes[rank]*Nx];
        param_setter(HxY_coefs,HyY_coefs,DzY_coefs,sigma_x, local_sigma_y,1,1,dt,pml_size,sizes[rank], Nx);

        
        double *IDz = new double[Nx * sizes[rank]];
        double *ICHx = new double[Nx * sizes[rank]];
        double *ICHy = new double[Nx * sizes[rank]];
        

        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};

        std::cout << "RANK: "<< rank << " ABOVE ME: "<< above_me << " UNDER ME: " << under_me << " I think my y size is :  " <<sizes[rank] << '\n';

        for (int t=1; t<iterations; ++t)
        {   

            MPI_Recv(&sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t, MPI_COMM_WORLD, &status);
            MPI_Send(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+iterations, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            
            update_H_pml(sim_space, Nx, sizes[rank], mux, muy, dx, dy, HxY_coefs, HyY_coefs,ICHx, ICHy);
            

            MPI_Recv(&sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t+(2*iterations), MPI_COMM_WORLD, &status);
            MPI_Send(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+(3*iterations), MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            update_E_pml(sim_space, Nx, sizes[rank], ep,dx, dy, dt, DzY_coefs,IDz);

            
            if(t == num_frame[which_frame])
            {
            SaveToFile(Nx*sizes[rank], &sim_space[Nx], "data/"+std::to_string(rank)+".txt");
            which_frame += 1;
            }

        }

            
        delete ep;
        delete mux ;
        delete muy;
        delete sim_space;

        std::fstream fs;
        fs.open("data/info.csv", std::fstream::in | std::fstream::out | std::fstream::trunc);
        fs << "Nx,"<<Nx<<'\n'
        << "Ny,"<<Ny<<'\n'
        << "Frames_saved,"<<frames_to_save<<'\n'
        << "Ranks,"<<size<<'\n'
        << "PML_size,"<< pml_size << '\n';
        
        for(int i=0; i<size; ++i)
        {
            fs << "RANK " << i << " size," << sizes[i] << '\n';
        }
        fs.close();



    }
    else if (rank < n_pml)
    {
        double *local_sigma_y = new double[sizes[rank]];
        for (int i=0; i<sizes[rank]; i++)
        {
            local_sigma_y[i] = sigma_x[pml_size-(i+c_sizes[rank-1]+1)];
            std::cout << "RANK: "<< rank  << " Sigma is " << local_sigma_y[i] << '\n';

        }
        
        Point *sim_space = new Point[(sizes[rank]+2)*Nx];

        double *ep = new double[sizes[rank]*Nx];
        double *mux = new double[sizes[rank]*Nx];
        double *muy = new double[sizes[rank]*Nx];


        for (int i=0; i<sizes[rank]; ++i)
        {
            for (int j=0; j<Nx; ++j)
            {
                ep[index(i,j,Nx)] = global_ep[index(i+c_sizes[rank-1],j,Nx)];
                mux[index(i,j,Nx)] = 1;
                muy[index(i,j,Nx)] = 1;
            }
        }
        
        PML_coefs *HxY_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *HyY_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *DzY_coefs = new PML_coefs[sizes[rank]*Nx];
        param_setter(HxY_coefs,HyY_coefs,DzY_coefs,sigma_x, local_sigma_y,1,1,dt,pml_size,sizes[rank], Nx);

        
        double *IDz = new double[Nx * sizes[rank]];
        double *ICHx = new double[Nx * sizes[rank]];
        double *ICHy = new double[Nx * sizes[rank]];
        
        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};

        std::cout << "RANK: "<< rank << " ABOVE ME: "<< above_me << " UNDER ME: " << under_me << " I think my y size is :  " <<sizes[rank] << '\n';

        for (int t=1; t<iterations; ++t)
        {   

            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t, &sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t,
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+iterations, &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+iterations,
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);
            
            update_H_pml(sim_space, Nx, sizes[rank], mux, muy, dx, dy, HxY_coefs, HyY_coefs,ICHx, ICHy);
            

            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t+(2*iterations), &sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t+(2*iterations),
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+(3*iterations), &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+(3*iterations),
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);

            update_E_pml(sim_space, Nx, sizes[rank], ep,dx, dy, dt, DzY_coefs,IDz);

            
            if(t == num_frame[which_frame])
            {
            SaveToFile(Nx*sizes[rank], &sim_space[Nx], "data/"+std::to_string(rank)+".txt");
            which_frame += 1;
            }


        }

            
            delete ep;
            delete mux ;
            delete muy;
            delete sim_space;

    }
    
    else if (rank == rank_source)
    {
        Point *sim_space = new Point[(sizes[rank]+2)*Nx];
        double *ep = new double[sizes[rank]*Nx];
        double *mux = new double[sizes[rank]*Nx];
        double *muy = new double[sizes[rank]*Nx];

        for (int j=0; j<Nx; ++j)
        {
            for (int i=0; i<sizes[rank]; ++i)
            {
                ep[index(i,j,Nx)] = global_ep[index(i+c_sizes[rank-1],j,Nx)];
                mux[index(i,j,Nx)] = c * dt;
                muy[index(i,j,Nx)] = c * dt;
            }
        }
        
        double *IHx = new double[2 * sizes[rank] * pml_size];
        double *IHy = new double[2 * sizes[rank] * pml_size];
        double *IDz = new double[2 * sizes[rank] * pml_size];

        double *ICHx = new double[2 * sizes[rank] * pml_size];
        double *ICHy = new double[2 * sizes[rank] * pml_size];
        double *ICDz = new double[2 * sizes[rank] * pml_size];

        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};

        std::cout << "RANK: "<< rank << " ABOVE ME: "<< above_me << " UNDER ME: " << under_me << " I think my y size is :  " <<sizes[rank] << '\n';
        
        double *sourceE = new double[iterations];
        double *sourceHx = new double[iterations];
        double *sourceHy = new double[iterations];
        inject_soft_source2(sourceE, sourceHx, sourceHy,iterations,dx,dy,dt, 1, 1, 1,1, source_sigma, t0);
        


        for (int t=1; t<iterations; ++t)
        {    
            for (int x=0; x<Nx ; ++x)
            {
            sim_space[index(sizes[rank] / 2, x, Nx)].InjectEz(sourceE[t]);
            sim_space[index(sizes[rank] / 2, x, Nx)].InjectDz(sourceE[t]);
            }
//             sim_space[index(sizes[rank] / 2, Nx/2, Nx)].InjectEz(sourceE[t]);
//             sim_space[index(sizes[rank] / 2, Nx/2, Nx)].InjectDz(sourceE[t]);

            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t, &sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t,
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+iterations, &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+iterations,
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);

            update_H_bulk(sim_space, Nx, sizes[rank], mux, muy, dx, dy, pml_size, HxX_coefs, HyX_coefs, IHx, IHy,ICHx, ICHy);
            
            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t+(2*iterations), &sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t+(2*iterations),
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+(3*iterations), &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+(3*iterations),
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);
            
            
            update_E_bulk(sim_space, Nx, sizes[rank], ep,dx, dy, dt, pml_size, DzX_coefs, IDz, ICDz); 
        
            if(t == num_frame[which_frame])
            {            
            SaveToFile(Nx*sizes[rank], &sim_space[Nx], "data/"+std::to_string(rank)+".txt");
            which_frame+=1;
            }

        }

        delete ep;
        delete mux ;
        delete muy;
        delete sim_space;
    }
    
    else if (rank == (size-1))
    {
        double *local_sigma_y = new double[sizes[rank]];

        if (rank == size - n_pml)
        {
            for (int i=0; i<sizes[rank]; i++)
            {
                local_sigma_y[i] = sigma_x[i];
                std::cout << "RANK: "<< rank  << " I tried to access element " << i << '\n';

            }
        }
        
        else
        {
            for (int i=0; i<sizes[rank]; i++)
            {
                local_sigma_y[i] = sigma_x[(c_sizes[rank] - c_sizes[size - n_pml]) +(i+1)];
            }

        }
        
        Point *sim_space = new Point[(sizes[rank]+2)*Nx];

        double *ep = new double[sizes[rank]*Nx];
        double *mux = new double[sizes[rank]*Nx];
        double *muy = new double[sizes[rank]*Nx];


        for (int i=0; i<sizes[rank]; ++i)
        {
            for (int j=0; j<Nx; ++j)
            {
                ep[index(i,j,Nx)] = global_ep[index(i+c_sizes[rank-1],j,Nx)];
                mux[index(i,j,Nx)] = 1;
                muy[index(i,j,Nx)] = 1;
            }
        }
        
        PML_coefs *HxY_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *HyY_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *DzY_coefs = new PML_coefs[sizes[rank]*Nx];
        param_setter(HxY_coefs,HyY_coefs,DzY_coefs,sigma_x, local_sigma_y,1,1,dt,pml_size,sizes[rank], Nx);

        
        double *IDz = new double[Nx * sizes[rank]];
        double *ICHx = new double[Nx * sizes[rank]];
        double *ICHy = new double[Nx * sizes[rank]];
        

        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};

        std::cout << "RANK: "<< rank << " ABOVE ME: "<< above_me << " UNDER ME: " << under_me << " I think my y size is :  " <<sizes[rank] << '\n';

        for (int t=1; t<iterations; ++t)
        {   

            MPI_Send(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t, MPI_COMM_WORLD);
            MPI_Recv( &sim_space[0], Nx, MPI_POINT, above_me,rank+t+iterations,
            MPI_COMM_WORLD, &status);
            MPI_Barrier(MPI_COMM_WORLD);
            
            update_H_pml(sim_space, Nx, sizes[rank], mux, muy, dx, dy, HxY_coefs, HyY_coefs,ICHx, ICHy);
            

            MPI_Send(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t+(2*iterations), MPI_COMM_WORLD);
            
            MPI_Recv(&sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+(3*iterations),
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);

            update_E_pml(sim_space, Nx, sizes[rank], ep,dx, dy, dt, DzY_coefs,IDz);

            
            if(t == num_frame[which_frame])
            {
            SaveToFile(Nx*sizes[rank], &sim_space[Nx], "data/"+std::to_string(rank)+".txt");
            which_frame += 1;
            }
        }
    }
    
    else if (rank >= (size - n_pml))
    {
        double *local_sigma_y = new double[sizes[rank]];

        if (rank == size - n_pml)
        {
            for (int i=0; i<sizes[rank]; i++)
            {
                local_sigma_y[i] = sigma_x[i];
            }
        }
        
        else
        {
            for (int i=0; i<sizes[rank]; i++)
            {
                local_sigma_y[i] = sigma_x[(c_sizes[rank] - c_sizes[size - n_pml]) +(i+1)];
            }
        }
        
        Point *sim_space = new Point[(sizes[rank]+2)*Nx];

        double *ep = new double[sizes[rank]*Nx];
        double *mux = new double[sizes[rank]*Nx];
        double *muy = new double[sizes[rank]*Nx];


        for (int i=0; i<sizes[rank]; ++i)
        {
            for (int j=0; j<Nx; ++j)
            {
                ep[index(i,j,Nx)] = global_ep[index(i+c_sizes[rank-1],j,Nx)];
                mux[index(i,j,Nx)] = 1;
                muy[index(i,j,Nx)] = 1;
            }
        }
        
        PML_coefs *HxY_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *HyY_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *DzY_coefs = new PML_coefs[sizes[rank]*Nx];
        param_setter(HxY_coefs,HyY_coefs,DzY_coefs,sigma_x, local_sigma_y,1,1,dt,pml_size,sizes[rank], Nx);

        
        double *IDz = new double[Nx * sizes[rank]];
        double *ICHx = new double[Nx * sizes[rank]];
        double *ICHy = new double[Nx * sizes[rank]];
        

        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};

        std::cout << "RANK: "<< rank << " ABOVE ME: "<< above_me << " UNDER ME: " << under_me << " I think my y size is :  " <<sizes[rank] << '\n';

        for (int t=1; t<iterations; ++t)
        {   

            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t, &sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t,
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+iterations, &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+iterations,
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);
            
            update_H_pml(sim_space, Nx, sizes[rank], mux, muy, dx, dy, HxY_coefs, HyY_coefs,ICHx, ICHy);
            

            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t+(2*iterations), &sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t+(2*iterations),
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+(3*iterations), &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+(3*iterations),
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);

            update_E_pml(sim_space, Nx, sizes[rank], ep,dx, dy, dt, DzY_coefs,IDz);

            
            if(t == num_frame[which_frame])
            {
            SaveToFile(Nx*sizes[rank], &sim_space[Nx], "data/"+std::to_string(rank)+".txt");
            which_frame += 1;
            }


        }

            
            delete ep;
            delete mux ;
            delete muy;
            delete sim_space;

        
    }

    else
    {
        Point *sim_space = new Point[(sizes[rank]+2)*Nx];
        double *ep = new double[sizes[rank]*Nx];
        double *mux = new double[sizes[rank]*Nx];
        double *muy = new double[sizes[rank]*Nx];

        for (int j=0; j<Nx; ++j)
        {
            for (int i=0; i<sizes[rank]; ++i)
            {
                ep[index(i,j,Nx)] = global_ep[index(i+c_sizes[rank-1],j,Nx)];
                mux[index(i,j,Nx)] = c * dt;
                muy[index(i,j,Nx)] = c * dt;
            }
        }
        
        double *IHx = new double[2 * sizes[rank] * pml_size];
        double *IHy = new double[2 * sizes[rank] * pml_size];
        double *IDz = new double[2 * sizes[rank] * pml_size];

        double *ICHx = new double[2 * sizes[rank] * pml_size];
        double *ICHy = new double[2 * sizes[rank] * pml_size];
        double *ICDz = new double[2 * sizes[rank] * pml_size];

        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};

        std::cout << "RANK: "<< rank << " ABOVE ME: "<< above_me << " UNDER ME: " << under_me << " I think my y size is :  " <<sizes[rank] << '\n';


        for (int t=1; t<iterations; ++t)
        {    

            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t, &sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t,
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+iterations, &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+iterations,
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);

            update_H_bulk(sim_space, Nx, sizes[rank], mux, muy, dx, dy, pml_size, HxX_coefs, HyX_coefs, IHx, IHy,ICHx, ICHy);
            
            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t+(2*iterations), &sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t+(2*iterations),
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+(3*iterations), &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+(3*iterations),
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);
            
            
            update_E_bulk(sim_space, Nx, sizes[rank], ep,dx, dy, dt, pml_size, DzX_coefs, IDz, ICDz); 
        
            if(t == num_frame[which_frame])
            {            
            SaveToFile(Nx*sizes[rank], &sim_space[Nx], "data/"+std::to_string(rank)+".txt");
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
    delete global_ep;
    
    MPI_Finalize();    
    return 0;
}
