#include <Ez_point.h>
#include <soft_source.h>
#include <PML_boundry.h>
#include <update_func_PML.h>
#include <objects.h>
#include<core_funcs.h>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include<stdio.h>

// #define PLANEWAVE //PLANEWAVE or POINTSOURCE
// #define SAVEFRAMES

int main(int argc, char *argv[])
{
    // Initilise and find 
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
    
////////////////////////////////////////////////////////////////////////////  
    
    //INITIAL PARAMS:                                                     
    double gridsize        {2.5e-6}; // Just gonna use square grids 
    double time            {100e-15};
    double lambda_min      {200e-9};                                            
    double pml_size_ratio  {10};                                            
    int pml_ratio          {3};                                                     
    int frames_to_save     {300};  
    
    
    double source_position[2]  {gridsize/2,gridsize/4.8};// {x,y}
    int objecttype             {2}; // 1=planar_convex_lens, 2=biconvel lens, 3 = double_slit
    double object_position[2]  {gridsize/2,3*gridsize/10}; // {x,y}
    
    // for lens:
    double rad_of_curvature   {2e-6};
    double lens_width         {1.5e-6};
    // for double slit:
    double slit_width         {50e-9};
    double slit_seperation    {500e-9};

////////////////////////////////////////////////////////////////////////////  
    
    
    // Calculated parameters:
    double dx {lambda_min / 20}; // at least 10 points per half wavelength
    double dy {lambda_min / 20};
    double dt {dx / (2*c)}; // best dispersion to resolution
    int iterations {950};    

    int Nx{gridsize / dx};
    int Ny{gridsize / dy};
    int pml_size{Nx / pml_size_ratio};
    
    // Creating the object:
    int center_x {object_position[0]/dx};
    int center_y {object_position[1]/dy};
    int radius_of_curvature_pix{rad_of_curvature/dx}; 
    int slit_seperation_pix{slit_seperation/dx};
    int width{};
    
    // makes an array with permitivity
    // each process chops from this 
    double *global_ep = new double[Nx*Ny];
    if (objecttype == 1)
    {
        width = lens_width / dx;
        lense_foc(global_ep,Nx,Ny,center_x, center_y, radius_of_curvature_pix, width);
    }
    else if (objecttype == 2)
    {
        width = lens_width / dx;
        lense_foc(global_ep,Nx,Ny,center_x, center_y, radius_of_curvature_pix, width);
    }
    else
    {
        width = slit_width / dx;
        double_slit(global_ep,Nx,Ny, center_y, center_x,slit_seperation_pix,width);
    }
    

    // creates graded sigma profile for PML boundry:
    double * sigma_x = new double[pml_size];
    double *sigma_y = new double[pml_size];
    sigma_setter(sigma_x,sigma_y, pml_size,dt);
        
    
    // for saving of frames: 
    int frame_interval = (frames_to_save + iterations - (iterations%frames_to_save))/frames_to_save;
    
    int which_frame{0};
    int num_frame[frames_to_save];
    for(int i=0;i<frames_to_save+1;++i)
    {
        num_frame[i] = i * frame_interval;
    }
    num_frame[0] = 1; // 0th iteration is a bit useless
    
   
    // Logic to figure out how to split grid between jobs: 
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
    

    int sizes[size]{}; // sizes of each rank
    int c_sizes[size]{}; // cumilative sizes 
    
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
    
    // cumilative sizes: 
    c_sizes[0] = sizes[0];
    for(int i=1; i<size; ++i)
    {
        c_sizes[i] = sizes[i] + c_sizes[i-1]; // pupulates cumilative size recursively
    }
    
    int source_position_x {source_position[0]/dx};
    int source_position_y {source_position[1]/dy};
    int rank_source =rank_index_y(source_position_y, Ny, sizes, size);

    
    if (rank == 0)
    {
        double start_time = MPI_Wtime();
        // creates local simulation space, with extra row above_me
        // and bellow for data from other ranks to be pushed into
        Point *sim_space = new Point[(sizes[rank]+2)*Nx];
        
        // creates conductivity profile for PML
        double *local_sigma_y = new double[sizes[rank]];
        for (int i=0; i<sizes[rank]; i++)
        {
            local_sigma_y[i] = sigma_x[pml_size-(i+1)];
        }
        
        // Chops out part of permitivity grid relavant to itself:
        double *ep = new double[sizes[rank]*Nx];
        for (int j=0; j<Nx; ++j)
        {
            for (int i=0; i<sizes[rank]; ++i)
            {
                ep[index(i,j,Nx)] = global_ep[index(i+c_sizes[rank-1],j,Nx)];
            }
        }


        // coficient matrix (taking into account ep, mu and pml:
        PML_coefs *Hx_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *Hy_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *Dz_coefs = new PML_coefs[sizes[rank]*Nx];
        param_setter(Hx_coefs,Hy_coefs,Dz_coefs,sigma_x, local_sigma_y,1,1,dt,pml_size,sizes[rank], Nx);
        
        // Holds integrated values:
        double *IDz = new double[Nx * sizes[rank]];
        double *ICHx = new double[Nx * sizes[rank]];
        double *ICHy = new double[Nx * sizes[rank]];
        
        // which ranks to communicate with:
        int under_me{(rank+1)%size};

        for (int t=1; t<iterations; ++t)
        {   

            MPI_Recv(&sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t, MPI_COMM_WORLD, &status);
            MPI_Send(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+iterations, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            
            update_H_pml(sim_space, Nx, sizes[rank], dx, dy, Hx_coefs, Hy_coefs,ICHx, ICHy);

            MPI_Recv(&sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t+(2*iterations), MPI_COMM_WORLD, &status);
            MPI_Send(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+(3*iterations), MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            update_E_pml(sim_space, Nx, sizes[rank], ep,dx, dy, dt, Dz_coefs,IDz);
            
            if(t == num_frame[which_frame])
            {
            SaveToFile(Nx*sizes[rank], &sim_space[Nx], "data/"+std::to_string(rank)+".txt");
            which_frame += 1;
            }

        }
        
        double endtime = MPI_Wtime();
        
        //save file with object map:
        std::fstream fs1;
        fs1.open("data/epmap.txt", std::fstream::in | std::fstream::out | std::fstream::trunc);
        for (int i=0; i<Nx*Ny; ++i)
        {
            fs1 << global_ep[i] << '\n';
        }
        fs1.close();
        
        // saves information on simulation parameters:
        std::fstream fs;
        fs.open("data/info.csv", std::fstream::in | std::fstream::out | std::fstream::trunc);
        fs << "Nx,"              << Nx<<'\n'
           << "Ny,"              << Ny<<'\n'
           << "Frames_saved,"    << frames_to_save<<'\n'
           << "Ranks,"           << size<<'\n'
           << "PML_size,"        << pml_size << '\n'
           << "Iterations,"      << iterations<<'\n'
           << "size,"            << Nx<<'\n'
           << "Time taken (s),"  << endtime - start_time << '\n';
        
        for(int i=0; i<size; ++i)
        {
            fs << "RANK " << i << " size," << sizes[i] << '\n';
        }
        fs.close();

        delete sim_space;
        delete local_sigma_y;
        delete ep;
        delete Hx_coefs;
        delete Hy_coefs;
        delete Dz_coefs;
        delete IDz;
        delete ICHx;
        delete ICHy;


    }
    else if (rank < n_pml)
    {
        // creates local simulation space, with extra row above_me
        // and bellow for data from other ranks to be pushed into
        Point *sim_space = new Point[(sizes[rank]+2)*Nx];
        
        // PML conductivity:
        double *local_sigma_y = new double[sizes[rank]];
        for (int i=0; i<sizes[rank]; i++)
        {
            local_sigma_y[i] = sigma_x[pml_size-(i+c_sizes[rank-1]+1)];
        }
        
        // Chops out part of permitivity grid relavant to itself:
        double *ep = new double[sizes[rank]*Nx];
        for (int j=0; j<Nx; ++j)
        {
            for (int i=0; i<sizes[rank]; ++i)
            {
                ep[index(i,j,Nx)] = global_ep[index(i+c_sizes[rank-1],j,Nx)];
            }
        }
                
        // coficient matrix (taking into account ep, mu and pml:
        PML_coefs *Hx_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *Hy_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *Dz_coefs = new PML_coefs[sizes[rank]*Nx];
        param_setter(Hx_coefs,Hy_coefs,Dz_coefs,sigma_x, local_sigma_y,1,1,dt,pml_size,sizes[rank], Nx);

        // Holds integrated values:
        double *IDz = new double[Nx * sizes[rank]];
        double *ICHx = new double[Nx * sizes[rank]];
        double *ICHy = new double[Nx * sizes[rank]];
        
        // which ranks to communicate with:
        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};

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
            
            update_H_pml(sim_space, Nx, sizes[rank], dx, dy, Hx_coefs, Hy_coefs,ICHx, ICHy);
            

            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t+(2*iterations), &sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t+(2*iterations),
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+(3*iterations), &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+(3*iterations),
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);

            update_E_pml(sim_space, Nx, sizes[rank], ep,dx, dy, dt, Dz_coefs,IDz);

            if(t == num_frame[which_frame])
            {
            SaveToFile(Nx*sizes[rank], &sim_space[Nx], "data/"+std::to_string(rank)+".txt");
            which_frame += 1;
            }
        }

        delete sim_space;
        delete local_sigma_y;
        delete ep;
        delete Hx_coefs;
        delete Hy_coefs;
        delete Dz_coefs;
        delete IDz;
        delete ICHx;
        delete ICHy;
    }
    
    else if (rank == rank_source) // asuming source is not literally inside the PML
    {
        // source parameters: 
        double source_sigma{10*dt};
        double t0 {6 * source_sigma};
        
        // creates local simulation space, with extra row above_me
        // and bellow for data from other ranks to be pushed into
        Point *sim_space = new Point[(sizes[rank]+2)*Nx];
        
        // Chops out part of permitivity grid relavant to itself: 
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
        
        PML_coefs *Hx_coefs = new PML_coefs[pml_size];
        PML_coefs *Hy_coefs = new PML_coefs[pml_size];
        PML_coefs *Dz_coefs = new PML_coefs[pml_size];
        param_setter_x(Hx_coefs, Hy_coefs, Dz_coefs,sigma_x, 1, 1, dt, pml_size);

        // Holds integrated values:
        double *IDz = new double[2 * sizes[rank] * pml_size];
        double *ICHx = new double[2 * sizes[rank] * pml_size];
        double *ICHy = new double[2 * sizes[rank] * pml_size];
        
        // which ranks to communicate with:
        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};
        
        //creates source array in advanced to save compute
        double *sourceE = new double[iterations];
        double *sourceHx = new double[iterations];
        double *sourceHy = new double[iterations];
        inject_soft_source(sourceE, sourceHx, sourceHy,iterations,dx,dy,dt, 1, 1, 1,1, source_sigma, t0);

        for (int t=1; t<iterations; ++t)
        {    

            for (int x=pml_size; x<Nx-pml_size ; ++x)
            {
            sim_space[index(source_position_y-c_sizes[rank-1]+1, x, Nx)].InjectEz(sourceE[t]);
            sim_space[index(source_position_y-c_sizes[rank-1]+1, x, Nx)].InjectDz(sourceE[t]);
            }
//             sim_space[index(source_position_y-c_sizes[rank-1]+1, source_position_x, Nx)].InjectEz(sourceE[t]);
//             sim_space[index(source_position_y-c_sizes[rank-1]+1, source_position_x, Nx)].InjectDz(sourceE[t]);


            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t, &sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t,
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+iterations, &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+iterations,
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);

            update_H_bulk(sim_space, Nx, sizes[rank], mux, muy, dx, dy, pml_size, Hx_coefs, Hy_coefs, ICHx, ICHy);

            
            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t+(2*iterations), &sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t+(2*iterations),
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+(3*iterations), &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+(3*iterations),
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);
            
            
            update_E_bulk(sim_space, Nx, sizes[rank], ep,dx, dy, dt, pml_size, Dz_coefs, IDz); 
            
            if(t == num_frame[which_frame])
            {
            SaveToFile(Nx*sizes[rank], &sim_space[Nx], "data/"+std::to_string(rank)+".txt");
            which_frame += 1;
            }

            
        }
        
        delete sim_space;
        delete ep;
        delete mux;
        delete muy;
        delete Hx_coefs;
        delete Hy_coefs;
        delete Dz_coefs;
        delete IDz;
        delete ICHx;
        delete ICHy;

    }
    
    else if (rank == (size-1))
    {
        // creates local simulation space, with extra row above_me
        // and bellow for data from other ranks to be pushed into
        Point *sim_space = new Point[(sizes[rank]+2)*Nx];
        
        // Creates PML sigma coficients
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
        
        // Chops out part of permitivity grid relavant to itself:
        double *ep = new double[sizes[rank]*Nx];
        for (int j=0; j<Nx; ++j)
        {
            for (int i=0; i<sizes[rank]; ++i)
            {
                ep[index(i,j,Nx)] = global_ep[index(i+c_sizes[rank-1],j,Nx)];
            }
        }
        
        PML_coefs *Hx_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *Hy_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *Dz_coefs = new PML_coefs[sizes[rank]*Nx];
        param_setter(Hx_coefs,Hy_coefs,Dz_coefs,sigma_x, local_sigma_y,1,1,dt,pml_size,sizes[rank], Nx);

        // Holds integrated values:
        double *IDz = new double[Nx * sizes[rank]];
        double *ICHx = new double[Nx * sizes[rank]];
        double *ICHy = new double[Nx * sizes[rank]];
        
        // which ranks to communicate with:
        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};

        for (int t=1; t<iterations; ++t)
        {   

            MPI_Send(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t, MPI_COMM_WORLD);
            MPI_Recv( &sim_space[0], Nx, MPI_POINT, above_me,rank+t+iterations,
            MPI_COMM_WORLD, &status);
            MPI_Barrier(MPI_COMM_WORLD);
            
            update_H_pml(sim_space, Nx, sizes[rank], dx, dy, Hx_coefs, Hy_coefs,ICHx, ICHy);
            

            MPI_Send(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t+(2*iterations), MPI_COMM_WORLD);
            
            MPI_Recv(&sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+(3*iterations),
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);

            update_E_pml(sim_space, Nx, sizes[rank], ep,dx, dy, dt, Dz_coefs,IDz);

            
            if(t == num_frame[which_frame])
            {
            SaveToFile(Nx*sizes[rank], &sim_space[Nx], "data/"+std::to_string(rank)+".txt");
            which_frame += 1;
            }
        }
        delete sim_space;
        delete local_sigma_y;
        delete ep;
        delete Hx_coefs;
        delete Hy_coefs;
        delete Dz_coefs;
        delete IDz;
        delete ICHx;
        delete ICHy;

    }
    
    else if (rank >= (size - n_pml))
    {
        // creates local simulation space, with extra row above_me
        // and bellow for data from other ranks to be pushed into
        Point *sim_space = new Point[(sizes[rank]+2)*Nx];

        // Creates PML sigma coficients
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
        
        // Chops out part of permitivity grid relavant to itself:
        double *ep = new double[sizes[rank]*Nx];
        for (int j=0; j<Nx; ++j)
        {
            for (int i=0; i<sizes[rank]; ++i)
            {
                ep[index(i,j,Nx)] = global_ep[index(i+c_sizes[rank-1],j,Nx)];
            }
        }

        // coficient matrix (taking into account ep, mu and pml:
        PML_coefs *Hx_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *Hy_coefs = new PML_coefs[sizes[rank]*Nx];
        PML_coefs *Dz_coefs = new PML_coefs[sizes[rank]*Nx];
        param_setter(Hx_coefs,Hy_coefs,Dz_coefs,sigma_x, local_sigma_y,1,1,dt,pml_size,sizes[rank], Nx);

        // Holds integrated values:
        double *IDz = new double[Nx * sizes[rank]];
        double *ICHx = new double[Nx * sizes[rank]];
        double *ICHy = new double[Nx * sizes[rank]];
        
        // which ranks to communicate with:
        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};

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
            
            update_H_pml(sim_space, Nx, sizes[rank], dx, dy, Hx_coefs, Hy_coefs,ICHx, ICHy);
            

            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t+(2*iterations), &sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t+(2*iterations),
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+(3*iterations), &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+(3*iterations),
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);

            update_E_pml(sim_space, Nx, sizes[rank], ep,dx, dy, dt, Dz_coefs,IDz);

            
            if(t == num_frame[which_frame])
            {
            SaveToFile(Nx*sizes[rank], &sim_space[Nx], "data/"+std::to_string(rank)+".txt");
            which_frame += 1;
            }
        }
    delete sim_space;
    delete local_sigma_y;
    delete ep;
    delete Hx_coefs;
    delete Hy_coefs;
    delete Dz_coefs;
    delete IDz;
    delete ICHx;
    delete ICHy;
    }

    else
    {
        // creates local simulation space, with extra row above_me
        // and bellow for data from other ranks to be pushed into
        Point *sim_space = new Point[(sizes[rank]+2)*Nx];
        
        // Chops out part of permitivity grid relavant to itself:
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
        
        // coficient matrix (taking into account ep, mu and pml:
        PML_coefs *Hx_coefs = new PML_coefs[pml_size];
        PML_coefs *Hy_coefs = new PML_coefs[pml_size];
        PML_coefs *Dz_coefs = new PML_coefs[pml_size];
        param_setter_x(Hx_coefs, Hy_coefs, Dz_coefs,sigma_x, 1, 1, dt, pml_size);

        // Holds integrated values:
        double *IDz = new double[2 * sizes[rank] * pml_size];
        double *ICHx = new double[2 * sizes[rank] * pml_size];
        double *ICHy = new double[2 * sizes[rank] * pml_size];

        // which ranks to communicate with:
        int above_me{(rank-1)%size};
        int under_me{(rank+1)%size};


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

            update_H_bulk(sim_space, Nx, sizes[rank], mux, muy, dx, dy, pml_size, Hx_coefs, Hy_coefs, ICHx, ICHy);
            
            MPI_Sendrecv(&sim_space[index(1,0,Nx)], Nx, MPI_POINT, above_me ,
            above_me+t+(2*iterations), &sim_space[index(sizes[rank]+1,0,Nx)] , Nx,
            MPI_POINT, under_me,rank+t+(2*iterations),
            MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&sim_space[index(sizes[rank],0,Nx)], Nx, MPI_POINT, under_me ,
            under_me+t+(3*iterations), &sim_space[0], Nx,
            MPI_POINT, above_me,rank+t+(3*iterations),
            MPI_COMM_WORLD, &status);

            MPI_Barrier(MPI_COMM_WORLD);
            
            
            update_E_bulk(sim_space, Nx, sizes[rank], ep,dx, dy, dt, pml_size, Dz_coefs, IDz); 
        
            if(t == num_frame[which_frame])
            {            
            SaveToFile(Nx*sizes[rank], &sim_space[Nx], "data/"+std::to_string(rank)+".txt");
            which_frame+=1;
            }

        }
        delete sim_space;
        delete ep;
        delete mux;
        delete muy;
        delete Hx_coefs;
        delete Hy_coefs;
        delete Dz_coefs;
        delete IDz;
        delete ICHx;
        delete ICHy;
    }

    delete sigma_x;
    delete sigma_y;
    delete global_ep;
    
    MPI_Finalize();    
    return 0;
}
