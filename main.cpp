#include "Ez_point.h"
#include "soft_source.h"
#include <update_func_metalic_mpi.h>
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

void spacer(std::string name)
{
    std::fstream fs;
    fs.open(name, std::fstream::in | std::fstream::out | std::fstream::app);
    fs << '\n' << "####" << '\n';
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
    MPI_Type_contiguous(6, MPI_DOUBLE, &MPI_POINT); //initilise
    MPI_Type_commit(&MPI_POINT);//activate
   
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

    
    //Point *sim_space = new Point[Ny*Nx];
    
    
    
//    for (int i=0; i<iterations; ++i)
 //   {
 //       std::cout << sourceE[i] << std::endl;
 //   }
    //sim_space[index(source_center,source_center, Nx)].SetEz(1);
    
    
    int worker_Ny_base = Ny / size ;
    std::cout<< "RANK " << rank << " : worker_Ny_Base is: " << worker_Ny_base<<std::endl;
    int worker_Ny_remainder = Ny % size;
    std::cout << "RANK " << rank << ": worker_Ny_Remainder is: " << worker_Ny_remainder<<std::endl;

    
    int source_center{50};
    int rank_occuring {50 / worker_Ny_base};
    std::cout << "RANK " << rank << ": Rank occuring is: " << rank_occuring<<std::endl;

    
    std::cout << rank <<std::endl;
    
    if (rank == 0)
    {
        
        int local_Ny = worker_Ny_base;
        int split_size{local_Ny / 2};
//         Point *sim_space_top = new Point[split_size*Nx];
//         Point *sim_space_bottom = new Point[(local_Ny-split_size)*Nx];
        
        double *ep = new double[local_Ny*Nx];
        double *mux = new double[local_Ny*Nx];
        double *muy = new double[local_Ny*Nx];
        
        std::cout << "RANK " << rank << ": OK SO I definately made my mu ep arrays"<<std::endl;

        
        for (int j=0; j<Nx; ++j)
        {
            for (int i=0; i<local_Ny; ++i)
            {
                std::cout << '\n' << "RANK " << rank << " i=" << i <<" j="<< j << " Index is:" << index(i,j,Nx) << "Array size: "<< local_Ny*Nx<< "Nx is " << Nx << "locaal ny " << local_Ny <<'\n' ;

                ep[index(i,j,Nx)] = 1.0;
                mux[index(i,j,Nx)] = c * dt;
                muy[index(i,j,Nx)] = c * dt;
            }
        }

        std::cout << "RANK " << rank << ": AND I managed to set them to something"<<std::endl;

        
//         Point *bottom = new Point[Nx];
//         Point *top = new Point[Nx];
//         Point *local_bottom = new Point[Nx];
//         Point *local_top = new Point[Nx];
        
        int above_me{size-1};
        int under_me{1};

        std::cout << "I am: "<< rank << " And I think above me lives: "<< above_me << " and under me lives" << under_me << std::endl ;
        
//         for (int t=1; t<=iterations; ++t)
//         {   
//             for (int k=0; k<Nx; ++k)
//             {
//                 top[k] = sim_space_bottom[k];
//             }
//             
//             for (int k=index(split_size-1, 0, Nx); k<index(split_size-1, Nx, Nx); ++k)
//             {
//                 bottom[k] = sim_space_top[k];
//             }
//             
//             
//             std::cout << "I'm process " << rank << "and im sending to" << above_me << '\n';
//             MPI_Send(top, Nx, MPI_POINT, above_me ,1, MPI_COMM_WORLD);
//             std::cout << "I'm process " << rank << "and im recieveing from" << above_me << '\n';
//             MPI_Recv(local_top, Nx,MPI_POINT,under_me,1, MPI_COMM_WORLD, &status);
// 
//             std::cout << "I'm process " << rank << "and im sending to" << under_me << '\n';
//             MPI_Send(bottom, Nx, MPI_POINT, under_me, 2, MPI_COMM_WORLD);
//             std::cout << "I'm process " << rank << "and im recieveing from" << under_me << '\n';
//             MPI_Recv(local_bottom, Nx,MPI_POINT,above_me,2, MPI_COMM_WORLD, &status);
//             
// //             MPI_Sendrecv(&top, Nx, MPI_POINT, above_me ,
// //                         1, &local_bottom, Nx,
// //                         MPI_POINT, rank,1,
// //                         MPI_COMM_WORLD, &status);
//             
//            
// 
//             //top
//             update_H_worker(sim_space_top, Nx, split_size, mux, muy, dx, dy,local_bottom);
//             update_E_master(sim_space_top, Nx, split_size, ep,dx, dy, dt);
//             
// //              MPI_Sendrecv(&bottom, Nx, MPI_POINT, under_me ,
// //                         1, &local_top, Nx,
// //                         MPI_POINT, rank,1,
// //                         MPI_COMM_WORLD, &status);
//              
//             update_H_master(sim_space_bottom, Nx, local_Ny-split_size, mux, muy, dx, dy);
//             update_E_worker(sim_space_bottom, Nx, local_Ny-split_size, ep,dx, dy, dt,local_top);
//             
// 
//             SaveToFile(Nx*local_Ny, sim_space_top, "data/"+std::to_string(rank)+"TOP.txt");
//             spacer("data/"+std::to_string(rank)+".txt");
//             SaveToFile(Nx*local_Ny, sim_space_bottom, "data/"+std::to_string(rank)+"BOTTOM.txt");
//             spacer("data/"+std::to_string(rank)+".txt");
//         }

         
//         delete ep;
//         delete mux ;
//         delete muy;
//         delete top;
//         delete bottom;
//         delete sim_space_top;
//         delete sim_space_bottom;
//         delete local_bottom;
//         delete local_top;
// 

    }
// /*
//     else if (rank == rank_occuring)
//     {
//             int local_Ny{};
//             if (rank <= worker_Ny_remainder)
//             {
//                 local_Ny = worker_Ny_base + 1;
//             }
//             else 
//             {
//                 local_Ny = worker_Ny_base;
//             }
//             
//             
//             double *ep = new double[local_Ny*Nx];
//             double *mux = new double[local_Ny*Nx];
//             double *muy = new double[local_Ny*Nx];
//             
//         std::cout << "RANK " << rank << ": OK SO I definately made my mu ep arrays"<<std::endl;
//         
//             for (int i=0; i<local_Ny; ++i)
//             {
//                 for (int j=0; j<local_Ny; ++j)
//                 {
//                     ep[index(i,j,Nx)] = 1.0;
//                     mux[index(i,j,Nx)] = c * dt;
//                     muy[index(i,j,Nx)] = c * dt;
//                 }
//             }
// 
//         std::cout << "RANK " << rank << ": AND I managed to set them to something"<<std::endl;
// 
//             
//             double *sourceE = new double[iterations];
//             double *sourceHx = new double[iterations];
//             double *sourceHy = new double[iterations];
//         
//             
//             int above_me{(rank-1)%size};
//             int under_me{(rank+1)%size};
// 
//             std::cout << "I am: "<< rank << " And I think above me lives: "<< above_me << " and under me lives" << under_me << std::endl ;
// 
//         
// //             inject_soft_source(sourceE,sourceHx,sourceHy,iterations,dx,dy,dt, 1, 1, 1,1, source_sigma, t0);
// 
//             Point *bottom = new Point[Nx];
//             Point *top = new Point[Nx];
//             Point *sim_space = new Point[local_Ny*Nx];
//             Point *local_bottom = new Point[Nx];
//             Point *local_top = new Point[Nx];
//             
//                     
//             int source_center_y {50 / worker_Ny_base};
//             
// //             for (int t=1; t<=iterations; ++t)
// //             {
// //                 
// //                     
// //             for (int k=0; k<Nx; ++k)
// //             {
// //                 top[k] = sim_space[k];
// //             }
// //             
// //             for (int k=index(local_Ny-1, 0, Nx); k<index(local_Ny-1, Nx, Nx); ++k)
// //             {
// //                 bottom[k] = sim_space[k];
// //             }
// //             
// //             std::cout << "I'm process " << rank << "and im sending to" << above_me << '\n';
// //             MPI_Send(top, Nx, MPI_POINT, above_me ,1, MPI_COMM_WORLD);
// //             std::cout << "I'm process " << rank << "and im recieveing from" << above_me << '\n';
// //             MPI_Recv(local_top, Nx,MPI_POINT,under_me,1, MPI_COMM_WORLD, &status);
// // 
// //             std::cout << "I'm process " << rank << "and im sending to" << under_me << '\n';
// //             MPI_Send(bottom, Nx, MPI_POINT, under_me, 2, MPI_COMM_WORLD);
// //             std::cout << "I'm process " << rank << "and im recieveing from" << under_me << '\n';
// //             MPI_Recv(local_bottom, Nx,MPI_POINT,above_me,2, MPI_COMM_WORLD, &status);
// //             
// // 
// //                 
// // //             MPI_Sendrecv(&top, Nx, MPI_POINT, above_me ,
// // //                             1, &local_bottom, Nx,
// // //                             MPI_POINT, rank,1,
// // //                             MPI_COMM_WORLD, &status);
// // //             
// // //             MPI_Sendrecv(&bottom, Nx, MPI_POINT, under_me ,
// // //                             1, &local_top, Nx,
// // //                             MPI_POINT, rank,1,
// // //                             MPI_COMM_WORLD, &status);
// //             
// //             sim_space[index(source_center_y,source_center, Nx)].InjectEz(sourceE[t]);
// //             sim_space[index(source_center_y,source_center, Nx)].InjectDz(sourceE[t]);
// //             sim_space[index(source_center_y,source_center, Nx)].InjectHx(sourceHx[t]);
// //             sim_space[index(source_center_y,source_center, Nx)].InjectHy(sourceHy[t]);
// //                 
// //             update_H_worker(sim_space, Nx, local_Ny, mux, muy, dx, dy, local_bottom);
// //             update_E_worker(sim_space, Nx, local_Ny, ep,dx, dy, dt, bottom);
// // 
// //             
// //             SaveToFile(Nx*local_Ny, sim_space, "data/"+std::to_string(rank)+".txt");
// //             spacer("data/"+std::to_string(rank)+".txt");
// //         
// // 
// // 
// //         }
//         
//         delete ep;
//         delete mux ;
//         delete muy;
//         delete top;
//         delete bottom;
//         delete sim_space;
//         delete local_bottom;
//         delete local_top;
// 
//     }
//     
//     else 
//     {
//         int local_Ny {};
//         if (rank <= worker_Ny_remainder)
//         {
//             local_Ny = worker_Ny_base + 1;
//         }
//         else 
//         {
//             local_Ny = worker_Ny_base;
//         }
//         
//         double *ep = new double[local_Ny*Nx];
//         double *mux = new double[local_Ny*Nx];
//         double *muy = new double[local_Ny*Nx];
//         
//             std::cout << "RANK " << rank << ": OK SO I definately made my mu ep arrays"<<std::endl;
// 
//         for (int i=0; i<Nx; ++i)
//         {
//             for (int j=0; j<local_Ny; ++j)
//             {
//                 ep[index(i,j,Nx)] = 1.0;
//                 mux[index(i,j,Nx)] = c * dt;
//                 muy[index(i,j,Nx)] = c * dt;
//             }
//         }
//         std::cout << "RANK " << rank << ": AND I managed to set them to something"<<std::endl;
// 
//         
//         Point *bottom = new Point[Nx];
//         Point *top = new Point[Nx];
//         Point *sim_space = new Point[local_Ny*Nx];
//         Point *local_bottom = new Point[Nx];
//         Point *local_top = new Point[Nx];
//         
// 
//         int above_me{(rank-1)%size};
//         int under_me{(rank+1)%size};
//         
//         std::cout << "I am: "<< rank << " And I think above me lives: "<< above_me << " and under me lives" << under_me << std::endl ;
// 
// /*        
//         for (int t=1; t<=iterations; ++t)
//         {       
//             
//         for (int k=0; k<Nx; ++k)
//         {
//             top[k] = sim_space[k];
//         }
//         
//         for (int k=index(local_Ny, 0, Nx); k<index(local_Ny-1, Nx, Nx); ++k)
//         {
//             bottom[k] = sim_space[k];
//         }
//         
//             std::cout << "I'm process " << rank << "and im sending to" << above_me << '\n';
//             MPI_Send(top, Nx, MPI_POINT, above_me ,1, MPI_COMM_WORLD);
//             std::cout << "I'm process " << rank << "and im recieveing from" << above_me << '\n';
//             MPI_Recv(local_top, Nx,MPI_POINT,under_me,1, MPI_COMM_WORLD, &status);
// 
//             std::cout << "I'm process " << rank << "and im sending to" << under_me << '\n';
//             MPI_Send(bottom, Nx, MPI_POINT, under_me, 2, MPI_COMM_WORLD);
//             std::cout << "I'm process " << rank << "and im recieveing from" << under_me << '\n';
//             MPI_Recv(local_bottom, Nx,MPI_POINT,above_me,2, MPI_COMM_WORLD, &status);
//             */
// 
//         
// //         MPI_Sendrecv(&top, Nx, MPI_POINT, above_me ,
// //                      1, &local_bottom, Nx,
// //                      MPI_POINT, rank,1,
// //                      MPI_COMM_WORLD, &status);
//         
// 
// /*
//         SaveToFile(Nx*local_Ny, sim_space, "data/"+std::to_string(rank)+".txt");
//         spacer("data/"+std::to_string(rank)+".txt");
//  
//         update_H_worker(sim_space, Nx, local_Ny, mux, muy, dx, dy, local_bottom);
//         
// //         MPI_Sendrecv(&bottom, Nx, MPI_POINT, under_me ,
// //                 1, &local_top, Nx,
// //                 MPI_POINT, rank,1,
// //                 MPI_COMM_WORLD, &status);
//         
//         update_E_worker(sim_space, Nx, local_Ny, ep,dx, dy, dt, bottom);
// 
//         */
// // //        }
//         
//         delete ep;
//         delete mux ;
//         delete muy;
//         delete top;
//         delete bottom;
//         delete sim_space;
//         delete local_bottom;
//         delete local_top;
// 
//     }*/
    
    MPI_Finalize();    
    return 0;
}
