#include <core_funcs.h>
#include <Ez_point.h>
#include <iostream>
#include <fstream>

/*
* Saves the whole 2D array (saved in a 1D array and manually indexed into a text file of a specified name.
* @param size Array size
* @param row array to save
* @param name save to name / savefile_path
*/
void SaveToFile(int size, Point row[], std::string name)
{
    std::fstream fs;
    fs.open(name, std::fstream::in | std::fstream::out | std::fstream::app);
    for(int i; i<size; ++i)
    {
        fs << row[i];
    }
}

/*
* Finds which rank a point on the x,y grid lives in.
* @param i y coordinate
* @param Ny Total y size
* @param *sizes Array holding sizes of each rank 
* @param rank_size How many ranks there are
* @return rank_occuring Rank where the specified y coordinate lives
*/
int rank_index_y(int i, int Ny, int *sizes, int rank_size)
{
    int rank_occuring {}; 
    for(int k=0; k<rank_size; ++k)
    {
        if (i-(sizes[k]-1) > 0)
        {
            i -= sizes[k]-1; // if it isnt in this rank, skip
        }
        else
        {
            rank_occuring = k;
            break;
        }
    }
    return rank_occuring;
}

/*
* Converts 2D array indicies into 1D array index
* @param i matrix row
* @param j matrix col
* @param Nx rowsize of matrix
* @return 1D index
*/
int index(int i,int j, int Nx){return ((i*Nx) + j);}
