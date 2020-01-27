#ifndef COREFUNC
#define COREFUNC
#include <Ez_point.h>
#include <iostream>
#include <fstream>

// Saves the whole 2D array (saved in a 1D array and manually indexed into a text file of a specified name.
void SaveToFile(int size, Point row[], std::string name);

// Finds which rank a point on the x,y grid lives in.
int rank_index_y(int i, int Ny, int *sizes, int rank_size);

// Converts 2D array indicies into 1D array index
int index(int i,int j, int Nx);

#endif 

