#include "Ez_point.h"
#include <iostream>
#include <fstream>
#include<string>

void SaveToFile(int size, Point row[], std::string name)
{
    std::fstream fs;
    fs.open(name, std::fstream::in | std::fstream::out | std::fstream::app);
    for(int i; i<size; ++i)
    {
        fs << row[i];
    }

}

int main()
{
    std::string savefile_path = "data";
    Point *array = new Point[100];
    std::string name = savefile_path + "/lol.txt";
    SaveToFile(100, array, name);

    
    return 0;
}
