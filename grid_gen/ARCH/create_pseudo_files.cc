#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <sstream>
#include "grid.h"

using namespace std;

// ------------------------------------------------------------------------------
// Write cell into cell partition file
// ------------------------------------------------------------------------------
template <int DIM>
void Grid<DIM>::write_cell (Cell<DIM>& cell, unsigned int part_tag, bool g_cell)
{
    string file = part_dir+"/cell_file_";
    stringstream ss;
    ss<<part_tag;
    file+=ss.str();
    ofstream c_file;
    c_file.open(file.c_str(), ios::app);
    for(unsigned int i=0; i<DIM+1; i++)
       c_file << cell.vertex[i] << " ";
    c_file << g_cell
           << endl;
    c_file.close();
}

template class Grid<2>;
template class Grid<3>;