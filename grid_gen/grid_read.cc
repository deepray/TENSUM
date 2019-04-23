#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <sstream>
#include "grid.h"

using namespace std;

//------------------------------------------------------------------------------
// Read grid from file
//------------------------------------------------------------------------------
template <int DIM>
void Grid<DIM>::read ()
{
   string commandline0 = "rm -rf " + part_dir;
   string commandline1 = "mkdir " + part_dir;
   system(commandline0.c_str());
   system(commandline1.c_str());

   read_gmsh ();	  
   create_partition_files();
   remove_pseudo_files();
}

template class Grid<2>;
template class Grid<3>;