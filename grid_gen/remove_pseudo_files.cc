#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include "grid.h"

using namespace std;



// ------------------------------------------------------------------------------
// Removes pseudo files cell.dat, face.dat, node.dat
// ------------------------------------------------------------------------------
template <int DIM>
void Grid<DIM>::remove_pseudo_files ()
{
    cout<<"  --- removing pseudo files ... "<<endl;
    
    string commandline = "rm " + part_dir + "/cell.dat";
    system(commandline.c_str());
    
    commandline = "rm " + part_dir + "/face.dat";
    system(commandline.c_str());
    
    commandline = "rm " + part_dir + "/node.dat";
    system(commandline.c_str());
    
    
    
}

template class Grid<2>;
template class Grid<3>;

