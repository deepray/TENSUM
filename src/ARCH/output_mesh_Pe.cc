#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <sstream>
#include "fv.h"
#include "parameter.h"
#include "mpi_utils.h"
using namespace std;
//------------------------------------------------------------------------------
// Writes Local mesh Peclet number for nodes
//------------------------------------------------------------------------------
void FiniteVolume::output_mesh_Pe(string DIR)
{   
   string filename_ind;
   stringstream ss;
   ss << counter <<"_" <<get_proc_loc_id();
   
   if     (counter <= 9)    filename_ind = "000";
   else if(counter <= 99)   filename_ind = "00";
   else if(counter <= 999)  filename_ind = "0";
   else if(counter <= 9999) filename_ind = "";
   else
      MPI_LOC_ERR("Peclet number output: counter is too large !!!");
   
   filename_ind+=ss.str();
   string filename = DIR + "/mesh_Pe_" + filename_ind + ".dat";
   
   ofstream PE_NO(filename.c_str());

   for(unsigned int i=0; i<grid.n_vertex; ++i)
   {
      PE_NO << setprecision(15) 
            << grid.vertex[i].coord.x<< " "
            << grid.vertex[i].coord.y<< " "
            << grid.dx_max[i] * primitive[i].velocity.norm()/param.material.viscosity(primitive[i].temperature)<<endl;
   }

   PE_NO.close();
   
}   