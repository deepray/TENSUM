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
// Writes heat flux across specified boundary surface. This has been shifted to 
// the FiniteVolume class due to the requirement of boundary conditions
//------------------------------------------------------------------------------
void FiniteVolume::output_surface_hf(string DIR)
{   
   if(param.write_surfaces.size() == 0) return;
   
   string surffilename_ind;
   stringstream ss;
   ss << counter <<"_" <<get_proc_loc_id();
   
   if     (counter <= 9)    surffilename_ind = "000";
   else if(counter <= 99)   surffilename_ind = "00";
   else if(counter <= 999)  surffilename_ind = "0";
   else if(counter <= 9999) surffilename_ind = "";
   else
      MPI_LOC_ERR("Heatflux output: counter is too large !!!");
   
   surffilename_ind+=ss.str();

   const int nsurf = param.write_surfaces.size();
   vector<ofstream*> ff (nsurf);
   map<int,int> type_to_idx;
   for(int i=0; i<nsurf; ++i)
   {
      stringstream ss;
      ss << param.write_surfaces[i];
      string ffilename = DIR + "/q" + surffilename_ind + "_" + ss.str() + ".dat";
      ff[i] = new ofstream(ffilename.c_str());
      MPI_LOC_ASSERT (ff[i]->is_open());
      type_to_idx.insert(pair<int,int>(param.write_surfaces[i], i));
   }
   for(unsigned int i=0; i<grid.bface.size(); ++i)
   {
      const int type = grid.bface[i].type;
      const unsigned int v0 = grid.bface[i].vertex[0];
      const unsigned int v1 = grid.bface[i].vertex[1];
      map<int,int>::iterator it;
      it = type_to_idx.find(type);
      if(it != type_to_idx.end())
      {   
         BoundaryCondition& bc = param.boundary_condition[type];
         vector<PrimVar> prim_state(2);
         vector<EntVar> ent_state(2);
         prim_state[0] = primitive[v0];
         prim_state[1] = primitive[v1];
         bc.apply(grid.vertex[v0].coord, int_step_time, grid.bface[i], prim_state[0]);
         bc.apply(grid.vertex[v1].coord, int_step_time, grid.bface[i], prim_state[1]);
         for(unsigned int j=0; j<2; ++j)
            ent_state[j] = param.material.prim2ent(prim_state[j]);
         
         // viscous force, using vertex gradients:
         double T = (prim_state[0].temperature +
                     prim_state[1].temperature)/2.0;
         double mu = param.material.viscosity(T);
         
         double k = mu * param.material.Cp / param.material.prandtl;
         
         Vector gradT0, gradT1, gradT;
         double fact0 = 1.0/(param.material.gas_const*pow(ent_state[0].entl,2.0));
         double fact1 = 1.0/(param.material.gas_const*pow(ent_state[1].entl,2.0));
         
         gradT0.equ(dE[v0][4],fact0);
         gradT1.equ(dE[v1][4],fact1);
         gradT.equ(gradT0,gradT1,0.5,0.5);
         
         Vector normal,center;
         normal.equ(grid.bface[i].normal,(1.0/ grid.bface[i].measure));
         double qn = gradT*normal;
         qn *= (-k);
         center.equ(grid.vertex[v0].coord,grid.vertex[v1].coord,0.5,0.5);
         
         const int j = it->second;
         *ff[j] << center.x << "  " << center.y << "  "
                << qn << "  "
                << endl;
      }
   }
   for(int i=0; i<nsurf; ++i)
   {
      ff[i]->close();
   }
   
}   