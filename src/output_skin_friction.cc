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
// Writes skin-friction on specified boundary surface. This has been shifted to 
// the FiniteVolume class due to the requirement of boundary conditions
//------------------------------------------------------------------------------
void FiniteVolume::output_surface_sf(string DIR)
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
      MPI_LOC_ERR("Skin-friction output: counter is too large !!!");
   
   surffilename_ind+=ss.str();

   const int nsurf = param.write_surfaces.size();
   vector<ofstream*> ff (nsurf);
   map<int,int> type_to_idx;
   for(int i=0; i<nsurf; ++i)
   {
      stringstream ss;
      ss << param.write_surfaces[i];
      string ffilename = DIR + "/f" + surffilename_ind + "_" + ss.str() + ".dat";
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
         
         Vector gradU0, gradU1, gradU, gradV0, gradV1, gradV;
         double fact0 = 1.0/ent_state[0].entl;
         double fact1 = 1.0/ent_state[1].entl;
         
         gradU0.equ(dE[v0][1],dE[v0][4],-fact0,fact0*ent_state[0].entv.x/ent_state[0].entl );
         gradU1.equ(dE[v1][1],dE[v1][4],-fact1,fact1*ent_state[1].entv.x/ent_state[1].entl );
         gradU.equ(gradU0,gradU1,0.5,0.5);
         
         gradV0.equ(dE[v0][2],dE[v0][4],-fact0,fact0*ent_state[0].entv.y/ent_state[0].entl );
         gradV1.equ(dE[v1][2],dE[v1][4],-fact1,fact1*ent_state[1].entv.y/ent_state[1].entl );
         gradV.equ(gradV0,gradV1,0.5,0.5);
         
         double div = gradU.x + gradV.y;
         double sxx = 2.0 * mu * (gradU.x - (1.0/3.0) * div);
         double syy = 2.0 * mu * (gradV.y - (1.0/3.0) * div);
         double sxy = mu * (gradU.y + gradV.x);
         Vector normal,center;
         normal.equ(grid.bface[i].normal,(1.0/ grid.bface[i].measure));
         double fx = (sxx * normal.x + sxy * normal.y);
         double fy = (sxy * normal.x + syy * normal.y);
         double cf = fx*normal.y - fy*normal.x;
         center.equ(grid.vertex[v0].coord,grid.vertex[v1].coord,0.5,0.5);
         
         const int j = it->second;
         *ff[j] << center.x << "  " << center.y << "  "
                << cf << "  "
                << endl;
      }
   }
   for(int i=0; i<nsurf; ++i)
   {
      ff[i]->close();
   }
   
}   