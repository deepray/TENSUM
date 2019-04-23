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
// Reduction of norm of mean and variance from neighbouring processors
//------------------------------------------------------------------------------
void FiniteVolume::reduce_stat_norm()
{   
    int ntime = mc_time.size();
    int nvars = NVAR + param.write_variables.size();
    double SBUF_MEAN[nvars*ntime], RBUF_MEAN[nvars*ntime];
    double SBUF_VAR[nvars*ntime], RBUF_VAR[nvars*ntime];
    
	
	for(unsigned int t = 0; t<ntime; ++t)
	{
	   SBUF_MEAN[t*nvars + 0] = prim_mean_l1[t].temperature;
	   SBUF_MEAN[t*nvars + 1] = prim_mean_l1[t].velocity.x;
	   SBUF_MEAN[t*nvars + 2] = prim_mean_l1[t].velocity.y;
	   SBUF_MEAN[t*nvars + 3] = prim_mean_l1[t].velocity.z;
	   SBUF_MEAN[t*nvars + 4] = prim_mean_l1[t].pressure;
	   
	   SBUF_VAR[t*nvars + 0] = prim_var_l1[t].temperature;
	   SBUF_VAR[t*nvars + 1] = prim_var_l1[t].velocity.x;
	   SBUF_VAR[t*nvars + 2] = prim_var_l1[t].velocity.y;
	   SBUF_VAR[t*nvars + 3] = prim_var_l1[t].velocity.z;
	   SBUF_VAR[t*nvars + 4] = prim_var_l1[t].pressure;
	   
	   int ind = 5; 
	   if(has_density)
	   {
	      SBUF_MEAN[t*nvars + ind]  = den_mean_l1[t];
	      SBUF_VAR[t*nvars + ind++] = den_var_l1[t];
	   }  
	   if(has_mach)
	   {
	      SBUF_MEAN[t*nvars + ind]  = mach_mean_l1[t];
	      SBUF_VAR[t*nvars + ind++] = mach_var_l1[t];
	   }
	   if(has_entropy)
	   {
	      SBUF_MEAN[t*nvars + ind]  = ent_mean_l1[t];
	      SBUF_VAR[t*nvars + ind++] = ent_var_l1[t];
	   }
	   if(has_vorticity)
	   {
	      SBUF_MEAN[t*nvars + ind]  = vor_mean_l1[t];
	      SBUF_VAR[t*nvars + ind++] = vor_var_l1[t];
	   } 
	      
	}   
    
    MPI_Allreduce(&SBUF_MEAN, &RBUF_MEAN, nvars*ntime, MPI_DOUBLE, MPI_SUM, grid.run_comm);
    MPI_Allreduce(&SBUF_VAR, &RBUF_VAR, nvars*ntime, MPI_DOUBLE, MPI_SUM, grid.run_comm);
    
    for(unsigned int t = 0; t<ntime; ++t)
	{
	   prim_mean_l1[t].temperature    = RBUF_MEAN[t*nvars + 0];
	   prim_mean_l1[t].velocity.x     = RBUF_MEAN[t*nvars + 1];
	   prim_mean_l1[t].velocity.y     = RBUF_MEAN[t*nvars + 2];
	   prim_mean_l1[t].velocity.z     = RBUF_MEAN[t*nvars + 3];
	   prim_mean_l1[t].pressure       = RBUF_MEAN[t*nvars + 4];
	   
	   prim_var_l1[t].temperature     = RBUF_VAR[t*nvars + 0];
	   prim_var_l1[t].velocity.x      = RBUF_VAR[t*nvars + 1];
	   prim_var_l1[t].velocity.y      = RBUF_VAR[t*nvars + 2];
	   prim_var_l1[t].velocity.z      = RBUF_VAR[t*nvars + 3];
	   prim_var_l1[t].pressure        = RBUF_VAR[t*nvars + 4];
	   
	   int ind = 5; 
	   if(has_density)
	   {
	      den_mean_l1[t]   = RBUF_MEAN[t*nvars + ind];
	      den_var_l1[t]    = RBUF_VAR[t*nvars + ind++];
	   }  
	   if(has_mach)
	   {
	      mach_mean_l1[t]  = RBUF_MEAN[t*nvars + ind];
	      mach_var_l1[t]   = RBUF_VAR[t*nvars + ind++];
	   }
	   if(has_entropy)
	   {
	      ent_mean_l1[t]   = RBUF_MEAN[t*nvars + ind];
	      ent_var_l1[t]    = RBUF_VAR[t*nvars + ind++];
	   }
	   if(has_vorticity)
	   {
	      vor_mean_l1[t]   = RBUF_MEAN[t*nvars + ind];
	      vor_var_l1[t]    = RBUF_VAR[t*nvars + ind++];
	   } 
	      
	}  
}