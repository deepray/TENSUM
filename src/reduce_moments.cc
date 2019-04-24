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
// Reduction of moments among processors handling a single type of mesh partition
//------------------------------------------------------------------------------
void FiniteVolume::reduce_moments()
{
	t_mom_red.start_time();
	
	int rank = get_proc_id();
	unsigned int NPROC = get_comm_size();
	unsigned int nvars = NVAR + param.write_variables.size();
	int tag;
	for(unsigned int t=0; t<time_instance.size(); ++t)
    {
	   if(rank < NPART)
	   {
		   double DIFF[grid.n_vertex*nvars], RBUF[grid.n_vertex*nvars], RBUF2[grid.n_vertex*nvars];
		   int n_old = param.sample_end_ind - param.sample_start_ind + 1;
		   int n_new = n_old;
		   int n_other;
		   for(unsigned int i = rank + NPART; i<NPROC; i+=NPART)
		   {
			  tag = i + t*NPROC;
			  MPI_Recv(&n_other,1,MPI_INT,i,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);	
			  	  
			  n_new = n_old + n_other;
			  double fact1 = double(n_other)/(double)n_new;
			  double fact2 = fact1*n_old;
		   
	          tag = i + (t+time_instance.size())*NPROC;
			  MPI_Recv(&DIFF, grid.n_vertex*nvars, MPI_DOUBLE, i, tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	
              tag = i + (t+2*time_instance.size())*NPROC;
			  MPI_Recv(&RBUF, grid.n_vertex*nvars, MPI_DOUBLE, i, tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			  
			  tag = i + (t+3*time_instance.size())*NPROC;
			  MPI_Recv(&RBUF2, grid.n_vertex*nvars, MPI_DOUBLE, i, tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			  for(unsigned int j=0; j<grid.n_vertex; ++j)
			  {
			   
				  DIFF[j*nvars+0] -= primitive_fm[t][j].pressure;
				  DIFF[j*nvars+1] -= primitive_fm[t][j].velocity.x;
				  DIFF[j*nvars+2] -= primitive_fm[t][j].velocity.y;
				  DIFF[j*nvars+3] -= primitive_fm[t][j].velocity.z;
				  DIFF[j*nvars+4] -= primitive_fm[t][j].temperature;
			   
				  primitive_fm[t][j].pressure    += DIFF[j*nvars+0]*fact1;
				  primitive_fm[t][j].velocity.x  += DIFF[j*nvars+1]*fact1;
				  primitive_fm[t][j].velocity.y  += DIFF[j*nvars+2]*fact1;
				  primitive_fm[t][j].velocity.z  += DIFF[j*nvars+3]*fact1;
				  primitive_fm[t][j].temperature += DIFF[j*nvars+4]*fact1;
				  
				  primitive_variance[t][j].pressure    += RBUF[j*nvars+0] + pow(DIFF[j*nvars+0],2.0)*fact2;
				  primitive_variance[t][j].velocity.x  += RBUF[j*nvars+1] + pow(DIFF[j*nvars+1],2.0)*fact2;
				  primitive_variance[t][j].velocity.y  += RBUF[j*nvars+2] + pow(DIFF[j*nvars+2],2.0)*fact2;
				  primitive_variance[t][j].velocity.z  += RBUF[j*nvars+3] + pow(DIFF[j*nvars+3],2.0)*fact2;
				  primitive_variance[t][j].temperature += RBUF[j*nvars+4] + pow(DIFF[j*nvars+4],2.0)*fact2;
				  
				  primitive_sm[t][j].pressure    += RBUF2[j*nvars+0];
				  primitive_sm[t][j].velocity.x  += RBUF2[j*nvars+1];
				  primitive_sm[t][j].velocity.y  += RBUF2[j*nvars+2];
				  primitive_sm[t][j].velocity.z  += RBUF2[j*nvars+3];
				  primitive_sm[t][j].temperature += RBUF2[j*nvars+4];
			   
			   
				  int ind = j*nvars+4;
			    
			   
				  if(has_density)
				  {
					 ind++;
					 DIFF[ind] -= density_fm[t][j];
					 density_fm[t][j] += DIFF[ind]*fact1;
					 density_variance[t][j] += RBUF[ind] + pow(DIFF[ind],2.0)*fact2;
					 density_sm[t][j] += RBUF2[ind];
				  }   
				  if(has_mach)
				  {
					 ind++;
					 DIFF[ind] -= mach_fm[t][j];
					 mach_fm[t][j] += DIFF[ind]*fact1;
					 mach_variance[t][j] += RBUF[ind] + pow(DIFF[ind],2.0)*fact2;
					 mach_sm[t][j] += RBUF2[ind];
				  }
				  if(has_entropy)
				  {
					 ind++;
					 DIFF[ind] -= entropy_fm[t][j];
					 entropy_fm[t][j] += DIFF[ind]*fact1;
					 entropy_variance[t][j] += RBUF[ind] + pow(DIFF[ind],2.0)*fact2;
					 entropy_sm[t][j] += RBUF2[ind];
				  }  
				  if(has_vorticity)
				  {
					 ind++;
					 DIFF[ind] -= vorticity_fm[t][j];
					 vorticity_fm[t][j] += DIFF[ind]*fact1;
					 vorticity_variance[t][j] += RBUF[ind] + pow(DIFF[ind],2.0)*fact2;
					 vorticity_sm[t][j] += RBUF2[ind];
				  }     
				  n_old = n_new;
			  }
				  
		   }
		   double fact3 = 1.0/(double(n_new) - 1.0);
		   double fact4 = 1.0/(double(n_new));
		   for(unsigned int j=0; j<grid.n_vertex; ++j)
		   {
			  primitive_variance[t][j] *= fact3; 
			  primitive_sm[t][j] *= fact4; 
			  if(has_density)
			  {
				 density_variance[t][j] *=fact3;
				 density_sm[t][j] *=fact4; 
			  }	 
			  if(has_mach)
			  {
				 mach_variance[t][j] *=fact3; 
				 mach_sm[t][j] *=fact4; 
			  }	 
			  if(has_vorticity)
			  {
				 vorticity_variance[t][j] *=fact3;  
				 vorticity_sm[t][j] *=fact4; 
			  }	 
			  if(has_entropy)
			  {
				 entropy_variance[t][j] *=fact3;   
				 entropy_sm[t][j] *=fact4; 
			  }	 
		   }
		
	   }
	   else
	   {
		   int n_loc_samples = param.sample_end_ind - param.sample_start_ind + 1;
		   int destination_proc_id = get_proc_id()%NPART;
		   int tag = rank + t*NPROC;
		   MPI_Ssend(&n_loc_samples, 1, MPI_INT, destination_proc_id,tag,MPI_COMM_WORLD);
		   
		   double SBUF[grid.n_vertex*nvars];
		   for(unsigned int i=0; i<grid.n_vertex; ++i)
		   {
			  SBUF[i*nvars+0] = primitive_fm[t][i].pressure;
			  SBUF[i*nvars+1] = primitive_fm[t][i].velocity.x;
			  SBUF[i*nvars+2] = primitive_fm[t][i].velocity.y;
			  SBUF[i*nvars+3] = primitive_fm[t][i].velocity.z;
			  SBUF[i*nvars+4] = primitive_fm[t][i].temperature;
		   
			  int ind = i*nvars+4;
		   
			  if(has_density)
				 SBUF[++ind] = density_fm[t][i];
			  if(has_mach)
				 SBUF[++ind] = mach_fm[t][i];
			  if(has_entropy)
				 SBUF[++ind] = entropy_fm[t][i];   
			  if(has_vorticity)
				 SBUF[++ind] = vorticity_fm[t][i];    
		   
		   }
		   
		   tag = rank + (t+time_instance.size())*NPROC;
		   MPI_Ssend(&SBUF, grid.n_vertex*nvars, MPI_DOUBLE, destination_proc_id,tag,MPI_COMM_WORLD);
		   
		   for(unsigned int i=0; i<grid.n_vertex; ++i)
		   {
			  SBUF[i*nvars+0] = primitive_variance[t][i].pressure;
			  SBUF[i*nvars+1] = primitive_variance[t][i].velocity.x;
			  SBUF[i*nvars+2] = primitive_variance[t][i].velocity.y;
			  SBUF[i*nvars+3] = primitive_variance[t][i].velocity.z;
			  SBUF[i*nvars+4] = primitive_variance[t][i].temperature;
			  int ind = i*nvars+4;
			  if(has_density)
				 SBUF[++ind] = density_variance[t][i];
			  if(has_mach)
				 SBUF[++ind] = mach_variance[t][i];
			  if(has_entropy)
				 SBUF[++ind] = entropy_variance[t][i];   
			  if(has_vorticity)
				 SBUF[++ind] = vorticity_variance[t][i];      
		   }
		   
		   tag = rank + (t+2*time_instance.size())*NPROC;
		   MPI_Ssend(&SBUF, grid.n_vertex*nvars, MPI_DOUBLE, destination_proc_id, tag,MPI_COMM_WORLD);
		   
		   for(unsigned int i=0; i<grid.n_vertex; ++i)
		   {
			  SBUF[i*nvars+0] = primitive_sm[t][i].pressure;
			  SBUF[i*nvars+1] = primitive_sm[t][i].velocity.x;
			  SBUF[i*nvars+2] = primitive_sm[t][i].velocity.y;
			  SBUF[i*nvars+3] = primitive_sm[t][i].velocity.z;
			  SBUF[i*nvars+4] = primitive_sm[t][i].temperature;
			  int ind = i*nvars+4;
			  if(has_density)
				 SBUF[++ind] = density_sm[t][i];
			  if(has_mach)
				 SBUF[++ind] = mach_sm[t][i];
			  if(has_entropy)
				 SBUF[++ind] = entropy_sm[t][i];   
			  if(has_vorticity)
				 SBUF[++ind] = vorticity_sm[t][i];      
		   }
		   
		   tag = rank + (t+3*time_instance.size())*NPROC;
		   MPI_Ssend(&SBUF, grid.n_vertex*nvars, MPI_DOUBLE, destination_proc_id, tag,MPI_COMM_WORLD);
	   }	
    }
    
    t_mom_red.add_time();
}