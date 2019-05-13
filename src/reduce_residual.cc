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
// Reduction of residual from neighbouring processors
//------------------------------------------------------------------------------
void FiniteVolume::reduce_residual()
{
    t_res_red.start_time();
    
    // STORING SHARED VERTEX RESIDUALS IN BUFFERS AND PERFORMING NON-BLOCKING
	// SENDS AND RECEIVES
	int gsize = grid.intp_mpi_groupings.size();
	int vsize[gsize];
	vector< vector <MPI_Request> > s_request, r_request;
	vector< VARARRAY > SBUF;
	vector< vector< VARARRAY > > RBUF;
	SBUF.resize(gsize);
	RBUF.resize(gsize);
	s_request.resize(gsize);
	r_request.resize(gsize);
	int buff_size[gsize];
   
	for(int i=0; i<gsize; ++i)
	{
	   int plen = grid.intp_mpi_groupings[i].proc_list.size();
	   RBUF[i].resize(plen - 1);
	   s_request[i].resize(plen - 1);
	   r_request[i].resize(plen - 1);
	}   
	
	for(int i=0; i<gsize; i++)
	{
	   vsize[i] = grid.intp_mpi_groupings[i].vertices.size();
	   buff_size[i] = vsize[i]*NVAR;
	   //RBUF[i].buf = new double[buff_size[i]];
	   SBUF[i].buf = new double[buff_size[i]];
	   int ind = 0;
	   for(int j = 0; j<vsize[i]; j++)
	   {
		  int locv=grid.g2l.find(grid.intp_mpi_groupings[i].vertices[j])->second;
		  SBUF[i].buf[ind++] = residual[locv].mass_flux;
		  SBUF[i].buf[ind++] = residual[locv].momentum_flux.x;
		  SBUF[i].buf[ind++] = residual[locv].momentum_flux.y;
		  SBUF[i].buf[ind++] = residual[locv].momentum_flux.z;
		  SBUF[i].buf[ind++] = residual[locv].energy_flux;
 
	   }
	   int sr_ind = 0;
	   set<int>::iterator it=grid.intp_mpi_groupings[i].proc_list.begin();
	   //for(int j=0; j<grid.intp_mpi_groupings[i].proc_list.size();++j)
	   for(int j=0; SafeLess(j,grid.intp_mpi_groupings[i].proc_list.size());++j)
	   {    
		  //if(*it!=get_proc_id())
		  if(SafeNeq(*it,get_proc_id()))
		  {
			  MPI_Isend(&SBUF[i].buf[0], buff_size[i], MPI_DOUBLE,*it,
						grid.intp_mpi_groupings[i].tag,MPI_COMM_WORLD,&s_request[i][sr_ind]);
			  RBUF[i][sr_ind].buf = new double[buff_size[i]];
			  MPI_Irecv(&RBUF[i][sr_ind].buf[0], buff_size[i], MPI_DOUBLE,*it,
						grid.intp_mpi_groupings[i].tag,MPI_COMM_WORLD,&r_request[i][sr_ind]);  
			  sr_ind++;          
		  }
		  ++it;                      
	   }   
	}
	
	// CHECKING IF NON-BLOCKING SENDS AND RECEIVES ARE COMPLETE AND ADDING CONTRIBUTIONS
	MPI_Status status;
	for(int i=0; i<gsize; i++)
	{
	   //for(int r=0; r<r_request[i].size(); ++r)
	   for(int r=0; SafeLess(r,r_request[i].size()); ++r)
	   {
		  MPI_Wait(&r_request[i][r],&status);   
		  int ind = 0;
		  for(int j = 0; j<vsize[i]; j++)
		  {
			 int locv=grid.g2l.find(grid.intp_mpi_groupings[i].vertices[j])->second;
			 residual[locv].mass_flux       += RBUF[i][r].buf[ind++];
			 residual[locv].momentum_flux.x += RBUF[i][r].buf[ind++];
			 residual[locv].momentum_flux.y += RBUF[i][r].buf[ind++];
			 residual[locv].momentum_flux.z += RBUF[i][r].buf[ind++];
			 residual[locv].energy_flux     += RBUF[i][r].buf[ind++]; 
		  } 
	   }
	   RBUF[i].clear();
	}
	RBUF.clear();
   
	for(int i=0; i<gsize; i++)
	   //for(int r=0; r<s_request[i].size(); ++r)
	   for(int r=0; SafeLess(r,s_request[i].size()); ++r)
		  MPI_Wait(&s_request[i][r],&status);   
	SBUF.clear();      
    
    t_res_red.add_time();     
}