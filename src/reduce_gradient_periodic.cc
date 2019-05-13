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
// Reduction of gradient from neighbouring periodic processors
//------------------------------------------------------------------------------
void FiniteVolume::reduce_gradient_periodic()
{
   for(unsigned int i=0; i<grid.n_periodic; ++i)
   {
	  vector<Vector> dE_sum;
      dE_sum.resize(NVAR);
      for(unsigned int k=0; k<NVAR;++k)
        dE_sum[k]=0.0;
      
	  for(unsigned int j=0; j<grid.periodic_lists[i].P_ij.size(); ++j)
	  {
		 int v = grid.periodic_lists[i].P_ij[j];
		 double w = 1.0/grid.mpi_node_share.find(v)->second;
		 
		 for(unsigned int k=0; k<NVAR;++k)
		    dE_sum[k].sadd(dE[v][k],w);
	  }
	  int p0 = grid.periodic_lists[i].P_ij[0];
	  dE[p0] = dE_sum;
   }
   
   if(NPART > 1)
   {
	  // STORING SHARED VERTEX GRADIENTS IN BUFFERS AND PERFORMING NON-BLOCKING
	  // SENDS AND RECEIVES
	  t_grad_red.start_time();
	  int gsize = grid.periodic_mpi_groupings.size();
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
		 int plen = grid.periodic_mpi_groupings[i].proc_list.size();
		 RBUF[i].resize(plen - 1);
		 s_request[i].resize(plen - 1);
		 r_request[i].resize(plen - 1);
	  }   
	
	  for(int i=0; i<gsize; i++)
	  {
		 vsize[i] = grid.periodic_mpi_groupings[i].vertices.size();
		 buff_size[i] = vsize[i]*NVAR*NDIR;
		 SBUF[i].buf = new double[buff_size[i]];
		 int ind = 0;
		 for(int j = 0; j<vsize[i]; j++)
		 {
			int v = grid.periodic_mpi_groupings[i].vertices[j];
			for(int k=0;k<NVAR;k++)
			{
			   SBUF[i].buf[ind++] = dE[v][k].x; 
			   SBUF[i].buf[ind++] = dE[v][k].y; 
			   SBUF[i].buf[ind++] = dE[v][k].z;
			}   
		 }
		 int sr_ind = 0;
		 for(set<int>::iterator it=grid.periodic_mpi_groupings[i].proc_list.begin();
			 it!=grid.periodic_mpi_groupings[i].proc_list.end(); ++it)
		 {    
			//if(*it!=get_proc_id())
			if(SafeNeq(*it,get_proc_id()))
			{
				MPI_Isend(&SBUF[i].buf[0], buff_size[i], MPI_DOUBLE,*it,
						  grid.periodic_mpi_groupings[i].tag,MPI_COMM_WORLD,&s_request[i][sr_ind]);
				RBUF[i][sr_ind].buf = new double[buff_size[i]];
				MPI_Irecv(&RBUF[i][sr_ind].buf[0], buff_size[i], MPI_DOUBLE,*it,
						  grid.periodic_mpi_groupings[i].tag,MPI_COMM_WORLD,&r_request[i][sr_ind]);          
				sr_ind++;  
			}                    
		 }   
	  }
   
	  // CHECKING IF NON-BLOCKING SENDS AND RECEIVES ARE COMPLETE AND ADDING CONTRIBUTIONS
	  t_grad_red.start_time();
	  MPI_Status status;
	  for(int i=0; i<gsize; ++i)
	  {
		 //for(int r=0; r<r_request[i].size(); ++r)
		 for(int r=0; SafeLess(r,r_request[i].size()); ++r)
		 {
			MPI_Wait(&r_request[i][r],&status);   
			int ind = 0;
			for(int j = 0; j<vsize[i]; j++)
			{
			   int v = grid.periodic_mpi_groupings[i].vertices[j];
			   for(int k=0;k<NVAR;k++)
			   {
				  dE[v][k].x += RBUF[i][r].buf[ind++];
				  dE[v][k].y += RBUF[i][r].buf[ind++]; 
				  dE[v][k].z += RBUF[i][r].buf[ind++];
			   }   
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
   
   }  
   
   for(unsigned int i=0; i<grid.n_periodic; ++i)
    {
       int p0 = grid.periodic_lists[i].P_ij[0];
       for(unsigned int j=1; j<grid.periodic_lists[i].P_ij.size(); ++j)
       {
          int v = grid.periodic_lists[i].P_ij[j];
          dE[v] = dE[p0];
       }
    }  
   
   t_grad_red.add_time(); 
}