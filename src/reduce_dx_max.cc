#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "grid.h"
#include "mpi_utils.h"

extern Dimension dim;
extern bool debug;
extern bool verbose;
extern bool PERIODIC;
using namespace std;

//------------------------------------------------------------------------------
// Reduce dx_max
//------------------------------------------------------------------------------
void Grid::reduce_dx_max ()
{   

    double dx_max_0;

    for(unsigned int i=0; i<n_periodic; ++i)
    {
       dx_max_0 = 0;
       for(unsigned int j=0; j<periodic_lists[i].P_ij.size(); ++j)
       {
          int v = periodic_lists[i].P_ij[j];
          dx_max_0 = max(dx_max_0,dx_max[v]);
       }
       int p0 = periodic_lists[i].P_ij[0];
       
       dx_max[p0] = dx_max_0;
    }
    
    if(NPART > 1)
    {
	   // STORING PERIODIC VERTEX AREAS IN BUFFERS AND PERFORMING NON-BLOCKING
	   // SENDS AND RECEIVES
	   int gsize = periodic_mpi_groupings.size();
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
		  int plen = periodic_mpi_groupings[i].proc_list.size();
		  RBUF[i].resize(plen - 1);
		  s_request[i].resize(plen - 1);
		  r_request[i].resize(plen - 1);
	   }   
	
	   for(int i=0; i<gsize; i++)
	   {
		  vsize[i] = periodic_mpi_groupings[i].vertices.size();
		  buff_size[i] = vsize[i];
		  //RBUF[i].buf = new double[buff_size[i]];
		  SBUF[i].buf = new double[buff_size[i]];
		  int ind = 0;
		  for(int j = 0; j<vsize[i]; j++)
		  {
			 int v = periodic_mpi_groupings[i].vertices[j];
			 SBUF[i].buf[ind++] = dx_max[v]; 
		  }
		  int sr_ind = 0;
		  set<int>::iterator it=periodic_mpi_groupings[i].proc_list.begin();
		  //for(int j=0; j<periodic_mpi_groupings[i].proc_list.size();++j)
		  for(int j=0; SafeLess(j,periodic_mpi_groupings[i].proc_list.size());++j)
		  {    
			 //if(*it!=get_proc_id())
			 if(SafeNeq(*it,get_proc_id()))
			 {
				 MPI_Isend(&SBUF[i].buf[0], buff_size[i], MPI_DOUBLE,*it,
						   periodic_mpi_groupings[i].tag,MPI_COMM_WORLD,&s_request[i][sr_ind]);
				 RBUF[i][sr_ind].buf = new double[buff_size[i]];
				 MPI_Irecv(&RBUF[i][sr_ind].buf[0], buff_size[i], MPI_DOUBLE,*it,
						   periodic_mpi_groupings[i].tag,MPI_COMM_WORLD,&r_request[i][sr_ind]);  
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
				int v = periodic_mpi_groupings[i].vertices[j];
				dx_max[v] = max(dx_max[v],RBUF[i][r].buf[ind++]);
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
	
	for(unsigned int i=0; i<n_periodic; ++i)
    {
       int p0 = periodic_lists[i].P_ij[0];
       for(unsigned int j=1; j<periodic_lists[i].P_ij.size(); ++j)
       {
          int v = periodic_lists[i].P_ij[j];
          dx_max[v] = max(dx_max[v],dx_max[p0]);
       }
    }

}
