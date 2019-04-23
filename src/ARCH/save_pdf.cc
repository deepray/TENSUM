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
// Saving pdf
//------------------------------------------------------------------------------
void FiniteVolume::save_pdf()
{
    int rank = get_proc_id();
	unsigned int nsmp = param.n_samples;
	unsigned int ntime = mc_time.size();
	unsigned int nvars = NVAR + param.write_variables.size();
    
    ofstream pdf_out_f;
    if(rank < NPART)
	{
        for(map<int,int>::iterator it = grid.l2g_mc_pt.begin(); it!=grid.l2g_mc_pt.end(); ++it)
     	{
			int globv = it->second + 1; // Back to gmsh indexing
			stringstream ss;
			ss << globv<<".dat";
			string file_path = PDF_DIR + "/pdf_n" + ss.str();
		
			pdf_out_f.open(file_path.c_str());
			MPI_LOC_ASSERT(pdf_out_f.is_open());
			pdf_out_f  << "NVAR"<<" "
					   <<  nvars<<" "
					   << "time_frames"<<" "
					   << ntime<<" "
					   << "samples"<<" "
					   << nsmp <<endl;
					
			pdf_out_f  << "pressure" << " "
					   << "u" << " " 
					   << "v" << " "
					   << "w" << " "
					   << "temperature"; 
			if(has_density)	        
			   pdf_out_f <<" "<< "density";
			if(has_mach)	        
			   pdf_out_f <<" "<< "mach";
			if(has_entropy)	        
			   pdf_out_f <<" "<< "entropy";
			if(has_vorticity)	        
			   pdf_out_f <<" "<< "vorticity";	
			pdf_out_f << endl;		        	        	        		        

            for(unsigned int i=0; i<N_MC_GROUPS; ++i)
            {
			   stringstream mc_id;
			   mc_id << i;
			   string file_path = PDF_DIR + "/" + mc_id.str()+"_pdf_n" + ss.str();
			   ifstream pdf_in_f;
			   pdf_in_f.open(file_path.c_str());
			   MPI_LOC_ASSERT(pdf_in_f.is_open());
			   string line;
			   while (getline (pdf_in_f,line) )
                  pdf_out_f << line <<endl;
			   pdf_in_f.close();
			}
			pdf_out_f.close();
		}   
	}

}