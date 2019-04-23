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
// Saving sample pdf for each MC group
//------------------------------------------------------------------------------
void FiniteVolume::save_sample_pdf()
{
    int color = get_proc_id() / NPART;
    stringstream mc_id;
    mc_id << color;
    
    ofstream pdf_out_f;
    for(map<int,int>::iterator it = grid.l2g_mc_pt.begin(); it!=grid.l2g_mc_pt.end(); ++it)
	{
		int globv = it->second + 1; // Back to gmsh indexing
		int locv  = it->first;
		stringstream ss;
		ss << globv<<".dat";
		string file_path = PDF_DIR + "/" + mc_id.str() + "_pdf_n" + ss.str();
		pdf_out_f.open(file_path.c_str(),ios::app);
		MPI_LOC_ASSERT(pdf_out_f.is_open());
		pdf_out_f  << primitive[locv].pressure<< " "
				   << primitive[locv].velocity.x<< " "
				   << primitive[locv].velocity.y<< " "
				   << primitive[locv].velocity.z<< " "
				   << primitive[locv].temperature;
		if(has_density)	        
			pdf_out_f <<" "<< param.material.Density(primitive[locv]);
		if(has_mach)	        
			pdf_out_f <<" "<< param.material.Mach(primitive[locv]);
		if(has_entropy)	        
			pdf_out_f <<" "<< param.material.Entropy(primitive[locv]);
		if(has_vorticity)	        
			pdf_out_f <<" "<< param.material.Vorticity(primitive[locv],dE[locv]);	
		pdf_out_f << endl;					
		pdf_out_f.close();
	}   

}