#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include <map>
#include "writer.h"
#include "mpi_utils.h"

extern Dimension dim;

using namespace std;

//------------------------------------------------------------------------------
// Add primitive variables defined at vertices
//------------------------------------------------------------------------------
void Writer::attach_data (vector<PrimVar>& data)
{
   MPI_LOC_ASSERT (!has_primitive);
   vertex_primitive = &data;
   has_primitive = true;
}

//------------------------------------------------------------------------------
// Add data defined at vertices
//------------------------------------------------------------------------------
void Writer::attach_data (vector<double>& data, std::string name)
{
   vertex_data.push_back (&data);
   vertex_data_name.push_back (name);
}

//------------------------------------------------------------------------------
// Specify which variables to write
//------------------------------------------------------------------------------
void Writer::attach_variables (const vector<string>& variables)
{
   if(variables.size() > 0)
      MPI_LOC_ASSERT (has_primitive);

   for(unsigned int i=0; i<variables.size(); ++i)
      if(variables[i]=="mach")
         write_mach = true;
      else if(variables[i]=="density")
         write_density = true;
      else if(variables[i]=="entropy")
         write_entropy = true;
      else if(variables[i]=="vorticity")
         write_vorticity = true; 
      else if(variables[i]=="mesh_Peclet")
         write_mesh_Peclet = true;      
      else
         MPI_LOC_ERR("Writer: unknown variable " << variables[i]);
}

//------------------------------------------------------------------------------
// Add entropy gradient values
//------------------------------------------------------------------------------
void Writer::attach_gradient (std::vector<std::vector<Vector> >& dE_)
{
   MPI_LOC_ASSERT (!has_gradient);
   dE = &dE_;
   has_gradient = true;
}
    

//------------------------------------------------------------------------------
// Call output function for saving solution to file
//------------------------------------------------------------------------------
void Writer::output (int part, int counter)
{
   string filename,filename_comm, filemaster;
   if     (counter <= 9)    filename = "sol000";
   else if(counter <= 99)   filename = "sol00";
   else if(counter <= 999)  filename = "sol0";
   else if(counter <= 9999) filename = "sol";
   else
      MPI_LOC_ERR("Writer::output: counter is too large !!!");
   
   filename_comm = filename;
   stringstream ss;
   ss << counter <<"_" <<part;
   filename += ss.str();

   if(format == "vtk")
   {
      filename += ".vtk";
      string file_path = SAMPLE_DIR + "/" +filename;
      output_vtk (file_path);
   }
   else
      MPI_LOC_ERR("Unrecongized output format "<<format<<". Only vtk output available at present");
      
   if(check_group_base())
   {
	   if(write_to_master)
	   {
		   for (int i=0; i<NPART; i++)
		   {
			  filemaster = filename_comm;
			  stringstream ss_comm;
			  ss_comm << counter <<"_"<<i<<"."<<format;
			  filemaster += ss_comm.str();  
			  master_file <<filemaster<<endl;
		   }    
	   }
   }

   // save solution at specified surfaces
   string surffilename_ind;
   if     (counter <= 9)    surffilename_ind = "000";
   else if(counter <= 99)   surffilename_ind = "00";
   else if(counter <= 999)  surffilename_ind = "0";
   else if(counter <= 9999) surffilename_ind = "";
   else
      MPI_LOC_ERR("Writer::output: counter is too large !!!");
      
   surffilename_ind += ss.str();
   output_surfaces (SAMPLE_DIR,surffilename_ind);      
   
}

//------------------------------------------------------------------------------
// Write data to vtk file
//------------------------------------------------------------------------------
void Writer::output_vtk (string filename)
{
   ofstream vtk;
   vtk.open (filename.c_str());

   vtk << "# vtk DataFile Version 3.0" << endl;
   vtk << "flo3d" << endl;
   vtk << "ASCII" << endl;
   if(time_mode == "unsteady")
   {
	   vtk << "DATASET UNSTRUCTURED_GRID" << endl;
	   vtk << "FIELD FieldData 1" <<endl;
	   vtk << "TIME 1 1 double"<<endl;
	   vtk << time <<endl;
   }
   vtk << "POINTS  " << grid->n_vertex << "  double" << endl;

   for(unsigned int i=0; i<grid->n_vertex; ++i)
   vtk << grid->vertex[i].coord.x << " " 
	   << grid->vertex[i].coord.y << " " 
	   << grid->vertex[i].coord.z << endl;	  

   vtk << "CELLS  " << grid->n_cell << " " << 4 * grid->n_cell << endl;
   for(unsigned int i=0; i<grid->n_cell; ++i)
   {
	   int v0 = grid->cell[i].vertex[0];
	   int v1 = grid->cell[i].vertex[1];
	   int v2 = grid->cell[i].vertex[2];

	   vtk << 3 << " " 
		   << v0<< " "
		   << v1 <<" "
		   << v2 << endl;
   }       

   vtk << "CELL_TYPES  " << grid->n_cell << endl;
   for(unsigned int i=0; i<grid->n_cell; ++i)
      vtk << 5 << endl;

   // Write vertex data
   if(vertex_data.size() > 0 || has_primitive) 
      vtk << "POINT_DATA  " << grid->n_vertex << endl;

   // If vertex primitive data is available, write to file
   if (has_primitive)
   {
      vtk << "SCALARS pressure double 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
         vtk << setprecision(12) <<(*vertex_primitive)[i].pressure << endl;

      vtk << "SCALARS temperature double 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
         vtk <<setprecision(12) << (*vertex_primitive)[i].temperature << endl;

      if(dim == axi)
      {
         vtk << "SCALARS Vtheta double 1" << endl;
         vtk << "LOOKUP_TABLE default" << endl;
         for(unsigned int i=0; i<grid->n_vertex; ++i)
            vtk <<setprecision(12) << (*vertex_primitive)[i].velocity.z << endl;
      }

      vtk << "VECTORS velocity double" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
		 vtk <<setprecision(12) 
		     << (*vertex_primitive)[i].velocity.x << " "
			 << (*vertex_primitive)[i].velocity.y << " "
			 << 0.0
			 << endl;	 
   }

   // Write mach number
   if(write_mach)
   {
      vtk << "SCALARS mach double 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
      {
		 double mach = material->Mach ( (*vertex_primitive)[i] );
		 vtk <<setprecision(12) << mach << endl;
		 vtk.flush();
      }
   }

   // Write density
   if(write_density)
   {
      vtk << "SCALARS density double 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
      {
         double density = material->Density ((*vertex_primitive)[i]);
		 vtk <<setprecision(12) << density << endl;
		 vtk.flush();	 
      }
   }
   
   // Write entropy
   if(write_entropy)
   {
      vtk << "SCALARS entropy double 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
      {
         double density = material->Entropy ((*vertex_primitive)[i]);
		 vtk <<setprecision(12) << density << endl;
		 vtk.flush();
      }
   }
   
   // write vorticity using entropy variable gradients
   if(write_vorticity)
   {
      // Check if gradient information is available
      MPI_LOC_ASSERT(has_gradient);
      vtk << "SCALARS vorticity double 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
      {
         double vorticity = material->gas_const*(*vertex_primitive)[i].temperature*(
							(*dE)[i][2].x + (*vertex_primitive)[i].velocity.y*(*dE)[i][4].x
							- (*dE)[i][1].y - (*vertex_primitive)[i].velocity.x*(*dE)[i][4].y);
		 vtk <<setprecision(12) << vorticity << endl;	 
      }
   }
   
   // write mesh Peclet number
   if(write_mesh_Peclet)
   {
      // Check if gradient information is available
      MPI_LOC_ASSERT(has_gradient);
      vtk << "SCALARS mesh_Peclet double 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
      {
         double Peclet = grid->dx_max[i]*material->Density ((*vertex_primitive)[i])*
                               (*vertex_primitive)[i].velocity.norm()/
                                material->viscosity((*vertex_primitive)[i].temperature);					
							
		 vtk <<setprecision(12) << Peclet << endl;	 
      }
   }

   // Write vertex data to file
   for(unsigned int d=0; d<vertex_data.size(); ++d)
   {
      vtk << "SCALARS  " << vertex_data_name[d] << "  double 1" << endl;
      vtk << "LOOKUP_TABLE default" << endl;
      for(unsigned int i=0; i<grid->n_vertex; ++i)
	  {
		  vtk <<setprecision(12) << (*vertex_data[d])[i] << endl;
		  vtk.flush();
	  }		 
   }
   
   // vtk << "VECTORS grad_V1 double" << endl;
//    for(unsigned int i=0; i<grid->n_vertex; ++i)
// 	  vtk <<setprecision(12) 
// 		  << (*dE)[i][0].x << " "
// 		  << (*dE)[i][0].y << " "
// 		  << (*dE)[i][0].z
// 		  << endl;
// 		  
//    vtk << "VECTORS grad_V2 double" << endl;
//    for(unsigned int i=0; i<grid->n_vertex; ++i)
// 	  vtk <<setprecision(12) 
// 		  << (*dE)[i][1].x << " "
// 		  << (*dE)[i][1].y << " "
// 		  << (*dE)[i][1].z
// 		  << endl;
// 		  
//    vtk << "VECTORS grad_V3 double" << endl;
//    for(unsigned int i=0; i<grid->n_vertex; ++i)
// 	  vtk <<setprecision(12) 
// 		  << (*dE)[i][2].x << " "
// 		  << (*dE)[i][2].y << " "
// 		  << (*dE)[i][2].z
// 		  << endl;
// 		  
//    vtk << "VECTORS grad_V4 double" << endl;
//    for(unsigned int i=0; i<grid->n_vertex; ++i)
// 	  vtk <<setprecision(12) 
// 		  << (*dE)[i][3].x << " "
// 		  << (*dE)[i][3].y << " "
// 		  << (*dE)[i][3].z
// 		  << endl;
// 		  
//    vtk << "VECTORS grad_V5 double" << endl;
//    for(unsigned int i=0; i<grid->n_vertex; ++i)
// 	  vtk <<setprecision(12) 
// 		  << (*dE)[i][4].x << " "
// 		  << (*dE)[i][4].y << " "
// 		  << (*dE)[i][4].z
// 		  << endl;		  		  		  		  
   
   
   vtk.close ();
}

//------------------------------------------------------------------------------
// Write solution at surfaces
// For each type, two files are created, one with data at vertices (pressure)
// and another at face centers (skin friction) [SECOND FILE WRITTEN BY A FUNCTION
// OF FiniteVolume AS BC APPLICATION WAS REQUIRED]
//------------------------------------------------------------------------------
void Writer::output_surfaces (string DIR_NAME, string index)
{
   if(surfaces.size() == 0) return;

   const int nsurf = surfaces.size();
   vector<ofstream*> fv (nsurf);
   map<int,int> type_to_idx;
   for(int i=0; i<nsurf; ++i)
   {
      stringstream ss;
      ss << surfaces[i];
      string vfilename = DIR_NAME + "/v" + index + "_" + ss.str() + ".dat";
      fv[i] = new ofstream(vfilename.c_str());
      MPI_LOC_ASSERT (fv[i]->is_open());
      type_to_idx.insert(pair<int,int>(surfaces[i], i));
   }
   for(unsigned int i=0; i<grid->bface.size(); ++i)
   {
      const int type = grid->bface[i].type;
      const unsigned int v0 = grid->bface[i].vertex[0];
      //const unsigned int v1 = grid->bface[i].vertex[1];
      map<int,int>::iterator it;
      it = type_to_idx.find(type);
      if(it != type_to_idx.end())
      {
         const int j = it->second;
         *fv[j] << grid->vertex[v0].coord.x << "  "
                << grid->vertex[v0].coord.y << "  "
                << (*vertex_primitive)[v0].pressure << "  "
                << (*vertex_primitive)[v0].temperature << "  "
                << (*vertex_primitive)[v0].velocity.x  << "  "
                << (*vertex_primitive)[v0].velocity.y  << "  "
                << (*vertex_primitive)[v0].velocity.z  << "  "
                << endl;
      }
   }
   for(int i=0; i<nsurf; ++i)
   {
      fv[i]->close();
   }
}


//------------------------------------------------------------------------------
// Write solution for restarting
//------------------------------------------------------------------------------
void Writer::output_restart (int iter, int counter, int head_tag, double time, double res0)
{
   MPI_LOC_ASSERT (has_primitive);

   // if(get_proc_loc_id() == 0)
//          cout << "Saving restart file(s) ...\n";
   ofstream fo;
   string restart_file_path = SAMPLE_DIR + "/restart_";
   stringstream ss;
   ss <<head_tag<<"_"<<get_proc_loc_id();
   restart_file_path += ss.str()+".dat";
   fo.open (restart_file_path.c_str());
   MPI_LOC_ASSERT (fo.is_open());

   for(unsigned int i=0; i<grid->n_vertex; ++i)
      fo << scientific << setprecision(15)
         << (*vertex_primitive)[i].temperature << "  "
         << (*vertex_primitive)[i].velocity.x  << "  "
         << (*vertex_primitive)[i].velocity.y  << "  "
         << (*vertex_primitive)[i].velocity.z  << "  "
         << (*vertex_primitive)[i].pressure    << endl;

   fo << iter << endl;
   fo << counter << endl;
   fo << time << endl;
   fo << res0 << endl;
   fo.close ();
}
