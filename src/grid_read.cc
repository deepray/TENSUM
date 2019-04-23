#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <sstream>
#include "grid.h"
#include "parameter.h"
#include "mpi_utils.h"

using namespace std;
extern bool verbose;
extern bool debug;

//------------------------------------------------------------------------------
// Read grid from file
//------------------------------------------------------------------------------
void Grid::read (const Parameter& param)
{
   cell_type = param.cell_type;
   unsigned int rank = get_proc_id();
   create_group_comm(run_comm);
   
   mpi_barrier(MPI_COMM_WORLD);
   read_part_mesh_file(param.part_dir,param.boundary_condition);
   mpi_barrier(run_comm);

   // At this stage, we have only boundary faces. We save this number.
   n_boundary_face = n_face;
   check_face_type (param.boundary_condition);
   if(rank == 0)
      cout << "\n  Creating mesh data structures... " << endl;
   preproc (param.boundary_condition);
  
   mpi_barrier(run_comm);
   
//    if(debug)
// 	  if(get_proc_id() == 0)
// 	  {
// 		 string commandline = "rm -rf probe_info.dat";
// 		 system(commandline.c_str());
// 		 ofstream probe_f;
// 		 probe_f.open("probe_info.dat");
// 		 probe_f << param.write_mc_points.size()<<endl;
// 	  }
   
   // if(debug)
// 	  for(unsigned int i=0; i<NPART; ++i)
// 	  {
// 		 if(get_proc_id() == i)
// 		 {
// 			ofstream probe_f;
// 			probe_f.open("probe_info.dat" , ios::app);
// 			for(map<int,int>::iterator it=l2g_mc_pt.begin(); it!=l2g_mc_pt.end(); it++)
// 			{
// 			   int locv = it->first;
// 			   int globv = it->second +1;
// 			   probe_f <<globv<<" "
// 					   <<i<<" "
// 					   <<locv<<" "
// 					   <<vertex[locv].coord.x<<" "
// 					   <<vertex[locv].coord.y<<" "<<endl;
// 			}
// 			probe_f.close();
// 		 }
// 		 mpi_barrier(run_comm);
// 	  }

   
   if(get_proc_id() == 0)
   {
	  string commandline = "rm -rf mesh_part_info.dat";
	  system(commandline.c_str());
   }  
   mpi_barrier(MPI_COMM_WORLD);
   for(unsigned int i=0; i<NPART; i++)
   {
      if(rank == i)
      {
         info ();
         fflush(stdout);
      }   
      MPI_Barrier(run_comm);
   }  
   fflush(stdout);
}

//------------------------------------------------------------------------------
// Print some grid information to screen
//------------------------------------------------------------------------------
void Grid::info ()
{
   double min_mcarea =  1.0e20;
   double max_mcarea = -1.0e20;
   double min_dcarea =  1.0e20;
   double max_dcarea = -1.0e20;
   double min_vradius =  1.0e20;
   double max_vradius = -1.0e20;

   for(unsigned int i=0; i<n_vertex; ++i)
   {
	  min_mcarea = min ( min_mcarea, mcarea[i] );
	  max_mcarea = max ( max_mcarea, mcarea[i] );
	  min_dcarea = min ( min_dcarea, dcarea[i] );
	  max_dcarea = max ( max_dcarea, dcarea[i] );
	  min_vradius = min ( min_vradius, vertex[i].radius );
	  max_vradius = max ( max_vradius, vertex[i].radius );
   }

   double min_face_length =  1.0e20;
   double max_face_length = -1.0e20;
   double min_fradius =  1.0e20;
   double max_fradius = -1.0e20;
   for(unsigned int i=0; i<n_face; ++i)
   {
      min_face_length = min ( min_face_length, face[i].normal.norm() );
      max_face_length = max ( max_face_length, face[i].normal.norm() );
      min_fradius = min ( min_fradius, face[i].radius );
      max_fradius = max ( max_fradius, face[i].radius );
   }
   
   double min_cell_area =  1.0e20;
   double max_cell_area = -1.0e20;
   double min_cradius =  1.0e20;
   double max_cradius = -1.0e20;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      min_cell_area = min ( min_cell_area, cell[i].area );
	  max_cell_area = max ( max_cell_area, cell[i].area );
	  min_cradius = min ( min_cradius, cell[i].radius );
	  max_cradius = max ( max_cradius, cell[i].radius );
   }
   
   string partion_info_file = "mesh_part_info.dat";
   ofstream part_f;
   part_f.open(partion_info_file.c_str() , ios::app);   	  
			 
   part_f << "\n------------------------------------------------------------" <<endl;
   part_f << "             Grid information for (local) processor "<<get_proc_id()<<":       "<<endl;
   part_f << "------------------------------------------------------------" <<endl;
   part_f << "   Number of vertices             = " << n_vertex_active << endl;
   part_f << "   Number of shared vertices      = " << n_vertex_mpi_bnd << endl;
   part_f << "   Number of triangles            = " << n_cell_active << endl;
   part_f << "   Number of interior faces       = " << n_face - n_boundary_face << endl;
   part_f << "   Number of boundary edges       = " << n_boundary_face << endl;
   part_f << setw(30) << "min" << setw(15) << "max" << endl;
   part_f << "  cell area    =  " << setw(15) << min_cell_area 
								<< setw(15) << max_cell_area << endl;
   part_f << "  median area  =  " << setw(15) << min_mcarea 
								<< setw(15) << max_mcarea << endl;
   part_f << "  dual area    =  " << setw(15) << min_dcarea
								<< setw(15) << max_dcarea << endl;
   part_f << "  face length  =  " << setw(15) << min_face_length
								<< setw(15) << max_face_length << endl;
   part_f << "  vertex radius = " << setw(15) << min_vradius 
								<< setw(15) << max_vradius << endl;
   part_f << "  face   radius = " << setw(15) << min_fradius
								<< setw(15) << max_fradius << endl;
   part_f << "  tri    radius = " << setw(15) << min_cradius
								<< setw(15) << max_cradius << endl;
   part_f << "------------------------------------------------------------" <<endl;   
  
   MPI_LOC_ASSERT (min_vradius > 0.0);
   MPI_LOC_ASSERT (min_fradius > 0.0);
   MPI_LOC_ASSERT (min_cradius > 0.0);
   fflush(stdout);
   part_f.close();                            
}

//------------------------------------------------------------------------------
// Check that all boundary faces have been assigned a bc type
//------------------------------------------------------------------------------
void Grid::check_face_type (const map<int,BoundaryCondition>& bc)
{
   for(unsigned int i=0; i<face.size(); ++i)
   {
      if(face[i].type == -1)
         MPI_LOC_ERR("Face = "<< i <<"in partition "<<get_proc_loc_id()<<" has not been assigned a tag");
      if(bc.find(face[i].type) == bc.end())
         MPI_LOC_ERR("check_face_type in partition "<<get_proc_loc_id()<<" :\n"
                 << "   No boundary condition specified for\n"
                 << "   face = " << i << " whose type = " << face[i].type << "\n"
                 << "   There may be more faces with similar problem.\n");
   }
}
