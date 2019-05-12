#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "fv.h"

using namespace std;
extern bool verbose;

//------------------------------------------------------------------------------
// For each force, create list of faces
//------------------------------------------------------------------------------
void FiniteVolume::create_force_face_list ()
{
   if(param.force_data.size() == 0)
      return;

   MPI_DISP("\n  Creating list of faces for force computation",verbose);   

   inv_force.resize (param.force_data.size());
   visc_force.resize (param.force_data.size());

   // For each force, count how many faces of a given type were found
   // This is meant to catch a mistake where a user specifies a face type
   // which does not exist in the grid.
   vector< map<int,int> > nface(param.force_data.size());
   for(unsigned int j=0; j<param.force_data.size(); ++j)
   {
      ForceData& force_data = param.force_data[j];
      for(unsigned int k=0; k<force_data.face_type.size(); ++k)
         nface[j].insert(pair<int,int>(force_data.face_type[k],0));
   }

   // Forces are computed only on boundary faces
   for(unsigned int i=0; i<grid.bface.size(); ++i)
   {
      int face_type = grid.bface[i].type;

      for(unsigned int j=0; j<param.force_data.size(); ++j)
      {
         ForceData& force_data = param.force_data[j];
         for(unsigned int k=0; k<force_data.face_type.size(); ++k)
            if(force_data.face_type[k] == face_type)
            {
               inv_force[j].face.push_back (i);
               visc_force[j].face.push_back (i);
               ++nface[j][face_type];
            }
      }
   }
   
   mpi_barrier(grid.run_comm);

   // Check for mistakes
   bool ok = true;

   // Check that all forces have faces
   for(unsigned int i=0; i<param.force_data.size(); ++i)
   {
      int fface_mpi;
      int fface_loc = inv_force[i].face.size();
      if(NPART > 1)
         MPI_Allreduce(&fface_loc,&fface_mpi,1,MPI_INT,MPI_SUM,grid.run_comm);
      else
         fface_mpi = fface_loc;   
      if(fface_mpi == 0)
      {
         if(get_proc_id() == 0)
         {
            cout << "Force " << param.force_data[i].name << " does not have any faces\n";
         }   
         ok = false;
      }
   }   

   // Check that face types actually were found in the grid
   for(unsigned int j=0; j<param.force_data.size(); ++j)
   {
      ForceData& force_data = param.force_data[j];
      for(unsigned int k=0; k<force_data.face_type.size(); ++k)
      {
         int fface_mpi;
         if(NPART > 1)
            MPI_Allreduce(&nface[j][force_data.face_type[k]],&fface_mpi,1,MPI_INT,MPI_SUM,
                       grid.run_comm);
         else 
            fface_mpi = nface[j][force_data.face_type[k]];               
                                     
         if(fface_mpi == 0)
         {
            if(get_proc_id() == 0)
            {
               cout << "Force: " << param.force_data[j].name << " has face type ";
               cout << force_data.face_type[k] << " but it is not present in the grid\n";
            }   
            ok = false;
         }
      }   
   }

   if(!ok)
      mpi_err_finalize();

   if(get_proc_id() == 0)
   {
      if(param.force_data.size() == 0)
         cout << "  --- No forces found\n";
      else
         cout << "  --- Found " << param.force_data.size() << " forces\n";
   }      
}

//------------------------------------------------------------------------------
// Compute forces
// TODO: Axisymmetric case
//------------------------------------------------------------------------------
void FiniteVolume::compute_forces (unsigned long int iter)
{
   // Do we have any forces to compute
   if(param.force_data.size() == 0) return;
   
   t_force_eval.start_time();

   // Recompute gradient needed for viscous force
   if(param.material.model == Material::ns) compute_gradients ();
   if(check_group_base())
   {
	   force_file << setw(6) << iter << " " << scientific << setw(15);
	   force_file << elapsed_time << setw(15);
   }

   for(unsigned int i=0; i<param.force_data.size(); ++i)
   {
      inv_force[i].value = 0.0;
      visc_force[i].value = 0.0;

      for(unsigned int j=0; j<inv_force[i].face.size(); ++j)
      {
         unsigned int face_no = inv_force[i].face[j];
         Vector& normal  = grid.bface[face_no].normal;

         unsigned int v0 = grid.bface[face_no].vertex[0];
         unsigned int v1 = grid.bface[face_no].vertex[1];
         
         int face_type = grid.bface[i].type;
         BoundaryCondition& bc = param.boundary_condition[face_type];
         
         vector<PrimVar> prim_state(2);
         vector<EntVar> ent_state(2);
         prim_state[0] = primitive[v0];
         prim_state[1] = primitive[v1];
         bc.apply(grid.vertex[v0].coord, int_step_time, grid.bface[i], prim_state[0]);
         bc.apply(grid.vertex[v1].coord, int_step_time, grid.bface[i], prim_state[1]);
         for(unsigned int j=0; j<2; ++j)
            ent_state[j] = param.material.prim2ent(prim_state[j]);
         
         // inviscid force   
         double pressure = 0.5 * (prim_state[0].pressure + prim_state[1].pressure);
         inv_force[i].value.sadd(normal,pressure);
         
         // viscous force, using vertex gradients:
         if(param.material.model == Material::ns)
         {
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
			visc_force[i].value.x += -(sxx * normal.x + sxy * normal.y);
			visc_force[i].value.y += -(sxy * normal.x + syy * normal.y);
         }
         


      }
      
      double inv_force_x, inv_force_y, visc_force_x, visc_force_y;
      if(NPART > 1)
      {
		  MPI_Allreduce(&inv_force[i].value.x,&inv_force_x,1,MPI_DOUBLE,MPI_SUM,grid.run_comm);
		  MPI_Allreduce(&inv_force[i].value.y,&inv_force_y,1,MPI_DOUBLE,MPI_SUM,grid.run_comm);
		  MPI_Allreduce(&visc_force[i].value.x,&visc_force_x,1,MPI_DOUBLE,MPI_SUM,grid.run_comm);
		  MPI_Allreduce(&visc_force[i].value.y,&visc_force_y,1,MPI_DOUBLE,MPI_SUM,grid.run_comm);
      }
      if(check_group_base())
      {
		  force_file << inv_force_x << "  " 
					 << inv_force_y << "  "
					 << visc_force_x << "  " 
					 << visc_force_y << "  " ;
      }           
   }
   if(check_group_base())
      force_file << endl;
   
   t_force_eval.add_time();
}
