#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "parameter.h"
#include "fv.h"
#include "writer.h"


extern Dimension dim;
extern bool restart;
extern bool preprocess;
extern bool bounds;
extern bool verbose;
extern std::vector<double> rnd_nos;

using namespace std;

bool check_for_stop_file ();  /*!<Checks for a stop file. If present then simulation ends.*/

//------------------------------------------------------------------------------
// Pre MC run intializations
//------------------------------------------------------------------------------
void FiniteVolume::initialize ()
{  
   t_init.start_time();
   
   if(get_proc_id() == 0)
      cout << "\n  Initializing memory for simulations ...\n";  
   primitive.resize (grid.n_vertex);
   conserved_old.resize (grid.n_vertex);
   residual.resize (grid.n_vertex);
   entropy_var.resize(grid.n_vertex);
   dt.resize (grid.n_vertex);
   phi.resize (grid.n_vertex);
   ssw.resize (grid.n_vertex);

   dE.resize(grid.n_vertex);
   
   dE_cell.resize(grid.n_cell);
   for(unsigned int i=0; i<grid.n_cell; ++i)
   {
      dE_cell[i].resize(NVAR);
      
   }  
   
   dE0.resize(NVAR);
   
   // Create time checkpoints to save variables
   if(param.time_mode == "unsteady" && param.n_time_stamps > 0)
   {
	   save_by_frequency = false;
	   mc_time.resize(param.n_time_stamps+1);
	   double mc_dt = param.final_time/(param.n_time_stamps);
	   for(unsigned int i =0; i<mc_time.size(); ++i)
	   {
		  mc_time[i] = i*mc_dt;
	   }
   }   
   else
   {   
	   save_by_frequency = true;
	   mc_time.resize(2);
	   mc_time[0] = 0.0;
   }
   
   
   if(param.time_scheme == "rk4")
   {
      residual2.resize(grid.n_vertex);
   }
   
   // Initializing MC variables
   if(param.online_stat)
   {
      primitive_fm.resize (mc_time.size());
      primitive_sm.resize (mc_time.size());
      primitive_variance.resize (mc_time.size());
      if(has_density)
      {
         density_fm.resize (mc_time.size());
         density_sm.resize (mc_time.size());
         density_variance.resize (mc_time.size());
      }
      if(has_mach)
      {
         mach_fm.resize (mc_time.size());
         mach_sm.resize (mc_time.size());
         mach_variance.resize (mc_time.size());
      }
      if(has_vorticity)
      {
         vorticity_fm.resize (mc_time.size());
         vorticity_sm.resize (mc_time.size());
         vorticity_variance.resize (mc_time.size());
      }
      if(has_entropy)
      {
         entropy_fm.resize (mc_time.size());
         entropy_sm.resize (mc_time.size());
         entropy_variance.resize (mc_time.size());
      }   
	  for(unsigned int i =0; i<mc_time.size(); ++i)
	  {
		 primitive_fm[i].resize (grid.n_vertex);
		 primitive_sm[i].resize (grid.n_vertex);
		 primitive_variance[i].resize (grid.n_vertex);
		 if(has_density)
		 {
			density_fm[i].resize (grid.n_vertex);
			density_sm[i].resize (grid.n_vertex);
			density_variance[i].resize (grid.n_vertex);
		 }
		 if(has_mach)
		 {
			mach_fm[i].resize (grid.n_vertex);
			mach_sm[i].resize (grid.n_vertex);
			mach_variance[i].resize (grid.n_vertex);
		 }
		 if(has_vorticity)
		 {
			vorticity_fm[i].resize (grid.n_vertex);
			vorticity_sm[i].resize (grid.n_vertex);
			vorticity_variance[i].resize (grid.n_vertex);
		 }
		 if(has_entropy)
		 {
			entropy_fm[i].resize (grid.n_vertex);
			entropy_sm[i].resize (grid.n_vertex);
			entropy_variance[i].resize (grid.n_vertex);
		 }
	  }
   
	  for(unsigned int i =0; i<mc_time.size(); ++i)
	  {
		 for(unsigned int j=0; j<grid.n_vertex; ++j)
		 { 
			primitive_fm[i][j] = 0;
			primitive_sm[i][j] = 0;
			primitive_variance[i][j] = 0;
			if(has_density)
			{
			   density_fm[i][j]   = 0;
			   density_sm[i][j]   = 0;
			   density_variance[i][j]   = 0;
			}
			if(has_mach)
			{
			   mach_fm[i][j]   = 0;
			   mach_sm[i][j]   = 0;
			   mach_variance[i][j]   = 0;
			}
			if(has_vorticity)
			{
			   vorticity_fm[i][j]   = 0;
			   vorticity_sm[i][j]   = 0;
			   vorticity_variance[i][j]   = 0;
			}
			if(has_entropy)
			{
			   entropy_fm[i][j]   = 0;
			   entropy_sm[i][j]   = 0;
			   entropy_variance[i][j]   = 0;
			}
		 }
	  }  
	  
	  prim_mean_l1.resize(mc_time.size());
	  prim_var_l1.resize(mc_time.size());	
   
	  if(has_density)
	  {
		  den_mean_l1.resize(mc_time.size(),0);
		  den_var_l1.resize(mc_time.size(),0);	
	  }	   
	  if(has_mach)
	  {
		  mach_mean_l1.resize(mc_time.size(),0);
		  mach_var_l1.resize(mc_time.size(),0);	
	  }
	  if(has_entropy)
	  {
		  ent_mean_l1.resize(mc_time.size(),0);
		  ent_var_l1.resize(mc_time.size(),0);	
	  }
	  if(has_vorticity)
	  {
		  vor_mean_l1.resize(mc_time.size(),0);
		  vor_var_l1.resize(mc_time.size(),0);	
	  }
   }	   
   
   t_init.add_time();

}

//------------------------------------------------------------------------------
// Intialize for sample run
//------------------------------------------------------------------------------
void FiniteVolume::initialize_sample_run (int sample_id)
{  
   t_sample_init.start_time();
   
   for(unsigned int i=0; i<grid.n_vertex; ++i)
   {
      dE[i].resize(NVAR);
      for(unsigned int j=0; j<NVAR; ++j)
         dE[i][j] = 0;
   }
   
   // Making sample directory
   SAMPLE_DIR = "SAMPLE_";
   stringstream ss;
   ss <<param.SAMPLE_LIST[sample_id];
   SAMPLE_DIR += ss.str();
   string res_file_path = SAMPLE_DIR + "/residual.dat";
   string iter_file_path = SAMPLE_DIR + "/stat_iters.dat";
   string force_file_path = SAMPLE_DIR + "/force.dat";
   string global_file_path = SAMPLE_DIR + "/global.dat";
   string tag_file = SAMPLE_DIR + "/LAST_RESTART_HEAD.txt";
   
   if(check_group_base())
   {
      if(!restart)
      {
		  string commandline0 = "rm -rf " + SAMPLE_DIR;
		  system(commandline0.c_str());
      }
      string commandline1 = "mkdir -p " + SAMPLE_DIR;
      system(commandline1.c_str());
   }   
   mpi_barrier(grid.run_comm);
   // If we are in preprocess mode, then dont open
   // these files. This avoids overwriting existing 
   // files.
   if(!preprocess)
   {
	    
		if(check_group_base())
		{ 	   
		   // MASTER FILE FOR VIEWING SOLUTION
		   if(param.write_soln)
		   {
              if(!restart)
              {
				 string master_file_path = SAMPLE_DIR + "/master_file.visit";
				 master_file.open (master_file_path.c_str());
				 master_file << "!NBLOCKS "<<NPART<<std::endl;
				 master_file.close ();
			  }
			  write_to_master = true;
			  
		   }	  
		   
		   if(bounds)
		   {
		      string bound_file_path = SAMPLE_DIR + "/bounds.dat";
              bound_file.open(bound_file_path.c_str(), ios::app);
              if(!restart)
              {
				  bound_file << left << setw(16*11) << setfill('-') <<""<<endl;
				  bound_file << left << setw(17) << setfill(' ') << "TIME"
							 << left << setw(17) << setfill(' ') << "TEMP_MIN"
							 << left << setw(17) << setfill(' ') <<"TEMP_MAX"
							 << left << setw(17) << setfill(' ') <<"P_MIN"
							 << left << setw(17) << setfill(' ') <<"P_MAX"
							 << left << setw(17) << setfill(' ') <<"U_MIN"
							 << left << setw(17) << setfill(' ') <<"U_MAX"
							 << left << setw(17) << setfill(' ') <<"V_MIN"
							 << left << setw(17) << setfill(' ') <<"V_MAX"
							 << left << setw(17) << setfill(' ') <<"W_MIN"
							 << left << setw(17) << setfill(' ') <<"W_MAX"<<endl;
				 bound_file << left << setw(16*11) << setfill('-') << ""<<endl;            
             }
		   }
		}
		mpi_barrier(grid.run_comm);
   }

 
   // For zero entropy gradients
   for (unsigned int i=0; i<NVAR; i++)
	   dE0[i] = 0.0;

   //set_seed(sample_id);
   
   if(param.rnd_per_sample > 0)
      extract_rnd_nos(sample_id);
      
   elapsed_time = 0.0;   
   
   ifstream fi_tag;
   bool restart_available = false;
   int head_tag;
   if(restart)
   {
	  fi_tag.open(tag_file.c_str());
	  if(fi_tag.is_open())
	  {
		 fi_tag>>head_tag;
		 fi_tag.close();
		 restart_available = true;
		 if(check_group_base())
            cout << "    Reading restart file(s) for "<<SAMPLE_DIR<<"  ...\n";
	  }
	  else
	  {
	     if(check_group_base())
            cout << "    Restart file(s) for "<<SAMPLE_DIR<<" does not exist."<<
            " Staring from T=0 ...\n";  
	  }
   }
   
      
   bool continue_solve = false;
   // If restart option specified, read previous solution from file
   if(restart && restart_available)
   {  
      ifstream fi;
      string restart_file_path = SAMPLE_DIR + "/restart_";
      stringstream ss;
      ss <<head_tag<<"_"<<get_proc_loc_id();
      restart_file_path += ss.str()+".dat";
      fi.open (restart_file_path.c_str());
      MPI_LOC_ASSERT (fi.is_open());
      for(unsigned int i=0; i<grid.n_vertex; ++i)
         fi >> primitive[i].temperature
            >> primitive[i].velocity.x
            >> primitive[i].velocity.y
            >> primitive[i].velocity.z
            >> primitive[i].pressure;
      fi >> last_iter;
      fi >> counter;
      fi >> elapsed_time;
      fi >> residual_norm_total0;
      fi.close ();
      if(check_group_base())
      {     
         if(elapsed_time >= param.final_time)
            cout<<"    SAMPLE "<<param.SAMPLE_LIST[sample_id]<<" has already reached final time!!"<<endl;
         else if(last_iter >= param.max_iter && param.time_mode == "steady")
            cout<<"    SAMPLE "<<param.SAMPLE_LIST[sample_id]<<" has already max iteration count!!"<<endl;   
         else if(elapsed_time > 0.0)
         {
            stringstream ss1;
            ss1 << last_iter;
            if(param.write_log)
            {
               string res_file_dummy_path = SAMPLE_DIR + "/residual_dummy.dat";
               string commandline0 = "awk '$1 <= " + ss1.str() + "' "
                                     + res_file_path.c_str() + " > " + res_file_dummy_path.c_str()
                                     + " && mv " + res_file_dummy_path.c_str() + " "
                                     + res_file_path.c_str();
               system(commandline0.c_str());                      
                                     
            }
            if(param.force_data.size() > 0)
            {
               string force_file_dummy_path = SAMPLE_DIR + "/force_dummy.dat";
               string commandline0 = "awk '$1 <= " + ss1.str() + "' "
                                     + force_file_path.c_str() + " > " + force_file_dummy_path.c_str()
                                     + " && mv " + force_file_dummy_path.c_str() + " "
                                     + force_file_path.c_str();
               system(commandline0.c_str()); 
            }
            if(param.has_global)
            {
               string global_file_dummy_path = SAMPLE_DIR + "/global_dummy.dat";
               string commandline0 = "awk '$1 <= " + ss1.str() + "' "
                                     + global_file_path.c_str() + " > " + global_file_dummy_path.c_str()
                                     + " && mv " + global_file_dummy_path.c_str() + " "
                                     + global_file_path.c_str();
               system(commandline0.c_str()); 
            }
            continue_solve = true;
    
         }     
      }        
   }
   else
   {
      for(unsigned int i=0; i<grid.n_vertex; ++i)
      {
         primitive[i] = param.initial_condition.value (grid.vertex[i].coord);
         MPI_LOC_ASSERT (primitive[i].temperature  > 0.0);
         MPI_LOC_ASSERT (primitive[i].pressure > 0.0);
      }
      last_iter = 0;
   }
   mpi_barrier(grid.run_comm);
   restart_head_tag = 0;
   for(unsigned int i=0; i<grid.n_vertex; ++i)
      entropy_var[i] = param.material.prim2ent(primitive[i]); 
      
      
   // RESIDUAL FILE
   if(check_group_base())
   {
	  if(param.write_log)
	  {
		 res_file.open(res_file_path.c_str(), ios::app);
		 if(continue_solve)
			res_file << "-------------------------------------------------------------"<<endl;  
	  }
   
	  // ITERATION LIST FOR SAMPLES USED TO EVAL MEAN AND VARIANCE
	  // if(param.write_soln)
// 	  {
// 		 iter_file.open(iter_file_path.c_str(), ios::app);
// 	  }

	  if(param.force_data.size() > 0)
	  {
		 force_file.open (force_file_path.c_str(), ios::app);		     
	  }	  

	  if(param.has_global == true)
	  {
		 global_file.open (global_file_path.c_str(), ios::app);
	  }  
	  
	  if(param.find_error)
	  {
		 density_exact = 0.0, entropy_exact = 0.0, mach_exact = 0.0,
		 density_abs_err = 0.0, entropy_abs_err = 0.0, mach_abs_err = 0.0,
		 density_err_L1 = 0.0, entropy_err_L1 = 0.0, mach_err_L1 = 0.0,
		 density_err_L2 = 0.0, entropy_err_L2 = 0.0, mach_err_L2 = 0.0,
		 density_err_Linf = 0.0, entropy_err_Linf = 0.0, mach_err_Linf = 0.0;
	  
		 string err_file_path = SAMPLE_DIR + "/error_norms.dat";
		 err_file.open (err_file_path.c_str());                   
	  }  
   } 
   mpi_barrier(grid.run_comm);
   
   t_sample_init.add_time();   
}

//------------------------------------------------------------------------------
// Set seed for perturbation
//------------------------------------------------------------------------------
// void FiniteVolume::set_seed(int sample_id)
// {  
//     seed = sample_id*time(0);
// }

//------------------------------------------------------------------------------
// Extracts random number from file
//------------------------------------------------------------------------------
// void FiniteVolume::extract_rnd_nos(int sample_id)
// {  
//     ifstream rnd_file;
//     rnd_file.open(param.rnd_file_loc.c_str());
//     int total_rnd;
//     rnd_file >> total_rnd;
//     int rel_hook = param.rnd_file_hook % total_rnd;
//     int start_pos = (rel_hook + param.SAMPLE_LIST[sample_id]*param.rnd_per_sample) % total_rnd;
//     int current_pos;
//     
//     double dummy;
//     for(current_pos=1; current_pos<start_pos; ++current_pos)
//     {
//        rnd_file >> dummy;
//     }   
//        
//     unsigned int rnd_ind = 0;
//     rnd_nos.clear();
//     while(rnd_ind < param.rnd_per_sample)
//     {
//        if(current_pos <= total_rnd)
//        {
//           rnd_file >> dummy;
//           rnd_nos.push_back(dummy);
//           rnd_ind++;
//        }  
//        else
//        {
//           rnd_file.clear();
//           rnd_file.seekg(0, ios::beg);
//           rnd_file >> dummy;
//           current_pos = 1;
//           rnd_file >> dummy;
//           rnd_nos.push_back(dummy);
//           rnd_ind++;
//        }
//     }
//     
//     rnd_file.close();
//     // cout<<"FOR SAMPLE :"<<sample_id<<endl;
// //     for(unsigned int i=0; i<rnd_nos.size(); ++i)
// //       cout<<rnd_nos[i]<<" ";
// //     cout<<std::endl; 
//     
// }

//------------------------------------------------------------------------------
// Extracts random number from file
//------------------------------------------------------------------------------
void FiniteVolume::extract_rnd_nos(int sample_id)
{  
    ifstream rnd_file;
    rnd_file.open(param.rnd_file_loc.c_str());
    int total_samples, sid;
    double dummy;
    rnd_file >> total_samples;
    bool found = false;
    rnd_nos.clear();
    rnd_nos.resize(param.rnd_per_sample,0);
    
//     if(check_group_base())
//     {
//         cout<<"FOR SAMPLE :"<<param.SAMPLE_LIST[sample_id]<<" "<<found<<endl;
//         for(unsigned int i=0; i<rnd_nos.size(); ++i)
//              cout<<rnd_nos[i]<<", ";
//         cout<<endl; 
//     }
    
    int count = 1;
    while(count<= total_samples && !found)
    {
        rnd_file >> sid;
        if(sid == param.SAMPLE_LIST[sample_id])
        {
            for(int i=0; i<param.rnd_per_sample; i++)
                rnd_file >> rnd_nos[i];
            found = true;  
        }      
        else
            for(int i=0; i<param.rnd_per_sample; i++)
                rnd_file >> dummy;
        count++;                  
    }
    
    rnd_file.close();
    
//     if(check_group_base())
//     {
//         cout<<"FOR SAMPLE :"<<param.SAMPLE_LIST[sample_id]<<" "<<found<<endl;
//         for(unsigned int i=0; i<rnd_nos.size(); ++i)
//              cout<<rnd_nos[i]<<", ";
//         cout<<endl; 
//     }
    
    if(!found)
        MPI_LOC_ERR("Sample "<<param.SAMPLE_LIST[sample_id]<<" not present available in "<<param.rnd_file_loc.c_str()<<"\n");

    
}

//------------------------------------------------------------------------------
// End sample run (closing associated files)
//------------------------------------------------------------------------------
void FiniteVolume::end_sample_run()
{  
   if(check_group_base())
   {
	   
	   if(param.write_log)
	      res_file.close();
	   // iter_file.close();
	   if(param.force_data.size() > 0)
		  force_file.close ();
	   if(global_file.is_open())
		  global_file.close ();
	   if(bound_file.is_open())
	      bound_file.close();	  
   } 
   write_sample_time_output();  
   t_sample_init.reset();
   t_sample_solve.reset();
   t_grad_eval.reset();
   t_grad_red.reset();
   t_dt_eval.reset();
   t_dt_red.reset();
   t_res_eval.reset();
   t_ires_eval.reset();
   t_vres_eval.reset();
   t_res_red.reset();
   t_res_norm.reset();
   t_res_norm_red.reset();
   t_sol_update.reset();
   t_sol_out.reset();
   t_force_eval.reset();
   t_global_eval.reset();
   t_log.reset();
}

//------------------------------------------------------------------------------
// Increment moments at given time 
//------------------------------------------------------------------------------
void FiniteVolume::increment_moments(const int sample_id)
{  
   t_mom_inc.start_time();
   // Adding moments
   for(unsigned int i=0; i<grid.n_vertex; ++i)
   {
      if(sample_id == param.sample_start_ind)
      {
         primitive_fm[mc_t_ind][i] = primitive[i];
         primitive_sm[mc_t_ind][i].cdot(primitive[i],primitive[i]);
         primitive_variance[mc_t_ind][i] = 0.0;
         if(has_density)
		 {
			double val = param.material.Density(primitive[i]);
			density_fm[mc_t_ind][i]   = val;
			density_sm[mc_t_ind][i] = val*val;
			density_variance[mc_t_ind][i]  = 0.0;
		 }
		 if(has_mach)
		 {
			double val = param.material.Mach(primitive[i]);
			mach_fm[mc_t_ind][i]   = val;
			mach_sm[mc_t_ind][i] = val*val;
			mach_variance[mc_t_ind][i]   = 0;
		 }
		 if(has_vorticity)
		 {
			double val = param.material.Vorticity(primitive[i],dE[i]);
			vorticity_fm[mc_t_ind][i]   = val;
			vorticity_sm[mc_t_ind][i] = val*val;
			vorticity_variance[mc_t_ind][i]   = 0;
		 }
		 if(has_entropy)
		 {
			double val = param.material.Entropy(primitive[i]);
			entropy_fm[mc_t_ind][i]   = val;
			entropy_sm[mc_t_ind][i] = val*val;
			entropy_variance[mc_t_ind][i]   = 0;
		 }
      }
      else
      {
         PrimVar Mk, Mkm1, dummy1, dummy2;
         double sMk, sMkm1, val;
         double k = sample_id-param.sample_start_ind+1;
         Mkm1 = primitive_fm[mc_t_ind][i];
    	 Mk.equ(primitive[i],Mkm1,1.0/k,(k-1.0)/k);
    	 primitive_fm[mc_t_ind][i] = Mk;  
    	 dummy1.equ(primitive[i],Mkm1,1.0,-1.0);
    	 dummy2.equ(primitive[i],Mk,1.0,-1.0);
    	 primitive_sm[mc_t_ind][i].addcdot(primitive[i],primitive[i]);
	     primitive_variance[mc_t_ind][i].addcdot(dummy1,dummy2);
	     
	     if(has_density)
	     {
	        val =  param.material.Density(primitive[i]);
	        sMkm1 = density_fm[mc_t_ind][i];
			sMk   = sMkm1 + (val - sMkm1)/k;
            density_fm[mc_t_ind][i] = sMk;
            density_sm[mc_t_ind][i]+=val*val;  
			density_variance[mc_t_ind][i] += (val - sMkm1)*(val - sMk);
	     }
	     if(has_mach)
	     {
	        val =  param.material.Mach(primitive[i]);
	        sMkm1 = mach_fm[mc_t_ind][i];
			sMk   = sMkm1 + (val - sMkm1)/k;
            mach_fm[mc_t_ind][i] = sMk;  
            mach_sm[mc_t_ind][i]+=val*val;
			mach_variance[mc_t_ind][i] += (val - sMkm1)*(val - sMk);
	     }
	     if(has_vorticity)
	     {
	        val =  param.material.Vorticity(primitive[i],dE[i]);
	        sMkm1 = vorticity_fm[mc_t_ind][i];
			sMk   = sMkm1 + (val - sMkm1)/k;
            vorticity_fm[mc_t_ind][i] = sMk;  
            vorticity_sm[mc_t_ind][i]+=val*val;
			vorticity_variance[mc_t_ind][i] += (val - sMkm1)*(val - sMk);
	     }
	     if(has_entropy)
	     {
	        val =  param.material.Entropy(primitive[i]);
	        sMkm1 = entropy_fm[mc_t_ind][i];
			sMk   = sMkm1 + (val - sMkm1)/k;
            entropy_fm[mc_t_ind][i] = sMk;  
            entropy_sm[mc_t_ind][i]+=val*val;
			entropy_variance[mc_t_ind][i] += (val - sMkm1)*(val - sMk);
	     }  
      }
   } 
   
   t_mom_inc.add_time();

}

//------------------------------------------------------------------------------
// Compute derivatives of entropy variables
//------------------------------------------------------------------------------
void FiniteVolume::compute_gradients ()
{
   t_grad_eval.start_time();
   
   vector<PrimVar> state(2);
   vector<EntVar> ent_state(2);

   for(unsigned int i=0; i<grid.n_vertex; ++i)
      for(unsigned int j=0; j< NVAR; ++j)
        dE[i][j] = 0.0;       
   
   // EVALUATION FOR CELLS WITH VERTICES SHARED BY OTHER PROCESSORS
   for(unsigned int it = 0; it<grid.cell_mpi_bnd.size(); ++it)
   {
      
      int i = grid.cell_mpi_bnd[it];
		  
	  unsigned int v0 = grid.cell[i].vertex[0];
	  unsigned int v1 = grid.cell[i].vertex[1];
	  unsigned int v2 = grid.cell[i].vertex[2];

	  Vector& normal0 = grid.cell[i].normal[0];
	  Vector& normal1 = grid.cell[i].normal[1];
	  Vector& normal2 = grid.cell[i].normal[2];

	  // E gradient
	  dE_cell[i][0].equ(normal0, normal1, normal2, 
						   entropy_var[v0].entf,
						   entropy_var[v1].entf,
						   entropy_var[v2].entf);
	  dE_cell[i][1].equ(normal0, normal1, normal2, 
						   entropy_var[v0].entv.x,
						   entropy_var[v1].entv.x,
						   entropy_var[v2].entv.x); 
	 
	  dE_cell[i][2].equ(normal0, normal1, normal2, 
						   entropy_var[v0].entv.y,
						   entropy_var[v1].entv.y,
						   entropy_var[v2].entv.y); 
											   
	  dE_cell[i][3].equ(normal0, normal1, normal2, 
						   entropy_var[v0].entv.z,
						   entropy_var[v1].entv.z,
						   entropy_var[v2].entv.z);                                                

	  dE_cell[i][4].equ(normal0, normal1, normal2, 
						   entropy_var[v0].entl,
						   entropy_var[v1].entl,
						   entropy_var[v2].entl);                        

   }
   
   // FIX THIS FOR GENERAL BOUNDARY TRIANGLES ACCORDING TO NOTES
   for(unsigned int it=0; it<grid.d_bnd_shared_face.size(); ++it) 
   {
      int i = grid.d_bnd_shared_face[it];
      
      MPI_LOC_ASSERT(grid.cell[grid.bface[i].lcell].mpi_bnd_active);
      int face_type = grid.bface[i].type;
      if(!grid.bface[i].periodic)
      {
		 BoundaryCondition& bc = param.boundary_condition[face_type];

		 unsigned int v0 = grid.bface[i].vertex[0];
		 unsigned int v1 = grid.bface[i].vertex[1];
		 unsigned int cl = grid.bface[i].lcell;

		 state[0] = primitive[v0];
		 state[1] = primitive[v1];
		 ent_state[0] = entropy_var[v0];
		 ent_state[1] = entropy_var[v1];

		 bc.apply(grid.vertex[v0].coord, int_step_time, grid.bface[i], state[0]);
		 bc.apply(grid.vertex[v1].coord, int_step_time, grid.bface[i], state[1]);
		   
		 for(unsigned int j=0; j<2; ++j)
			ent_state[j] = param.material.prim2ent(state[j]);


		 dE_cell[cl][0].sadd(grid.bface[i].normal,
							 (ent_state[0].entf + ent_state[1].entf
							  - entropy_var[v0].entf - entropy_var[v1].entf)/2.0);
		 dE_cell[cl][1].sadd(grid.bface[i].normal,
							 (ent_state[0].entv.x + ent_state[1].entv.x
							  - entropy_var[v0].entv.x - entropy_var[v1].entv.x)/2.0);
		 dE_cell[cl][2].sadd(grid.bface[i].normal,
							 (ent_state[0].entv.y + ent_state[1].entv.y
							  - entropy_var[v0].entv.y - entropy_var[v1].entv.y)/2.0);
		 dE_cell[cl][3].sadd(grid.bface[i].normal, 
							 (ent_state[0].entv.z + ent_state[1].entv.z
							  - entropy_var[v0].entv.z - entropy_var[v1].entv.z)/2.0);               
		 dE_cell[cl][4].sadd(grid.bface[i].normal, 
							 (ent_state[0].entl + ent_state[1].entl
							  - entropy_var[v0].entl - entropy_var[v1].entl)/2.0);  
      }                            
   }
   
   // ADDING CONTRIBUTIONS TO SHARED VERTICES
   for(unsigned int it = 0; it<grid.cell_mpi_bnd.size(); ++it)
   {
      
      int i = grid.cell_mpi_bnd[it];
	  unsigned int v0 = grid.cell[i].vertex[0];
	  unsigned int v1 = grid.cell[i].vertex[1];
	  unsigned int v2 = grid.cell[i].vertex[2];

	  dE[v0][0] += dE_cell[i][0];
	  dE[v1][0] += dE_cell[i][0];
	  dE[v2][0] += dE_cell[i][0];

	  dE[v0][1] += dE_cell[i][1];
	  dE[v1][1] += dE_cell[i][1];
	  dE[v2][1] += dE_cell[i][1];
  
	  dE[v0][2] += dE_cell[i][2];
	  dE[v1][2] += dE_cell[i][2];
	  dE[v2][2] += dE_cell[i][2];
  
	  dE[v0][3] += dE_cell[i][3];
	  dE[v1][3] += dE_cell[i][3];
	  dE[v2][3] += dE_cell[i][3];
  
	  dE[v0][4] += dE_cell[i][4];
	  dE[v1][4] += dE_cell[i][4];
	  dE[v2][4] += dE_cell[i][4];

   }
   
   // STORING SHARED VERTEX GRADIENTS IN BUFFERS AND PERFORMING NON-BLOCKING
   // SENDS AND RECEIVES
   t_grad_red.start_time();
   int gsize = grid.intp_mpi_groupings.size();
   int vsize[gsize];
   vector< vector <MPI_Request> > s_request, r_request;
   vector< VARARRAY > SBUF;
   vector< vector< VARARRAY > > RBUF;
   int buff_size[gsize];
   if(NPART > 1)
   {
	  SBUF.resize(gsize);
	  RBUF.resize(gsize);
	  s_request.resize(gsize);
	  r_request.resize(gsize);
   
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
		 buff_size[i] = vsize[i]*NVAR*NDIR;
		 SBUF[i].buf = new double[buff_size[i]];
		 int ind = 0;
		 for(int j = 0; j<vsize[i]; j++)
		 {
			int locv=grid.g2l.find(grid.intp_mpi_groupings[i].vertices[j])->second;
			for(int k=0;k<NVAR;k++)
			{
			   SBUF[i].buf[ind++] = dE[locv][k].x; 
			   SBUF[i].buf[ind++] = dE[locv][k].y; 
			   SBUF[i].buf[ind++] = dE[locv][k].z;
			}   
		 }
		 int sr_ind = 0;
		 for(set<int>::iterator it=grid.intp_mpi_groupings[i].proc_list.begin();
			 it!=grid.intp_mpi_groupings[i].proc_list.end(); ++it)
		 {    
			if(*it!=get_proc_id())
			{
				MPI_Isend(&SBUF[i].buf[0], buff_size[i], MPI_DOUBLE,*it,
						  grid.intp_mpi_groupings[i].tag,MPI_COMM_WORLD,&s_request[i][sr_ind]);
				RBUF[i][sr_ind].buf = new double[buff_size[i]];
				MPI_Irecv(&RBUF[i][sr_ind].buf[0], buff_size[i], MPI_DOUBLE,*it,
						  grid.intp_mpi_groupings[i].tag,MPI_COMM_WORLD,&r_request[i][sr_ind]);          
				sr_ind++;  
			}                    
		 }   
	  }
	  t_grad_red.add_time();
   }
   
   // EVALUATION FOR INTERIOR CELLS
   for(unsigned int it = 0; it<grid.cell_mpi_int.size(); ++it)
   {
      int i = grid.cell_mpi_int[it];
		  
	  unsigned int v0 = grid.cell[i].vertex[0];
	  unsigned int v1 = grid.cell[i].vertex[1];
	  unsigned int v2 = grid.cell[i].vertex[2];

	  Vector& normal0 = grid.cell[i].normal[0];
	  Vector& normal1 = grid.cell[i].normal[1];
	  Vector& normal2 = grid.cell[i].normal[2];

	  // E gradient
	  dE_cell[i][0].equ(normal0, normal1, normal2, 
						   entropy_var[v0].entf,
						   entropy_var[v1].entf,
						   entropy_var[v2].entf);
	  dE_cell[i][1].equ(normal0, normal1, normal2, 
						   entropy_var[v0].entv.x,
						   entropy_var[v1].entv.x,
						   entropy_var[v2].entv.x); 
	 
	  dE_cell[i][2].equ(normal0, normal1, normal2, 
						   entropy_var[v0].entv.y,
						   entropy_var[v1].entv.y,
						   entropy_var[v2].entv.y); 
											   
	  dE_cell[i][3].equ(normal0, normal1, normal2, 
						   entropy_var[v0].entv.z,
						   entropy_var[v1].entv.z,
						   entropy_var[v2].entv.z);                                                

	  dE_cell[i][4].equ(normal0, normal1, normal2, 
						   entropy_var[v0].entl,
						   entropy_var[v1].entl,
						   entropy_var[v2].entl);                        

   }
      
   for(unsigned int it=0; it<grid.d_bnd_unshared_face.size(); ++it)
   {
      int i = grid.d_bnd_unshared_face[it];
      int face_type = grid.bface[i].type;
      if(!grid.bface[i].periodic)
      {
		 BoundaryCondition& bc = param.boundary_condition[face_type];

		 unsigned int v0 = grid.bface[i].vertex[0];
		 unsigned int v1 = grid.bface[i].vertex[1];     
		 unsigned int cl = grid.bface[i].lcell;

		 state[0] = primitive[v0];
		 state[1] = primitive[v1];
		 ent_state[0] = entropy_var[v0];
		 ent_state[1] = entropy_var[v1];

		 bc.apply(grid.vertex[v0].coord, int_step_time, grid.bface[i], state[0]);
		 bc.apply(grid.vertex[v1].coord, int_step_time, grid.bface[i], state[1]);
		   
		 for(unsigned int j=0; j<2; ++j)
			ent_state[j] = param.material.prim2ent(state[j]);


		 dE_cell[cl][0].sadd(grid.bface[i].normal,
							 (ent_state[0].entf + ent_state[1].entf
							  - entropy_var[v0].entf - entropy_var[v1].entf)/2.0);
		 dE_cell[cl][1].sadd(grid.bface[i].normal,
							 (ent_state[0].entv.x + ent_state[1].entv.x
							  - entropy_var[v0].entv.x - entropy_var[v1].entv.x)/2.0);
		 dE_cell[cl][2].sadd(grid.bface[i].normal,
							 (ent_state[0].entv.y + ent_state[1].entv.y
							  - entropy_var[v0].entv.y - entropy_var[v1].entv.y)/2.0);
		 dE_cell[cl][3].sadd(grid.bface[i].normal, 
							 (ent_state[0].entv.z + ent_state[1].entv.z
							  - entropy_var[v0].entv.z - entropy_var[v1].entv.z)/2.0);               
		 dE_cell[cl][4].sadd(grid.bface[i].normal, 
							 (ent_state[0].entl + ent_state[1].entl
							  - entropy_var[v0].entl - entropy_var[v1].entl)/2.0);  
      }                            
   }
   
   // ADDING CONTRIBUTIONS TO INTERIOR VERTICES
   for(unsigned int it = 0; it<grid.cell_mpi_int.size(); ++it)
   {
      
      int i = grid.cell_mpi_int[it];
      unsigned int v0 = grid.cell[i].vertex[0];
	  unsigned int v1 = grid.cell[i].vertex[1];
	  unsigned int v2 = grid.cell[i].vertex[2];
	  
	  dE[v0][0] += dE_cell[i][0];
	  dE[v1][0] += dE_cell[i][0];
	  dE[v2][0] += dE_cell[i][0];

	  dE[v0][1] += dE_cell[i][1];
	  dE[v1][1] += dE_cell[i][1];
	  dE[v2][1] += dE_cell[i][1];
  
	  dE[v0][2] += dE_cell[i][2];
	  dE[v1][2] += dE_cell[i][2];
	  dE[v2][2] += dE_cell[i][2];
  
	  dE[v0][3] += dE_cell[i][3];
	  dE[v1][3] += dE_cell[i][3];
	  dE[v2][3] += dE_cell[i][3];
  
	  dE[v0][4] += dE_cell[i][4];
	  dE[v1][4] += dE_cell[i][4];
	  dE[v2][4] += dE_cell[i][4];
   }


   if(NPART > 1)
   {
	  // CHECKING IF NON-BLOCKING SENDS AND RECEIVES ARE COMPLETE AND ADDING CONTRIBUTIONS
	  t_grad_red.start_time();
	  MPI_Status status;
	  for(int i=0; i<gsize; ++i)
	  {
		 for(int r=0; r<r_request[i].size(); ++r)
		 {
			MPI_Wait(&r_request[i][r],&status);   
			int ind = 0;
			for(int j = 0; j<vsize[i]; j++)
			{
			   int locv=grid.g2l.find(grid.intp_mpi_groupings[i].vertices[j])->second;
			   for(int k=0;k<NVAR;k++)
			   {
				  dE[locv][k].x += RBUF[i][r].buf[ind++];
				  dE[locv][k].y += RBUF[i][r].buf[ind++]; 
				  dE[locv][k].z += RBUF[i][r].buf[ind++];
			   }   
			} 
		 }
		 RBUF[i].clear();
	  }
	  RBUF.clear();
   
	  for(int i=0; i<gsize; i++)
		 for(int r=0; r<s_request[i].size(); ++r)
			MPI_Wait(&s_request[i][r],&status);   
	  SBUF.clear();      
	  t_grad_red.add_time();
   }

   if(PERIODIC)
      reduce_gradient_periodic();

   // SCALING VERTEX GRADIENTS
   for(unsigned int i = 0; i<grid.n_vertex; ++i)
   {
      MPI_LOC_ASSERT(grid.mcarea[i]>0);
	  double fact = 1.0/(6.0 * grid.mcarea[i]);
	  dE[i][0] *= fact;
	  dE[i][1] *= fact;
	  dE[i][2] *= fact; 
	  dE[i][3] *= fact;
	  dE[i][4] *= fact;
   }

   // SCALING CELL GRADIENTS
   for(unsigned int i = 0; i<grid.n_cell; ++i)
   {
	  double fact = 1.0/(2.0 * grid.cell[i].area);
	  dE_cell[i][0] *= fact;
	  dE_cell[i][1] *= fact;
	  dE_cell[i][2] *= fact;
	  dE_cell[i][3] *= fact;
	  dE_cell[i][4] *= fact;
   }
 
   t_grad_eval.add_time();
}



//------------------------------------------------------------------------------
// Compute inviscid residual for each cell
//------------------------------------------------------------------------------
void FiniteVolume::compute_inviscid_residual ()
{
   t_ires_eval.start_time();
   
   // Loop over interior faces and accumulate flux   
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      vector<PrimVar> state(2);

      unsigned int cl = grid.face[i].vertex[0];
      unsigned int cr = grid.face[i].vertex[1];
      MPI_LOC_ASSERT(grid.vertex[cl].active);
      MPI_LOC_ASSERT(grid.vertex[cr].active);

      FluxData data;
      data.ssw = ssw[cl]+ssw[cr];      

      Flux flux;

      param.material.num_flux ( grid.vertex[cl].coord, grid.vertex[cr].coord,
                                   primitive[cl], primitive[cr], 
                                   state[0], state[1], 
                                   dE[cl], dE[cr],
                                   grid.face[i].normal, data, flux);
          
      
      double ip_bc = 1;
      if(grid.face[i].type == -2)
         ip_bc = 0.5;
      residual[cl].sadd(flux,grid.face[i].radius*ip_bc);
      residual[cr].sadd(flux,-grid.face[i].radius*ip_bc);
   }
   
   // Loop over boundary faces and accumulate flux
   for(unsigned int i=0; i<grid.bface.size(); ++i)
   {
      if(!grid.bface[i].periodic)
      {
		 vector<PrimVar> state(2);
		 Flux flux;

		 int face_type = grid.bface[i].type;
		 BoundaryCondition& bc = param.boundary_condition[face_type];
		 //cout<<"CHECK"<<endl;

		 unsigned int cl = grid.bface[i].vertex[0];
		 unsigned int cr = grid.bface[i].vertex[1];
		 MPI_LOC_ASSERT(grid.vertex[cl].active);
		 MPI_LOC_ASSERT(grid.vertex[cr].active);

		 FluxData data;
		 data.ssw = 0;

		 state[0] = primitive[cl];
		 state[1] = primitive[cl];
	  
		 Vector vel = state[0].velocity;
		 double pre = state[0].pressure;
		 double ent = -param.material.Density(state[0])*param.material.Entropy(state[0])/(param.material.gamma - 1.0);
	  
		 //t_inv_bc_apply.start_time();
		 bc.apply (grid.vertex[cl].coord, int_step_time, grid.bface[i], state);
		 //t_inv_bc_apply.add_time();
	  
		 // OLD METHOD---------------------------------------------------------------------------------
		 if(bc.type == BC::farfield)
			param.material.num_flux ( grid.vertex[cl].coord, grid.vertex[cl].coord,
									  state[0], state[1], 
									  state[0], state[1], 
									  dE0,dE0,
									  grid.bface[i].normal, data, flux);
		 else
		 {
			param.material.inviscid_bflux(state[0],state[1],grid.bface[i].normal,flux);
		 }  
		 //---------------------------------------------------------------------------------------------
	  
		 //NEW METHOD-----------------------------------------------------------------------------------
		 //Still needs to be tested
		 //bc.inviscid_states(grid.bface[i],state,vel,pre,ent);
		 //param.material.inviscid_bflux(primitive[cl],vel,pre,ent,grid.bface[i].normal,flux);
		 //---------------------------------------------------------------------------------------------  
   
	  
		 //residual[cl] += flux * (0.5 * grid.vertex[cl].radius);
		 residual[cl].sadd(flux,(0.5 * grid.vertex[cl].radius));

		 state[0] = primitive[cr];
		 state[1] = primitive[cr];
		 vel = state[0].velocity;
		 pre = state[0].pressure;
		 ent = -param.material.Density(state[0])*param.material.Entropy(state[0])/(param.material.gamma - 1.0);
	  
		 //t_inv_bc_apply.start_time();
		 bc.apply (grid.vertex[cr].coord, int_step_time, grid.bface[i], state);
		 //t_inv_bc_apply.add_time();
	  
		 //t_inv_bflux.start_time();
	  
		 // OLD METHOD---------------------------------------------------------------------------------
		 if(bc.type == BC::farfield)
			param.material.num_flux ( grid.vertex[cr].coord, grid.vertex[cr].coord,
									  state[0], state[1], 
									  state[0], state[1], 
									  dE0,dE0,
									  grid.bface[i].normal, data, flux);
		 else
		 {
			param.material.inviscid_bflux(state[0],state[1],grid.bface[i].normal,flux);
		 }   
		 //---------------------------------------------------------------------------------------------  
	   
		 //NEW METHOD-----------------------------------------------------------------------------------
		 //Still needs to be tested 
		 //bc.inviscid_states(grid.bface[i],state,vel,pre,ent);
		 //param.material.inviscid_bflux(primitive[cr],vel,pre,ent,grid.bface[i].normal,flux); 
		 //---------------------------------------------------------------------------------------------   
	   
		 //t_inv_bflux.add_time();     
		 
		 //residual[cr] += flux * (0.5 * grid.vertex[cr].radius);
		 residual[cr].sadd(flux,(0.5 * grid.vertex[cr].radius));
      }
   }
   
   t_ires_eval.add_time();
}

//------------------------------------------------------------------------------
// Compute viscous residual for each cell
//------------------------------------------------------------------------------
void FiniteVolume::compute_viscous_residual ()
{
   t_vres_eval.start_time();
   
   for(unsigned int i = 0; i<grid.n_cell; ++i)
   {
      unsigned int v0 = grid.cell[i].vertex[0];
      unsigned int v1 = grid.cell[i].vertex[1];
      unsigned int v2 = grid.cell[i].vertex[2];
      
	  if(grid.cell[i].nbr_bface.size() == 0)//(grid.cell[i].type == 0)// if it is an interior triangle
	  										                          // or boundary triangles with no  
	  										                          // non-periodic bface
	  {   
		 EntVar ent_avg;
		 double fact = 1.0/3.0;
		 ent_avg.equ(entropy_var[v0],entropy_var[v1],entropy_var[v2],fact,fact,fact);
		 Flux flux0, flux1, flux2;
		 Vector& normal0 = grid.cell[i].normal[0]; 
		 Vector& normal1 = grid.cell[i].normal[1];
		 Vector& normal2 = grid.cell[i].normal[2];

		 param.material.viscous_flux (ent_avg, 
									  dE_cell[i], 
									  normal0, flux0,
									  normal1, flux1,
									  normal2, flux2);                                  
						   
		 // negative sign as the inward normals are used                          
		 residual[v0].sadd(flux0,-(0.5 * grid.cell[i].radius));
		 residual[v1].sadd(flux1,-(0.5 * grid.cell[i].radius));
		 residual[v2].sadd(flux2,-(0.5 * grid.cell[i].radius));
	  }
   }
   
   
   // Diffusive flux from boundary triangles
   for(unsigned int i=0; i<grid.bface.size(); ++i)
   {
      
      //int face_type = grid.bface[i].type;
      unsigned int cl = grid.bface[i].lcell;
      if(grid.cell[cl].nbr_bface.size() > 0 && !grid.bface[i].periodic) // Only boundary triangles with no periodic face
      {
		 EntVar ent_avg;
		 ent_avg = 0.0;
		 for(unsigned j=0;j<grid.cell[cl].nbr_bface.size();++j)
		 {
			int ind = grid.cell[cl].nbr_bface[j];
			int f = grid.cell[cl].face[ind];
			int f_type = grid.face[f].type;
			BoundaryCondition& bc = param.boundary_condition[f_type];
			PrimVar state;
			EntVar ent_state;
			for(unsigned int k=0;k<2;++k)
			{
			   unsigned int vb = grid.face[f].vertex[k];
			   state = primitive[vb];
			   bc.apply(grid.vertex[vb].coord, int_step_time, grid.face[f], state);
			   ent_state = param.material.prim2ent(state);
			   ent_avg+=ent_state;
			}
	     }
	     ent_avg*=(0.5*1.0/grid.cell[cl].nbr_bface.size());
	     
	     unsigned int v0 = grid.cell[cl].vertex[0];
		 unsigned int v1 = grid.cell[cl].vertex[1];
		 unsigned int v2 = grid.cell[cl].vertex[2];
		 unsigned int v0b = grid.bface[i].vertex[0];
		 unsigned int v1b = grid.bface[i].vertex[1];
		 
		 Flux flux0, flux1, flux2,fluxb;
		 Vector& normal0 = grid.cell[cl].normal[0]; 
		 Vector& normal1 = grid.cell[cl].normal[1];
		 Vector& normal2 = grid.cell[cl].normal[2];
		 Vector& normalb = grid.bface[i].normal;

		 param.material.viscous_flux (ent_avg, 
									  dE_cell[cl], 
									  normal0, flux0,
									  normal1, flux1,
									  normal2, flux2);                                  
		 // negative sign as the inward normals are used                          
		 residual[v0].sadd(flux0,-(0.5 * grid.cell[cl].radius));
		 residual[v1].sadd(flux1,-(0.5 * grid.cell[cl].radius));
		 residual[v2].sadd(flux2,-(0.5 * grid.cell[cl].radius));
		 
		 // if((grid.vertex[v0].coord.x < 0.2 && grid.vertex[v0].coord.y > 1.99) 
// 		     || (grid.vertex[v0].coord.x > 3.7 && grid.vertex[v0].coord.y > 1.99)) 
//             cout<<i<<" "<<grid.vertex[v0].coord.x<<" "<<grid.vertex[v0].coord.y<<" : "
//                       <<residual[v0].mass_flux<<" "
//                       <<residual[v0].momentum_flux.x<<" "
//                       <<residual[v0].momentum_flux.y<<" "
//                       <<residual[v0].energy_flux<<endl;
//                       
//          if((grid.vertex[v1].coord.x < 0.2 && grid.vertex[v1].coord.y > 1.99) 
// 		     || (grid.vertex[v1].coord.x > 3.7 && grid.vertex[v1].coord.y > 1.99)) 
//             cout<<i<<" "<<grid.vertex[v1].coord.x<<" "<<grid.vertex[v1].coord.y<<" : "
//                       <<residual[v1].mass_flux<<" "
//                       <<residual[v1].momentum_flux.x<<" "
//                       <<residual[v1].momentum_flux.y<<" "
//                       <<residual[v1].energy_flux<<endl;             
// 		 if((grid.vertex[v2].coord.x < 0.2 && grid.vertex[v2].coord.y > 1.99) 
// 		     || (grid.vertex[v2].coord.x > 3.7 && grid.vertex[v2].coord.y > 1.99)) 
//             cout<<i<<" "<<grid.vertex[v2].coord.x<<" "<<grid.vertex[v2].coord.y<<" : "
//                       <<residual[v2].mass_flux<<" "
//                       <<residual[v2].momentum_flux.x<<" "
//                       <<residual[v2].momentum_flux.y<<" "
//                       <<residual[v2].energy_flux<<endl;
		 
		 // Viscous flux contribution from boundary triangle across boundary face
		 int f_type = grid.bface[i].type;
		 BoundaryCondition& bc = param.boundary_condition[f_type];
		 param.material.viscous_bflux (ent_avg, 
									   bc.adiabatic,
									   dE_cell[cl], 
									   normalb, fluxb);   
		 //t_visc_bflux.add_time();                                    
		 residual[v0b].sadd(fluxb,(0.5 * grid.vertex[v0b].radius));
		 residual[v1b].sadd(fluxb,(0.5 * grid.vertex[v1b].radius));
		 
		 // if((grid.vertex[v0b].coord.x < 0.2 && grid.vertex[v0b].coord.y > 1.99) 
// 		     || (grid.vertex[v0b].coord.x > 3.7 && grid.vertex[v0b].coord.y > 1.99)) 
//             cout<<i<<" "<<grid.vertex[v0b].coord.x<<" "<<grid.vertex[v0b].coord.y<<" : "
//                       <<residual[v0b].mass_flux<<" "
//                       <<residual[v0b].momentum_flux.x<<" "
//                       <<residual[v0b].momentum_flux.y<<" "
//                       <<residual[v0b].energy_flux<<endl;
//                       
//          if((grid.vertex[v1b].coord.x < 0.2 && grid.vertex[v1b].coord.y > 1.99) 
// 		     || (grid.vertex[v1b].coord.x > 3.7 && grid.vertex[v1b].coord.y > 1.99)) 
//             cout<<i<<" "<<grid.vertex[v1b].coord.x<<" "<<grid.vertex[v1b].coord.y<<" : "
//                       <<residual[v1b].mass_flux<<" "
//                       <<residual[v1b].momentum_flux.x<<" "
//                       <<residual[v1b].momentum_flux.y<<" "
//                       <<residual[v1b].energy_flux<<endl;    
		 	
// 			   
// 
// 	 
// 		 }
// 		 
// 		 BoundaryCondition& bc = param.boundary_condition[face_type];
// 
// 		 // cell adjacent to boundary face
// 		 // compute average state on this cell
// 		 //unsigned int cl = grid.bface[i].lcell; 
// 		 EntVar ent_avg;
// 		 // for(unsigned j=0;j<3;++j)
//    // 	  {
//    // 		 int f = grid.cell[i].face[j];
//    // 		 int f_type = grid.face[f].type;
//    // 		 if(f_type != -1 && f_.type != -2)
//    // 		 {
//    // 		    BoundaryCondition& bc = param.boundary_condition[f_type];
//    // 		    
//    // 		 }
//    // 	 
//    // 	  }
// 		 if(grid.bface[i].periodic)
// 		 {
// 			unsigned int v0 = grid.cell[i].vertex[0];
// 			unsigned int v1 = grid.cell[i].vertex[1];
// 			unsigned int v2 = grid.cell[i].vertex[2];
// 			double fact = 1.0/3.0;
// 			ent_avg.equ(entropy_var[v0],entropy_var[v1],entropy_var[v2],fact,fact,fact);
// 			Flux flux0, flux1, flux2, fluxb;
// 			Vector& normal0 = grid.cell[i].normal[0]; 
// 			Vector& normal1 = grid.cell[i].normal[1];
// 			Vector& normal2 = grid.cell[i].normal[2];
// 
// 			// Viscous flux contribution from boundary triangle across inner faces
// 			param.material.viscous_flux (ent_avg, 
// 										 dE_cell[i], 
// 										 normal0, flux0,
// 										 normal1, flux1,
// 										 normal2, flux2);                                  
// 						   
// 			// negative sign as the inward normals are used                          
// 			residual[v0].sadd(flux0,-(0.5 * grid.cell[i].radius));
// 			residual[v1].sadd(flux1,-(0.5 * grid.cell[i].radius));
// 			residual[v2].sadd(flux2,-(0.5 * grid.cell[i].radius));
// 			
// 			
// 		 }
// 		 else
// 		 {		 
// 			double factor = 1; // 1 for cell.type 1, 0.5 for cell.type 2
// 			// if(grid.cell[cl].type == 2)
//    // 			factor = 0.5;
// 			 
// 			unsigned int v0b = grid.bface[i].vertex[0];
// 			unsigned int v1b = grid.bface[i].vertex[1];
// 			unsigned int v0, v1, v2; //v2 will correspond to the inner node
// 			v0 = v1 = v2 = 0;
// 			Vector normal0, normal1, normal2;
// 	 
// 		   // t_visc_normal_find.start_time();
// 			for(unsigned int j=0;j<3;++j)
// 			{
// 			   if(grid.cell[cl].vertex[j] == v0b)
// 			   {
// 				  v0 = v0b;
// 				  normal0 = grid.cell[cl].normal[j];
// 			   }
// 			   else if(grid.cell[cl].vertex[j] == v1b)
// 			   {
// 				  v1 = v1b;
// 				  normal1 = grid.cell[cl].normal[j];
// 			   }   
// 			   else
// 			   {
// 				  v2 = grid.cell[cl].vertex[j];
// 				  normal2 = grid.cell[cl].normal[j];
// 			   }
// 			}   
// 			//t_visc_normal_find.add_time();
// 		   
// 			vector<PrimVar> state(2);
// 			vector<EntVar> ent_state(2);
// 			state[0] = primitive[v0];
// 			state[1] = primitive[v1];
// 	 
// 			//t_visc_bc_apply.start_time();
// 			bc.apply(grid.vertex[v0].coord, int_step_time, seed, grid.bface[i], state[0]);
// 			bc.apply(grid.vertex[v1].coord, int_step_time, seed, grid.bface[i], state[1]);
// 			//t_visc_bc_apply.add_time();
// 	 
// 			for(unsigned int j=0; j<2; ++j)
// 			   ent_state[j] = param.material.prim2ent(state[j]);
// 	 
// 			ent_avg.equ(ent_state[0],ent_state[1],0.5,0.5);
// 			//EntVar ent_avg = (ent_state[0] + ent_state[1] + entropy_var[v2])/3.0;
// 	 
// 			Flux flux0, flux1, flux2, fluxb;
// 	 
// 			// Viscous flux contribution from boundary triangle across inner faces
// 			//t_visc_flux.start_time();  
// 			param.material.viscous_flux (ent_avg, 
// 										 dE_cell[cl], 
// 										 normal0, flux0,
// 										 normal1, flux1,
// 										 normal2, flux2);   
// 			//t_visc_flux.add_time();                               
// 			// negative sign as the inward normals are used                                    
// 			residual[v0].sadd(flux0,-(0.5 * grid.cell[cl].radius)*factor);
// 			residual[v1].sadd(flux1,-(0.5 * grid.cell[cl].radius)*factor);
// 			residual[v2].sadd(flux2,-(0.5 * grid.cell[cl].radius)*factor);
//    
// 			// Viscous flux contribution from boundary triangle across boundary face
// 			//t_visc_bflux.start_time();  
// 			param.material.viscous_bflux (ent_avg, 
// 										  bc.adiabatic,
// 										  dE_cell[cl], 
// 										  normal2, fluxb);   
// 			//t_visc_bflux.add_time();                                    
// 			residual[v0].sadd(fluxb,-(0.5 * grid.vertex[v0].radius));
// 			residual[v1].sadd(fluxb,-(0.5 * grid.vertex[v1].radius));
// 		 }
      }
      
   } 
   
   t_vres_eval.add_time();  
}


//------------------------------------------------------------------------------
// Compute residual for each cell
//------------------------------------------------------------------------------
void FiniteVolume::compute_residual ()
{
   t_res_eval.start_time();
   
   // Set residual vector to zero
   for(unsigned int i=0; i<grid.n_vertex; ++i)
      residual[i].zero ();  

   // Compute cell and vertex gradients
   compute_gradients ();

   // Inviscid fluxes
   compute_inviscid_residual ();
   
   
   // Viscous fluxes
   if(param.material.model == Material::ns)
      compute_viscous_residual ();
   
   if(NPART > 1)  
	  reduce_residual();  
   
   if(PERIODIC)
	  reduce_residual_periodic();
    
   t_res_eval.add_time();
    
}

//------------------------------------------------------------------------------
// Compute time step
//------------------------------------------------------------------------------
void FiniteVolume::compute_dt ()
{
   t_dt_eval.start_time();
   
   if (param.time_step > 0.0)
   {
      dt_global = param.time_step;
      for(unsigned int i=0; i<grid.n_vertex; ++i)
         dt[i] = dt_global;
      return;
   }

   for(unsigned int i=0; i<grid.n_vertex; ++i)
      dt[i] = 0.0;
      
   for(unsigned int i=0; i<grid.n_face; ++i)
   {
      double area = grid.face[i].measure;
      
      double ip_bc = 1.0;
      if(grid.face[i].type == -2)
         ip_bc = 0.5;

      unsigned int cl = grid.face[i].vertex[0];
      double vel_normal_left = primitive[cl].velocity * grid.face[i].normal;
      double c_left = param.material.sound_speed (primitive[cl]);
      dt[cl] += (fabs(vel_normal_left) + c_left * area)*ip_bc;

      unsigned int cr = grid.face[i].vertex[1];
      double vel_normal_right = primitive[cr].velocity * grid.face[i].normal;
      double c_right = param.material.sound_speed (primitive[cr]);
      dt[cr] += (fabs(vel_normal_right) + c_right * area)*ip_bc;
   }
   
   
   //mpi_barrier();

   // Boundary faces
   for(unsigned int i=0; i<grid.bface.size(); ++i)
   {
      if(!grid.bface[i].periodic)
      {
		 double area = grid.bface[i].measure;

		 unsigned int cl = grid.bface[i].vertex[0];
		 double vel_normal_left = primitive[cl].velocity * grid.bface[i].normal;
		 double c_left = param.material.sound_speed (primitive[cl]);
		 dt[cl] += 0.5 * (fabs(vel_normal_left) + c_left * area);

		 unsigned int cr = grid.bface[i].vertex[1];
		 double vel_normal_right = primitive[cr].velocity * grid.bface[i].normal;
		 double c_right = param.material.sound_speed (primitive[cr]);
		 dt[cr] += 0.5 * (fabs(vel_normal_right) + c_right * area);
      }
   }
   
   // Add viscous time step restriction
   if(param.material.model == Material::ns)
   {
      // Interior faces
      for(unsigned int i=0; i<grid.n_face; ++i)
      {
         double area = grid.face[i].measure;
         
         double ip_bc = 1.0;
         if(grid.face[i].type == -2)
            ip_bc = 0.5;

         unsigned int cl = grid.face[i].vertex[0];
         unsigned int cr = grid.face[i].vertex[1];
         Vector ds;
         ds.equ(grid.vertex[cr].coord,grid.vertex[cl].coord,1.0,-1.0);
         double dlength = ds.norm();
         double mu_l = param.material.viscosity (primitive[cl].temperature);
         double mu_r = param.material.viscosity (primitive[cr].temperature);
         double rho_l = param.material.Density(primitive[cl]);
         double rho_r = param.material.Density(primitive[cr]);

         dt[cl] += (mu_l * area / dlength / rho_l)*ip_bc;
         dt[cr] += (mu_r * area / dlength / rho_r)*ip_bc;
      }
   }

   // Compute global time step
   dt_global = 1.0e20;   

   // REDUCE POINTWISE dt
   if(NPART > 1)
	  reduce_dt();
   
   if(PERIODIC)
	  reduce_dt_periodic();

   for(unsigned int i=0; i<grid.n_vertex; ++i)
   {
	  dt[i] = param.cfl * grid.dcarea[i] / dt[i];
	  dt_global = min( dt_global, dt[i] );
   }
   
   if(NPART > 1)
   {
	  // REDUCE GLOBAL dt
	  double dt_red;
	  MPI_Allreduce(&dt_global,&dt_red,1,MPI_DOUBLE,MPI_MIN,grid.run_comm);
	  dt_global = dt_red;
   }
   

   // For unsteady simulation, use global time step
   if(param.time_mode == "unsteady")
   {
      // Adjust time step so that final time is exactly reached
      if(elapsed_time + dt_global > mc_time[mc_t_ind])
      {
         dt_global = mc_time[mc_t_ind] - elapsed_time;
         mc_time_reached = true;
      }   
      for(unsigned int i=0; i<grid.n_vertex; ++i)
         dt[i] = dt_global;
   }
   
   t_dt_eval.add_time();
}

//------------------------------------------------------------------------------
// Store old conserved variables for multi-stage RK
//------------------------------------------------------------------------------
void FiniteVolume::store_conserved_old ()
{
   for(unsigned int i=0; i<grid.n_vertex; ++i)
       conserved_old[i] = param.material.prim2con (primitive[i]);
}

//------------------------------------------------------------------------------
// Update solution to new time level
//------------------------------------------------------------------------------
void FiniteVolume::update_solution (const unsigned int r)
{
   t_sol_update.start_time();
   double factor;
   ConVar conserved;

   if(param.time_scheme == "ssprk3")
   {
      //unsigned int grid.n_vertex_active = grid.active_vertex.size();
      for(unsigned int i = 0; i<grid.n_vertex; ++i)
      {
		 factor      = dt[i] / (grid.vertex[i].radius * grid.dcarea[i]);
		 conserved   = param.material.prim2con (primitive[i]);
		 conserved *= b_rk[r];
		 conserved.sadd(conserved_old[i],residual[i],a_rk[r],-factor*b_rk[r]);
		 primitive[i]= param.material.con2prim (conserved);
		 entropy_var[i] = param.material.prim2ent(primitive[i]);		 
      }
   }
//    else if(param.time_scheme == "rk1" ||
//            param.time_scheme == "rk3")
//    {
//       double f = 1.0/(param.n_rks - r);
//       for(unsigned int i=0; i<grid.n_vertex; ++i)
//       {
//           
//          factor      = f * dt[i] / (grid.vertex[i].radius * grid.dcarea[i]);
// 		 conserved.equ(conserved_old[i],residual[i],1.0, -factor);
// 		 primitive[i]= param.material.con2prim (conserved);
// 		 entropy_var[i] = param.material.prim2ent(primitive[i]);	 
//       }
//    }
   else if(param.time_scheme == "rk4")
   {
      for(unsigned int i=0; i<grid.n_vertex; ++i)
      {
         factor      = dt[i] / (grid.vertex[i].radius * grid.dcarea[i]);
         residual2[i].sadd(residual[i],1.0/rk4_rks[r]);
         if(r<3)
            conserved.equ(conserved_old[i],residual[i],1.0, -factor*rk4_rks[r+1]);
         else   
		    conserved.equ(conserved_old[i],residual2[i],1.0, -factor/6.0);
		 primitive[i]= param.material.con2prim (conserved);
		 entropy_var[i] = param.material.prim2ent(primitive[i]);	 
      }
   }
   else if(param.time_scheme == "heuns")
   {
      for(unsigned int i=0; i<grid.n_vertex; ++i)
      {
         factor = dt[i] / (grid.vertex[i].radius * grid.dcarea[i]);
		 conserved = param.material.prim2con (primitive[i]);
		 conserved *= (1.0 - r * 0.5);
		 conserved.equ(conserved_old[i],residual[i],r * 0.5,-factor*(1.0 - r * 0.5)); 
		 primitive[i]= param.material.con2prim (conserved);
		 entropy_var[i] = param.material.prim2ent(primitive[i]);
      }
   }
   
   // Apply strong bc
   if (param.bc_scheme == Parameter::strong)
      for(unsigned int i=0; i<grid.bface.size(); ++i)
      {
         int face_type = grid.bface[i].type;
         BoundaryCondition& bc = param.boundary_condition[face_type];
         unsigned int v0 = grid.bface[i].vertex[0];
         unsigned int v1 = grid.bface[i].vertex[1];
         bc.apply (grid.vertex[v0].coord, int_step_time, grid.bface[i], primitive[v0]);
         bc.apply (grid.vertex[v1].coord, int_step_time, grid.bface[i], primitive[v1]);
      }
      
   t_sol_update.add_time();   
}

//------------------------------------------------------------------------------
// Compute L1, L2 and Linfty norms of solution errors if exact solution is available
//------------------------------------------------------------------------------
void FiniteVolume::compute_error_norm ()
{
   map<int,double>::iterator it;
   double w;
   for(unsigned int i = 0; i<grid.n_vertex; ++i)
   {
      double area = grid.dcarea[i];
	  double nbs = grid.mpi_node_share.find(i)->second;
	  it = grid.periodic_node_weight.find(i);
	  if(it!=grid.periodic_node_weight.end())
	    w = it->second;
	  else
	    w = 1.0/nbs;
	  
	  prim_exact    = param.exact_solution.value (grid.vertex[i].coord,param.final_time);  
	  density_exact = param.material.Density(prim_exact);
	  entropy_exact = param.material.Entropy(prim_exact);
	  mach_exact    = param.material.Mach(prim_exact);
	  
	  prim_abs_err.absdiff(prim_exact,primitive[i]);
	  density_abs_err = fabs(density_exact - param.material.Density(primitive[i]));   
	  mach_abs_err    = fabs(mach_exact - param.material.Mach(primitive[i])); 
	  entropy_abs_err = fabs(entropy_exact - param.material.Entropy(primitive[i])); 
	  
	  prim_err_Linf.max(prim_abs_err);
	  prim_err_L1.sadd(prim_abs_err,area*w);
	  prim_err_L2.sqsadd(prim_abs_err,area*w);
	  
	  density_err_Linf = max(density_err_Linf,density_abs_err);
	  density_err_L1  += density_abs_err*area*w;
	  density_err_L2  += pow(density_abs_err,2.0)*area*w;
	  
	  mach_err_Linf = max(mach_err_Linf,mach_abs_err);
	  mach_err_L1  += mach_abs_err*area*w;
	  mach_err_L2  += pow(mach_abs_err,2.0)*area*w;
	  
	  entropy_err_Linf = max(entropy_err_Linf,entropy_abs_err);
	  entropy_err_L1  += entropy_abs_err*area*w;
	  entropy_err_L2  += pow(entropy_abs_err,2.0)*area*w;  
   }

   if(NPART > 1)
      reduce_err_norm();
   
   if(check_group_base())
   {
      err_file << "N = " << grid.n_vertex_total << "\n"
               << scientific
               << setprecision(8)
               << "TEMPERATURE (L1, L2, Linf): "
               << prim_err_L1.temperature << "   "
               << sqrt(prim_err_L2.temperature) << "   "
               << prim_err_Linf.temperature << "\n"
               << "PRESSURE (L1, L2, Linf): "
               << prim_err_L1.pressure << "   "
               << sqrt(prim_err_L2.pressure) << "   "
               << prim_err_Linf.pressure << "\n"
               << "VEL_X (L1, L2, Linf): "
               << prim_err_L1.velocity.x << "   "
               << sqrt(prim_err_L2.velocity.x) << "   "
               << prim_err_Linf.velocity.x << "\n"
               << "VEL_Y (L1, L2, Linf): "
               << prim_err_L1.velocity.y << "   "
               << sqrt(prim_err_L2.velocity.y) << "   "
               << prim_err_Linf.velocity.y << "\n"
               << "VEL_Z (L1, L2, Linf): "
               << prim_err_L1.velocity.z << "   "
               << sqrt(prim_err_L2.velocity.z) << "   "
               << prim_err_Linf.velocity.z << "\n"
               << "DENSITY (L1, L2, Linf): "
               << density_err_L1 << "   "
               << sqrt(density_err_L2) << "   "
               << density_err_Linf << "\n"
               << "MACH (L1, L2, Linf): "
               << mach_err_L1 << "   "
               << sqrt(mach_err_L2) << "   "
               << mach_err_Linf << "\n"
               << "ENTROPY (L1, L2, Linf): "
               << entropy_err_L1 << "   "
               << sqrt(entropy_err_L2) << "   "
               << entropy_err_Linf << "\n";
       err_file.close();        
   }
}


//------------------------------------------------------------------------------
// Compute L2 norm of mass, momentum and energy residuals
//------------------------------------------------------------------------------
void FiniteVolume::compute_residual_norm (const unsigned int iter)
{
   t_res_norm.start_time();
   
   // If using strong bc, set momentum residual to zero on noslip
   // surfaces.
   if (param.bc_scheme == Parameter::strong)
      for(unsigned int i=0; i<grid.bface.size(); ++i)
      {
         int face_type = grid.bface[i].type;
         BoundaryCondition& bc = param.boundary_condition[face_type];
         if(bc.type == BC::noslip)
         {
            unsigned int v0 = grid.bface[i].vertex[0];
            unsigned int v1 = grid.bface[i].vertex[1];
            residual[v0].momentum_flux = 0;
            residual[v1].momentum_flux = 0;
         }
      }

   residual_norm.mass_flux     = 0.0;
   residual_norm.momentum_flux = 0.0;
   residual_norm.energy_flux   = 0.0;

   // Sum of squares for each component
   map<int,double>::iterator it;
   double w;
   for(unsigned int i = 0; i<grid.n_vertex; ++i)
   {
      double area = grid.dcarea[i];
	  double nbs = grid.mpi_node_share.find(i)->second;
	  it = grid.periodic_node_weight.find(i);
	  if(it!=grid.periodic_node_weight.end())
	    w = it->second;
	  else
	    w = 1.0/nbs;
	  residual_norm.mass_flux       += pow(residual[i].mass_flux       / area, 2)*w;
	  residual_norm.momentum_flux.x += pow(residual[i].momentum_flux.x / area, 2)*w;
	  residual_norm.momentum_flux.y += pow(residual[i].momentum_flux.y / area, 2)*w;
	  residual_norm.momentum_flux.z += pow(residual[i].momentum_flux.z / area, 2)*w;
	  residual_norm.energy_flux     += pow(residual[i].energy_flux     / area, 2)*w;	  
   }

   if(NPART>1)
	  reduce_res_norm();
   

   // Take square root and normalize by n_vertex
   residual_norm.mass_flux       = sqrt (residual_norm.mass_flux)       / grid.n_vertex_total;
   residual_norm.momentum_flux.x = sqrt (residual_norm.momentum_flux.x) / grid.n_vertex_total;
   residual_norm.momentum_flux.y = sqrt (residual_norm.momentum_flux.y) / grid.n_vertex_total;
   residual_norm.momentum_flux.z = sqrt (residual_norm.momentum_flux.z) / grid.n_vertex_total;
   residual_norm.energy_flux     = sqrt (residual_norm.energy_flux)     / grid.n_vertex_total;

   // Total residual of all components
   residual_norm_total = pow(residual_norm.mass_flux, 2) +
                         residual_norm.momentum_flux.square () +
                         pow(residual_norm.energy_flux, 2);
   residual_norm_total = sqrt (residual_norm_total);

   // Copy residual in first iteration for normalization
   if(iter == 0)
   {
      residual_norm_total0 = residual_norm_total;
      if(check_group_base())
      {
		  if(residual_norm_total0 == 0.0)
		  {
			 cout << "  WARNING: Initial residual is zero !!!\n";
			 cout << "  WARNING: Setting it to 1.0\n";
			 residual_norm_total0 = 1.0;
		  }
      }		  
   }

   residual_norm_total /= residual_norm_total0;
   
   t_res_norm.add_time();
}

//------------------------------------------------------------------------------
// Compute L1 norm of mean and variance (locally)
//------------------------------------------------------------------------------
void FiniteVolume::compute_stat_norm ()
{
   t_l1_norm.start_time();
   
   for(unsigned int t = 0; t<mc_time.size(); ++t)
   {
	  for(unsigned int i = 0; i<grid.n_vertex; ++i)
	  {
		 double area = grid.dcarea[i];
		 double nbs = grid.mpi_node_share.find(i)->second;
		 prim_mean_l1[t].temperature  += area*abs(primitive_fm[t][i].temperature)/nbs;
		 prim_mean_l1[t].velocity.x   += area*abs(primitive_fm[t][i].velocity.x)/nbs;
		 prim_mean_l1[t].velocity.y   += area*abs(primitive_fm[t][i].velocity.y)/nbs;
		 prim_mean_l1[t].velocity.z   += area*abs(primitive_fm[t][i].velocity.z)/nbs;
		 prim_mean_l1[t].pressure     += area*abs(primitive_fm[t][i].pressure)/nbs;
		 
		 prim_var_l1[t].temperature   += area*abs(primitive_sm[t][i].temperature)/nbs;
		 prim_var_l1[t].velocity.x    += area*abs(primitive_sm[t][i].velocity.x)/nbs;
		 prim_var_l1[t].velocity.y    += area*abs(primitive_sm[t][i].velocity.y)/nbs;
		 prim_var_l1[t].velocity.z    += area*abs(primitive_sm[t][i].velocity.z)/nbs;
		 prim_var_l1[t].pressure      += area*abs(primitive_sm[t][i].pressure)/nbs;
		 
		 if(has_density)
		 {
		    den_mean_l1[t]            += area*abs(density_fm[t][i])/nbs;
		    den_var_l1[t]             += area*abs(density_sm[t][i])/nbs;
		 } 
		 
		 if(has_mach)
		 {
		    mach_mean_l1[t]            += area*abs(mach_fm[t][i])/nbs;
		    mach_var_l1[t]             += area*abs(mach_sm[t][i])/nbs;
		 }
		 
		 if(has_entropy)
		 {
		    ent_mean_l1[t]            += area*abs(entropy_fm[t][i])/nbs;
		    ent_var_l1[t]             += area*abs(entropy_sm[t][i])/nbs;
		 }
		 
		 if(has_vorticity)
		 {
		    vor_mean_l1[t]            += area*abs(vorticity_fm[t][i])/nbs;
		    vor_var_l1[t]             += area*abs(vorticity_sm[t][i])/nbs;
		 }  
	  }
   }	
   
   t_l1_norm.add_time();  

}

//------------------------------------------------------------------------------
// Log messages to screen and file
//------------------------------------------------------------------------------
void FiniteVolume::log_messages (const unsigned int iter)
{
   t_log.start_time();
   if(check_group_base())
   {
	   if(param.time_mode == "steady")
	   {
		  // File output
		  res_file  << setw(8) << iter << "  " 
					<< scientific
					<< setprecision (4) 
					<< dt_global << "  " 
					<< residual_norm_total << "  "
					<< residual_norm.mass_flux << "  "
					<< residual_norm.momentum_flux.x << "  "
					<< residual_norm.momentum_flux.y << "  "
					<< residual_norm.momentum_flux.z << "  "
					<< residual_norm.energy_flux
					<< endl;
	   }
	   else
	   {
		  // File output
		  res_file  << setw(8) << iter << "  " 
					<< scientific
					<< setprecision (12) 
					<< dt_global << "  " 
					<< elapsed_time << "  "
					<< residual_norm_total
					<< endl;
	   }
   }
   if(bounds)
      compute_bounds (iter);
   t_log.add_time();   
}

//------------------------------------------------------------------------------
// Save solution to file for visualization
//------------------------------------------------------------------------------
void FiniteVolume::output (const unsigned int iter, bool write_variables)
{
   t_sol_out.start_time();

   Writer writer (grid, param.material, param.write_format, param.write_surfaces, 
                  elapsed_time, SAMPLE_DIR);
   writer.attach_data (primitive);
   writer.attach_gradient (dE);
   if(param.write_variables.size() > 0 && write_variables == true)
      writer.attach_variables (param.write_variables);
   
   writer.master_write_check(write_to_master);     
   writer.output (get_proc_loc_id(),counter, elapsed_time);
   output_surface_sf(SAMPLE_DIR);
   output_surface_hf(SAMPLE_DIR);
   if(param.save_mesh_Pe)
      output_mesh_Pe(SAMPLE_DIR);
   if(param.time_mode == "unsteady") 
      ++counter;
   else
      write_to_master = false;   
   
   t_sol_out.add_time();
}

//------------------------------------------------------------------------------
// Save solution to file for restart
//------------------------------------------------------------------------------
void FiniteVolume::output_restart (int iter)
{
   Writer writer (grid,SAMPLE_DIR);
   writer.attach_data (primitive);
   writer.output_restart (iter,counter,restart_head_tag,elapsed_time,residual_norm_total0);
   mpi_barrier(grid.run_comm);
   if(check_group_base())
   {
      string fil = SAMPLE_DIR + "/LAST_RESTART_HEAD.txt";
      ofstream fo;
      fo.open (fil.c_str());
      fo << restart_head_tag <<endl;
      fo.close();  
   }   
   restart_head_tag = (restart_head_tag+1) % 2;   
}

//------------------------------------------------------------------------------
// Save ensemble mean to file for visualization
//------------------------------------------------------------------------------
void FiniteVolume::output_mean()
{
   t_mean_out.start_time();
   
   Writer writer (grid, param.write_format, MEAN_DIR);
   writer.attach_data (primitive_fm[mc_t_ind]);
   if(has_density)
      writer.attach_data (density_fm[mc_t_ind],"density");
   if(has_mach)
      writer.attach_data (mach_fm[mc_t_ind],"mach");
   if(has_vorticity)
      writer.attach_data (vorticity_fm[mc_t_ind],"vorticity");
   if(has_entropy)
      writer.attach_data (entropy_fm[mc_t_ind],"entropy");         
   writer.master_write_check(write_to_master_mean);     
   writer.output (get_proc_id(),mc_t_ind, mc_time[mc_t_ind]);
 
   t_mean_out.add_time();
}

//------------------------------------------------------------------------------
// Save ensemble second moments to file for visualization
//------------------------------------------------------------------------------
void FiniteVolume::output_sec_mom()
{
   t_mean_out.start_time();
   
   Writer writer (grid, param.write_format, SEC_MOM_DIR);
   writer.attach_data (primitive_sm[mc_t_ind]);
   if(has_density)
      writer.attach_data (density_sm[mc_t_ind],"density");
   if(has_mach)
      writer.attach_data (mach_sm[mc_t_ind],"mach");
   if(has_vorticity)
      writer.attach_data (vorticity_sm[mc_t_ind],"vorticity");
   if(has_entropy)
      writer.attach_data (entropy_sm[mc_t_ind],"entropy");         
   writer.master_write_check(write_to_master_sec_mom);     
   writer.output (get_proc_id(),mc_t_ind, mc_time[mc_t_ind]);
 
   t_mean_out.add_time();
}

//------------------------------------------------------------------------------
// Save ensemble var to file for visualization
//------------------------------------------------------------------------------
void FiniteVolume::output_var()
{
   t_var_out.start_time();
   
   Writer writer (grid, param.write_format, VAR_DIR);
   writer.attach_data (primitive_variance[mc_t_ind]);
   if(has_density)
      writer.attach_data (density_variance[mc_t_ind],"density");
   if(has_mach)
      writer.attach_data (mach_variance[mc_t_ind],"mach");
   if(has_vorticity)
      writer.attach_data (vorticity_variance[mc_t_ind],"vorticity");
   if(has_entropy)
      writer.attach_data (entropy_variance[mc_t_ind],"entropy");
   writer.master_write_check(write_to_master_var);     
   writer.output (get_proc_id(),mc_t_ind, mc_time[mc_t_ind]);

   t_var_out.add_time();
}

//------------------------------------------------------------------------------
// Find minimum and maximum values in the solution
//------------------------------------------------------------------------------
void FiniteVolume::compute_bounds (const unsigned int iter)
{
   PrimVar prim_min;
   PrimVar prim_max;

   prim_min.temperature=  1.0e20;
   prim_min.velocity.x =  1.0e20;
   prim_min.velocity.y =  1.0e20;
   prim_min.velocity.z =  1.0e20;
   prim_min.pressure   =  1.0e20;

   prim_max.temperature= -1.0e20;
   prim_max.velocity.x = -1.0e20;
   prim_max.velocity.y = -1.0e20;
   prim_max.velocity.z = -1.0e20;
   prim_max.pressure   = -1.0e20;

   for(unsigned int i=0; i<grid.n_vertex; ++i)
   {
      prim_min.min(primitive[i]);
      prim_max.max(primitive[i]);
   }

   if(NPART > 1)
   {
      //Reducing bounds
	  double min_red,max_red;
	  MPI_Allreduce(&prim_min.temperature,&min_red,1,MPI_DOUBLE,MPI_MIN,grid.run_comm);
	  MPI_Allreduce(&prim_max.temperature,&max_red,1,MPI_DOUBLE,MPI_MAX,grid.run_comm);
	  prim_min.temperature = min_red;
	  prim_max.temperature = max_red;
   
	  MPI_Allreduce(&prim_min.velocity.x,&min_red,1,MPI_DOUBLE,MPI_MIN,grid.run_comm);
	  MPI_Allreduce(&prim_max.velocity.x,&max_red,1,MPI_DOUBLE,MPI_MAX,grid.run_comm);
	  prim_min.velocity.x = min_red;
	  prim_max.velocity.x = max_red;
   
	  MPI_Allreduce(&prim_min.velocity.y,&min_red,1,MPI_DOUBLE,MPI_MIN,grid.run_comm);
	  MPI_Allreduce(&prim_max.velocity.y,&max_red,1,MPI_DOUBLE,MPI_MAX,grid.run_comm);
	  prim_min.velocity.y = min_red;
	  prim_max.velocity.y = max_red;
   
	  MPI_Allreduce(&prim_min.velocity.z,&min_red,1,MPI_DOUBLE,MPI_MIN,grid.run_comm);
	  MPI_Allreduce(&prim_max.velocity.z,&max_red,1,MPI_DOUBLE,MPI_MAX,grid.run_comm);
	  prim_min.velocity.z = min_red;
	  prim_max.velocity.z = max_red;
   
	  MPI_Allreduce(&prim_min.pressure,&min_red,1,MPI_DOUBLE,MPI_MIN,grid.run_comm);
	  MPI_Allreduce(&prim_max.pressure,&max_red,1,MPI_DOUBLE,MPI_MAX,grid.run_comm);
	  prim_min.pressure = min_red;
	  prim_max.pressure = max_red;
   }
   
   bound_file << left << setw(17) << setfill(' ') << setprecision(10) <<elapsed_time
			  << left << setw(17) << setfill(' ') << setprecision(10) <<prim_min.temperature
			  << left << setw(17) << setfill(' ') << setprecision(10) <<prim_max.temperature
			  << left << setw(17) << setfill(' ') << setprecision(10) <<prim_min.pressure 
			  << left << setw(17) << setfill(' ') << setprecision(10) <<prim_max.pressure 
			  << left << setw(17) << setfill(' ') << setprecision(10) <<prim_min.velocity.x
			  << left << setw(17) << setfill(' ') << setprecision(10) <<prim_max.velocity.x
			  << left << setw(17) << setfill(' ') << setprecision(10) <<prim_min.velocity.y
			  << left << setw(17) << setfill(' ') << setprecision(10) <<prim_max.velocity.y
			  << left << setw(17) << setfill(' ') << setprecision(10) <<prim_min.velocity.z
			  << left << setw(17) << setfill(' ') << setprecision(10) <<prim_max.velocity.z<<endl;

   if (prim_min.temperature < 0.0 ||
       prim_min.pressure    < 0.0)
   {
         output (iter, false);
         MPI_LOC_ERR("Error:: Observed negative temperature/pressure in solution");
   }
}

//------------------------------------------------------------------------------
// Compute some quantities like global KE
//------------------------------------------------------------------------------
void FiniteVolume::compute_global (unsigned int iter)
{
   if (!param.has_global) return;

   t_global_eval.start_time();
   if(get_proc_loc_id() == 0)
   {
      global_file << iter << "  " << elapsed_time;
   }   
   
   double global_KE = 0;
   double global_entropy = 0;
   map<int,double>::iterator it;
   double w;
   for(unsigned int i=0; i<grid.n_vertex; ++i)
   {
	  double nbs = grid.mpi_node_share.find(i)->second;
	  it = grid.periodic_node_weight.find(i);
	  if(it!=grid.periodic_node_weight.end())
	    w = it->second;
	  else
	    w = 1.0/nbs;
	  global_KE += (0.5 * 
				   param.material.Density (primitive[i]) * 
				   primitive[i].velocity.square() *
				   grid.dcarea[i])*w; 
	  global_entropy += (- param.material.Density (primitive[i])*param.material.Entropy (primitive[i])/
	                       (param.material.gamma - 1.0))*grid.dcarea[i]*w;			               
   }
   
   if(NPART > 1)
   {
	  double global_KE_MPI;
	  MPI_Allreduce(&global_KE,&global_KE_MPI,1,MPI_DOUBLE,MPI_SUM,grid.run_comm);
	  global_KE = global_KE_MPI;
   
	  double global_entropy_MPI;
	  MPI_Allreduce(&global_entropy,&global_entropy_MPI,1,MPI_DOUBLE,MPI_SUM,grid.run_comm);
	  global_entropy = global_entropy_MPI;
   }
   
   if(get_proc_loc_id() == 0)
   {
	  global_file << "  " << global_KE
	              << "  " << global_entropy
	              << endl;
   }   
   
   t_global_eval.add_time();

}

//------------------------------------------------------------------------------
// Perform time marching iterations
//------------------------------------------------------------------------------
void FiniteVolume::solve (const int sample_id)
{
   t_sample_solve.start_time();
   
   unsigned int iter = last_iter;
   residual_norm_total = 1.0e20;
   unsigned int last_output_iter = 0;
   
   compute_gradients ();
   
   if(param.write_soln && !restart)
   {
	  output (0);
	  // if(check_group_base())
// 		 iter_file << 0 << endl;
   }   
   if(param.online_stat)
   {
	  increment_moments(sample_id);
	  //save_sample_pdf();
   }  
   
   if(param.time_mode == "steady")
      counter++;  // Solns from here on saved with index 1 for steady flow
   
   for(unsigned int i =0; i<mc_time.size()-1; ++i)
   {
	  if(elapsed_time >= mc_time[i] && elapsed_time < mc_time[i+1])
	     mc_t_ind = i+1;
   }

   found_stop_file = check_for_stop_file ();
   
   unsigned int pseudo_iter = 0;
   if(param.time_mode == "steady")
      pseudo_iter = iter;  
   
   while (residual_norm_total > param.min_residue &&
          pseudo_iter < param.max_iter && 
          elapsed_time < param.final_time && !found_stop_file)
   {
      store_conserved_old ();
      compute_dt ();
      
      if(param.liou_fix)
      {
         compute_ssw(); 
      }   
      
      if(param.time_scheme == "rk4")
      {
         for(unsigned int i=0; i<grid.n_vertex; ++i)
         {
            residual2[i].zero();
         }
      }
      
      for(unsigned int r=0; r<param.n_rks; ++r)
      {
         if(param.time_scheme == "heuns")
            int_step_time = elapsed_time + r*dt_global;
            
         compute_residual();

         if(r == param.n_rks-1)
            compute_residual_norm (iter);
         update_solution (r);
      }

      ++iter;   
	  if(param.time_mode == "steady")
		 pseudo_iter = iter;
		 
      elapsed_time += dt_global;
      if(param.write_log)
         log_messages (iter);

      compute_forces (iter);
      compute_global (iter);
      
      if(iter % param.write_frequency == 0 && save_by_frequency && param.write_soln)
      {
         output (iter);
         last_output_iter = iter;
         found_stop_file = check_for_stop_file ();
      }
      else if(mc_time_reached && param.time_mode == "unsteady")
      {
         if(param.write_soln)
		 {
			output (iter);
			last_output_iter = iter;
            found_stop_file = check_for_stop_file ();
		 }
		 if(param.online_stat)
		 {
			increment_moments(sample_id);
			//save_sample_pdf();
		 }
		 // if(param.write_soln && check_group_base())
// 			iter_file << counter << endl;
		 mc_t_ind++; 
		 mc_time_reached = false;	 
      }
// 	  if(iter % param.write_frequency == 0 || mc_time_reached) 
// 	  {
// 		 if(mc_time_reached )
// 		 {
// 			if(param.online_stat)
// 			{
// 			   increment_moments(sample_id);
// 			   save_sample_pdf();
// 			}
// 			if(param.write_soln && check_group_base())
// 			   iter_file << counter << endl;
// 			mc_t_ind++; 
// 			mc_time_reached = false;
// 		 }
// 		 if(param.write_soln)
// 		 {
// 			output (iter);
// 		 }   
// 		 last_output_iter = iter;
// 		 found_stop_file = check_for_stop_file ();
// 	  }
      if(param.write_restart && (iter % param.restart_write_frequency == 0 || elapsed_time==param.final_time)) 
         output_restart (iter);

   }
   
   if(param.time_mode == "steady")  // Saving final output for steady simulations
   {
      if(param.write_soln)
         output (iter);
      if(param.online_stat)
	  {
		 increment_moments(sample_id);
		 //save_sample_pdf();
	  }
   }
   
   if(param.find_error)
      compute_error_norm();
   // Save final solution
   // if(iter != last_output_iter && param.write_soln && elapsed_time != param.final_time)
//       output (iter);   
   
   t_sample_solve.add_time();
}


//------------------------------------------------------------------------------
// This is where the real work starts
//------------------------------------------------------------------------------
void FiniteVolume::run ()
{
   t_total.start_time();
   
   // Read grid from file
   t_grid.start_time();
   grid.read (param);
   t_grid.add_time();
   
   mpi_barrier(MPI_COMM_WORLD);

   //If -p flag given on command line, then we stop
   if(preprocess)
	  return;
   
   if(param.n_samples > 1)
   {
      MPI_DISP("\n  Initiated "<<N_MC_GROUPS<< " Monte-Carlo processor group with "<<NPART<<" processors each.",verbose);	  
   }   
   else
   {
      MPI_DISP("\n  Initiated solver with "<<NPART<<" processors.",verbose);	 
   }       
   
   create_force_face_list ();
   
   // Set initial condition
   initialize ();  
   
   //cout<<get_proc_id()<<endl;
   if(get_proc_id()==0)
   {
      cout<<"  Starting simulations ...\n";
      cout<<"    Samples completed:\n";
   
   }   
   
   mpi_barrier(MPI_COMM_WORLD);
   
   for(int ind = param.sample_start_ind; ind <= param.sample_end_ind; ++ind)
   {
	   int sample_no = param.SAMPLE_LIST[ind];
	      
	   counter  = 0; // For file indexing
	   mc_t_ind = 0;  // Index for mc time snap shots
	   
	   initialize_sample_run(ind); 
	   
	   // Solve the sample
	   if((elapsed_time < param.final_time && param.time_mode == "unsteady") 
	       || (last_iter < param.max_iter && param.time_mode == "steady"))
          solve (ind);
	   
	   end_sample_run(); 

	   if(check_group_base())
	      cout<<"     --- "<<sample_no<<" ... "<<endl;
	   
   }	     
   mpi_barrier(MPI_COMM_WORLD);
   
   if(param.online_stat && !found_stop_file)
   {
	   if(get_proc_id() == 0)
		  cout<<"  Reducing moments ...\n";   
	   reduce_moments();   
	   
	   mpi_barrier(MPI_COMM_WORLD);
	   if(get_proc_id() == 0)
		  cout<<"  Saving sample mean and variance and second moments ...\n";
	   if(get_proc_id() < NPART)
	   {
		  for(mc_t_ind = 0; mc_t_ind < mc_time.size(); ++mc_t_ind)
		  {
			 output_mean();
			 output_sec_mom();
			 output_var();
			 
		  }   
	   }      
	   if(get_proc_id() == 0)
		  cout<<"  Saving L1 norms of mean and variance ...\n";   
	   compute_stat_norm(); 
	   mpi_barrier(MPI_COMM_WORLD);
	   reduce_stat_norm();
	   mpi_barrier(MPI_COMM_WORLD);
	   save_stat_norm();
	   
	   // if(get_proc_id() == 0)
// 		  cout<<"  Saving pdf of probe points...\n";  
// 	   mpi_barrier(MPI_COMM_WORLD);	  
// 	   if(get_proc_id() == 0)
// 		  cout<<"  Building and saving final pdf ...\n";  
// 	   //build_pdf();
// 	   save_pdf();
	   mpi_barrier(MPI_COMM_WORLD);
   }
   
   t_total.add_time();
   write_time_output();     
}
