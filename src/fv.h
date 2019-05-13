#ifndef __FV_H__
#define __FV_H__

#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "parameter.h"
#include "material.h"
#include "grid.h"
#include "mpi_utils.h"
#include "timer.h"

extern bool preprocess;
extern bool restart;

//------------------------------------------------------------------------------
//! Finite Volume class 
/*!
  This is the main class for the finite volume solver
*/
//------------------------------------------------------------------------------
class FiniteVolume
{
   public:
      FiniteVolume (char* file, char* loc) 
      :
      has_density (false),
      has_mach (false),
      has_entropy (false),
      has_vorticity (false)
      { 
         param.file = file;
         param.part_dir = loc;
         param.part_dir += "/PARTITION";
         param.read ();
         if(restart && param.online_stat)
             MPI_ERR("error in input.param: restart facility currently unavailable with online statistics.")
          
         if(param.online_stat && !preprocess)
         {
			 // MASTER FILE FOR VIEWING SOLUTION STATISTICS AND OTHER STAT FOLDERS
			 MEAN_DIR    = "MEAN_ONLINE";
			 SEC_MOM_DIR = "SECOND_MOMENT_ONLINE";
			 VAR_DIR     = "VARIANCE_ONLINE";
// 			 PDF_DIR     = "PDF_ONLINE";
// 			 L1_DIR      = "STAT_L1_ONLINE";
			 
			 for(unsigned int i=0; i<param.write_variables.size(); ++i)
                if(param.write_variables[i]=="mach")
                   has_mach = true;
                else if(param.write_variables[i]=="density")
                   has_density = true;
                else if(param.write_variables[i]=="entropy")
                   has_entropy = true;
                else if(param.write_variables[i]=="vorticity")
                   has_vorticity = true;   
                else
                   MPI_ERR("FV: unknown variable " << param.write_variables[i])
		 
			 if(get_proc_id() == 0)
			 {
				std::string commandline0 = "rm -rf " + MEAN_DIR;
				std::string commandline1 = "mkdir " + MEAN_DIR;
				system(commandline0.c_str());
				system(commandline1.c_str());
				
				std::string master_file_path = MEAN_DIR + "/master_file.visit";
				master_file_mean.open (master_file_path.c_str());
			    master_file_mean << "!NBLOCKS "<<NPART<<std::endl;
			    master_file_mean.close ();
			
				commandline0 = "rm -rf " + VAR_DIR;
				commandline1 = "mkdir " + VAR_DIR;
				system(commandline0.c_str());
				system(commandline1.c_str());
				
				master_file_path = VAR_DIR + "/master_file.visit";
				master_file_var.open (master_file_path.c_str());
			    master_file_var << "!NBLOCKS "<<NPART<<std::endl;
			    master_file_var.close ();
			    
			    commandline0 = "rm -rf " + SEC_MOM_DIR;
				commandline1 = "mkdir " + SEC_MOM_DIR;
				system(commandline0.c_str());
				system(commandline1.c_str());
				
				master_file_path = SEC_MOM_DIR + "/master_file.visit";
				master_file_sec_mom.open (master_file_path.c_str());
			    master_file_sec_mom << "!NBLOCKS "<<NPART<<std::endl;
			    master_file_sec_mom.close ();
			    
			    // commandline0 = "rm -rf " + PDF_DIR;
// 				commandline1 = "mkdir " + PDF_DIR;
// 				system(commandline0.c_str());
// 				system(commandline1.c_str());
				
				commandline0 = "rm -rf " + L1_DIR;
				commandline1 = "mkdir " + L1_DIR;
				system(commandline0.c_str());
				system(commandline1.c_str());
			 }
			 write_to_master_mean = true;
			 write_to_master_sec_mom = true;
			 write_to_master_var = true;
		 }
		 mpi_barrier(MPI_COMM_WORLD);
		 time_instance_reached = false;
         
      };
      ~FiniteVolume () 
      {
		 
         dE0.resize(0);
            
      };
      void run ();

   private:
      std::ofstream res_file;           /*!< Output stream for saving residual file */
      std::ofstream force_file;         /*!< Output stream for saving force file */
      std::ofstream global_file;        /*!< Output stream for saving global quantities file */
      std::ofstream bound_file;         /*!< Output stream for solution bounds */
      std::ofstream err_file;           /*!< Output stream for solution error norms */
      std::ofstream master_file;        /*!< Output stream for saving masterfile for 
                                             partitioned solution files */
      std::ofstream master_file_mean;   /*!< Output stream for saving masterfile for 
                                             partitioned sample mean files */
      std::ofstream master_file_sec_mom;/*!< Output stream for saving masterfile for 
                                             partitioned sample second moment files */                                       
      std::ofstream master_file_var;    /*!< Output stream for saving masterfile for 
                                             partitioned sample variance files */
      
      Parameter param;                  
      Grid      grid;

      std::vector<PrimVar> primitive;          /*!< vector of primitive variables at nodes*/
      std::vector<ConVar>  conserved_old;      /*!< vector of conserved (old) variables at nodes.
                                                    This is needed for RK updates */
      std::vector<EntVar>  entropy_var;        /*!< vector of entropy variables at nodes*/
      std::vector<Flux>    residual;           /*!< vector of residuals at nodes*/
      std::vector<Flux>    residual2;          /*!< needed for RK4*/
      std::vector<Vector>  dE0;                /*!< vector to store gradients. For each
                                                    component of the (entropy) variable 
                                                    it stores a Vector corresponding to
                                                    its partial derivatives in each 
                                                    coordinate direction */
      std::vector<std::vector<Vector> >  dE;   /*!< The first
            										vector corresponds to nodes, the second
            										vector to variable components. The final
            										inner Vector corresponds to the partial
            										derivatives in each coordinate direction */
      std::vector<std::vector<Vector> > dE_cell;/*!<The first
            										vector corresponds to cells, the second
            										vector to variable components. The final
            										inner Vector corresponds to the partial
            										derivatives in each coordinate direction */
      
      std::vector<PrimVar> phi;                 /*!<used in Barth-Jesperson limiter*/
      std::vector<double>  ssw;  				/*!<used in Liou's fix*/
      Flux                 residual_norm;       /*!<residual norm for each component*/
      double               residual_norm_total; /*!<total residual norm*/
      double               residual_norm_total0;/*!<initial total residual norm*/
      std::vector<double>  dt;        			/*!<local time step for each node*/
      double               dt_global;           /*!<global time step*/
      std::vector<Force>   inv_force;           /*!<vector of inviscid force value for chosen faces*/
      std::vector<Force>   visc_force;          /*!<vector of viscous force value for chosen faces*/
      double               elapsed_time;        /*!<current elapsed time in simulation*/
      double               int_step_time;       /*!<intermediate time step. This is required
                                                    for Heun's method for time integration*/
      bool                 found_stop_file;                                              
      unsigned long int    last_iter;           /*!<this corresponds to the last iteration 
                                                    in the restart file. If no restart file 
                                                    is used, then it is set to 0*/
      // bool                 save_by_frequency;    /*!<if this is true, then the solution files
//       												are saved based on write_frequency*/
      bool                 write_to_master;     /*!<data is written to solution master file
            										only if this is true. This ensures that
            										the master file data is compatible with
            										both steady and unsteady simulations.*/      										
      bool                 write_to_master_mean;/*!<data is written to mean master file
            										only if this is true. This ensures that
            										the master file data is compatible with
            										both steady and unsteady simulations.*/
      bool                 write_to_master_sec_mom;/*!<data is written to second moment master file
            										only if this is true. This ensures that
            										the master file data is compatible with
            										both steady and unsteady simulations.*/      										
      bool                 write_to_master_var; /*!<data is written to variance master file
            										only if this is true. This ensures that
            										the master file data is compatible with
            										both steady and unsteady simulations.*/
      std::string          SAMPLE_DIR;          /*!<path to directory where solution files 
      												for a sample are stored*/
      std::string          MEAN_DIR;            /*!<path to directory where solution files 
      												for a sample mean are stored*/
      std::string          SEC_MOM_DIR;         /*!<path to directory where solution files 
      												for a sample second moments are stored*/												
      std::string          VAR_DIR;				/*!<path to directory where solution files 
      												for a variance are stored*/
      // std::string          PDF_DIR;             /*!<path to directory where pdf at probe nodes
//        											    are stored*/
      std::string          L1_DIR;              /*!<path to directory L1 norm of mean and variance
       											    are stored*/ 										    											    												
      int                  counter;             /*!<counter required to index solution files.
      												This is incremented only for unsteady
      												simulations.*/
      bool                 has_density,         /*!<If true, data structures for density 
      												statistics are constructed*/
                           has_mach, 			/*!<If true, data structures for mach number 
      												statistics are constructed*/
                           has_entropy,         /*!<If true, data structures for physical entropy 
      												statistics are constructed*/
                           has_vorticity;       /*!<If true, data structures for vorticity 
      												statistics are constructed*/
      
      int   restart_head_tag;					/*!<0 or 1. Used for saving restart files in
      												alternating sets*/							
      
      // For checking run times
      Timer t_total,         /*!<Total simulation time */
            t_grid,          /*!<Time for grid generation */
            t_init,          /*!<Time for pre-sample initialization */
            t_sample_init,   /*!<Time for sample initialization */
            t_sample_solve,  /*!<Time for sample solve (each sample) */
            t_grad_eval,     /*!<Time for gradient evaluation (each sample) */
            t_grad_red,      /*!<Time for gradient reduction (each sample) */
            t_mom_inc,       /*!<Time for moment increment (each sample) */
            t_mom_red,       /*!<Time for reduce moments */
//             t_pdf_eval,      /*!<Time for pdf evaluation (each sample) */
            t_dt_eval,       /*!<Time for dt evaluation (each sample) */
            t_dt_red,        /*!<Time for dt reduction (each sample) */
            t_res_eval,      /*!<Time for residual reduction (each sample) */
            t_ires_eval,     /*!<Time for inviscid residual evaluation (each sample) */
            t_vres_eval,     /*!<Time for viscous residual evaluation (each sample) */
            t_res_red,       /*!<Time for residual reduction (each sample) */
            t_res_norm,      /*!<Time for residual norm evaluation (each sample) */
            t_res_norm_red,  /*!<Time for residual norm reduction (each sample) */
            t_sol_update,    /*!<Time for solution update (each sample) */
            t_sol_out,       /*!<Time for solution output (each sample) */
            t_force_eval,    /*!<Time for face force evaluation (each sample) */
            t_global_eval,   /*!<Time for face force evaluation (each sample) */
            t_log,           /*!<Time for log output (each sample) */
//             t_build_pdf,     /*!<Time for building pdf */
            t_l1_norm,       /*!<Time for evaluation of L1 norm of statistics */
            t_mean_out,      /*!<Time for mean output */
            t_var_out;       /*!<Time for var output */
            
      // Monte-Carlo variable
      std::vector<std::vector<PrimVar> > primitive_fm;   /*!<First moment of primitive variables. 
       														 The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/
      std::vector<std::vector<double> > density_fm;      /*!<First moment of density. 
       														 The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/
      std::vector<std::vector<double> > mach_fm;         /*!<First moment of mach number. 
       														 The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/
      std::vector<std::vector<double> > entropy_fm;      /*!<First moment of physical entropy. 
       														 The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/
      std::vector<std::vector<double> > vorticity_fm;    /*!<First moment of vorticity. 
       														 The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/
      std::vector<std::vector<PrimVar> >primitive_sm;    /*!<Second moment primitive  
       														 variable. The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/
      std::vector<std::vector<double> >density_sm;   	 /*!<Second moment density.
      														 The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/
      std::vector<std::vector<double> > mach_sm;		 /*!<Second moment mach 
      														 number. The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/
      std::vector<std::vector<double> > entropy_sm;      /*!<Second moment physical 
      														 entropy. The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/
      std::vector<std::vector<double> > vorticity_sm;    /*!<Second moment vorticity. 
      														 The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/
      std::vector<std::vector<PrimVar> >primitive_variance;    /*!<Variance of primitive  
       														 variable. The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/
      std::vector<std::vector<double> >density_variance;   	 /*!<Variance of density.
      														 The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/
      std::vector<std::vector<double> > mach_variance;		 /*!<Variance of mach 
      														 number. The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/
      std::vector<std::vector<double> > entropy_variance;      /*!<Variance of physical 
      														 entropy. The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/
      std::vector<std::vector<double> > vorticity_variance;    /*!<Variance of vorticity. 
      														 The first vector corresponds
       														 to the number of time frames,
       														 the second one corresponds to
       														 the mesh nodes*/ 														 
      std::vector<PrimVar> prim_mean_l1;   /*!< l1 norm of mean of primary variable*/
      std::vector<PrimVar> prim_var_l1;	   /*!< l1 norm of variance of primary variable*/
      std::vector<double> den_mean_l1;     /*!< l1 norm of mean of density*/
      std::vector<double>  den_var_l1;     /*!< l1 norm of variance of density*/
      std::vector<double>  mach_mean_l1;   /*!< l1 norm of mean of mach number*/
      std::vector<double>  mach_var_l1;    /*!< l1 norm of variance of density*/
      std::vector<double> ent_mean_l1;     /*!< l1 norm of mean of physical entropy*/
      std::vector<double> ent_var_l1;      /*!< l1 norm of variance of density*/
      std::vector<double> vor_mean_l1;     /*!< l1 norm of mean of vorticity*/
      std::vector<double> vor_var_l1; 	   /*!< l1 norm of variance of density*/												 
      
      std::vector<double>  time_instance;  /*!<Vector of times at which soln/statistics are evaluated.
      									       This is useful only for unsteady flows.  */
//       unsigned int  time_instance_ind;        /*!<This varies from 0 to total time frame - 1.
//        					                          When a sample run is finished, this is reset to 0*/
      bool time_instance_reached;            /*!<Indicates if a time frame of soln/statistical evaluation
      									         has been reached for unsteady flows.*/
      
      // Errors when exact solution evailable									 
      PrimVar prim_exact,prim_abs_err, prim_err_L1, prim_err_L2,prim_err_Linf;
      double density_exact , entropy_exact , mach_exact ,
             density_abs_err , entropy_abs_err , mach_abs_err ,
             density_err_L1 , entropy_err_L1 , mach_err_L1 ,
             density_err_L2 , entropy_err_L2 , mach_err_L2 ,
             density_err_Linf , entropy_err_Linf , mach_err_Linf ;									 
      									 
      void compute_ssw();                          /*!<Evaluates shock indicators for Liou's fix*/
      void initialize ();                          /*!<Initializes data structures for solver*/
      void initialize_sample_run (int sample_id);  /*!<Initializes data structures for a sample*/
      void end_sample_run();                       /*!<Cleans up after a sample run*/
  
      void compute_gradients ();                   /*!<Computes cell and nodal gradients of 
      												   entropy variables*/
      void store_conserved_old ();				   /*!<Stores old conserved variables for 
      												   multi-level RK*/
      void compute_inviscid_residual ();		   /*!<Computes inviscid residual*/
      void compute_viscous_residual ();            /*!<Computes viscous residual*/
      void compute_residual ();				       /*!<Computes total residual*/
      void compute_dt ();						   /*!<Computes local time step*/
      
      void compute_residual_norm (const unsigned long int iter);  /*!<Computes residual norm*/
      void log_messages (const unsigned long int iter);           /*!<Prints residual norm into 
      													  		 log file*/
      void update_solution (const unsigned int r);			 /*!<Updates solution in RK
      													  		 step*/
      void solve (const int sample_id);						 /*!<Main solve for a sample*/
      void compute_bounds ();        					     /*!<Computes solution bounds*/
      void output (bool write_variables = true);			 /*!<Writes output files*/	
      void output_restart (unsigned long int iter);						 /*!<Writes restart files*/   
      void output_surface_sf (std::string DIR);              /*!<Writes surface skin-friction*/
      void output_surface_hf (std::string DIR);             /*!<Writes surface heat flux*/
      void output_mesh_Pe(std::string DIR);                  /*!<Writes mesh-Pe number*/
      void create_force_face_list ();                        /*!<Creates force-face listing*/
      void compute_forces (unsigned long int iter);				 /*!<Computes forces*/
      void compute_global (unsigned long int iter);				 /*!<Computes global quantities*/
      void compute_error_norm ();				             /*!<Computes error norms*/
      
      void reduce_gradient_periodic();           /*!<Reduces gradients across periodic nodes*/
      void reduce_residual();			         /*!<Reduces residual across shared nodes*/
      void reduce_residual_periodic();			 /*!<Reduces residual across periodic nodes*/
      void reduce_res_norm();			 /*!<Reduces residual norm across partitions*/
      void reduce_err_norm();			 /*!<Reduces error norm across partitions*/
      void reduce_dt();                  /*!<Reduces time step across shared nodes*/
      void reduce_dt_periodic();                  /*!<Reduces time step across periodic nodes*/
      void write_time_output();          /*!<Writes simulation timing info*/
      void write_sample_time_output();   /*!<Writes sample simulation timing info*/
      
      void increment_moments(const int sample_id);   /*!<Reduces moments across shared nodes*/
      // void find_pdf (const int sample_id);           /*!<Evaluates sample pdf contribution*/
//       void save_sample_pdf();					     /*!<Saves sample pdfs from each MC group*/
//       void save_pdf();                               /*!<Saves final pdfs*/
      void reduce_moments ();                        /*!<Reduces moments*/
      void output_mean();                            /*!<Saves sample mean*/
      void output_sec_mom();                         /*!<Saves sample second moments*/
      void output_var();                             /*!<Saves sample variance*/
      void extract_rnd_nos(int sample_id);           /*!<Extracts random from file */													 
      void compute_stat_norm();                      /*!<Computes L1 norm of mean and variance*/
      void reduce_stat_norm();						 /*!<Reduces L1 norm of mean and variance*/	
      void save_stat_norm();		                 /*!<Saves L1 norm of mean and variance*/						 
      
};

#endif
