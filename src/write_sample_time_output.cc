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
// Timing output for sample
//------------------------------------------------------------------------------
void FiniteVolume::write_sample_time_output()
{
      // Time reductions
      t_sample_init.reduce(grid.run_comm);
      t_sample_solve.reduce(grid.run_comm);
      t_grad_eval.reduce(grid.run_comm);
      t_grad_red.reduce(grid.run_comm);
      t_dt_eval.reduce(grid.run_comm);
      t_dt_red.reduce(grid.run_comm);
      t_res_eval.reduce(grid.run_comm);
      t_ires_eval.reduce(grid.run_comm);
      t_vres_eval.reduce(grid.run_comm);
      t_res_red.reduce(grid.run_comm);
      t_res_norm.reduce(grid.run_comm);
      t_res_norm_red.reduce(grid.run_comm);
      t_sol_update.reduce(grid.run_comm);
      t_sol_out.reduce(grid.run_comm);
      t_force_eval.reduce(grid.run_comm);
      t_global_eval.reduce(grid.run_comm);
      t_log.reduce(grid.run_comm);
      
      mpi_barrier(grid.run_comm);
      
      if(check_group_base())
      {
         ofstream time_log;
         string t_file_name = SAMPLE_DIR + "/time_log.log";
         time_log.open(t_file_name.c_str()); 
         time_log <<"Time log for sample (minutes)"<<endl<<endl; 
         time_log << left << setw(25) << setfill(' ') << "" 
                  << left << setw(17) << setfill(' ') << "T_MIN"
                  << left << setw(17) << setfill(' ') << "T_MAX"
                  << left << setw(17) << setfill(' ') << "T_AVG" <<endl;
         time_log << left << setw(25) << setfill(' ') << "SAMPLE_INITIALIZE" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_sample_init.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_sample_init.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_sample_init.t_avg << endl;
         time_log << left << setw(25) << setfill(' ') << "SAMPLE_SOLVE" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_sample_solve.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_sample_solve.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_sample_solve.t_avg << endl;
         time_log << left << setw(25) << setfill(' ') << "GRADIENT_EVALUATION" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_grad_eval.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_grad_eval.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_grad_eval.t_avg << endl;              
         time_log << left << setw(25) << setfill(' ') << "GRADIENT_REDUCE" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_grad_red.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_grad_red.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_grad_red.t_avg << endl; 
         time_log << left << setw(25) << setfill(' ') << "TIME_STEP_EVALUATION" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_dt_eval.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_dt_eval.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_dt_eval.t_avg << endl; 
         time_log << left << setw(25) << setfill(' ') << "TIME_STEP_REDUCE" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_dt_red.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_dt_red.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_dt_red.t_avg << endl; 
         time_log << left << setw(25) << setfill(' ') << "RESIDUAL_EVALUATION" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_res_eval.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_res_eval.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_res_eval.t_avg << endl; 
         time_log << left << setw(25) << setfill(' ') << "INVISCID_RESIDUAL" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_ires_eval.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_ires_eval.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_ires_eval.t_avg << endl; 
         time_log << left << setw(25) << setfill(' ') << "VISCOUS_RESIDUAL" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_vres_eval.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_vres_eval.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_vres_eval.t_avg << endl;                                                
         time_log << left << setw(25) << setfill(' ') << "RESIDUAL_REDUCTION" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_res_red.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_res_red.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_res_red.t_avg << endl; 
         time_log << left << setw(25) << setfill(' ') << "RESIDUAL_NORM" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_res_norm.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_res_norm.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_res_norm.t_avg << endl; 
         time_log << left << setw(25) << setfill(' ') << "RESIDUAL_NORM_REDUCE" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_res_norm_red.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_res_norm_red.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_res_norm_red.t_avg << endl;    
         time_log << left << setw(25) << setfill(' ') << "SOLUTION_UPDATE" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_sol_update.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_sol_update.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_sol_update.t_avg << endl; 
         time_log << left << setw(25) << setfill(' ') << "SOLUTION_WRITE" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_sol_out.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_sol_out.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_sol_out.t_avg << endl; 
         time_log << left << setw(25) << setfill(' ') << "FORCE_EVALUATION" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_force_eval.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_force_eval.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_force_eval.t_avg << endl; 
         time_log << left << setw(25) << setfill(' ') << "GLOBAL_EVALUATION" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_global_eval.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_global_eval.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_global_eval.t_avg << endl; 
         time_log << left << setw(25) << setfill(' ') << "LOG_WRITE" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_log.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_log.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_log.t_avg << endl;                                                      
         time_log.close();
      }                          
      mpi_barrier(grid.run_comm);
}