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
// Timing output
//------------------------------------------------------------------------------
void FiniteVolume::write_time_output()
{

      string remove_old_log = "rm -rf time_log.log";
      system(remove_old_log.c_str());
      mpi_barrier(MPI_COMM_WORLD);
      
      // Time reductions
      t_total.reduce(MPI_COMM_WORLD);
      t_grid.reduce(MPI_COMM_WORLD);
      t_init.reduce(MPI_COMM_WORLD);
      t_mom_inc.reduce(MPI_COMM_WORLD);
      t_mom_red.reduce(MPI_COMM_WORLD);
      // t_pdf_eval.reduce(MPI_COMM_WORLD);
//       t_build_pdf.reduce(MPI_COMM_WORLD);
      t_l1_norm.reduce(MPI_COMM_WORLD);
      t_mean_out.reduce(MPI_COMM_WORLD);
      t_var_out.reduce(MPI_COMM_WORLD);
      
      mpi_barrier(MPI_COMM_WORLD);
      
      if(get_proc_id() == 0)
      {
         ofstream time_log;
         time_log.open("time_log.log"); 
         time_log <<"Time log for full simulation (minutes)"<<endl<<endl; 
         time_log << left << setw(25) << setfill(' ') << "" 
                  << left << setw(17) << setfill(' ') << "T_MIN"
                  << left << setw(17) << setfill(' ') << "T_MAX"
                  << left << setw(17) << setfill(' ') << "T_AVG" <<endl;
         time_log << left << setw(25) << setfill(' ') << "TOTAL_TIME" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_total.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_total.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_total.t_avg << endl;
         time_log << left << setw(25) << setfill(' ') << "INITIALIZE" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_init.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_init.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_init.t_avg << endl;
         time_log << left << setw(25) << setfill(' ') << "MESH_DATA_CREATION" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_grid.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_grid.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_grid.t_avg << endl;
         time_log << left << setw(25) << setfill(' ') << "MOMENT_INCREMENT" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_mom_inc.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_mom_inc.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_mom_inc.t_avg << endl;              
         time_log << left << setw(25) << setfill(' ') << "MOMENT_REDUCE" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_mom_red.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_mom_red.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_mom_red.t_avg << endl; 
         // time_log << left << setw(25) << setfill(' ') << "PDF_EVALUATION" 
//                   << left << setw(17) << setfill(' ') << setprecision(3)<< t_pdf_eval.t_min
//                   << left << setw(17) << setfill(' ') << setprecision(3)<< t_pdf_eval.t_max
//                   << left << setw(17) << setfill(' ') << setprecision(3)<< t_pdf_eval.t_avg << endl; 
//          time_log << left << setw(25) << setfill(' ') << "BUILD_PDF" 
//                   << left << setw(17) << setfill(' ') << setprecision(3)<< t_build_pdf.t_min
//                   << left << setw(17) << setfill(' ') << setprecision(3)<< t_build_pdf.t_max
//                   << left << setw(17) << setfill(' ') << setprecision(3)<< t_build_pdf.t_avg << endl; 
         time_log << left << setw(25) << setfill(' ') << "L1_NORM_OF_STATS" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_l1_norm.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_l1_norm.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_l1_norm.t_avg << endl; 
         time_log << left << setw(25) << setfill(' ') << "WRITE_MEAN" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_mean_out.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_mean_out.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_mean_out.t_avg << endl; 
         time_log << left << setw(25) << setfill(' ') << "WRITE_VARIANCE" 
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_var_out.t_min
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_var_out.t_max
                  << left << setw(17) << setfill(' ') << setprecision(3)<< t_var_out.t_avg << endl;                                                                                                     
         time_log.close();
      }                      
}