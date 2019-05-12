#ifndef __PARAMETER_H__
#define __PARAMETER_H__

#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <cmath>
#include <set>
#include <algorithm>
#include "reader.h"
#include "material.h"
#include "force.h"
#include "ic.h"
#include "exact.h"
#include "bc.h"
#include "constants.h"

//!  Parameter class. 
/*!
  Reads solver parameter file and sets various simulation parameters
*/
class Parameter
{
   public:
      Parameter ()
      {};
      
      char* file;               /*!< Name of parameter file. */

      enum BCScheme
      {
         strong,     /*!< Strong implementation of boundary condition at boundary nodes */
         weak        /*!< Weak implementation of boundary condtions at boundary nodes */
      };
      
      // Grid section   
      CellType    cell_type; 
      std::string part_dir;     /*!< Relative path of directory containing mesh partition 
                                     files*/
      
      // Numerics sections
      std::string time_mode;        /*!< Options: steady, unsteady */
      std::string time_scheme;      /*!< Options: heuns, ssprk3, rk4 */
      unsigned int n_rks;           /*!< Number of rk steps */
      unsigned long int max_iter;   /*!< Maximum number of allowed iterations.
                                         Useful for steady simulations. Note that
                                         this does not effect unsteady simulations.*/
      double time_step;             /*!< Specifies global time_step to be used */
      double cfl;                   /*!< Specifies CFL condition. NOTE: cannot assign positive
                                     values for CFL and time_step simultaneously */
      double final_time;            /*!< Final time. Useful for unsteady simulations. NOTE: 
                                     steady simulations, final_time is set to 1.0e20 */
      double min_residue;           /*!< Stopping criteria for steady simulations */         
      double Cpen;                  /*!< Penalty term for Nitsche penalty method  */
      BCScheme bc_scheme;           /*!< Set boundary condition type  */
      bool liou_fix;                /*!< Use liou fix for handling Carbuncle */
      std::vector<int> SAMPLE_LIST; /*!< List of sample indices for indexing Monte-Carlo
                                         samples */
      unsigned int n_samples;       /*!< Total number of samples */  
      int  sample_start_ind,        /*!< Sample start index for an MC processor group */
           sample_end_ind;          /*!< Sample end index for an MC processor group */                      

      unsigned int rnd_per_sample;  /*!< Number of random numbers needed per sample */
      std::string rnd_file_loc ;    /*!< Path of file containing list of random number*/
      // unsigned int rnd_file_hook;   /*!< Position in random number file from where extraction
//                                          should commence */
      
      // Material section
      Material material;
      
      // IC section
      InitialCondition initial_condition;
      
      // Exact solution section
      ExactSolution exact_solution;
      bool exact_available;

      // BC section
      std::map<int,BoundaryCondition> boundary_condition; /*!< Maps boundary face to 
                                                               BC type */ 
      
      // Integral section
      std::vector<ForceData> force_data;         /*!< Stores face IDs for force evaluation */ 

      // Output section
      std::string  write_format;                 /*!< Options: vtk (currently the only 
       												  format available */
      unsigned int write_frequency;              /*!< Number of iterations after which
                                                      a solution file is saved (only used for
                                                      steady simulations */
      unsigned int restart_write_frequency;      /*!< Number of iterations after which
                                                      a restart solution file is saved */
                                                      
      std::vector<std::string> write_variables;  /*!< Additional variables to be written
                                                      in output file */
      int n_time_stamps;                         /*!< Number of time stamps at which 
                                                      samples are saved (excluding the
                                                      initial solution). NOTE: 
                                                      i) For steady deterministic flows,
                                                      this parameter is ignored and solution
                                                      files are saved on the basis of the 
                                                      write_frequency parameter. Only the 
                                                      initial and final solution files
                                                      are saved.
                                                      ii) For steady MC simulations, solution
                                                      is saved based on the write_frequency,
                                                      but the statistics will be evaluated
                                                      only for the initial state and the final
                                                      state of the solution. Thus, n_time_stamps
                                                      is ignored.
                                                      iii) For unsteady deterministic flows, 
                                                      if n_time_stamps > 0, then
                                                      write_frequency is ignored, otherwise
                                                      files are saved on the basis of the 
                                                      write_frequency parameter. 
                                                      iv) For unsteady MC simulations, we need
                                                      to choose n_time_stamps > 0 as 
                                                      write_frequency need not lead to 
                                                      statistical sampling at the same time
                                                      across all samples. In other words,
                                                      the write_frequency is ignored.*/   
      bool online_stat;                          /*!< If true, then online statistics 
                                                      of MC samples is performed */                                              
      std::vector<int> write_surfaces;           /*!< Tags of boundaries at which solution
                                                      is saved (also skin friction) */
      // std::set<int> write_mc_points;             /*!< Node ids (as in mesh file) at which 
//                                                       solution is stored to evaluate the PDF */
      bool write_restart;                        /*!< If true, a restart file is written */
      bool write_log;                            /*!< If true, log files are written/displayed */
      bool write_soln;                           /*!< If true, solution files are written */
      bool has_global;                           /*!< If true, global quantities are evaluated
                                                      and saved */
      // bool save_mesh_Pe;                         /*!< If true, then mesh Peclet number is evaluated 
//       												  and saved */  
      bool find_error;                            /*!< If true, then the Lp norms of the errors (provided 
      													the exact solution is available) */  												                                                

      void read ();                              /*!< Function to read parameter file */

   private:
      void read_constants (Reader&);           /*!< Reads constants section of parameter file */
      void read_grid (Reader&); 			   /*!< Reads grid section of parameter file */
      void read_numeric (Reader&);             /*!< Reads numeric section of parameter file */
      void read_material (Reader&);            /*!< Reads material section of parameter file */
      void read_initial_condition (Reader&);   /*!< Reads IC section of parameter file */
      void read_exact (Reader&);               /*!< Reads exact solution section of parameter file */
      void read_boundary (Reader&);            /*!< Reads BC section of parameter file */
      void read_integrals (Reader&);           /*!< Reads integral section of parameter file */
      void read_output (Reader&);              /*!< Reads output section of parameter file */
      void set_sample_ids();                   /*!< Sets sample start and end IDs for 
                                                    MC group */
      void rearrange_samples();                /*!< Rearranges ordering of samples to ensure 
                       								samples are simulated in an ascending
													order by various MC groups */                                              
};

#endif
