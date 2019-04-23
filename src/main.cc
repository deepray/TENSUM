#include <iostream>
#include <cstdlib>
#include <string>
#include <map>
#include "fv.h"
#include "mpi_utils.h"

using namespace std;

Dimension dim; 
bool debug;                   /*!< Extern variable. Switched on by the -d option. If true,
                                   then primary and dual cells are printed and saved. */
bool restart;                 /*!< Extern variable. Switched on by the -r option. If true,
                                   then initial solution set from the restart file 
                                   (only if restart file is available). */
bool preprocess;              /*!< Extern variable. Switched on by the -p option. If true,
                                   then the solver does everything except the actual 
                                   solve. */
bool bounds;                  /*!< Extern variable. Switched on by the -b option. If true,
                                   then solution bounds are shown. */
bool verbose;                 /*!< Extern variable. Switched on by the -v option. If true,
                                   then verbose mode is active. */    
bool PERIODIC;                /*!< Extern variable. It is set to true if periodic boundary 
								   conditions are specified in solver parameter file. */                                                                      
map<string,double> constants; /*!< Extern variable. Map of constants specified in the 
                                   solver parameter file. */
unsigned int NPART;           /*!< Extern variable. Number of mesh partitions. */
unsigned int N_MC_GROUPS;     /*!< Extern variable. Number of Monte-Carlo groups. */

vector<double> rnd_nos;       /*!< Extern variable. Stores random numbers needed for each 
								   Monte-Carlo sample. */

void process_command_line (int argc, char* argv[], int& ifile, int& iloc, unsigned int& NPART);

int main(int argc, char* argv[])
{
   mpi_init(argc,argv); 
   int rank = get_proc_id();
   
   if(rank == 0)
   {
      cout << "==============================================================================\n";
      cout << "   Starting TENSUM: 2D entropy stable unstructured mesh solver\n";   
      cout << "   ---  Author: Deep Ray, TIFR-CAM, Bangalore\n";
      cout << "   ---  Date  : 30th October, 2015\n";
      cout << "   TENSUM is parallelized and can be used for Monte-Carlo simulations.\n";
      cout << "   (The code is based on the TAXIS solver by Praveen C., TIFR-CAM Bangalore.)\n";
      cout << "==============================================================================\n\n";
   }   
   int ifile, iloc;
   process_command_line (argc, argv, ifile, iloc, NPART);
   PERIODIC = false; // Default
   
   if(get_comm_size() % NPART !=0)
      MPI_ERR("The number of mesh partitions must divide the number of processors called "
              <<"by mpi EXACTLY!! Please change this and re-run the code.")
   
   N_MC_GROUPS = get_comm_size()/NPART;
   
   FiniteVolume problem (argv[ifile],argv[iloc]);
   problem.run ();
   
   if(rank == 0)
   {
      cout << "\n\n---------------Solver has finished!!----------------------\n";   
   }
   mpi_finalize();
   
}
