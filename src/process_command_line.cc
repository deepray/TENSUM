#include <iostream>
#include <cstring>
#include <cstdlib>
#include "constants.h"
#include "mpi_utils.h"

extern Dimension dim;
extern bool debug;
extern bool restart;
extern bool preprocess;
extern bool bounds;
extern bool verbose;
using namespace std;

void show_options ();

//------------------------------------------------------------------------------
// Get command line flags and input file
//------------------------------------------------------------------------------ 
/*!
  Function processing the command line statement and run options/parameters.
*/
void process_command_line (int   argc,
                           char* argv[],
                           int&  ifile,
                           int&  iloc,
                           unsigned int&  NPART)
{ 
   if(argc < 6)
      show_options ();

   // By default, all are off
   debug          = false;
   restart        = false;
   preprocess     = false;
   bounds         = false;
   verbose        = false;

   // Default is 2d
   dim            = two;

   bool found_input_file = false;
   bool found_npart      = false;
   bool found_loc        = false;
   
   int i = 1;

   while (i < argc)
   {
      if(strcmp(argv[i],"-axi")==0)
      {
         dim = axi;
      }
      else if(strcmp(argv[i],"-d")==0)
      {
         debug = true;
      }
      else if(strcmp(argv[i],"-r")==0)
      {
         restart = true;
      }
      else if(strcmp(argv[i],"-p")==0)
      {
         preprocess = true;
      }
      else if(strcmp(argv[i],"-b")==0)
      {
         bounds = true;
      }
      else if(strcmp(argv[i],"-v")==0)
      {
         verbose = true;
      }
      else if(strcmp(argv[i],"-i")==0)
      {
         MPI_ASSERT (i+1 < argc); // check that there is another argument
         ifile = i+1;
         ++i;
         found_input_file = true;
      }
      else if(strcmp(argv[i],"-npart")==0)
      {
         MPI_ASSERT (i+1 < argc); // check that there is another argument
         NPART = atoi(argv[i+1]);
         ++i;
         found_npart = true;
      }
      else if(strcmp(argv[i],"-L")==0)
      {
         MPI_ASSERT(i+1 < argc); // check that there is another argument
         iloc = i+1;
         ++i;
         found_loc = true;
      }
      else
      {
         if(get_proc_id() == 0)
            cout << "Unknown command line flag: " << argv[i] << endl;
         show_options ();
      }

      ++i;
   }

   if(!found_input_file || !found_npart || !found_loc)
      show_options ();
}

//------------------------------------------------------------------------------
// Print command line options available
//------------------------------------------------------------------------------
/*!
  Function to display command line options
*/
void show_options ()
{
   int rank = get_proc_id();
   if(rank == 0)
   {
	   cout << "Valid flags are:\n";
	   cout << "   -i filename   Specify input file name (required)\n";
	   cout << "   -npart N      Specify number of mesh partitions N (required)\n";
	   cout << "   -L path       Location of mesh partition folder (required)\n";
	   cout << "   -axi          Axisymmetric flow (optional)\n";
	   cout << "   -d            Enable debug mode (optional)\n";
	   cout << "   -r            Read restart file for initial condition (optional)\n";
	   cout << "   -p            Do everything but do not solve (optional)\n";
	   cout << "   -b            Compute min/max range of solution (optional)\n";
	   cout << "   -v            Verbose (optional)\n";
   }
   mpi_err_finalize();
}
