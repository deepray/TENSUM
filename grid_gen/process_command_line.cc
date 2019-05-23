#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cassert>

//extern bool periodic;

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
                           int&  NPART,
                           int& DIM)
{ 
   if(argc < 9)
      show_options ();

   bool found_file  = false;
   bool found_loc   = false;
   bool found_dim   = false;
   bool found_npart = false;
   //periodic = false;
   
   
   int i = 1;

   while (i < argc)
   {
      if(strcmp(argv[i],"-D")==0)
      {
         assert(i+1 < argc); // check that there is another argument
         DIM = atoi(argv[i+1]);
         ++i;
         found_dim = true;
      }
      else if(strcmp(argv[i],"-I")==0)
      {
         assert(i+1 < argc); // check that there is another argument
         ifile = i+1;
         ++i;
         found_file = true;
      }
      else if(strcmp(argv[i],"-P")==0)
      {
         assert(i+1 < argc); // check that there is another argument
         NPART = atoi(argv[i+1]);
         ++i;
         found_npart = true;
      }
      else if(strcmp(argv[i],"-L")==0)
      {
         assert(i+1 < argc); // check that there is another argument
         iloc = i+1;
         ++i;
         found_loc = true;
      }
      else
      {
         std::cout << "Unknown command line flag: " << argv[i] << std::endl;
         show_options ();
      }

      ++i;
   }

   if ( (!found_file || !found_dim || !found_loc || !found_npart) ||
        (DIM != 2 && DIM !=3) ||
        NPART < 1
      )  
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
   std::cout << "Valid flags are:\n";
   std::cout << "   -D dim      dim=2,3 (required)\n";
   std::cout << "   -I file     Specify mesh file name (required)\n";
   std::cout << "   -P npart    Specify number of mesh partitions (required)\n";
   std::cout << "   -L loc      Location where partition directory is to be created (required)\n";   
}
