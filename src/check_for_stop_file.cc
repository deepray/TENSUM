#include <iostream>
#include <fstream>
#include "mpi_utils.h"

using namespace std;

bool check_for_stop_file ()
{
   ifstream ifile("STOP");
   if(ifile)
      if(get_proc_id() == 0)
         cout << "Found STOP file !!!\n";
   return ifile;
}
