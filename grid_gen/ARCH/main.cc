#include <iostream>
#include <cstdlib>
#include <string>
#include <map>
#include "grid.h"

// 
// template struct Face<2>;
// template struct Face<3>;
// 
// template struct Cell<2>;
// template struct Cell<3>;

using namespace std;

//bool periodic; 
                          
void process_command_line (int argc, char* argv[], int& ifile, int& iloc, int& DIM);

int main(int argc, char* argv[])
{
   cout << "======================================================\n";
   cout << "   Starting mesh preprocessor for partitioning\n";   
   cout << "   ---  Author: Deep Ray, TIFR-CAM, Bangalore\n";
   cout << "   ---  Date  : 11th December, 2015\n";
   cout << "======================================================\n\n";  
   
   int ifile, iloc, DIM;
   process_command_line (argc, argv, ifile, iloc, DIM);
   
   if(DIM == 2)
   {
      Grid<2> grid(argv[ifile], argv[iloc]);
      grid.read();
   }   
   else if(DIM == 3)
   {
      cout<<"WARNING: 3D face matching not implemented yet. Need to create edges as well"<<endl;
	  exit(0);
      Grid<3> grid(argv[ifile], argv[iloc]);
      grid.read();
   }   
   
   
   cout << "\n\n---------------Partitioned files created!!----------------------\n";      
   
}


