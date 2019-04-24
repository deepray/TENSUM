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
// Saving L1 norm of mean and variance
//------------------------------------------------------------------------------
void FiniteVolume::save_stat_norm()
{   
   if(get_proc_id() == 0)
   {
       ofstream mean_l1_f;
       
       string mean_file = L1_DIR + "/mean_l1_norm.dat";
       
       mean_l1_f.open(mean_file.c_str());
       mean_l1_f << " TIME"
                 << " TEMPERATURE"
                 << " U"
                 << " V"
                 << " W"
                 << " PRESSURE";
       if(has_density)
          mean_l1_f << " DENSITY";
       if(has_mach)
          mean_l1_f << " MACH";
       if(has_vorticity)
          mean_l1_f << " VORTICITY"; 
       if(has_entropy)
          mean_l1_f << " ENTROPY";     
       mean_l1_f << "\n";
       
       for(unsigned int t = 0; t<time_instance.size(); ++t)
       {
           mean_l1_f << t
                     <<" "<<prim_mean_l1[t].temperature
                     <<" "<<prim_mean_l1[t].velocity.x
                     <<" "<<prim_mean_l1[t].velocity.y
                     <<" "<<prim_mean_l1[t].velocity.z
                     <<" "<<prim_mean_l1[t].pressure;      
           if(has_density)
			  mean_l1_f << " "<<den_mean_l1[t];
		   if(has_mach)
			  mean_l1_f << " "<<mach_mean_l1[t];
		   if(has_vorticity)
			  mean_l1_f << " "<<vor_mean_l1[t];
		   if(has_entropy)
			  mean_l1_f << " "<<ent_mean_l1[t];     
		   mean_l1_f << "\n";           
       }                 
       mean_l1_f.close();
       
       
       ofstream var_l1_f;
       
       string var_file = L1_DIR + "/var_l1_norm.dat";
       
       var_l1_f.open(var_file.c_str());
       var_l1_f << " TIME"
                 << " TEMPERATURE"
                 << " U"
                 << " V"
                 << " W"
                 << " PRESSURE";
       if(has_density)
          var_l1_f << " DENSITY";
       if(has_mach)
          var_l1_f << " MACH";
       if(has_vorticity)
          var_l1_f << " VORTICITY"; 
       if(has_entropy)
          var_l1_f << " ENTROPY";     
       var_l1_f << "\n";
       
       for(unsigned int t = 0; t<time_instance.size(); ++t)
       {
           var_l1_f << t
                     <<" "<<prim_var_l1[t].temperature
                     <<" "<<prim_var_l1[t].velocity.x
                     <<" "<<prim_var_l1[t].velocity.y
                     <<" "<<prim_var_l1[t].velocity.z
                     <<" "<<prim_var_l1[t].pressure;          
           if(has_density)
			  var_l1_f << " "<<den_var_l1[t];
		   if(has_mach)
			  var_l1_f << " "<<mach_var_l1[t];
		   if(has_vorticity)
			  var_l1_f << " "<<vor_var_l1[t];
		   if(has_entropy)
			  var_l1_f << " "<<ent_var_l1[t];     
		   var_l1_f << "\n";          
       }                 
       var_l1_f.close();
   }
}   