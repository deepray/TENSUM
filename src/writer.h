#ifndef __WRITER_H__
#define __WRITER_H__

#include <vector>
#include <string>
#include <fstream>
#include "grid.h"
#include "material.h"
#include <fstream>
#include "mpi_utils.h"


//------------------------------------------------------------------------------
//! Writing class 
/*!
  This class is responsible of writing the solutions into output files
*/
//------------------------------------------------------------------------------
class Writer
{
   public:
      Writer (const Grid&     grid, std::string SAMPLE_DIR)
         : 
         grid (&grid),
         has_primitive (false),
         has_gradient (false),
         write_mach (false),
         write_density (false),
         write_entropy (false),
         write_vorticity (false),
         SAMPLE_DIR (SAMPLE_DIR)
         {};
      Writer (const Grid&      grid,
              const Material&  material,
              std::string      format,
              std::vector<int> surfaces,
              const double time,
              std::string time_mode,
              std::string SAMPLE_DIR) 
         : 
         grid (&grid),
         material (&material),
         format (format),
         surfaces (surfaces),
         has_primitive (false),
         has_gradient (false),
         write_mach (false),
         write_density (false),
         write_entropy (false),
         write_vorticity (false),
         SAMPLE_DIR(SAMPLE_DIR),
         time_mode(time_mode),
         time (time)
         {
            if(check_group_base())
            {
               std::string master_file_path = SAMPLE_DIR + "/master_file.visit";
		       master_file.open (master_file_path.c_str(),std::fstream::app);
		    }
         };
         
      Writer (const Grid&      grid,
              std::string      format,
              const double time,
              std::string time_mode,
              std::string DIR) 
         : 
         grid (&grid),
         format (format),
         has_primitive (false),
         has_gradient (false),
         write_mach (false),
         write_density (false),
         write_entropy (false),
         write_vorticity (false),
         SAMPLE_DIR(DIR),
         time_mode(time_mode),
         time (time)
         {
            if(get_proc_id() == 0)
            {
               std::string master_file_path = DIR + "/master_file.visit";
		       master_file.open (master_file_path.c_str(),std::fstream::app);
		    }
         };       
      void master_write_check(bool val){ write_to_master = val;}   
      void attach_data (std::vector<PrimVar>& data);
      void attach_data (std::vector<double>& data, std::string name);
      void attach_variables (const std::vector<std::string>& variables);
      void attach_gradient (std::vector<std::vector<Vector> >& dE);
      void output (int part, int counter);
      void output_vtk (std::string filename);
      void output_surfaces (std::string DIR_NAME, std::string index);
      void output_level (std::string levelfilename);
      void output_restart (unsigned long int iter, int counter,int head_tag, double time, double res0);
      
      ~Writer()
      {
         master_file.close();
      };

   private:

      const Grid*      grid;
      const Material*  material;
      std::string      format;
      std::vector<int> surfaces;

      std::vector< std::vector<double>* > vertex_data;
      std::vector<std::string> vertex_data_name;

      std::vector<PrimVar>* vertex_primitive;
      std::vector<std::vector<Vector> >* dE;
      bool has_primitive;
      bool has_gradient;
      bool write_mach;
      bool write_density;
      bool write_entropy;
      bool write_vorticity;
      bool write_mesh_Peclet;
      std::ofstream master_file;
      bool write_to_master;
      std::string SAMPLE_DIR;
      std::string time_mode;
      double time;

};

#endif
