#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <cstdlib>
#include "parameter.h"
#include "mpi_utils.h"

using namespace std;
extern map<string,double> constants;
extern bool verbose;

//------------------------------------------------------------------------------
// Read parameters from file
//------------------------------------------------------------------------------
void Parameter::read ()
{
   if(get_proc_id() == 0)
      cout << "  Reading input file " << file << endl;
   Reader fin(file);

   read_grid (fin);
   read_numeric (fin);
   read_material (fin);
   read_constants (fin);
   read_initial_condition (fin);
   read_boundary (fin);
   read_exact (fin);
   read_integrals (fin);
   read_output (fin);
   set_sample_ids();
   rearrange_samples();
}

//------------------------------------------------------------------------------
// Read grid section
//------------------------------------------------------------------------------
void Parameter::read_grid (Reader &fin)
{
   MPI_DISP("  --- Reading grid section",verbose);

   string input;

   fin.begin_section ("grid");

   fin.entry ("cell");
   fin >> input;
   if(input=="median")
      cell_type = median;
   else if(input=="voronoi")
      cell_type = voronoi;
   else
      MPI_ERR("read_grid: unknown cell type " << input
      		  << ". Possible options: \n"
      		  <<" median \n voronoi \n");
   
   fin.end_section ();
}

//------------------------------------------------------------------------------
// Read numeric section
//------------------------------------------------------------------------------
void Parameter::read_numeric (Reader &fin)
{
   MPI_DISP("  --- Reading numeric section",verbose);

   string input;

   fin.begin_section ("numeric");

   fin.entry ("time_mode");
   fin >> time_mode;
   if (time_mode != "steady" && time_mode != "unsteady")
   {
      MPI_ERR("read_numeric: unknown time_mode " << time_mode
               << ". Possible options:\n"
               << " steady\n unsteady\n");  
   }

   fin.entry ("time_scheme");
   fin >> time_scheme;
   if (time_scheme != "heuns"  &&
       time_scheme != "ssprk3" &&
       time_scheme != "rk4")
   {        
       MPI_ERR("read_numeric: unknown time_scheme " << time_scheme
               << ". Possible options:\n"
               << " heuns\n ssprk3\n rk4\n");
   }        
   if(time_scheme=="heuns")  n_rks = 2; 
   if(time_scheme=="ssprk3") n_rks = 3;
   if(time_scheme=="rk4")    n_rks = 4;

   fin.entry ("time_step");
   fin >> time_step;
   MPI_ASSERT(time_step >= 0.0);

   fin.entry ("cfl");
   fin >> cfl;
   MPI_ASSERT(cfl >= 0.0); 

   fin.entry ("max_iter");
   fin >> max_iter;
   MPI_ASSERT(max_iter > 0);

   fin.entry ("final_time");
   fin >> final_time;
   MPI_ASSERT(final_time >= 0.0)

   fin.entry ("min_residue");
   fin >> min_residue;
   MPI_ASSERT(min_residue >= 0.0)

   fin.begin_section("reconstruct");
   fin.entry("method");
   fin >> input;
   if(input == "first")
         material.reconstruct_scheme = Material::first;
   else if(input == "second")
	  material.reconstruct_scheme = Material::second;   
   else if(input == "minmod")
   {
	  material.reconstruct_scheme = Material::tvd_minmod;
   }   
   else if(input == "tvb_minmod")
   {
	  material.reconstruct_scheme = Material::tvb_minmod;  
	  fin.entry("M_value"); 
	  fin >> material.tvb_minmod_coef;
	  MPI_ASSERT(material.tvb_minmod_coef > 0.0);
   } 
   else if(input == "eno2")
	  material.reconstruct_scheme = Material::eno2;  
   else if(input == "van_albada")
   {
	  material.reconstruct_scheme = Material::van_albada;
   }      
   else
	  MPI_ERR("read_numeric: unknown reconstruction scheme " << input
			   << ". Possible options:\n"
			   << " first \n second \n minmod \n eno2 \n"
			   <<" tvb_minmod <tvb_minmod_coeff> \n van_albada\n");
   fin.end_section();

   
   fin.entry ("bc_scheme");
   fin >> input;
   if(input == "strong")
      bc_scheme = Parameter::strong;
   else if(input == "weak")
      bc_scheme = Parameter::weak;
   else
      MPI_ERR("read_numeric: unknown bc scheme " << input
              << ". Possible options:\n"
			  << " strong \n weak \n");

   // fin.entry ("nitsche_pen");
//    fin >> Cpen;
//    MPI_ASSERT(Cpen >= 0.0);
   
   fin.entry ("liou_fix");
   fin >> input;
   if(input == "yes")
      liou_fix = true;
   else if(input == "no")
      liou_fix = false;
   else
      MPI_ERR("read_numeric: unknown option for liou_fix -> " << input
              << ". Possible options:\n"
			  << " yes \n no \n");
        
   fin.begin_section("sample_list");
   fin.begin_section("groups");
   while (!fin.eos())
   {
      unsigned int s_start, s_end;
      fin >> s_start;
      fin >> s_end;
      for(unsigned int i = s_start; i<=s_end; ++i)
         SAMPLE_LIST.push_back(i);
   } 
   fin.begin_section("free_list");
   while (!fin.eos())
   { 
      int sample_no;
      fin >> sample_no;
      SAMPLE_LIST.push_back(sample_no);
   }
   fin.end_section();

   
   fin.entry("rnd_file_loc");
   fin >> rnd_file_loc;
   
   fin.entry("rnd_per_sample");
   fin >> rnd_per_sample;
   MPI_ASSERT(rnd_per_sample >= 0);
   
   
   fin.end_section ();
   
   if(SAMPLE_LIST.empty())
      SAMPLE_LIST.push_back(0);
   
   sort( SAMPLE_LIST.begin(), SAMPLE_LIST.end() );
   SAMPLE_LIST.erase( unique( SAMPLE_LIST.begin(), SAMPLE_LIST.end() ), SAMPLE_LIST.end() );
   n_samples = SAMPLE_LIST.size();
   if(n_samples < N_MC_GROUPS)
      MPI_ERR("read_numeric: number of samples must be at least (no. of proc)/(no. of mesh partitions)"
              << " = "<<N_MC_GROUPS);

    
   
   if(get_proc_id() == 0)
   {
	   cout<<"  Requested samples: ";
	   int sample_block_start_id = SAMPLE_LIST[0];
	   int sample_block_end_id = SAMPLE_LIST[0];
	   unsigned int ind = 1;
	   while(ind<SAMPLE_LIST.size())
	   {
		   if(SAMPLE_LIST[ind] == sample_block_end_id + 1)
		   {
			  sample_block_end_id++;
			  ind++;
		   }   
		   else
		   {
			  if(sample_block_start_id == sample_block_end_id)
			     cout<<" "<<sample_block_start_id<<", ";
			  else   
			     cout<<" "<<sample_block_start_id<<"-"<<sample_block_end_id<<", ";
			  sample_block_start_id = SAMPLE_LIST[ind];
			  sample_block_end_id = SAMPLE_LIST[ind]; 
			  ind++;  
		   }   
	   }
	   if(sample_block_start_id == sample_block_end_id)
		  cout<<" "<<sample_block_start_id<<endl;
	   else   
		  cout<<" "<<sample_block_start_id<<"-"<<sample_block_end_id<<endl;
   }

   if(time_step == 0.0 && cfl == 0.0)
      MPI_ERR("Either time_step or cfl must be non-zero\n");
   if(time_step > 0.0 && cfl > 0.0)
      MPI_ERR("You have specified both time_step and cfl. "
           << "Specify only one of them, and set the other to zero\n");

   // For steady flow, no need for final time
   if(time_mode == "steady")
      final_time = 1.0e20;
   else
   {
      if (final_time <= 0.0)
         MPI_ERR("read_numeric: final_time must be positive for unsteady flow ");
      // For unsteady flow, final time must be specified, so dont put any limit
      // on max_iter
   }
   
   if(rnd_per_sample > 0)
   {
      ifstream ifile(rnd_file_loc);
      if(!ifile)
      {
         MPI_ERR("read_numeric: the random parameter file " << rnd_file_loc
                 <<" doesn't exist!!");
      }
   }   

}

//------------------------------------------------------------------------------
// Read material section
//------------------------------------------------------------------------------
void Parameter::read_material (Reader &fin)
{
   MPI_DISP("  --- Reading material section",verbose);

   string input;

   fin.begin_section ("material");

   fin.entry ("gamma");
   fin >> material.gamma;
   MPI_ASSERT(material.gamma > 1.0);

   fin.entry ("gas_const");
   fin >> material.gas_const;
   MPI_ASSERT(material.gas_const > 0.0);
   
   fin.begin_section ("viscosity");
   fin.entry ("model");
   fin >> input;
   if(input == "constant")
   {
      material.mu_model = Material::mu_constant;

      fin.entry ("mu_ref");
      fin >> material.mu_ref;
      MPI_ASSERT(material.mu_ref >= 0.0);
   }
   else if(input == "sutherland")
   {
      material.mu_model = Material::mu_sutherland;
      fin.entry ("mu_ref");
      fin >> material.mu_ref;
      fin.entry ("T_ref");
      fin >> material.T_ref;
      // fin.entry ("T_0");
//       fin >> material.T_0;
      
      MPI_ASSERT(material.mu_ref > 0.0);
      MPI_ASSERT(material.T_ref > 0.0);
      //MPI_ASSERT(material.T_0 > 0.0);
   }
   else if(input == "power")
   {
      material.mu_model = Material::mu_power;
      fin.entry ("mu_ref");
      fin >> material.mu_ref;
      fin.entry ("T_ref");
      fin >> material.T_ref;
      fin.entry ("omega");
      fin >> material.omega;
      
      MPI_ASSERT(material.mu_ref > 0.0);
      MPI_ASSERT(material.T_ref > 0.0);
      MPI_ASSERT(material.omega > 0.0);
   }
   else
      MPI_ERR("read_material: unknown viscosity type "<<input
               << ". Possible options:\n"
			  << " constant \n sutherland \n");
   fin.end_section();

   fin.entry ("prandtl");
   fin >> material.prandtl;
   MPI_ASSERT(material.prandtl > 0.0);

   fin.entry ("model");
   fin >> input;
   if(input == "euler")
      material.model = Material::euler;
   else if(input == "ns")
      material.model = Material::ns;
   else
      MPI_ERR("read_material: unknown flow model " << input
              << ". Possible options:\n"
			  << " euler \n ns \n");


   fin.entry ("flux");
   fin >> input;
   if(input == "kepes_tecno_roe")
      material.flux_scheme = Material::kepes_tecno_roe;  
   else if(input == "kepes_tecno_rusanov")
      material.flux_scheme = Material::kepes_tecno_rusanov;   
   else if(input == "kepec")
      material.flux_scheme = Material::kepec;
   else if(input == "kep")
      material.flux_scheme = Material::kep;
   else if(input == "roe_ec")
      material.flux_scheme = Material::roe_ec; 
   else if(input == "simple_avg")
      material.flux_scheme = Material::avg;         
   else   
      MPI_ERR("read_material: unknown flux scheme: " << input
              << ". Possible options:\n"
			  << " kepes_tecno_roe \n kepes_tecno_rusanov \n" 
			  << " kepec \n kep \n roe_ec \n simple_avg \n");

   
   fin.begin_section("flux_par");
   fin.entry("alpha");
   fin >>  material.ev_mod_alpha;
   MPI_ASSERT(material.ev_mod_alpha >= 0);
   fin.entry("beta");
   fin >> material.ev_mod_beta;
   MPI_ASSERT(material.ev_mod_beta >= 0);
   fin.end_section();
   
   // fin.entry ("flux_correction");
//    fin >> input;    
//    if(input=="yes")
//       material.order_correction= true;
//    else if(input=="no")
//       material.order_correction= false;
//    else
//       MPI_ERR("reader_material:: unknown input " << input<<" for flux_correction. "
//               << "Possible options:\n"
// 			  << " yes \n no \n");   
   
   // To test the balancing effects of physical viscosity, heat flux and numerical diss
   fin.begin_section("balance");
   fin.entry ("switch_active");
   fin >> input;
   if(input == "yes")
   {
      if(get_proc_id() == 0)
         cout <<"  --- Balancing switches active"<<endl;
      
      if(material.model != Material::ns)
         MPI_ERR("read_material: balance switch can be active only with ns model");

      fin.entry ("ND");
      fin >> input;
      if (input == "yes")
		  material.nd_active = true;
	  else if(input == "no")
		  material.nd_active = false;  
	  else
	      MPI_ERR("read_material: unknown option input " << input<<" for ND. "
              << "Possible options:\n"
			  << " yes \n no \n");
	
	  fin.entry ("HEAT");
	  fin >> input;
      if (input == "yes")
		  material.heat_active = true;
	  else if(input == "no")
		  material.heat_active = false;
	  else
	     MPI_ERR("read_material: unknown option input " << input<<" for HEAT. "
              << "Possible options:\n"
			  << " yes \n no \n");
	  
	  fin.entry ("VISC");
	  fin >> input;
      if (input == "yes")
		  material.visc_active = true;
	  else if(input == "no")
		  material.visc_active = false;
	  else
	      MPI_ERR("read_material: unknown option input " << input<<" for VISC. "
              << "Possible options:\n"
			  << " yes \n no \n");
   }
   else if(input == "no")
   {
      // visc and heat switches active is balancing switch is off
      fin >> input;
      fin >> input;
      fin >> input;
      fin >> input;
      fin >> input;
      fin >> input;
      material.visc_active = true;
      material.heat_active = true;
   }
   else
      MPI_ERR("read_material: unknown option input " << input<<" for switch_active. "
              << "Possible options:\n"
			  << " yes \n no \n");
	  
   fin.end_section ();
   
   fin.end_section ();

   material.initialize ();
}

//------------------------------------------------------------------------------
// Read constants section
//------------------------------------------------------------------------------
void Parameter::read_constants (Reader &fin)
{
   MPI_DISP("  --- Reading constants section",verbose);

   string input;
   double value;

   fin.begin_section ("constants");

   while (!fin.eos())
   {
      fin >> input;
      fin >> value;
      if(get_proc_id() == 0)
         cout << setw(16) << input << setw(16) << value << endl;
      constants.insert ( pair<string,double>(input, value) );
   }

}

//------------------------------------------------------------------------------
// Read initial condition
//------------------------------------------------------------------------------
void Parameter::read_initial_condition (Reader &fin)
{
   MPI_DISP("  --- Reading initial condition section",verbose);

   string input;

   fin.begin_section ("initial_condition");

   fin >> input;

   if(input == "library")
   {
      fin >> input;
      initial_condition.add (input);
   }
   else
   {
      MPI_ASSERT(input == "temperature");
      fin.getline (input);
      initial_condition.add ("temperature", input);

      fin.entry ("xvelocity");
      fin.getline (input);
      initial_condition.add ("xvelocity", input);

      fin.entry ("yvelocity");
      fin.getline (input);
      initial_condition.add ("yvelocity", input);

      fin.entry ("zvelocity");
      fin.getline (input);
      initial_condition.add ("zvelocity", input);

      fin.entry ("pressure");
      fin.getline (input);
      initial_condition.add ("pressure", input);
   }

   fin.end_section ();
}

//------------------------------------------------------------------------------
// Read boundary conditions
//------------------------------------------------------------------------------
void Parameter::read_boundary (Reader &fin)
{
   MPI_DISP("  --- Reading boundary section",verbose);

   string input;

   fin.begin_section ("boundary");

   while (!fin.eos())
   {
      vector<int> f_type;
      while(!fin.bos())
      {
         int tmp;
         fin >> tmp;
         f_type.push_back (tmp);
      }

      fin.entry ("type");
      string bc_type;
      fin >> bc_type;

      vector<string> variable, function;
      while(!fin.eos())
      {
         fin >> input;
         variable.push_back (input);
         fin.getline(input);
         function.push_back (input);
      }
      BoundaryCondition bc(material, bc_type, variable, function);
      for(unsigned int i=0; i<f_type.size(); ++i)
         boundary_condition.insert (pair<int,BoundaryCondition>(f_type[i], bc));
   }
}

//------------------------------------------------------------------------------
// Read exact solution
//------------------------------------------------------------------------------
void Parameter::read_exact (Reader &fin)
{
   MPI_DISP("  --- Reading exact solution section",verbose);

   string input;

   fin.begin_section ("exact_soln");
   
   fin.entry("available");   
   fin >> input;
   if(input == "yes")
      exact_available = true;
   else if(input == "no")
      exact_available = false;   
   else   
      MPI_ERR("read_eaxct: unknown option input " << input<<" for available. "
              << "Possible options:\n"
			  << " yes \n no \n");
   
   if(exact_available)
   {   
	  fin.entry("temperature");
	  fin.getline (input);
	  exact_solution.add ("temperature", input);

	  fin.entry ("xvelocity");
	  fin.getline (input);
	  exact_solution.add ("xvelocity", input);
  
	  fin.entry ("yvelocity");
	  fin.getline (input);
	  exact_solution.add ("yvelocity", input);
	  
	  fin.entry ("zvelocity");
	  fin.getline (input);
	  exact_solution.add ("zvelocity", input);

	  fin.entry ("pressure");
	  fin.getline (input);
	  exact_solution.add ("pressure", input);
	  
	  fin.end_section ();
   }
   else
   {
      while(!fin.eos())
        fin >> input;
   }     
}

//------------------------------------------------------------------------------
// Read integrals section
//------------------------------------------------------------------------------
void Parameter::read_integrals (Reader &fin)
{
   MPI_DISP("  --- Reading integrals section",verbose);

   string input;

   fin.begin_section ("integrals");

   while (!fin.eos())
   {
      fin >> input;
      if(input=="force") // integral type is force
      {
         fin.entry ("{");

         force_data.resize( force_data.size() + 1 );
         fin >> force_data.back().name; // name to identify force

         while (!fin.eos())
         {
            int face_type;
            fin >> face_type;
            force_data.back().face_type.push_back(face_type);
         }

         // Check that there was atleast one face type given
         MPI_ASSERT(force_data.back().face_type.size() > 0);
      }
   }

}

//------------------------------------------------------------------------------
// Read output section
//------------------------------------------------------------------------------
void Parameter::read_output (Reader &fin)
{
   MPI_DISP("  --- Reading output section",verbose);

   string input;

   fin.begin_section ("output");

   fin.entry ("format");
   fin >> write_format;
   MPI_ASSERT (write_format == "vtk");

   fin.entry ("frequency");
   fin >> write_frequency;
   MPI_ASSERT (write_frequency > 0);
   
   fin.entry ("time_stamps");
   fin >> n_time_stamps;
   MPI_ASSERT (n_time_stamps >= 0);
   
   fin.entry ("use_online_stat_tool");
   fin >> input;
   if (input == "yes")
	   online_stat = true;
   else if(input == "no")
	   online_stat = false;
   else
      MPI_ERR("reader_output:: unknown option "<<input<<" for use_online_stat_tool. "
              << "Possible options:\n"
			  << " yes \n no \n");
   
   fin.begin_section ("variables");
   while (!fin.eos())
   {
      fin >> input;
      MPI_ASSERT (input=="mach" || input=="density" || input=="entropy" || 
                  input=="vorticity" || "mesh_Peclet");
      write_variables.push_back (input);
      
      if(input=="mesh_Peclet") // Only allowed for Navier-Stokes with non-zero viscosity
      {
         MPI_ASSERT(material.model == Material::ns && material.mu_ref > 0.0);
      }
   }

   fin.begin_section ("surfaces");
   while (!fin.eos())
   {
      int id;
      fin >> id;
      write_surfaces.push_back (id);
   }

//    fin.begin_section ("pdf_points");
//    while (!fin.eos())
//    {
//       int id;
//       fin >> id;
//       write_mc_points.insert (id);
//    }


   fin.entry ("restart");
   fin >> input;
   if(input=="no")
      write_restart = false;
   else if(input=="yes")
      write_restart = true;
   else
      MPI_ERR("reader_output:: unknown option "<<input<<" for restart. "
              << "Possible options:\n"
			  << " yes \n no \n");

   fin.entry ("restart_frequency");
      fin >> restart_write_frequency;
   if(write_restart)   
      MPI_ASSERT(restart_write_frequency>0);   
   
   fin.entry ("log_output");
   fin >> input;
   if(input=="no")
      write_log = false;
   else if(input=="yes")
      write_log = true;
   else
      MPI_ERR("reader_output:: unknown option "<<input<<" for log_output. "
              << "Possible options:\n"
			  << " yes \n no \n");
   
   fin.entry ("soln_output");
   fin >> input;
   if(input=="no")
      write_soln = false;
   else if(input=="yes")
      write_soln = true;
   else
      MPI_ERR("reader_output:: unknown option "<<input<<" for soln_output. "
              << "Possible options:\n"
			  << " yes \n no \n");

   has_global = false;

   fin.entry ("save_global");
   fin >> input;
   if(input=="no")
      has_global= false;
   else if(input=="yes")
   {
      has_global= true;
   }
   else
      MPI_ERR("reader_output:: unknown option "<<input<<" for save_global. "
              << "Possible options:\n"
			  << " yes \n no \n");
      
   // fin.entry ("save_mesh_Pe");
//    fin >> input;
//    if(input=="no")
//       save_mesh_Pe= false;
//    else if(input=="yes")
//    {
//       save_mesh_Pe= true;
//    }
//    else
//       MPI_ERR("reader_output:: unknown option "<<input<<" for save_mesh_Pe. "
//               << "Possible options:\n"
// 			  << " yes \n no \n");  
      
   fin.entry ("find_error");
   fin >> input;
   if(input=="no")
      find_error= false;
   else if(input=="yes")
   {
      find_error= true;
   }
   else
      MPI_ERR("reader_output:: unknown option "<<input<<" for find_error. "
              << "Possible options:\n"
			  << " yes \n no \n");      

   fin.end_section ();
      
   if(time_mode == "unsteady" && n_samples > 1 && n_time_stamps == 0)
      MPI_ERR("reader_output:: time_stamps must be greater than 0 for unsteady MC "
              << " simulations.");   
       
   if(n_samples == 1 && online_stat == true)
      MPI_ERR("reader_output:: online statistics cannot be evaluated with only one sample.");
   
   if(find_error && !exact_available)
   {
      // if(get_proc_id() == 0)
//       {
//          cout<<"      WARNING! Exact solution has not been specified.";
//          cout<<" Turning off error evaluator"<<endl;
//       }   
//       find_error = false;
     MPI_ERR("reader_output:: Error evaluator activated in the absence of exact solution.");
   }
}

//------------------------------------------------------------------------------
// Set MC sample ids to be handled by the particular MC group
//------------------------------------------------------------------------------
void Parameter::set_sample_ids ()
{
   int sample_overflow = n_samples % N_MC_GROUPS;
   int sample_normal   = n_samples/N_MC_GROUPS;
   int group_id        = get_proc_id()/NPART;
   if( group_id < sample_overflow)
   {
      sample_start_ind = group_id*(sample_normal + 1);
      sample_end_ind   = sample_start_ind + sample_normal;
   }
   else
   {
      sample_start_ind = group_id*(sample_normal) + sample_overflow;
      sample_end_ind   = sample_start_ind + sample_normal - 1;
   }

}


//------------------------------------------------------------------------------
// Rearrange ordering of samples to ensure samples are simulated in an ascending
// order by various MC groups
//------------------------------------------------------------------------------
void Parameter::rearrange_samples ()
{
   vector<int> DUMMY_LIST;
   DUMMY_LIST.resize(n_samples);
   
   int ind_base  = 0;
   int next_ind  = ind_base;
   for (int i=0; i < n_samples; ++i)
   {
      DUMMY_LIST[i] = SAMPLE_LIST[next_ind];
      next_ind += N_MC_GROUPS;
      if(next_ind >= n_samples)
      {
         ind_base++;
         next_ind = ind_base;
      }   
   }
   
   SAMPLE_LIST = DUMMY_LIST;
   DUMMY_LIST.clear();  

}
