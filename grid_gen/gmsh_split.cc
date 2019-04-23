#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include "grid.h"

using namespace std;

//------------------------------------------------------------------------------
// Read grid from file
//------------------------------------------------------------------------------
template <int DIM>
void Grid<DIM>::gmsh_split ()
{
   unsigned int count, elem_type, ntags, n_elem, dummy;
   string line;
   cout << "\n  Reading and splitting gmsh grid file " << msh_file << endl;

   ifstream file;
   file.open (msh_file.c_str());
   assert( file.is_open() );
   
   file >> line;
   while (line !="$Nodes")
      file >> line;
   
   ofstream vfile;
   string vfile_name = part_dir+"/node.dat";
   cout << "\n  --- creating node file " << vfile_name.c_str() << endl;
   vfile.open(vfile_name.c_str());
   assert(vfile.is_open() );

   file >> n_vertex;
   assert(n_vertex > 0);
   vfile << n_vertex <<endl;
   
   Vertex dummy_vertex;
   for(unsigned int i=0; i<n_vertex; ++i)
   {
	  file >> count 
		   >> dummy_vertex.coord.x
		   >> dummy_vertex.coord.y
		   >> dummy_vertex.coord.z;
	  vfile << setprecision(15) << dummy_vertex.coord.x << " "
	        << setprecision(15) << dummy_vertex.coord.y << " "
	        << setprecision(15) << dummy_vertex.coord.z << endl;	   
   }
   vfile.close();		   
   
   file >> line;
   file >> line;
   
   ofstream ffile,cfile;
   string ffile_name = part_dir+"/face.dat";
   string cfile_name = part_dir+"/cell.dat";
   cout << "\n  --- creating face file " << ffile_name.c_str() 
        << " and cell file "<< cfile_name.c_str() << endl;
   ffile.open(ffile_name.c_str());
   assert(ffile.is_open() );
   cfile.open(cfile_name.c_str());
   assert(cfile.is_open() );

   // Read cells
   file >> n_elem;
   assert(n_elem > 0);

   for(unsigned int i=0; i<n_elem; ++i)
   {
      file >> count
           >> elem_type
           >> ntags;
      int tags[ntags];     
      if(elem_type == 1 && DIM == 3) // Line element in 3D is not required
      {
         for(unsigned int p=0; p<ntags; ++p)
            file >> dummy;
         file >> dummy >> dummy;     // Nodes of line element
      }
      else if(elem_type == geo_type[0]) // Face
      {
         n_face++;
         for(unsigned int p=1; p<=ntags; ++p)
            file >> tags[p]; 
         int vertex_id[DIM];   
         for(unsigned int j=0; j<DIM; ++j)
            file >> vertex_id[j];
         
         ffile << ntags << " ";
         for(unsigned int p=1; p<=ntags; ++p)
            ffile << tags[p] <<" ";   
         for(unsigned int j=0; j<DIM; ++j)
            ffile << vertex_id[j] << " ";
         ffile <<endl;          
                  
            
      }
      else if(elem_type == geo_type[1]) // Cell
      {  
         n_cell++;
         
         for(unsigned int p=1; p<=ntags; ++p)
            file >> tags[p]; 
              
         int vertex_id[DIM+1];   
         for(unsigned int j=0; j<DIM+1; ++j)
            file >> vertex_id[j];
         
         cfile << ntags << " ";
         for(unsigned int p=1; p<=ntags; ++p)
            cfile << tags[p] <<" ";   
         for(unsigned int j=0; j<DIM+1; ++j)
            cfile << vertex_id[j] << " ";
         cfile <<endl;      
      }
      else
      {
         cout<<"Unknown element type " << elem_type<<endl;
         exit(0);
      }   
   }
   
   file >> line;
   
   if(!file.eof())
   {
      file >> line;
      if(line == "$Periodic")
      {
         periodic_available = true;
         cout << "\n  --- periodic data available --> creating connections."<<endl;
         int p_elem, e_dim, n_count, master, slave;
         file >> p_elem;
         assert(p_elem > 0);
         map<int, vector<int> >::iterator itm1;
         vector<vector<int> >::iterator itvv1;
	     vector<int>::iterator itv1;
         for(unsigned int p=0; p<p_elem; ++p)
         {
            file >> e_dim >> master >> slave;
            if(DIM == 3)
            {
               getline(file,line); // Dummy Affine line for periodic dat in 3D
            }
            if(e_dim == 0 || (e_dim == 1 && DIM == 3))
            {
               file >> dummy >> master >> slave;
               //In C++ node indexing starts from 0
               master--;
               slave--;
   
			   vector<int> inner;
			   inner.push_back(master);
			   inner.push_back(slave);
			   periodic_node_sets.push_back(inner);
	  
	  
            }
            else if(e_dim == geo_type[0])  
            {
               periodic_face_tag.push_back(master);
               periodic_face_tag.push_back(slave);
               file >> n_count;
               for(int n=0; n<n_count; ++n)
               {
                  file >> master >> slave;
                  //In C++ node indexing starts from 0
				  master--;
				  slave--;
				  vector<int> inner;
			      inner.push_back(master);
			      inner.push_back(slave);
			      periodic_node_sets.push_back(inner);
                  
               }   
            }
         } 
         
         complete_connections();
         
         
      }
   }
   
   // Closing file
   ffile.close();
   cfile.close();
   file.close();

}

template class Grid<2>;
template class Grid<3>;