#include <iostream>
#include <fstream>
#include <string.h>
#include <cstdlib>
#include <sstream>
#include <set>
#include <algorithm>
#include "grid.h"

using namespace std;

//------------------------------------------------------------------------------
// Read a grid in gmsh format
//------------------------------------------------------------------------------
template <int DIM>
void Grid<DIM>::read_gmsh ()
{
   
   // Temperary/dummy variables and structures
   unsigned int ntags, tag, geo_tag, shared, dummy;
   string line;
   Cell<DIM> cell_tmp;
   
   // Split mesh file into different components
   gmsh_split();        
   
   ifstream vfile;
   string vfile_name = part_dir+"/node.dat";
   vfile.open (vfile_name.c_str());
   assert( vfile.is_open() );
   cout<<"\n  --- reading vertices... "<<endl;
   vertex.resize(n_vertex);                  
   vfile >> dummy;
   assert(n_vertex==dummy);
   for(unsigned int i=0; i<n_vertex; ++i)
   {
      vfile >> vertex[i].coord.x >> vertex[i].coord.y  >> vertex[i].coord.z;      
   }
   vfile.close();

   ifstream ffile;
   string ffile_name = part_dir+"/face.dat";
   ffile.open (ffile_name.c_str());
   assert( ffile.is_open() );
   
   cout<<"\n  --- reading faces... "<<endl;
   face.resize(n_face);
   for(unsigned int i=0; i<n_face; ++i)
   {
      ffile >> ntags;
      ffile >> face[i].type; // First tag is face type
      ffile >> geo_tag;
      if(find(periodic_face_tag.begin(),periodic_face_tag.end(),geo_tag)!=periodic_face_tag.end())
         face[i].periodic = true;
      else   
         face[i].periodic = false;
         
      if(ntags > 2)
      {
         // Dummy tags
		 // for(unsigned int p=1; p<ntags-2; ++p)
		 // 		 ffile >> tag; 
         
         // Number of partitions sharing this face
         // This should be 1 for 2D. Verify for 3D
         ffile >> tag;
         assert(tag == 1);
         
         // Processor tag
         ffile >> tag;
         face[i].part = tag - 1;
      }
      else
         face[i].part = 0;    
		 
           
         
      map<int,int>::iterator it = part_nface.find(face[i].part);
      if(it!= part_nface.end())
		 it->second++;
	  else
		 part_nface.insert( pair<int,int>(face[i].part,1));   
         
	  // Face data
	  for(unsigned int j=0; j<DIM; ++j)
		 ffile >> face[i].vertex[j];
         
	  // Node number must start from 0
	  // but gmsh starts from 1.
	  for(unsigned int j=0; j<DIM; ++j)
		  face[i].vertex[j]--;  
             
	  // Mark nodes as domain boundary nodes    
	  for(unsigned int j=0; j<DIM; ++j)
		  vertex[face[i].vertex[j]].d_bnd = true;           
   }
   ffile.close();
   
   ifstream cfile;
   string cfile_name = part_dir+"/cell.dat";
   cfile.open (cfile_name.c_str());
   assert( cfile.is_open() );
   cout<<"  --- reading and partitioning cells... "<<endl;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      cfile >> ntags;   
	  cfile >> tag 
		    >> tag; // Dummy tags
              
	  // Partition number
	  if(ntags > 2)
		 cfile >> shared;
      else
         shared = 1;
         		 
	  bool g_cell[shared];
	  int part_tag[shared];
	  for(unsigned int p=0; p<shared; ++p)   
	  {
		  if(ntags > 2)
		     cfile >> part_tag[p]; 
		  else
		     part_tag[p] = 1;    
		  partitions.insert(abs(part_tag[p])-1);
		  if(part_tag[p] < 0)
			 g_cell[p] = true;
		  else
			 g_cell[p] = false; 
		  
		  map<int,int>::iterator it = part_ncell.find(abs(part_tag[p])-1);
		  if(it!= part_ncell.end())
			 it->second++;
		  else
			 part_ncell.insert( pair<int,int>(abs(part_tag[p])-1,1));     
	  }    

	  // Node number must start from 0
	  // but gmsh starts from 1.
	  for(unsigned int j=0; j<DIM+1; ++j)
	  {
		  cfile >> cell_tmp.vertex[j]; 
		  cell_tmp.vertex[j]--; 
	  }    
				   
		   
	  for(unsigned int p=0; p<shared; ++p)   
	  {
		  write_cell(cell_tmp,abs(part_tag[p])-1,g_cell[p]);   
		  for(unsigned int q = 0; q < DIM+1; ++q)
		  {
			  int globv = cell_tmp.vertex[q]; 
			  
			  if(part_tag[p] > 0)
					vertex[globv].active_shared.insert(part_tag[p]-1);
			  
			  vertex[globv].all_shared.insert(abs(part_tag[p])-1);
						
		  }
		  
	  }	    
   }
   cfile.close ();
   assert(partitions.size() == NPART);
   
   for(unsigned int i=0;i<n_vertex;++i)
   {
      for(set<int>::iterator it=vertex[i].all_shared.begin(); 
	                 it!=vertex[i].all_shared.end();
	                 ++it)
	  {               
		 map<int,int>::iterator itt = part_nvertex.find(*it);
		 if(itt!= part_nvertex.end())
			itt->second++;
		 else
		    part_nvertex.insert( pair<int,int>(*it,1));  
	  }			
   }
   
   
   cout<<"  --- creating periodic listings (if available)... "<<endl;
   vector<vector<int> >::iterator itvv1;
   vector<int>::iterator itv1;
   periodic_set_part_list.resize(periodic_node_sets.size());
   for(int i=0;i<NPART;++i)
      part_periodic_set.insert(pair<int,int>(i,0));
   
   
   for(int i=0; i<periodic_node_sets.size();++i)
   {
      for(itv1 = periodic_node_sets[i].begin();itv1 != periodic_node_sets[i].end();++itv1)
      {
         int v = *itv1;
         for(set<int>::iterator it=vertex[v].active_shared.begin(); 
	                 it!=vertex[v].active_shared.end();
	                 ++it)
		 {               
			periodic_set_part_list[i].insert(*it);
		 }
      }
      for(int j=0;j<NPART;++j)
         if(periodic_set_part_list[i].count(j)!=0)
            part_periodic_set.find(j)->second++;
   }      
}

template class Grid<2>;
template class Grid<3>;