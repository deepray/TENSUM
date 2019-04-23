#include <string>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include "grid.h"

using namespace std;

// ------------------------------------------------------------------------------
// Intialize main partition files
// ------------------------------------------------------------------------------
// template <int DIM>
// void Grid<DIM>::initialize_partition_files ()
// {
//     unsigned int n_part = partitions.size();
//     
//     for(unsigned int i=0; i<n_part; ++i)
//     {
// 	   string file = part_dir+"/part_file_";
// 	   stringstream ss;
// 	   ss << i;
// 	   file += ss.str()+".dat";
// 	   ofstream p_file;
// 	   p_file.open(file.c_str());
// 	   p_file << n_part <<endl;
// 	   p_file << part_nvertex.find(i)->second<<" OF "<<n_vertex<<endl;
// 	   p_file.close();
//     }
// }



// ------------------------------------------------------------------------------
// Creates main partition files
// ------------------------------------------------------------------------------
template <int DIM>
void Grid<DIM>::create_partition_files ()
{
    cout<<"  --- writing final partition files into "<<part_dir<<"... "<<endl;
    const int n_part = partitions.size();
	vector<ofstream*> p_file (n_part);
	for(int i=0; i<n_part; ++i)
	{
	   stringstream ss;
	   ss << i;
	   string filename = part_dir+"/part_file_"+ ss.str() + ".dat";
	   p_file[i] = new ofstream(filename.c_str());
	   assert(p_file[i]->is_open());
	   *p_file[i] << "DIMENSION" << " " << DIM <<endl;
	   *p_file[i] << "PARTITIONS" << " " << n_part <<endl;
	   *p_file[i] << "NODES "<<part_nvertex.find(i)->second<<" OF "<<n_vertex<<endl;
	}
    
    
    // Writing vertices
    
    // FORMAT FOR VERTEX DATA:
    //------------------------------------------------------------------------------------
    // GLOB_ID, COORD_X, COORD_Y, COORD_Z, NUM_ACTIVE_SHARED_PART, ACTIVE_SHARED_PART_IDS,
    // ON_DOMAIN_BOUNDARY
    //------------------------------------------------------------------------------------
    // If NUM_ACTIVE_SHARED_PART == 0, then NUM_ACTIVE_SHARED_PART is set to 1, followed
    // by the ID of the partition in which it is a ghost.
    // If the vertex is on the domain boundary, then ON_DOMAIN_BOUNDARY == 1, else 0
    for(unsigned int i=0; i<n_vertex; ++i)
    {
	   stringstream active_list ;
	   for(set<int>::iterator it=vertex[i].active_shared.begin(); 
	                 it!=vertex[i].active_shared.end();
	                 ++it)
	      active_list<<*it<<" ";
	      
	   for(set<int>::iterator it=vertex[i].all_shared.begin(); 
	                 it!=vertex[i].all_shared.end();
	                 ++it)
	   {              
	      *p_file[*it] << i << " "
	                  << setprecision(16) << vertex[i].coord.x << " "
	                  << setprecision(16) << vertex[i].coord.y << " "
	                  << setprecision(16) << vertex[i].coord.z << " ";
	      if(!vertex[i].active_shared.empty())      
	      { 
	         if(vertex[i].active_shared.count(*it)!=0)
	            *p_file[*it] << vertex[i].active_shared.size()<< " "<< active_list.str() << " ";  
	         else
	            *p_file[*it] << 1 << " " << *it << " ";    
	      }
	      else
	      {
	         *p_file[*it] << 1 << " " << *it << " "; 
	      }    

	      *p_file[*it] << vertex[i].d_bnd <<endl;                
	   }  
    }
    
    // FORMAT FOR VERTEX PERIODIC DATA:
    //------------------------------------------------------
    // ID of P_i, #(P_ij), P_ij, #(L_i), L_i 
    //------------------------------------------------------
    // ID of P_i -> ranges from 0 to periodic_node_sets.size()-1
    // P_ij -> list of periodic points in partition j
    // L_i -> list of partitions associated with with points in P_i
    // NOTE: These lists only include active nodes
    for(unsigned int i=0; i<n_part; ++i)
    {
       *p_file[i] << "PERIODIC "<<part_periodic_set.find(i)->second <<endl;
       for(unsigned int j=0;j<periodic_node_sets.size();++j)
	   {
		   int cp = periodic_set_part_list[j].count(i);
		   if(cp!=0)
		   {
		      *p_file[i] << j << " ";
		      int count = 0;
		      stringstream P_ij ;
		      for(vector<int>::iterator itv1 = periodic_node_sets[j].begin();
		          itv1 != periodic_node_sets[j].end(); ++itv1)
		      {    
				 if(vertex[*itv1].active_shared.count(i)!=0)
				 {
                    P_ij<<*itv1<<" ";		   
                    count++;
                 }   
		      }   
		      *p_file[i] << count <<" "<< P_ij.str()<<" ";
		      *p_file[i] << periodic_set_part_list[j].size()<<" ";
		      for(set<int>::iterator it2=periodic_set_part_list[j].begin();
		          it2!=periodic_set_part_list[j].end();++it2)
		         *p_file[i]<<*it2<<" ";
		      *p_file[i]<<endl;  
		   }    
	   }   
	}      
	
    vertex.clear();
           
    
    // FORMAT FOR FACE DATA:
    //------------------------------------------------------------------------------------
    // FACE_TYPE, VERTEX_LIST (size DIM), PERIODIC
    //------------------------------------------------------------------------------------
    // If the face does not have a periodic partner, then PERIODIC == 0, else
    // PERIODIC == 1.
    

    for(unsigned int i=0; i<n_part; ++i)
    {
       //Writing faces
	   if(part_nface.find(i)!=part_nface.end())
	   {
		   *p_file[i] << "FACES "<<part_nface.find(i)->second <<" OF "<<n_face<<endl;
	   }
	   else
	       *p_file[i] << "FACES "<< 0 <<" OF "<<n_face<<endl;
	}
	
	for(unsigned int i=0; i<n_face; ++i)
	{  
	   int part = face[i].part;
	   *p_file[part] << face[i].type << " ";
	   
	   for(unsigned int j=0; j<DIM; ++j)
	      *p_file[part] << face[i].vertex[j] << " ";   
	      
	   *p_file[part] << face[i].periodic << endl;       
    }
    
    
    //Writing cells
    for(unsigned int i=0; i<n_part; ++i)
    {
	   
	   *p_file[i] << "CELLS "<<part_ncell.find(i)->second <<" OF "<<n_cell<<endl;
	   string c_file_name = part_dir+"/cell_file_";
	   stringstream ss;
	   ss << i;
	   c_file_name+=ss.str();
	   ifstream c_file;
	   c_file.open(c_file_name.c_str());
	   assert(c_file.is_open() );
	   bool g_cell;
	   for(int j=0; j<part_ncell.find(i)->second; ++j)
	   {
	       Cell<DIM> c;
	       for(int k=0; k<DIM+1; ++k)
		   {
			  c_file >> c.vertex[k];
			  *p_file[i] << c.vertex[k] << " ";
		   }
		   c_file >> g_cell;
	       *p_file[i] << g_cell<<endl;      
	   }
	   c_file.close();
	   string commandline = "rm " + c_file_name;
	   system(commandline.c_str());       
    }  
    
    for(int i=0; i<n_part; ++i)
    {
       p_file[i]->close();
    }
    
}

template class Grid<2>;
template class Grid<3>;

