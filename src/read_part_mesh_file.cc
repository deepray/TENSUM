#include <iostream>
#include <fstream>
#include <string.h>
#include <cstdlib>
#include <sstream>
#include <set>
#include <algorithm>
#include "grid.h"

using namespace std;
extern bool verbose;
extern bool PERIODIC;

//------------------------------------------------------------------------------
// Read a grid in gmsh format
//------------------------------------------------------------------------------
void Grid::read_part_mesh_file(string PART_DIR,const map<int,BoundaryCondition>& bc)
{
   
   // Temp/dummy variables and structures
   unsigned int npart,tag,int_dummy,part_no;
   string line;
   vector< set<int> > vertex_shared_list; // For each vertex, stores set of processors in 
                                          // which it is active
   vector< VERTEX_MPI_NBR > nbr_list;     // For each vertex in a partition,
                                          // stores number of adjoining ghost cells and 
                                          // active cells. 
   
   int rank = get_proc_id();
   //int loc_rank = get_proc_loc_id();
   
   if(rank == 0)
       cout << "\n  Reading partitioned grid files" << endl;

   ifstream file;
   string file_name = PART_DIR+"/part_file_";
   stringstream ss;
   ss<<get_proc_loc_id();
   file_name+=ss.str()+".dat";
   file.open (file_name.c_str());
   MPI_LOC_ASSERT (file.is_open() );

   file >> line;
   MPI_LOC_ASSERT(line=="DIMENSION")
   file >> int_dummy;
   if(int_dummy!=2)
      MPI_LOC_ERR("ERROR: Dimension of partitioned file is "<<int_dummy<<". The dimension must be 2!!")
   
   file >> line;
   MPI_LOC_ASSERT(line=="PARTITIONS")
   file >> npart;
   if(npart != NPART)
      MPI_LOC_ERR("ERROR: Number of mesh partitions in parameter file do not match the number of partitioned file!!")

   file >> line;
   MPI_LOC_ASSERT(line=="NODES")
   MPI_DISP("  --- reading vertices... ",verbose);
   // Read vertices
   file >> n_vertex;              
   MPI_LOC_ASSERT (n_vertex> 0);
   vertex.resize (n_vertex);
   vertex_shared_list.resize(n_vertex);
   file >> line >> n_vertex_total;
   MPI_LOC_ASSERT (n_vertex_total> 0);
   MPI_LOC_ASSERT (n_vertex <= n_vertex_total);

   for(unsigned int i=0; i<n_vertex; ++i)
   {
	      unsigned int shared_list_size;
	      int part, globv;
	      file >> globv
           >> vertex[i].coord.x
           >> vertex[i].coord.y
           >> vertex[i].coord.z;
           vertex[i].active = false;
           vertex[i].mpi_bnd_active = false;    
           g2l.insert ( pair<int,int>(globv,i) );
           l2g.insert ( pair<int,int>(i,globv) );
           
           file >> shared_list_size;
           for(unsigned int j=0;j<shared_list_size;++j)
           {
              file >> part;
              if(part >= NPART)
              {
                 cout<<"Error: Partition ID = "<<part<<" detected in mesh files.";
                 cout<<" Number of mesh partition mentioned in parameter file = "<<NPART<<endl;
                 MPI_Abort(MPI_COMM_WORLD,1);
              }
              vertex_shared_list[i].insert(get_proc_glob_id(part));
           }   
           file >> vertex[i].d_bnd;
   }  

   // int mc_pt[mc_pt_list.size()];
//    for(unsigned int i=0; i<mc_pt_list.size(); ++i)
//       mc_pt[i] = 0;
//    int ind = 0;   
//    for(set<int>::iterator it=mc_pt_list.begin();it!=mc_pt_list.end();++it )
//    {
//       if(g2l.count(*it-1) == 1)
//          mc_pt[ind] = 1;   
//       ind++;   
//    } 
//    int buff[mc_pt_list.size()];
//    MPI_Allreduce(&mc_pt,&buff,mc_pt_list.size(),MPI_INT,MPI_MAX,MPI_COMM_WORLD);

   // Checking if point of pdf evaluation is valid 
   // set<int>::iterator itp=mc_pt_list.begin();
//    bool okay = true;
//    stringstream ss2;
//    for(unsigned int i=0; i<mc_pt_list.size(); ++i)
//    {
// 	  if(buff[i] == 0)
// 	  {
// 		  ss2 << *itp <<"; ";
// 		  okay = false;
// 	  }
// 	  itp++;
//    }
//    if(!okay)      
// 	  MPI_LOC_ERR("Error: The following points for pdf evaluation do not exist in the mesh:\n"<<ss2.str());

   // Read periodic data if available
   file >> line;
   MPI_LOC_ASSERT(line=="PERIODIC")
   file >> n_periodic;
   MPI_LOC_ASSERT(n_periodic >= 0);
   if(PERIODIC)
   {
      periodic_lists.resize(n_periodic);
   }   
   for(unsigned int i=0; i<n_periodic; ++i)
   {
      file>>int_dummy;
      
      if(PERIODIC)
         periodic_lists[i].glob_id = int_dummy;
         
      file >> int_dummy;
      int v;
      if(PERIODIC)
         periodic_lists[i].P_ij.resize(int_dummy);
      for(unsigned int j=0; j<int_dummy; ++j)
      {
         file>>v;
         if(PERIODIC)
            periodic_lists[i].P_ij[j] = g2l.find(v)->second;;
      }   
      
      file >> int_dummy;
      for(unsigned int j=0; j<int_dummy; ++j)
      {
         file>>part_no;  
         if(PERIODIC) 
            periodic_lists[i].L_i.insert(get_proc_glob_id(part_no));
      }    
   }
   
   // Read faces
   MPI_DISP("  --- reading faces... ",verbose); 
   file >> line;
   MPI_LOC_ASSERT(line=="FACES")  
   file >> n_face; 
   MPI_LOC_ASSERT (n_face >= 0);
   if(n_face > 0)
      face.resize(n_face);
   file >> line >> n_face_total;
   MPI_LOC_ASSERT (n_face_total > 0);
   MPI_LOC_ASSERT (n_face_total >= n_face);
   for(unsigned int i=0; i<n_face; ++i)
   {
      unsigned int globv[2];
      file >> face[i].type
           >> globv[0]
           >> globv[1]
           >> face[i].periodic;          
      if((bc.find(face[i].type)->second).type!=BC::periodic)
         face[i].periodic = false;    
      else
         MPI_LOC_ASSERT(face[i].periodic);   
      // changing vertex ids to local ids
      face[i].vertex[0] = g2l.find(globv[0])->second;
      face[i].vertex[1] = g2l.find(globv[1])->second;
   }
   
   // Read cells
   MPI_DISP("  --- reading cells... ",verbose);
   file >> line;
   MPI_LOC_ASSERT(line=="CELLS") 
   file >> n_cell;
   MPI_LOC_ASSERT (n_cell > 0);
   cell.resize(n_cell);
   file >> line >> n_cell_total;
   MPI_LOC_ASSERT (n_cell_total > 0);
   MPI_LOC_ASSERT (n_cell_total >= n_cell);
   for(unsigned int i=0; i<n_cell; ++i)
   {
      unsigned int globv[3];
      file >> globv[0]
           >> globv[1]
           >> globv[2]
           >> cell[i].g_cell;
      // changing vertex ids to local ids
      cell[i].vertex[0] = g2l.find(globv[0])->second;
      cell[i].vertex[1] = g2l.find(globv[1])->second;
      cell[i].vertex[2] = g2l.find(globv[2])->second;
      
      if(!cell[i].g_cell) 
        n_cell_active++;
   }  
   
   file.close ();
   
   nbr_list.resize(vertex.size());
   //cout<<"check3:"<<rank<<endl;
   // CREATING NBR_LIST FOR VERTICES AND LOCAL LIST OF NON-GHOST CELLS
   MPI_DISP("  --- localising cell vertex ids and creating ghost cell list... ",verbose);
   int ind = 0;
   for(unsigned int j = 0; j < n_cell; ++j)
   {
	   for(int q=0;q<3;q++)
	   {
		   int v = cell[j].vertex[q];
		   if(!cell[j].g_cell)
		   {
			  if(!vertex[v].active)
			  {   
				 vertex[v].active = true;
				 n_vertex_active++; 
			  }
			  nbr_list[v].ngc_nbr++;
		   }
		   else
			  nbr_list[v].gc_nbr++;
		   
		   if(nbr_list[v].ngc_nbr > 0 && nbr_list[v].gc_nbr > 0)
		   {	  
		       if(!vertex[v].mpi_bnd_active)
		       {
		          vertex[v].mpi_bnd_active = true; 
		          ++n_vertex_mpi_bnd;
		       } 
		   }    
	   }
	   if(!cell[j].g_cell)
	      active_cell_ind.insert(pair<int,int>(j,ind++) );
	   
   }
   nbr_list.clear();
   
   // MAKING LIST OF MPI BOUNDARY CELLS AND INTERIOR CELLS
   // ALSO MAKING LIST OF GHOST CELLS WITH FACE(S) ON DOMAIN BOUNDARY
   MPI_DISP("  --- making MPI cell lists... ",verbose);
   for(unsigned int j = 0; j < n_cell; ++j)
   {
       cell[j].mpi_bnd_active = false;
	   if(!cell[j].g_cell)
	   {
	      for(int q=0;q<3;q++)
	      {
		      int v = cell[j].vertex[q];
			  if(vertex[v].mpi_bnd_active && !cell[j].mpi_bnd_active)
			  {
			     cell[j].mpi_bnd_active = true;
			     cell_mpi_bnd.push_back(j);
			  }  
		  }	  
		  
		  if(!cell[j].mpi_bnd_active)
		      cell_mpi_int.push_back(j);    
	   }
	   else
	   {
	      int v0 = cell[j].vertex[0];
	      int v1 = cell[j].vertex[1];
	      int v2 = cell[j].vertex[2];
	      
	      int bface_count = 0;
	      BOUNDARY_GCELL dummy_bgcell;
	      
	      if(vertex[v0].d_bnd && vertex[v1].d_bnd)
	      {
	         dummy_bgcell.bface_vertices.push_back(std::make_pair(v0,v1));
	         bface_count++;
	      }
	      if(vertex[v1].d_bnd && vertex[v2].d_bnd)
	      {
	         dummy_bgcell.bface_vertices.push_back(std::make_pair(v1,v2));
	         bface_count++;
	      }
	      if(vertex[v2].d_bnd && vertex[v0].d_bnd)
	      {
	         dummy_bgcell.bface_vertices.push_back(std::make_pair(v2,v0));
	         bface_count++;
	      }
	      
	      if(bface_count > 0)
	      {
	         dummy_bgcell.cell_id = j;
	         b_gcell.push_back(dummy_bgcell);
	      }
	   }
   }
   n_cell_mpi_bnd = cell_mpi_bnd.size();
   n_cell_mpi_int = cell_mpi_int.size();
   
   // CREATING INDEX MAPPING FOR ACTIVE VERTICES
   MPI_DISP("  --- listing active vertices... ",verbose);
   ind = 0;
   for(unsigned int j=0; j<n_vertex; j++)
   {
       if(vertex[j].active)
       {
          active_node_ind.insert(pair<int,int>(j,ind++) );
          int nv = vertex_shared_list[j].size();
          mpi_node_share.insert(pair<int,int>(j,nv) );   
       }
   }
   
   //cout<<"check4:"<<rank<<endl;
   // CHANGING MC NODE ID TO LOCAL ID AND MAKING L2G MAP
   // MPI_DISP("  --- creating local-global-local vertex maps for shared vertices... ",verbose);
//    for(set<int>::iterator it=mc_pt_list.begin(); it!=mc_pt_list.end(); it++)
//    {
//        int globv = *it - 1;
//        map<int,int>::iterator it1 = g2l.find(globv);
//        // Assigning mc pt data to the first processor in vertex shared shared list for the node
//        if(it1!=g2l.end())
//        {
//           int locv = it1->second;
//           set<int>::iterator it2 = vertex_shared_list[locv].begin();
//           if(*it2 == get_proc_id())
//              l2g_mc_pt.insert(pair<int,int>(locv,globv)); 
//        }   
//           
//    } 
   
   //REMOVING ALL BUT ACTIVE BOUNDARY VERTICES FROM LOCAL GLOBAL MAPPINGS
   map<int,int>::iterator it = g2l.begin();
   while(it!=g2l.end())
   {
	   int locv = it->second;
	   if(!vertex[locv].mpi_bnd_active)
	   {
		  g2l.erase(it++);
		  l2g.erase(l2g.find(locv));
	   }   
	   else
		  it++;   
   }
   
//    for(int i=0;i<get_comm_size();++i)
//    {
//       cout<<"PROC "<<i<<endl;
//       for(it=l2g.begin();it!=l2g.end();it++)
//          cout<<it->second<<" ";
//       cout<<endl;        
//       fflush(stdout);
//    }
//    
   //cout<<"check5:"<<rank<<endl;
   //MAKING PROCESSOR GROUPINGS (for inter-partition sharing of vertices)
   MPI_DISP("  --- making processor groupings... ",verbose);
   vector< SHARED_LIST > Local_vertex_shared_list;
   for(it=l2g.begin();it!=l2g.end();it++)
   {
	   int locv = it->first;
	   SHARED_LIST dummy;
	   dummy.vertex = locv;
	   dummy.proc_list = vertex_shared_list[locv];
	   Local_vertex_shared_list.push_back(dummy);
   }
   vertex_shared_list.clear();
   while(Local_vertex_shared_list.size()!=0)
   {
	  INTP_VERTEX_GROUPING_LIST dummy;
	  dummy.proc_list = Local_vertex_shared_list[0].proc_list;
	  int locv = Local_vertex_shared_list[0].vertex;
	  int globv = l2g.find(locv)->second;
	  dummy.vertices.push_back(globv);
	  Local_vertex_shared_list.erase(Local_vertex_shared_list.begin());
	  vector< SHARED_LIST >::iterator its = Local_vertex_shared_list.begin();
	  while(its!=Local_vertex_shared_list.end())
	  {
		  if((*its).proc_list == dummy.proc_list)
		  {
			 locv = (*its).vertex;
			 globv = l2g.find(locv)->second;
			 dummy.vertices.push_back(globv);
			 Local_vertex_shared_list.erase(its);
		  }
		  else
		  { 
			 its++;
		  }   
	  }
	  sort(dummy.vertices.begin(), dummy.vertices.end());
	  dummy.tag = -1;
	  intp_mpi_groupings.push_back(dummy);
	  if(rank < NPART)
	     all_intp_groupings.insert(dummy.proc_list);
   }
   // PERIODIC LISTINGS
   if(PERIODIC)
   {
	  vector< SHARED_LIST > Local_vertex_periodic_list;
	  for(unsigned int i=0;i<n_periodic;++i)
	  {
		  SHARED_LIST dummy;
		  dummy.vertex = periodic_lists[i].P_ij[0];
		  dummy.proc_list = periodic_lists[i].L_i;
		  if(dummy.proc_list.size() > 1)
		     Local_vertex_periodic_list.push_back(dummy);
	  }
	  while(Local_vertex_periodic_list.size()!=0)
	  {
		 INTP_VERTEX_GROUPING_LIST dummy;
		 dummy.proc_list = Local_vertex_periodic_list[0].proc_list;
		 dummy.vertices.push_back(Local_vertex_periodic_list[0].vertex);
		 Local_vertex_periodic_list.erase(Local_vertex_periodic_list.begin());
		 vector< SHARED_LIST >::iterator its = Local_vertex_periodic_list.begin();
		 while(its!=Local_vertex_periodic_list.end())
		 {
			 if((*its).proc_list == dummy.proc_list)
			 {
				dummy.vertices.push_back((*its).vertex);
				Local_vertex_periodic_list.erase(its);
			 }
			 else
			 { 
				its++;
			 }   
		 }
		 dummy.tag = -1;
		 periodic_mpi_groupings.push_back(dummy);
		 if(rank < NPART)
			all_periodic_groupings.insert(dummy.proc_list);
	  }
   }
   
   mpi_barrier(MPI_COMM_WORLD);
   MPI_DISP("  --- sharing processor list... ",verbose);
   // Sharing all groupings (inter partition)
   if(rank == 0)
   {
      for(int i = 1; i<NPART; i++)
      {
          get_all_groupings(all_intp_groupings); 
      }      
      for(int i = 1; i<NPART; i++)
      {
          tag = i;
          send_all_groupings(0,i,all_intp_groupings,tag);
      }    
   }
   else if(rank < NPART)
   {
      tag = rank;
      send_all_groupings(rank,0,all_intp_groupings,tag);
      all_intp_groupings.clear();
      get_all_groupings(all_intp_groupings);
   }   

   mpi_barrier(MPI_COMM_WORLD);
   
   // Sharing all groupings (periodic)
   if(PERIODIC)
   {
	  if(rank == 0)
	  {
		 for(int i = 1; i<NPART; i++)
		 {
			 get_all_groupings(all_periodic_groupings); 
		 }     
		 for(int i = 1; i<NPART; i++)
		 {
			 tag = i;
			 send_all_groupings(0,i,all_periodic_groupings,tag);  
		 }    
	  }
	  else if(rank < NPART)
	  {
		 tag = rank;
		 send_all_groupings(rank,0,all_periodic_groupings,tag);
		 all_periodic_groupings.clear();
		 get_all_groupings(all_periodic_groupings);
	  } 
	  mpi_barrier(MPI_COMM_WORLD);
   }
   
   
   MPI_DISP("  --- creating processor tags... ",verbose);
   
   // Setting group tags
   tag = 0;
   int total_share_tags = all_intp_groupings.size()*N_MC_GROUPS;
   int NPROC = get_comm_size();
   if(rank < NPART)
   {
      for(set<set<int> >::iterator its=all_intp_groupings.begin();
          its!=all_intp_groupings.end();++its) 
      {
         set<int> list = *its;
         for(unsigned int i=0; i<intp_mpi_groupings.size(); ++i)
            if(intp_mpi_groupings[i].proc_list == list)
               intp_mpi_groupings[i].tag = tag;
         tag++;
      }   
      int n_groups = all_intp_groupings.size();
      all_intp_groupings.clear();
      
      if(NPROC >= NPART)
      {
		 int buf[intp_mpi_groupings.size()];
		 for(unsigned int i=0; i<intp_mpi_groupings.size(); ++i)
			buf[i] = intp_mpi_groupings[i].tag;
		 
		 int tag_inc = n_groups;
		 for(unsigned int i=1; i<N_MC_GROUPS; ++i)  
		 {
			MPI_Send(&tag_inc, 1, MPI_INT, rank+NPART*i,rank,MPI_COMM_WORLD);
			MPI_Send(&buf, intp_mpi_groupings.size(), MPI_INT, rank+NPART*i,NPROC+rank,MPI_COMM_WORLD); 
			tag_inc+=n_groups;
	     }		
      }   
   }
   else
   {
      int source = rank%NPART;
      int tag = source;
      int tag_inc;
      int buf[intp_mpi_groupings.size()];
      MPI_Recv(&tag_inc, 1, MPI_INT, source,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv(&buf, intp_mpi_groupings.size(), MPI_INT, source,NPROC+tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
     
      for(unsigned int i=0; i<intp_mpi_groupings.size(); ++i)
      {
         intp_mpi_groupings[i].tag = buf[i]+tag_inc;  
      }   
   }
    
   mpi_barrier(MPI_COMM_WORLD);
   
   if(PERIODIC)
   {
      tag = total_share_tags;
      if(rank < NPART)
	  {
		 for(set<set<int> >::iterator its=all_periodic_groupings.begin();
			 its!=all_periodic_groupings.end();++its) 
		 {
			set<int> list = *its;
			for(unsigned int i=0; i<periodic_mpi_groupings.size(); ++i)
			   if(periodic_mpi_groupings[i].proc_list == list)
				  periodic_mpi_groupings[i].tag = tag;
			tag++;
		 }   
		 int n_groups = all_periodic_groupings.size();
		 all_periodic_groupings.clear();
	  
		 if(NPROC >= NPART)
		 {
			int buf[periodic_mpi_groupings.size()];
			for(unsigned int i=0; i<periodic_mpi_groupings.size(); ++i)
			   buf[i] = periodic_mpi_groupings[i].tag;
		 
			int tag_inc = n_groups;
			for(unsigned int i=1; i<N_MC_GROUPS; ++i)  
			{
			   MPI_Send(&tag_inc, 1, MPI_INT, rank+NPART*i,rank,MPI_COMM_WORLD);
			   MPI_Send(&buf, periodic_mpi_groupings.size(), MPI_INT, rank+NPART*i,NPROC+rank,MPI_COMM_WORLD); 
			   tag_inc+=n_groups;
			}		
		 }   
	  }
	  else
	  {
		 int source = rank%NPART;
		 int tag = source;
		 int tag_inc;
		 int buf[periodic_mpi_groupings.size()];
		 MPI_Recv(&tag_inc, 1, MPI_INT, source,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		 MPI_Recv(&buf, periodic_mpi_groupings.size(), MPI_INT, source,NPROC+tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
	 
		 for(unsigned int i=0; i<periodic_mpi_groupings.size(); ++i)
		 {
			periodic_mpi_groupings[i].tag = buf[i]+tag_inc;  
	     }		
	  }
	  mpi_barrier(MPI_COMM_WORLD);
   }   
   
   // mpi_barrier(MPI_COMM_WORLD);
//    for(int n =0 ; n<NPROC; ++n)
//    {
//       if(get_proc_id() == n)
//       {
//          cout<<"PROC :"<<n<<endl;
//          cout<<"SHARED LIST ::"<<endl;
//          for(unsigned int i=0; i<intp_mpi_groupings.size(); ++i)
//          {
//              for(set<int>::iterator it=intp_mpi_groupings[i].proc_list.begin();
//                  it!=intp_mpi_groupings[i].proc_list.end();++it)
//                  cout<<*it<<" ";
// 			 cout<<intp_mpi_groupings[i].tag<<endl;
// 	     }
// 	     cout<<"PERIODIC LIST ::"<<endl;
//          for(unsigned int i=0; i<periodic_mpi_groupings.size(); ++i)
//          {
//              for(set<int>::iterator it=periodic_mpi_groupings[i].proc_list.begin();
//                  it!=periodic_mpi_groupings[i].proc_list.end();++it)
//                  cout<<*it<<" ";
// 			 cout<<periodic_mpi_groupings[i].tag<<endl;
// 	     }		 
//      }
//      mpi_barrier(MPI_COMM_WORLD);
//    }   
      
}

