#ifndef __GRID_H__
#define __GRID_H__

#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <iostream>
#include <cassert>
#include "vec.h"
#define NDIR 3

//extern bool periodic;

//------------------------------------------------------------------------------
//!  Cell structure 
/*!
  Describes the data structures associated with a mesh cell (primary cell)
*/
//------------------------------------------------------------------------------
template <int DIM>
struct Cell
{
      unsigned int vertex[DIM+1];       /*!< ids of vertices */ 
};

//------------------------------------------------------------------------------
//!  Vertex structure 
/*!
  Describes the data structures associated with a vertex/node
*/
//------------------------------------------------------------------------------
struct Vertex
{
   Vector coord;						    /*!< node coordinates */ 
   std::set<int> all_shared;                /*!< stores set of processors in the vertex 
									             is a present */
   std::set<int> active_shared;                  /*!< stores set of processors in the vertex 
   												 is active */
   bool d_bnd;						        /*!< true if the vertex lies on domain boundary */													 
};   

//------------------------------------------------------------------------------
//!  Face structure 
/*!
  Describes the data structures associated with a face
*/
//------------------------------------------------------------------------------
template <int DIM>
struct Face
{
   unsigned int vertex[DIM+1];        /*!< two vertices of the face */
   int          type; 			      /*!< face tag from mesh file*/
   int          part;                 /*!< partition the face belongs to*/
   bool         periodic;             /*!< true if face lies on a periodic boundary*/
};

//------------------------------------------------------------------------------
//!  Grid class 
/*!
  Deals with the grid data structures and mesh partitioning
*/
//------------------------------------------------------------------------------
template <int DIM>
class Grid
{
   public:
      Grid(char *file, char *loc);
      
      unsigned int geo_type[2];   /*!< geometry type of elements being read */ 				 
      unsigned int n_vertex;			/*!< number of vertices in the mesh */
      unsigned int n_face;	            /*!< number of total boundary faces in the mesh */
      unsigned int n_cell;	            /*!< number of total cells in the mesh */
      std::vector<Vertex> vertex;		/*!< list of vertices */
      std::vector<Face<DIM> > face;		    /*!< list of faces */
      std::map<int,std::vector<int> > periodic_nodes;   /*!< maps between periodic node IDs.
       														A node may have more than one periodic
       														counterpart, thus a vector is needed.
                                                            Needed for handling periodic mesh */                                     
      std::vector<int> periodic_face_tag;   /*!< geo_tags of periodic face IDs.
                                             Needed for handling periodic mesh */                                                                                       
      std::map<int,int> part_nface;        /*!< maps between partition ID and number of boundary
      											faces */
      std::map<int,int> part_ncell;        /*!< maps between partition ID and number of cells */
      std::map<int,int> part_nvertex;      /*!< maps between partition ID and number of vertices */											
      std::set<int> partitions;            /*!< Ids for all partitions */
      std::string msh_file;                /*!< mesh file obtained from gmsh*/
      std::string part_dir;                /*!< directory where partition files are to be stores*/
      bool periodic_available;              /*!< True if periodic data is available in mesh file*/
      std::vector<std::vector<int> > periodic_node_sets;/*!< vector of periodic (vector) node sets */
      std::vector<std::set<int> > periodic_set_part_list;/*!< vector of sets listing partitions corresponding to periodic sets */
      std::map<int,int> part_periodic_set;/*!< maps between partition ID and number of periodic sets */	
      /*!
	  * This is the main function of this class
	  */
      void read ();  


   private:
      /*!
	  * Reads gmsh files and creates pseudo partition files
	  */
      void read_gmsh ();
      
      /*!
	  * Splits gmsh files into nodes, faces, cells and periodic data files
	  */
      void gmsh_split ();
      
      /*!
	  * Completes periodic links between nodes. This is used to complete links between
	  * corner nodes and edge nodes (3D).
	  */
      void complete_connections();
      
      /*!
	  * Checks if faces are equal. Returns 0 if the faces do not match.
	  * Returns 1 if the faces match and have the same orientation
	  * Returns -1 if the faces match and have opposite orientations
	  * @param[in] master and slave Face object to be compared
	  */
      //int compare_faces (Face<DIM> &master, Face<DIM> &slave);
      
      /*!
	  * Writes cell data into its corresponding cell partition file
	  * @param[in] cell Cell object to be written
	  * @param[in] part_tag Partition the cell belongs to (C++ indexing)
	  * @param[in] g_cell True if it is a ghost cell in the partition
	  */
      void write_cell (Cell<DIM> &cell, unsigned int part_tag, bool g_cell);
      
      /*!
	  * Intializes main mesh partition files
	  */
      //void initialize_partition_files();
      
      /*!
	  * Creates main mesh partition files
	  */
      void create_partition_files();

};

template <int DIM>
inline
Grid<DIM>::Grid(char *file, char *loc)
{
   n_vertex = n_face = n_cell = 0;
   periodic_available = false;
   msh_file  = file;
   part_dir  = loc;
   part_dir +=  "/PARTITION";
   if(DIM == 2)
   {
	  geo_type[0] = 1; geo_type[1] = 2;
   }
   else if(DIM == 3)
   {
	  geo_type[0] = 2; geo_type[1] = 4;
   }
};

#endif
