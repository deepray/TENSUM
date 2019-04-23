#ifndef __GRID_H__
#define __GRID_H__

#include <vector>
#include <set>
#include "parameter.h"
#include "vec.h"
#include "face.h"
#include "mpi.h"
#include "mpi_utils.h"
#define NDIR 3

//------------------------------------------------------------------------------
//!  Cell class 
/*!
  Describes the data structures associated with a mesh cell (primary cell)
*/
//------------------------------------------------------------------------------
class Cell
{
   public:
      Vector       centroid;        /*!< cell centroid coordinates */ 
      unsigned int vertex[3];       /*!< ids of vertices */ 
      int          face[3];			/*!< ids of faces */ 
      double       area;			/*!< cell area (primary) */ 
      double       radius;			/*!< used for 3D axi-symmetric flows */ 
      std::vector<Vector> normal; 	/*!< inward face normals */ 
      // int          type; 			/*!<  0 : interior cell
//        									  
//        									  1 : cell with one boundary edge (not including periodic)
//        									  
//        									  2 : cell with two boundary edges (not including periodic)*/ 
      bool         g_cell;          /*!< true if the cell is a ghost cell in a 
                                         given MPI partition */ 
      bool         mpi_bnd_active;  /*!< true if the cell is an active cell in the given
                                         partition and has an shared node */ 
      std::vector<int> nbr_bface;	    /*!<  indices of faces which are boundary faces (non periodic)
                                              : these can take the values 0,1 or 2*/                                   
};


//------------------------------------------------------------------------------
//!  Vertex class 
/*!
  Describes the data structures associated with a vertex/node
*/
//------------------------------------------------------------------------------
struct Vertex
{
   Vector coord;						    /*!< node coordinates */ 
   double radius;							/*!< used for 3D axi-symmetric flows */ 
   bool active;								/*!< true if the vertex is an active vertex in
   												 the given partition */ 
   bool mpi_bnd_active;						/*!< true if the vertex is an active vertex in
   												 the given partition and is also a
   												 shared node */ 
   bool d_bnd;						        /*!< true if the vertex lies on domain boundary */												 
};

//------------------------------------------------------------------------------
/*!
  Structure describing the set of processors sharing a vertex
*/
//------------------------------------------------------------------------------
struct SHARED_LIST
{
   int vertex;                  /*!< vertex id */ 
   std::set<int> proc_list;     /*!< set of processors sharing this vertex */ 
};

//------------------------------------------------------------------------------
/*!
  Structure describing the the ghost and active cells neighbouring a vertex
*/
//------------------------------------------------------------------------------
struct VERTEX_MPI_NBR
{
   int    gc_nbr;    /*!< number of adjoining ghost cells */ 
   int    ngc_nbr;   /*!< number of adjoining active cells */  
};

//------------------------------------------------------------------------------
/*!
  Structure describing the grouping for processors with shared vertices
*/
//------------------------------------------------------------------------------
struct INTP_VERTEX_GROUPING_LIST
{
   std::set<int> proc_list;   /*!< set of processors sharing a group of vertices */
   std::vector<int> vertices; /*!< shared vertices*/
   int tag;                   /*!< tag for group */
};

//------------------------------------------------------------------------------
/*!
  Structure for boundary ghost cells
*/
//------------------------------------------------------------------------------
struct BOUNDARY_GCELL
{
   int cell_id;                                /*!< local id of ghost cell */
   std::vector<std::pair<int, int> > bface_vertices; /*!< pair of vertices forming domain boundary face */
};

//------------------------------------------------------------------------------
/*!
  Periodic vertex data
*/
//------------------------------------------------------------------------------
struct PERIODIC_LISTS
{
   int glob_id          ; /*!< global ID of P_i*/
   std::vector<int> P_ij; /*!< subset of P_i of periodic nodes in partition j */
   std::set<int> L_i;  /*!< partition sharing list for list P_i */
};

//------------------------------------------------------------------------------
//!  Grid class 
/*!
  Deals with the grid data structures and mesh partitioning
*/
//------------------------------------------------------------------------------
class Grid
{
   public:
      Grid () 
      : 
      	 n_vertex (0),
      	 n_periodic (0),
      	 n_cell (0),
      	 n_face (0),
      	 n_boundary_face (0),
      	 n_vertex_active (0),
      	 n_cell_active (0),
      	 n_cell_mpi_bnd (0),
      	 n_cell_mpi_int (0),
      	 n_vertex_total (0), 
      	 n_cell_total (0),
      	 n_vertex_mpi_bnd (0)
      {MESH_DATA_DIR = "MESH_DATA";};
      CellType cell_type;				/*!< median or voronoi */
      unsigned int n_vertex;			/*!< number of vertices in a partition */
      unsigned int n_periodic;			/*!< number of periodic vertex subsets P_ij */
      unsigned int n_cell;				/*!< number of cells in a partition */	
      unsigned int n_face;				/*!< number of faces in a partition */
      unsigned int n_boundary_face;		/*!< number of boundary (mesh) faces in a partition */
      std::vector<Vertex> vertex;		/*!< list of vertices in a partition */
      std::vector<Cell>   cell;			/*!< list of cells in a partition */
      std::vector<Face>   face;			/*!< list of faces in a partition */
      std::vector<Face>   bface;		/*!< list of boundary faces in a partition */	
      std::vector<BOUNDARY_GCELL> b_gcell; /*!< list of boundary ghost cells (needed for voronoi cells) */	
      std::vector<double> mcarea;		/*!< list of median cell areas (for median mesh) */	
      std::vector<double> dcarea;		/*!< list of dual cell areas */
      std::vector<double> dx_max;		/*!< max of distance of node from neighbouring nodes */	
      std::map <int,int> l2g, 			/*!< local-2-global index maps for shared nodes of the 
                                             partition */
      				     g2l;           /*!< global-2-local index maps for shared nodes of the 
                                             partition */
  
      std::vector< INTP_VERTEX_GROUPING_LIST > intp_mpi_groupings; /*!< information about 
      																	various processor 
      																	groups the current 
      																	partition belongs 
      																	to, based on shared 
      																	vertices. The vertex IDS
      																	are global */
      std::vector< INTP_VERTEX_GROUPING_LIST > periodic_mpi_groupings; /*!< information about 
      																	    various processor 
      																	    groups the current 
      																	    partition belongs 
      																	    to via, 
      																	    based on periodic 
      																	    vertices.The vertex IDS
      																	    are local */																	
      std::vector< PERIODIC_LISTS > periodic_lists; /*!< Lists of periodic subsets and their
        												 partition sharing information */
      std::set< std::set <int> > all_intp_groupings;    /*!< list of all existing processor
      														 groupings with shared nodes */
      std::set< std::set <int> > all_periodic_groupings;    /*!< list of all existing processor
      														     groupings with periodic nodes */														 
      MPI_Comm run_comm;  /*!< the sample run group communicator for the particular processor */
      
      unsigned int n_vertex_active, /*!< number of active vertices in a partition */
                   n_cell_active,   /*!< number of active cells in a partition */
                   n_cell_mpi_bnd,  /*!< number of active cells in a partition with at 
                   						 least one shared node */
                   n_cell_mpi_int;  /*!< number of active cells in a partition with no 
                   						 shared nodes */
      int n_mpi_face;               /*!< number of inter-partition faces in the given 
                                         partition */
      
      std::map<int,int> active_node_ind;   // list of pairs (i,j) where i is the local index of an active node
							               // while j is the local index of this node in the collection of all
                                           // nodes. Used in writer.cc while writing solution to file
      std::map<int,int> active_cell_ind;   // list of pairs (i,j) where i is the local index of an active cell
										   // while j is the local index of this cell in the collection of all
										   // active cells.
      std::map<int,int> mpi_node_share;    /*!< list of pairs (i,j) where i is the local index 
      											of an active node while j is the number of 
      											processors sharing it. This is used while evaluating 
      											the residual norm and global
      											quantities in fv.cc */  
                                                
      std::map<int,double> periodic_node_weight; /*!< list of pairs (i,j) where i is the local index 
      											of an active periodic node while j is the inverse of
      											the number of other active periodic nodes in the partition. 
      											This is used while evaluating the residual norm and global
      											quantities in fv.cc */                                                                                 
      
      unsigned int n_vertex_total,         /*!< number of nodes in the full mesh */
                   n_cell_total, 		   /*!< number of cells in the full mesh */
                   n_face_total; 		   /*!< number of boundary faces in the full mesh */
                   
      unsigned int n_vertex_mpi_bnd;       /*!< number of vertices which lie on an 
      											inter-partition boundary for the given 
      											partition */
      std::vector<int> cell_mpi_bnd;       /*!< list of local indices of active cells that 
      											have at least one shared node */
      std::vector<int> cell_mpi_int;       /*!< list of local indices of active cells that 
      											have no shared nodes */ 
      std::vector<int> d_bnd_shared_face;  /*!< list of mesh boundary faces with a shared 
      										    node */
      std::vector<int> d_bnd_unshared_face;/*!< list of mesh boundary faces with no shared 
      										    nodes */
      std::map<int,int> l2g_mc_pt;         /*!< local to global id map for point at which 
      											pdf is evaluated */
      
      std::string MESH_DATA_DIR;           /*!< directory where mesh details such as cells
      											and area data files are stored */
      
      /*!
	  * This is the main function of this class
	  */
      void read (const Parameter& param);  


   private:
      /*!
	  * Reads partitioned mesh files and creates basic grid mappings
	  * @param[in] PART_DIR directory containing partitioned mesh files
	  */
      void read_part_mesh_file(std::string PART_DIR,
                               const std::map<int,BoundaryCondition>& bc);
      
      /*!
	  * Check that all boundary faces have been assigned a bc type
	  * @param[in] bc mapping between face tags and boundary conditions
	  */
      void check_face_type (const std::map<int,BoundaryCondition>& bc);
      
      /*!
	  * Creates all the necessary mesh data structures
	  */
      void preproc (const std::map<int,BoundaryCondition>& bc);
      
      /*!
	  * Finds cell centroid
	  */
      void compute_cell_centroid ();
      
      /*!
	  * Finds cell face centroid
	  */
      void compute_face_centroid ();
      
      /*!
	  * Finds cell area
	  */
      void compute_cell_area ();
      
      /*!
	  * Reduce dual cell areas at periodic points
	  */
      void reduce_dual_area ();
      
      /*!
	  * Finds face normals
	  */
      void compute_face_normal_and_area ();
      
      /*!
	  * Initially, only boundary faces are listed in the mesh file. This function generates
	  * the remaining (interior) faces
	  */
      void add_face (const Face& new_face);
      
      /*!
	  * Finds the vertices (two for interior faces) opposite the faces
	  */
      void find_vertex_opposite_face ();
      
      /*!
	  * This function makes all the interior faces. It calls the function add_face()
	  */
      void make_faces ();
      
      /*!
	  * Finds faces for all primary cells
	  */
      void find_cell_faces ();
      
      /*!
	  * Prints mesh information
	  */
      void info ();
      
      /*!
	  * Removes faces with with zero area
	  */
      void remove_empty_faces ();
      
      /*!
	  * Prints primary and dual cell into a file
	  */
      void print_cells();
      
      /*!
	  * Prints dual cell areas into file
	  */
      void print_areas();
      
      /*!
	  * Finds radius for axisymmetric case
	  */
      void compute_radius ();
      
      /*!
	  * Finds dx_max for Peclet number
	  */
      void compute_dx_max ();
      
      /*!
	  * Mark boundary edges of cells on the boundary
	  */
      void mark_cell_bface(const std::map<int,BoundaryCondition>& bc);
      
      /*!
	  * Removes ghost cells and nodes from data structures
	  */
      void remove_ghosts();
      
      /*!
	  * Creates list of weights for active periodic nodes
	  */
      void create_periodic_weights ();
      
      /*!
	  * Checks if grid data structures have correct labels/tags
	  */
      void check_data();

      std::vector< std::vector<unsigned int> > node_face;  /*!< Required for making faces */

};

#endif
