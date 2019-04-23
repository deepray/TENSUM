#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "grid.h"
#include "mpi_utils.h"

extern Dimension dim;
extern bool debug;
extern bool verbose;
extern bool PERIODIC;
using namespace std;

//------------------------------------------------------------------------------
// Compute cell Centroid
//------------------------------------------------------------------------------
void Grid::compute_cell_centroid ()
{
   MPI_DISP("  --- finding cell centroid ... ",verbose);

   if(cell_type == median)
   {
      MPI_DISP("       ------> median cell",verbose);
      double fact = 1.0/3.0;   
      for(unsigned int i=0; i<n_cell; ++i)
      { 
         unsigned int v0, v1, v2;
         v0 = cell[i].vertex[0];
         v1 = cell[i].vertex[1];
         v2 = cell[i].vertex[2];
         cell[i].centroid.equ(vertex[v0].coord, 
                              vertex[v1].coord, 
                              vertex[v2].coord,
                              fact,fact,fact);
      }   
   }
   else if(cell_type == voronoi)
   {
      MPI_DISP("       ... voronoi cell",verbose);
      for(unsigned int i=0; i<n_cell; ++i)
      { 
         unsigned int n1 = cell[i].vertex[0];
         unsigned int n2 = cell[i].vertex[1];
         unsigned int n3 = cell[i].vertex[2];

         double x1 = vertex[n1].coord.x;
         double y1 = vertex[n1].coord.y;
         double x2 = vertex[n2].coord.x;
         double y2 = vertex[n2].coord.y;
         double x3 = vertex[n3].coord.x;
         double y3 = vertex[n3].coord.y;

         double l1  = pow(x2 - x3, 2) + pow(y2 - y3, 2);
         double l2  = pow(x3 - x1, 2) + pow(y3 - y1, 2);
         double l3  = pow(x1 - x2, 2) + pow(y1 - y2, 2);

         double beta1 = max(0.0, l2 + l3 - l1);
         double beta2 = max(0.0, l3 + l1 - l2);
         double beta3 = max(0.0, l1 + l2 - l3);

         // This fix is supposed to remove very small cv faces.
         // I am not totally happy with this one.
         if(beta1 < beta2/2.0 && beta1 < beta3/2.0) beta1=0.0;
         if(beta2 < beta3/2.0 && beta2 < beta1/2.0) beta2=0.0;
         if(beta3 < beta1/2.0 && beta3 < beta2/2.0) beta3=0.0;

         double det = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);
         double b1  = 0.5*( (x2 - x1)*(x2 + x1) + (y2 - y1)*(y2 + y1) );
         double b2  = 0.5*( (x3 - x1)*(x3 + x1) + (y3 - y1)*(y3 + y1) );
         double xc  = ( (y3-y1)*b1 - (y2-y1)*b2)/det;
         double yc  = (-(x3-x1)*b1 + (x2-x1)*b2)/det;

         if(beta1 == 0.0)
         {
            cell[i].centroid.x = 0.5*(x2 + x3);
            cell[i].centroid.y = 0.5*(y2 + y3);
         }
         else if(beta2 == 0.0)
         {
            cell[i].centroid.x = 0.5*(x3 + x1);
            cell[i].centroid.y = 0.5*(y3 + y1);
         }
         else if(beta3 == 0.0)
         {
            cell[i].centroid.x = 0.5*(x1 + x2);
            cell[i].centroid.y = 0.5*(y1 + y2);
         }
         else
         {
            cell[i].centroid.x = xc;
            cell[i].centroid.y = yc;
         }

      }
      // check if centroid (circumcenter) is on boundary face
      // If yes, then set centroid to geometric centroid
      // Otherwise, boundary vertices become decoupled as there is no flux
      // between the two vertices
      unsigned int count = 0;
      double fact = 1.0/3.0;
      for(unsigned int i=0; i<bface.size(); ++i)
      {
         int lcell = bface[i].lcell;
         Vector dr;
         dr.equ(bface[i].centroid,cell[lcell].centroid,1.0,-1.0);
         if(dr.norm() < 1.0e-14)
         {
            unsigned int v0 = cell[lcell].vertex[0];
            unsigned int v1 = cell[lcell].vertex[1];
            unsigned int v2 = cell[lcell].vertex[2];
            cell[lcell].centroid.equ(vertex[v0].coord,
                                     vertex[v1].coord,
                                     vertex[v2].coord,
                                     fact,fact,fact);
            ++count;
         }
      }
      
      //Correction for boundary ghost cells
      for(unsigned int i=0; i<b_gcell.size(); ++i)
      {
         int lcell = b_gcell[i].cell_id;
         double dr_norm = 1.0e10;
         int v0,v1,v2;
         for(unsigned int j=0; j<b_gcell[i].bface_vertices.size(); ++j)
         {
            Vector dr, face_centroid;
            v0 = b_gcell[i].bface_vertices[j].first;
            v1 = b_gcell[i].bface_vertices[j].second;
            face_centroid.equ(vertex[v0].coord,vertex[v1].coord,0.5,0.5);
            dr.equ(face_centroid,cell[lcell].centroid,1.0,-1.0);
            dr_norm=min(dr_norm,dr.norm());
         }
         if(dr_norm < 1.0e-14)
         {
            v0 = cell[lcell].vertex[0];
            v1 = cell[lcell].vertex[1];
            v2 = cell[lcell].vertex[2];
            cell[lcell].centroid.equ(vertex[v0].coord,
                                     vertex[v1].coord,
                                     vertex[v2].coord,
                                     fact,fact,fact);        
            ++count;
         }
      }
      // if(count > 0)
//          cout << "  " << count << " boundary faces had circumcenters in partition "<<get_proc_id()<< endl;
   }
}

//------------------------------------------------------------------------------
// Compute face Centroid
//------------------------------------------------------------------------------
void Grid::compute_face_centroid ()
{
   MPI_DISP("  --- computing face centers ...",verbose);   
   for(unsigned int i=0; i<n_face; ++i)
   { 
      unsigned int v0 = face[i].vertex[0];
      unsigned int v1 = face[i].vertex[1];
      face[i].centroid.equ(vertex[v0].coord,vertex[v1].coord,0.5,0.5);
   }   

   for(unsigned int i=0; i<bface.size(); ++i)
   { 
      unsigned int v0 = bface[i].vertex[0];
      unsigned int v1 = bface[i].vertex[1];
      bface[i].centroid.equ(vertex[v0].coord,vertex[v1].coord,0.5,0.5);
   }   
}
//------------------------------------------------------------------------------
// Compute cell areas
//------------------------------------------------------------------------------
void Grid::compute_cell_area ()
{
   MPI_DISP("  --- computing cell areas ...",verbose);

   mcarea.resize (n_vertex);
   dcarea.resize (n_vertex);

   for(unsigned int i=0; i<n_vertex; ++i)
   {
      mcarea[i] = 0.0;
      dcarea[i] = 0.0;
   }

   for(unsigned int i=0; i<n_cell; ++i)
   {
      unsigned int v0 = cell[i].vertex[0];
      unsigned int v1 = cell[i].vertex[1];
      unsigned int v2 = cell[i].vertex[2];

      // compute area as vector cross product
      Vector l1,l2,area;
      l1.equ(vertex[v1].coord,vertex[v0].coord,1.0,-1.0);
      l2.equ(vertex[v2].coord,vertex[v0].coord,1.0,-1.0);
      area.cprod(l1,l2);
      cell[i].area = 0.5 * area.z;

      MPI_LOC_ASSERT ( cell[i].area > 0.0 );

      // contribution to dual cell
      mcarea[v0] += cell[i].area / 3.0;
      mcarea[v1] += cell[i].area / 3.0;
      mcarea[v2] += cell[i].area / 3.0;
   }

   if(cell_type == median)
      dcarea = mcarea;
   else
   {
      for(unsigned int i=0; i<n_cell; ++i)
      {
         unsigned int n0, n1, n2;
         for(unsigned int j=0; j<3; ++j)
         {
            n0 = cell[i].vertex[j];
            if(j==0)
               n1 = cell[i].vertex[2];
            else
               n1 = cell[i].vertex[j-1];
            if(j==2)
               n2 = cell[i].vertex[0];
            else
               n2 = cell[i].vertex[j+1];
            vector<Vector> point(5);
            point[0] = vertex[n0].coord;
            point[1].equ(vertex[n0].coord,vertex[n2].coord,0.5,0.5);
            point[2] = cell[i].centroid;
            point[3].equ(vertex[n0].coord,vertex[n1].coord,0.5,0.5);
            point[4] = point[0];

            double area = 0;
            for(unsigned int k=0; k<4; ++k)
               area += point[k].x * point[k+1].y - point[k+1].x * point[k].y;
            if(area < 0.0)
               MPI_LOC_ERR("Dual cell area is non-positive\n"
                                    <<"   Area     = " << area << "\n"
                                    <<"   Triangle = " << i << "\n");
            
            area *= 0.5;
            dcarea[n0] += area;
         }
      }
   }
   if(PERIODIC)
   {
      MPI_DISP("  --- reducing dual cell areas at periodic nodes ...",verbose);
      reduce_dual_area();
   }   
   mpi_barrier(MPI_COMM_WORLD);

}

//------------------------------------------------------------------------------
// Compute face normals
//------------------------------------------------------------------------------
void Grid::compute_face_normal_and_area ()
{

   MPI_DISP("  --- computing dual face data ...",verbose);
   // Check orientation of interior faces
   for(unsigned int i=0; i<n_face; ++i)
   {
      unsigned int v0 = face[i].vertex[0];
      unsigned int v1 = face[i].vertex[1];
      //Vector dr = vertex[v1].coord - vertex[v0].coord;
      Vector dr;
      dr.equ(vertex[v1].coord,vertex[v0].coord,1.0,-1.0);
      
      Vector normal;
      normal.x = +dr.y;
      normal.y = -dr.x;
      normal.z =  0.0;
            
      // Check orintation of boundary face
      unsigned int cl = face[i].lcell;
      Vector dcf;
      dcf.equ(face[i].centroid,cell[cl].centroid,1.0,-1.0);
      double d = dcf * normal;
      if(d < 0.0)
      {
         face[i].vertex[0] = v1;
         face[i].vertex[1] = v0;
      }
   }
   
   // interior faces
   for(unsigned int i=0; i<n_face; ++i)
   {
      unsigned int cl = face[i].lcell;
      Vector dr;
      dr.equ(face[i].centroid,cell[cl].centroid,1.0,-1.0);

      if(face[i].type == -1 || face[i].type == -2) // interior edge, has right cell also
      {
         unsigned int cr = face[i].rcell;
         dr.sadd(cell[cr].centroid,face[i].centroid,1.0,-1.0);
      }

      face[i].normal.x = -dr.y;
      face[i].normal.y =  dr.x;
      face[i].normal.z =  0.0;

      face[i].measure = face[i].normal.norm();
  }

   // boundary faces
   for(unsigned int i=0; i<bface.size(); ++i)
   {
      unsigned int v0 = bface[i].vertex[0];
      unsigned int v1 = bface[i].vertex[1];
      Vector dr;
      dr.equ(vertex[v1].coord,vertex[v0].coord,1.0,-1.0);

      bface[i].normal.x =  dr.y;
      bface[i].normal.y = -dr.x;
      bface[i].normal.z =  0.0;

      bface[i].measure = bface[i].normal.norm();
      
      // Check orintation of boundary face
      unsigned int cl = bface[i].lcell;
      Vector dcf;
      dcf.equ(bface[i].centroid,cell[cl].centroid,1.0,-1.0);
      double d = dcf * bface[i].normal;
      if(d < 0.0)
      {
         bface[i].vertex[0] = v1;
         bface[i].vertex[1] = v0;
         bface[i].normal *= -1.0;
      }
   }

   // Gradient of P1 basis function on each triangle
   // You need to divide by 2*area to get gradient
   Vector dr01, dr12, dr20;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      cell[i].normal.resize(3);

      unsigned int v0 = cell[i].vertex[0];
      unsigned int v1 = cell[i].vertex[1];
      unsigned int v2 = cell[i].vertex[2];

      dr01.x =  (vertex[v0].coord.y - vertex[v1].coord.y);
      dr12.x =  (vertex[v1].coord.y - vertex[v2].coord.y);
      dr20.x =  (vertex[v2].coord.y - vertex[v0].coord.y);

      dr01.y = -(vertex[v0].coord.x - vertex[v1].coord.x);
      dr12.y = -(vertex[v1].coord.x - vertex[v2].coord.x);
      dr20.y = -(vertex[v2].coord.x - vertex[v0].coord.x);

      dr01.z =  0.0;
      dr12.z =  0.0;
      dr20.z =  0.0;

      cell[i].normal[0] = dr12;
      cell[i].normal[1] = dr20;
      cell[i].normal[2] = dr01;
   }


}

//------------------------------------------------------------------------------
// Add new_face to face list
//------------------------------------------------------------------------------
void Grid::add_face (const Face& new_face)
{
   bool found = false;
   unsigned int n = 0;

   // Take any vertex of this face
   unsigned int v = new_face.vertex[0];

   // Loop over all existing faces of vertex v
   while (n<node_face[v].size() && !found)
   {
      unsigned int f = node_face[v][n];

      if(face[f].lcell==-1) // Boundary face not filled yet
      {
         if(face[f] == new_face)
         {
            face[f].lcell   = new_face.lcell;
            found = true;
         }
      }
      else if(face[f].rcell == -1) // Boundary or interior face
      {
         if(face[f] == new_face)
         {
            face[f].rcell   = new_face.lcell;
            found = true;
         }
      }

      ++n;
   }

   if(!found) // This is a new face
   {
      face.resize (n_face+1);
      face[n_face].type = -1; 
      face[n_face].vertex [0] = new_face.vertex [0];
      face[n_face].vertex [1] = new_face.vertex [1];
      face[n_face].lcell      = new_face.lcell; // left cell
      face[n_face].rcell      = -1; // right cell to be found
      face[n_face].periodic   = false;

      // Add this face to its two vertices
      for(unsigned int i=0; i<2; ++i)
      {
         v = new_face.vertex[i];
         node_face[v].push_back (n_face);
      }

      ++n_face;
   }

}

//------------------------------------------------------------------------------
// Find vertex opposite each face; for boundary only left vertex present
//------------------------------------------------------------------------------
void Grid::find_vertex_opposite_face ()
{
   for(unsigned int i=0; i<n_face; ++i)
   {
	  int vl = -1;
	  unsigned int v0 = face[i].vertex[0];
	  unsigned int v1 = face[i].vertex[1];
	  int cl = face[i].lcell;
	  for(unsigned int j=0; j<3; ++j)
		 if(cell[cl].vertex[j] != v0 && cell[cl].vertex[j] != v1)
			vl = cell[cl].vertex[j];
	  MPI_LOC_ASSERT (vl != -1);
	  face[i].lvertex = vl;

	  // interior face and mpi shared face: find right vertex also
	  if(face[i].type == -1 || face[i].type == -2)
	  {
		 int vr = -1;
		 int cr = face[i].rcell;
		 for(unsigned int j=0; j<3; ++j)
			if(cell[cr].vertex[j] != v0 && cell[cr].vertex[j] != v1)
			   vr = cell[cr].vertex[j];
		 MPI_LOC_ASSERT (vr != -1);
		 face[i].rvertex = vr;
	  }
	  else
		 face[i].rvertex = -1;
   }
}

//------------------------------------------------------------------------------
// Create interior faces and connectivity data
//------------------------------------------------------------------------------
void Grid::make_faces ()
{
   MPI_DISP("  --- creating faces ...",verbose);
   unsigned int i;

   node_face.resize (n_vertex);

   // Existing boundary faces
   for(i=0; i<n_face; ++i)
   {
      face[i].lcell   = -1;
      face[i].rcell   = -1;

      // Add this face to the two vertices
      for(unsigned int j=0; j<2; ++j)
      {
         unsigned int v = face[i].vertex[j];
         node_face[v].push_back(i);
      }
   }

   Face new_face;

   for(i=0; i<n_cell; ++i)
   {
      // first face
      new_face.vertex[0] = cell[i].vertex[0];
      new_face.vertex[1] = cell[i].vertex[1];
      new_face.lcell     = i;
      new_face.periodic  = false;
      add_face (new_face);

      // second face
      new_face.vertex[0] = cell[i].vertex[1];
      new_face.vertex[1] = cell[i].vertex[2];
      new_face.lcell     = i;
      new_face.periodic  = false;
      add_face (new_face);

      // third face
      new_face.vertex[0] = cell[i].vertex[2];
      new_face.vertex[1] = cell[i].vertex[0];
      new_face.lcell     = i;
      new_face.periodic  = false;
      add_face (new_face);
   }

   // Removing ghost faces and changing tag of inter-partition faces
   n_mpi_face = 0;
   vector<Face>::iterator it = face.begin();
   while(it!=face.end())
   {
      int lcell = (*it).lcell;
      int rcell = (*it).rcell;
	  
	  (*it).d_bnd_and_shared = false;
	  if((*it).type == -1) 
	  {
	     if((*it).rcell == -1) // GHOST DOMAIN BOUNDARY FACES
		    face.erase(it); 
		 else if(cell[lcell].g_cell && cell[rcell].g_cell) // GHOST INTERIOR FACES
		    face.erase(it);    
		 else if((cell[lcell].g_cell && !cell[rcell].g_cell) ||
		         (!cell[lcell].g_cell && cell[rcell].g_cell)) // INTER_PART FACES
		 // NOTE: A ghost face can have both its nodes as mpi boundary active. So can't
		 // use that as a check. The same can be true for an interior face
		 {
		    (*it).type = -2; // Type for inter-processor faces
		    n_mpi_face++; 
		    it++;
	     }   
	     else
	        it++; 
	  }
	  else 
	  {
		 // If domain boundary face has a shared vertex, then set d_bnd_and_active = true
		 if(vertex[(*it).vertex[0]].mpi_bnd_active || vertex[(*it).vertex[1]].mpi_bnd_active)
		 {
		    (*it).d_bnd_and_shared = true;
		 }    
		 it++;
	  }
	  
   }
   n_face = face.size();
   
   MPI_DISP("  --- checking face data ...",verbose);

   // Now check that face data is complete
   for(i=0; i<n_face; ++i)
   {
      // Check left cell exists
      MPI_LOC_ASSERT(face[i].lcell != -1);
      if(face[i].type == -1 || face[i].type == -2) // Interior face and inter-partition face
      {                                             //, check right cell
            MPI_LOC_ASSERT(face[i].rcell != -1);
      }   
   }

   find_vertex_opposite_face ();

   // Copy boundary faces into bface
   for(i=0; i<n_face; ++i)
      if(face[i].type != -1 && face[i].type != -2)
      {
         bface.push_back(face[i]);  
         if(face[i].d_bnd_and_shared)
            d_bnd_shared_face.push_back(i);
         else 
            d_bnd_unshared_face.push_back(i); 
      }
   
   // for(i=0; i<bface.size(); ++i)
//    { 
//       int v0 = bface[i].vertex[0];
//       int v1 = bface[i].vertex[1];
//       cout<<bface[i].type<<" "
//            <<vertex[v0].coord.x<<" "
//             <<vertex[v0].coord.y<<" "
//             <<vertex[v1].coord.x<<" "
//             <<vertex[v1].coord.y<<" "
//             <<bface[i].periodic<<endl;
//    }         
   // Free memory of node_face since we dont need it any more
   for(i=0; i<n_vertex; ++i)
      node_face[i].resize (0);
   node_face.resize (0);
   
}

//------------------------------------------------------------------------------
// Find faces surrounding a cell
//------------------------------------------------------------------------------
void Grid::find_cell_faces ()
{
   MPI_DISP("  --- finding faces for each cell ...",verbose);
   
   unsigned int i, j;
   int lcell, rcell;
   
   // First set all faces to -1
   for(i=0; i<n_cell; ++i)
   {
      cell[i].face[0] = -1;
      cell[i].face[1] = -1;
      cell[i].face[2] = -1;
   }
   
   for(i=0; i<n_face; ++i)
   {
      lcell = face[i].lcell;
      j = 0;
      while(cell[lcell].face[j] != -1)
         ++j;
      cell[lcell].face[j] = i;
            
      rcell = face[i].rcell;
      if(rcell != -1)
      { 
         j = 0;
         while(cell[rcell].face[j] != -1)
            ++j;
         cell[rcell].face[j] = i;
      }
    }
   
}

//------------------------------------------------------------------------------
// Find points connected to a point
//------------------------------------------------------------------------------
void Grid::remove_empty_faces()
{
   MPI_DISP("  --- removing faces with zero area ...",verbose);
   vector<Face> tmp;
   for(unsigned int i=0; i<n_face; ++i)
   {
      if(face[i].measure > 0.0) tmp.push_back (face[i]);
   }
   face.resize (0);
   face.resize (tmp.size());
   face = tmp;
   tmp.resize(0);

   n_face = face.size ();
}

//------------------------------------------------------------------------------
// Print triangles and dual cells into file
// Can be visualized with gnuplot
//------------------------------------------------------------------------------
void Grid::print_cells ()
{
   string file_tri = MESH_DATA_DIR+"/tri_", file_dual = MESH_DATA_DIR+"/dual_";
   stringstream ss;
   ss <<get_proc_loc_id();
   file_tri += ss.str() + ".dat"; 
   file_dual += ss.str() + ".dat";
   
   
   if(get_proc_id() == 0)
   {
      cout << "  --- writing triangles into "<<MESH_DATA_DIR<<"/tri*.dat ...\n";
      cout << "  --- writing dual cells into "<<MESH_DATA_DIR<<"/dual*.dat ...\n";
   }   

   ofstream tri(file_tri.c_str());
   ofstream dual(file_dual.c_str());

   for(unsigned int i=0; i<n_cell; ++i)
   {
	  unsigned int n0 = cell[i].vertex[0];
	  unsigned int n1 = cell[i].vertex[1];
	  unsigned int n2 = cell[i].vertex[2];
  
	  tri << vertex[n0].coord.x << "  " << vertex[n0].coord.y << "\n";
	  tri << vertex[n1].coord.x << "  " << vertex[n1].coord.y << "\n";
	  tri << vertex[n2].coord.x << "  " << vertex[n2].coord.y << "\n";
	  tri << vertex[n0].coord.x << "  " << vertex[n0].coord.y << "\n";
	  tri << "\n";
  
	  Vector p01,p12,p20;
	  
	
	  p01.equ(vertex[n0].coord,vertex[n1].coord,0.5,0.5);
	  p12.equ(vertex[n1].coord,vertex[n2].coord,0.5,0.5);
	  p20.equ(vertex[n0].coord,vertex[n2].coord,0.5,0.5);

	  dual << p01.x << "  " << p01.y << "\n";
	  dual << cell[i].centroid.x << "  " << cell[i].centroid.y << "\n";
	  dual << "\n";

	  dual << p12.x << "  " << p12.y << "\n";
	  dual << cell[i].centroid.x << "  " << cell[i].centroid.y << "\n";
	  dual << "\n";

	  dual << p20.x << "  " << p20.y << "\n";
	  dual << cell[i].centroid.x << "  " << cell[i].centroid.y << "\n";
	  dual << "\n";
   }
   for(unsigned int i=0; i<bface.size(); ++i)
   {
      unsigned int v0 = bface[i].vertex[0];
      unsigned int v1 = bface[i].vertex[1];
      dual << vertex[v0].coord.x << "  " << vertex[v0].coord.y << endl;
      dual << vertex[v1].coord.x << "  " << vertex[v1].coord.y << endl;
      dual << endl;
   }
   tri.close();
   dual.close();
}

//------------------------------------------------------------------------------
// Print dual cells areas. Needed for external quadrature
//------------------------------------------------------------------------------
void Grid::print_areas ()
{
   
   string file_mesh = MESH_DATA_DIR+"/area_meshp_";
   stringstream ss;
   ss <<get_proc_id();
   file_mesh += ss.str() + ".dat";
   
   string commandline = "rm -rf " + file_mesh;
   system(commandline.c_str());
   
   mpi_barrier(run_comm);
   
   if(get_proc_id() == 0)
   {
      cout << "  --- writing mesh area (dual) into "<<MESH_DATA_DIR<<"/area_meshp_*.dat ...\n";
   }   

   ofstream meshlist(file_mesh.c_str());
   
   meshlist<<n_vertex<<endl;
   for(unsigned int i=0; i<n_vertex; ++i)
   {
      double nbs = mpi_node_share.find(i)->second;
      meshlist << setprecision(15) << dcarea[i]<< " "<<nbs<<"\n";
   }

   meshlist.close();
}

//------------------------------------------------------------------------------
// Radius for axisymmetric case.
//------------------------------------------------------------------------------
void Grid::compute_radius ()
{
   // If 2d, set all radii to one
   if(dim == two)
   {
      for(unsigned int i=0; i<n_vertex; ++i)
         vertex[i].radius = 1.0;
      for(unsigned int i=0; i<n_face; ++i)
         face[i].radius = 1.0;
      for(unsigned int i=0; i<n_cell; ++i)
         cell[i].radius = 1.0;
      return;
   }

   MPI_DISP("  --- axisymmetric case: computing radii ...",verbose);

   // Radius for vertices
   for(unsigned int i=0; i<n_vertex; ++i)
      vertex[i].radius = vertex[i].coord.x;

   // Radius for faces
   for(unsigned int i=0; i<n_face; ++i)
      face[i].radius = face[i].centroid.x;

   // Radius for cells
   for(unsigned int i=0; i<n_cell; ++i)
      cell[i].radius = cell[i].centroid.x;   
}

//------------------------------------------------------------------------------
// Find dx_max for each vertex
//------------------------------------------------------------------------------
void Grid::compute_dx_max ()
{
   // Initialize dx_max
   dx_max.resize(n_vertex,0.0);
   
   for(unsigned int i=0; i<n_face; ++i) 
   {
     int v0 = face[i].vertex[0];
     int v1 = face[i].vertex[1];
     Vector dr;
     dr.equ(vertex[v1].coord,vertex[v0].coord,1.0,-1.0);
     double len = dr.norm();
     dx_max[v0] = max(dx_max[v0],len);
     dx_max[v1] = max(dx_max[v1],len);
   }  
}

//------------------------------------------------------------------------------
// Mark boundary faces of boundary cells
//------------------------------------------------------------------------------
void Grid::mark_cell_bface (const std::map<int,BoundaryCondition>& bc)
{
   MPI_DISP("  --- marking boundary faces of boundary cells ...",verbose);
   
   for(unsigned int i=0; i<n_cell; ++i)
   {
      
	   for(unsigned j=0;j<3;++j)
	   {
		  int f = cell[i].face[j];
		  if(face[f].type != -1 && face[f].type != -2)
		  {  
			if( (bc.find(face[f].type)->second).type!=BC::periodic)
			{
			   cell[i].nbr_bface.push_back(j);
			}	  
		  }  
	  
	   }     
   }
}

//------------------------------------------------------------------------------
// Remove all ghost cells and nodes (no longer needed for the simulation)
//------------------------------------------------------------------------------
void Grid::remove_ghosts ()
{  
   MPI_DISP("  --- removing ghost data ...",verbose);
   // Deleting areas corresponding to ghost nodes:
   vector<double>::iterator itma = mcarea.begin();
   vector<double>::iterator itda = dcarea.begin();
   for(unsigned int i = 0 ;i<n_vertex; ++i)
   {
      if(vertex[i].active) 
	  {
	  	  itma++;
	  	  itda++;
	  }	   
      else
      {
          mcarea.erase(itma);   
          dcarea.erase(itda);
      }      
   }
   MPI_LOC_ASSERT(mcarea.size()==n_vertex_active);
   MPI_LOC_ASSERT(dcarea.size()==n_vertex_active);
   
   // Changing local vertex indexing in l2g and g2l
   map<int,int>::iterator itl = g2l.begin();
   while(itl!=g2l.end())
   {
	   int loc_v = itl->second;
	   itl->second = active_node_ind.find(loc_v)->second;
	   itl++;   
   }
   
   // Cannot change the key of a map directly
   itl = l2g.begin();
   map <int,int> l2g_new;
   while(itl!=l2g.end())
   {
	   int loc_v = itl->first;
	   int glob_v = itl->second;
	   loc_v = active_node_ind.find(loc_v)->second;
	   l2g_new.insert( pair<int,int>(loc_v,glob_v));
	   itl++;   
   }
   l2g.clear();
   l2g = l2g_new;
   l2g_new.clear();
   MPI_LOC_ASSERT(l2g.size()==g2l.size());

   itl = l2g_mc_pt.begin();
   while(itl!=l2g_mc_pt.end())
   {
	   int loc_v = itl->first;
	   int glob_v = itl->second;
	   loc_v = active_node_ind.find(loc_v)->second;
	   //cout<<itl->first<<"-->"<<loc_v<<endl;
	   l2g_new.insert( pair<int,int>(loc_v,glob_v));
	   itl++;   
   }
   l2g_mc_pt.clear();
   l2g_mc_pt = l2g_new;
   l2g_new.clear();

   // Changing local vertex indexing in mpi_node_share
   itl = mpi_node_share.begin();
   map <int,int> mpi_node_share_new;
   while(itl!=mpi_node_share.end())
   {
	   int loc_v = itl->first;
	   int nbc = itl->second;
	   loc_v = active_node_ind.find(loc_v)->second;
	   mpi_node_share_new.insert( pair<int,int>(loc_v,nbc));
	   itl++;   
   }
   mpi_node_share.clear();
   mpi_node_share = mpi_node_share_new;
   mpi_node_share_new.clear();
   
   // Changing local vertex indexing in periodic list
   if(PERIODIC)
   {
	  for(unsigned int i=0; i<n_periodic; ++i)
	  {
		 for(unsigned int j=0; j<periodic_lists[i].P_ij.size(); ++j)
		 {
			int loc_v = periodic_lists[i].P_ij[j];
			periodic_lists[i].P_ij[j] = active_node_ind.find(loc_v)->second;
		 }
	  }
	  for(unsigned int i=0; i<periodic_mpi_groupings.size(); ++i)
	  {
		 for(unsigned int j=0; j<periodic_mpi_groupings[i].vertices.size(); ++j)
		 {
			int loc_v = periodic_mpi_groupings[i].vertices[j];
			periodic_mpi_groupings[i].vertices[j] = active_node_ind.find(loc_v)->second;
		 }
	  }
   }
   
   // Changing local cell indexing in cell_mpi_bnd and cell_mpi_int
   for(unsigned int i = 0; i < cell_mpi_bnd.size(); ++i)
   {
	   int ind = cell_mpi_bnd[i];
	   cell_mpi_bnd[i] = active_cell_ind.find(ind)->second;   
   }
   for(unsigned int i = 0; i < cell_mpi_int.size(); ++i)
   {
	   int ind = cell_mpi_int[i];
	   cell_mpi_int[i] = active_cell_ind.find(ind)->second;   
   }
   MPI_LOC_ASSERT(cell_mpi_bnd.size()== n_cell_mpi_bnd);
   MPI_LOC_ASSERT(cell_mpi_int.size()== n_cell_mpi_int);
   
   
   // Deleting ghost nodes
   vector<Vertex>::iterator itv = vertex.begin();
   while(itv!=vertex.end())
   {
      if((*itv).active)    
		  itv++;
      else
          vertex.erase(itv);      
   }
   MPI_LOC_ASSERT(vertex.size() == n_vertex_active);
   n_vertex = n_vertex_active;
   
   // Changing node/cell index associated with each face to active node/cell index.
   // NOTE: rcell for interpartition faces of type = -2 now meaningless
   vector<Face>::iterator itf = face.begin();
   while(itf!=face.end())
   {   
      unsigned int v0 = (*itf).vertex[0];
      unsigned int v1 = (*itf).vertex[1];
      (*itf).vertex[0] = active_node_ind.find(v0)->second;
      (*itf).vertex[1] = active_node_ind.find(v1)->second;	
      
      unsigned int lcell = (*itf).lcell;
      unsigned int rcell = (*itf).rcell;
      
      if((*itf).type != -2) // lcell is active for interior and domain bnd faces
         (*itf).lcell = active_cell_ind.find(lcell)->second;    
      else // for mpi boundary faces, it is not known which one of lcell or rcell is active. 
      {   
		 if(!cell[lcell].g_cell)
			(*itf).lcell = active_cell_ind.find(lcell)->second;
		 if(!cell[rcell].g_cell)
			(*itf).rcell = active_cell_ind.find(rcell)->second;    
      }  
      
      if((*itf).type == -1) // rcell is active for interior
         (*itf).rcell = active_cell_ind.find(rcell)->second; 
	  
	  itf++;      
   }
   
   // Changing node index associated with each cell to active cell index. 
   // Also deleting ghost cells
   vector<Cell>::iterator itc = cell.begin();
   while(itc!=cell.end())
   {
      if(!(*itc).g_cell)
      {
		  for(unsigned j=0;j<3;++j)
		  {
			 unsigned int v = (*itc).vertex[j];
			 (*itc).vertex[j] = active_node_ind.find(v)->second;     
		  }     
		  itc++;
      }
      else
          cell.erase(itc);      
   }
   MPI_LOC_ASSERT(cell.size()==n_cell_active);
   n_cell = n_cell_active;
   
   itf = bface.begin();
   while(itf!=bface.end())
   {   
      unsigned int v0 = (*itf).vertex[0];
      unsigned int v1 = (*itf).vertex[1];
      (*itf).vertex[0] = active_node_ind.find(v0)->second;
      (*itf).vertex[1] = active_node_ind.find(v1)->second;	
      
      unsigned int lcell = (*itf).lcell;
      (*itf).lcell = active_cell_ind.find(lcell)->second;
	  
	  itf++;      
   }

   
}

//------------------------------------------------------------------------------
// Creating periodic node weights
//------------------------------------------------------------------------------
void Grid::create_periodic_weights ()
{
   MPI_DISP("  --- creating periodic node weights ...",verbose);
   for(unsigned int i=0; i<n_periodic; ++i)
   {
	  double psize = periodic_lists[i].P_ij.size();
	  double Lsize = periodic_lists[i].L_i.size();
	  double w = 1.0/(psize*Lsize);
	  for(unsigned int j=0; j<psize; ++j)
	  {
		 int v = periodic_lists[i].P_ij[j];
		 periodic_node_weight.insert ( pair<int,double>(v,w) );
	  }
   }
}


//------------------------------------------------------------------------------
// Checking if all structures have been created properly
//------------------------------------------------------------------------------
void Grid::check_data ()
{
   MPI_DISP("  --- verifying data structures ...",verbose);
   // CELL DATA
   for(unsigned int i = 0; i<n_cell; ++i)
   {
        MPI_LOC_ASSERT(!cell[i].g_cell);
        for(unsigned int q = 0; q<3; q++)
        {
           unsigned int v = cell[i].vertex[q];
           MPI_LOC_ASSERT(vertex[v].active);
           if(vertex[v].mpi_bnd_active)
              MPI_LOC_ASSERT(cell[i].mpi_bnd_active);  
        }  
        
        for(unsigned int q = 0; q<3; q++)
        {
           unsigned int f = cell[i].face[q];
           if(face[f].d_bnd_and_shared)
              MPI_LOC_ASSERT(cell[i].mpi_bnd_active); 
        }
   }

   // VERTEX DATA
   for(unsigned int i = 0; i<n_vertex; ++i)
        MPI_LOC_ASSERT(vertex[i].active);
   
   // FACE DATA
   for(unsigned int i = 0; i<n_face; ++i)
   {
        unsigned int v0 = face[i].vertex[0];
        unsigned int v1 = face[i].vertex[1];
        MPI_LOC_ASSERT(vertex[v0].active && vertex[v1].active);
        unsigned int lcell = face[i].lcell;
        unsigned int rcell = face[i].rcell;
        
        if(face[i].type == -1)
        {
           MPI_LOC_ASSERT(!cell[rcell].g_cell);
        }   
        else if(face[i].type == -2)
           MPI_LOC_ASSERT(vertex[v0].mpi_bnd_active && vertex[v1].mpi_bnd_active);
        
        if(face[i].d_bnd_and_shared)
        {
           MPI_LOC_ASSERT(cell[lcell].mpi_bnd_active);
           MPI_LOC_ASSERT(vertex[v0].mpi_bnd_active || vertex[v1].mpi_bnd_active);
        }   
   }
   
   for(unsigned int i = 0; i<bface.size(); ++i)
   {
        MPI_LOC_ASSERT(face[i].type!=-1 && face[i].type!=-2);
        unsigned int v0 = face[i].vertex[0];
        unsigned int v1 = face[i].vertex[1];
        unsigned int lcell = face[i].lcell;
        //cout<<get_proc_id()<<" "<<bface[i].normal.x<<" "<<bface[i].normal.y<<endl;
        //cout<<get_proc_id()<<" ("<<vertex[v0].coord.x<<","<<vertex[v0].coord.y<<
        //                    ") ("<<vertex[v1].coord.x<<","<<vertex[v1].coord.y<<")"<<endl;
        MPI_LOC_ASSERT(vertex[v0].active && vertex[v1].active);
        if(vertex[v0].mpi_bnd_active || vertex[v1].mpi_bnd_active)
        {   
           MPI_LOC_ASSERT(face[i].d_bnd_and_shared); 
           MPI_LOC_ASSERT(cell[lcell].mpi_bnd_active);  
        }
        
        if(face[i].d_bnd_and_shared)
        {
           MPI_LOC_ASSERT(cell[lcell].mpi_bnd_active);
           MPI_LOC_ASSERT(vertex[v0].mpi_bnd_active || vertex[v1].mpi_bnd_active);
        }  
         
        // cout<<vertex[v0].coord.x<<" "<<vertex[v0].coord.y<< " ------"
//             <<vertex[v1].coord.x<<" "<<vertex[v1].coord.y<<endl;
   }
   
   // Emptying data structures no longer needed
   active_node_ind.clear();
   active_cell_ind.clear();
   
}



//------------------------------------------------------------------------------
// Preprocess the grid
//------------------------------------------------------------------------------
void Grid::preproc (const std::map<int,BoundaryCondition>& bc)
{
   make_faces ();
   find_cell_faces ();
   compute_face_centroid ();
   compute_cell_centroid ();
   compute_cell_area ();
   compute_face_normal_and_area ();
   if(cell_type == voronoi)
      remove_empty_faces ();
   remove_ghosts();
   if(debug)
   {
      if(get_proc_id() == 0)
      {
         string commandline0 = "rm -rf " + MESH_DATA_DIR;
		 string commandline1 = "mkdir " + MESH_DATA_DIR;
		 system(commandline0.c_str());
		 system(commandline1.c_str());
      }
      mpi_barrier(MPI_COMM_WORLD);
      if(get_proc_id() < NPART)
      {
         print_cells();
         print_areas();
      }   
   }      
   mark_cell_bface(bc);
   compute_dx_max();
   compute_radius ();
   if(PERIODIC)
      create_periodic_weights ();
   check_data();
   // HACK
   // string point_res = "area";
//    stringstream sss;
//    sss <<get_proc_loc_id();
//    point_res += sss.str() + ".dat";
//    ofstream point_res_file;
//    point_res_file.open (point_res.c_str());
//    for(unsigned int i=0; i<n_vertex; ++i)
//    {
// 	   point_res_file << vertex[i].coord.x << "  "
// 					  << vertex[i].coord.y << "  "
// 					  << dcarea[i] <<endl;
//    }
//    point_res_file.close ();
}
