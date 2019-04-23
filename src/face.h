#ifndef __FACE_H__
#define __FACE_H__

#include "vec.h"

//------------------------------------------------------------------------------
//! Face class 
/*!
  This describes the face data structures
*/
//------------------------------------------------------------------------------
class Face
{
   public:
      Face () {};
      Face (const Face&);            /*!< copy constructor */
      unsigned int vertex[2];        /*!< two vertices of the face */
      int          lcell, 			 /*!< primary cell to the left of the face (always 
                                          exists) */
                   rcell;            /*!< primary cell to the right of the face (exists 
                                          for interior faces) */
      int          lvertex,          /*!< vertex (left) opposite the face (exists 
                                          for interior faces) */
                   rvertex; 		 /*!< vertex (right) opposite the face (exists 
                                          for interior faces) */
      Vector       normal;           /*!< dual face normal */
      Vector       centroid; 		 /*!< face centroid */
      int          type; 			 /*!<  -2 : for inter processor face
       									  
       									   -1 : for interior face
       									  
       									    tag from mesh file for boundary face */ 
     
      double       measure;			 /*!< face length */
      double       radius;		     /*!< used for 3D axi-symmetric flows */
      bool         d_bnd_and_shared; /*!< true if face is domain boundary and shared */
      bool         periodic;         /*!< true if the face is periodic*/

      bool operator== (const Face& face) const; /*!< checks if two faces are the same
      												 by matching their verices */
};

// copy constructor
inline
Face::Face (const Face& face)
{
   vertex[0] = face.vertex[0];
   vertex[1] = face.vertex[1];
   lcell     = face.lcell;
   rcell     = face.rcell;
   lvertex   = face.lvertex;
   rvertex   = face.rvertex;
   normal    = face.normal;
   centroid  = face.centroid;
   type      = face.type;
   measure   = face.measure;
   radius    = face.radius;
   d_bnd_and_shared = face.d_bnd_and_shared;
   periodic = face.periodic;
}

// Check if two faces are same by checking their vertices
inline
bool Face::operator== (const Face& face) const
{

   if ( vertex[0]==face.vertex[0] &&
        vertex[1]==face.vertex[1]) return true;

   if ( vertex[0]==face.vertex[1] &&
        vertex[1]==face.vertex[0]) return true;

   return false;
}

#endif
