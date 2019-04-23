#ifndef __VEC_H__
#define __VEC_H__

#include <cmath>

//------------------------------------------------------------------------------
//! Vector class 
/*!
  This is describes 3D vcetor as an object, along with various operators
  and functions
*/
//------------------------------------------------------------------------------
class Vector
{
   public:
      double x, y, z;
      
      /*!
	  * Operator assigning one vector to the current one
	  * @param[in] rhs
	  */
      Vector& operator=  (const Vector& rhs);
      
      /*!
	  * Operator assigning a scalar value all components of the current vector
	  * @param[in] scalar
	  */
      Vector& operator=  (const double& scalar);
      
      /*!
	  * Operator adding a vector to the current one
	  * @param[in] rhs
	  */
      Vector& operator+= (const Vector& rhs);
      
      /*!
	  * Operator subtracting a vector to the current one
	  * @param[in] rhs
	  */
      Vector& operator-= (const Vector& rhs);
      
      /*!
	  * Operator multiplying the current vector by a contstant
	  * @param[in] scalar
	  */
      Vector& operator*= (const double& scalar);
      
      /*!
	  * Operator dividing the current vector by a contstant
	  * @param[in] scalar
	  */
      Vector& operator/= (const double& scalar);
      
      /*!
	  * Function returning the dot product of a vector and the current one
	  * @param[in] vec
	  */
      double  operator*  (const Vector& vec) const; 
      
      /*!
	  * Function assigning the cross product of two vectors to the current one
	  * @param[in] v1
	  * @param[in] v2
	  */
      Vector& cprod (const Vector& v1,const Vector& v2); 
      
      /*!
	  * Function assigning a scaled vector to the current one
	  * @param[in] v1
	  * @param[in] c1
	  */
      Vector& equ(const Vector& v1, const double& c1);
      
      /*!
	  * Function assigning a linear combination of two vector to the current one
	  * @param[in] v1
	  * @param[in] v2
	  * @param[in] c1
	  * @param[in] c2
	  */
      Vector& equ(const Vector& v1,const Vector& v2, const double& c1,const double& c2); 
      
      /*!
	  * Function assigning a linear combination of three vector to the current one
	  * @param[in] v1
	  * @param[in] v2
	  * @param[in] v3
	  * @param[in] c1
	  * @param[in] c2
	  * @param[in] c3
	  */
      Vector& equ(const Vector& v1,const Vector& v2, const Vector& v3,
                   const double& c1,const double& c2, const double& c3);
      
      /*!
	  * Function adding a scaled vector to the current one
	  * @param[in] v1
	  * @param[in] c1
	  */
      Vector& sadd(const Vector& v1, const double& c1);   
      
      /*!
	  * Function adding the square (component-wise) of a vector to the current one
	  * @param[in] v
	  */
      Vector& sqadd(const Vector& v);  
      
      /*!
	  * Function adding a linear combination of two vector to the current one
	  * @param[in] v1
	  * @param[in] v2
	  * @param[in] c1
	  * @param[in] c2
	  */        
      Vector& sadd(const Vector& v1,const Vector& v2, const double& c1,const double& c2); 
      
      /*!
	  * Function adding a linear combination of three vector to the current one
	  * @param[in] v1
	  * @param[in] v2
	  * @param[in] v3
	  * @param[in] c1
	  * @param[in] c2
	  * @param[in] c3
	  */
      Vector& sadd(const Vector& v1,const Vector& v2, const Vector& v3,
                   const double& c1,const double& c2, const double& c3);
      
      /*!
	  * Function returning squared norm of the current vector
	  */
      double  square () const;
      
      /*!
	  * Function returning norm of the current vector
	  */
      double  norm () const;
};

// Assign one vector to another
inline
Vector& Vector::operator= (const Vector& rhs){
   x = rhs.x;
   y = rhs.y;
   z = rhs.z;

   return *this;
}

// Assign one vector to another
inline
Vector& Vector::operator= (const double& scalar){
   x = scalar;
   y = scalar;
   z = scalar;

   return *this;
}

// Add vector to given vector: this = this + rhs
inline
Vector& Vector::operator+= (const Vector& rhs){
   x += rhs.x;
   y += rhs.y;
   z += rhs.z;

   return *this;
}

// Subtract vector from given vector: this = this - rhs
inline
Vector& Vector::operator-= (const Vector& rhs){
   x -= rhs.x;
   y -= rhs.y;
   z -= rhs.z;

   return *this;
}

// Multiply vector by scalar and copy result to same vector
inline
Vector& Vector::operator*= (const double& scalar){
   x *= scalar;
   y *= scalar;
   z *= scalar;

   return *this;
}

// Divide vector by scalar and copy result to same vector
inline
Vector& Vector::operator/= (const double& scalar){
   x /= scalar;
   y /= scalar;
   z /= scalar;

   return *this;
}

// L2 norm square of vector
inline
double Vector::square () const
{
   return x*x + y*y + z*z;
}

// L2 norm of vector
inline
double Vector::norm () const
{
   return std::sqrt(x*x + y*y + z*z);
}

// Dot product of two vectors
inline
double Vector::operator* (const Vector& vec) const
{
   return x * vec.x + y * vec.y + z * vec.z;
}


// Cross product of two vectors
inline
Vector& Vector::cprod (const Vector& v1,const Vector& v2)
{

   x = v1.y * v2.z - v1.z * v2.y;
   y = v1.z * v2.x - v1.x * v2.z;
   z = v1.x * v2.y - v1.y * v2.x;

   return *this;
}

// Scaled addition of a vector and storage directly into another vector
inline
Vector& Vector::equ(const Vector& v1, const double& c1){
   x = v1.x*c1;
   y = v1.y*c1;
   z = v1.z*c1;

   return *this;
}


// Scaled addition of two vectors and storage directly into another vector
inline
Vector& Vector::equ(const Vector& v1,const Vector& v2, const double& c1,const double& c2){
   x = v1.x*c1 + v2.x*c2;
   y = v1.y*c1 + v2.y*c2;
   z = v1.z*c1 + v2.z*c2;

   return *this;
}

// Scaled addition of three vectors and storage directly into another vector
inline
Vector& Vector::equ(const Vector& v1,const Vector& v2, const Vector& v3,
                     const double& c1,const double& c2, const double& c3)
{                     
   x = v1.x*c1 + v2.x*c2 + v3.x*c3;
   y = v1.y*c1 + v2.y*c2 + v3.y*c3;
   z = v1.z*c1 + v2.z*c2 + v3.z*c3;

   return *this;
}

// Scaled and add a vector to an existing vector
inline
Vector& Vector::sadd(const Vector& v1, const double& c1){
   x += v1.x*c1;
   y += v1.y*c1;
   z += v1.z*c1;;

   return *this;
}


// Scaled addition of two vectors to an existing vector
inline
Vector& Vector::sadd(const Vector& v1,const Vector& v2, const double& c1,const double& c2){
   x += v1.x*c1 + v2.x*c2;
   y += v1.y*c1 + v2.y*c2;
   z += v1.z*c1 + v2.z*c2;

   return *this;
}

// Scaled addition of three vectors to an existing vector
inline
Vector& Vector::sadd(const Vector& v1,const Vector& v2, const Vector& v3,
                     const double& c1,const double& c2, const double& c3)
{                     
   x += v1.x*c1 + v2.x*c2 + v3.x*c3;
   y += v1.y*c1 + v2.y*c2 + v3.y*c3;
   z += v1.z*c1 + v2.z*c2 + v3.z*c3;

   return *this;
}

// Squared (component-wise) of three vectors to an existing vector
inline
Vector& Vector::sqadd(const Vector& v)
{                     
   x += v.x*v.x;
   y += v.y*v.y;
   z += v.z*v.z;

   return *this;
}

#endif
