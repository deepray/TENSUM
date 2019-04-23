#ifndef __PRIMVAR_H__
#define __PRIMVAR_H__

#include "vec.h"

//------------------------------------------------------------------------------
//! Primitive Variable class 
/*!
  This is describes primitive variable as an object, along with various operators
  and functions
*/
//------------------------------------------------------------------------------
class PrimVar
{
   public:
      PrimVar(){temperature=0; pressure=0; velocity=0;};
      double temperature, pressure;
      Vector velocity;

      /*!
	  * Function saving the component-wise product of two primitive variables into the 
	    current one
	  * @param[in] prim_var1
	  * @param[in] prim_var2
	  */
      PrimVar&  cdot(const PrimVar& prim_var1, const PrimVar& prim_var2); // componentwise multi
      
      /*!
	  * Function adding the component-wise product of two primitive variables into the 
	    current one
	  * @param[in] prim_var1
	  * @param[in] prim_var2
	  */
      PrimVar&  addcdot(const PrimVar& prim_var1, const PrimVar& prim_var2); // componentwise multi
      
      /*!
	  * Operator multiplying the current primitive variable by a scalar
	  * @param[in] scalar 
	  */
      PrimVar& operator*= (const double& scalar);
      
      /*!
	  * Operator assigning a scalar value to all primitive variable components
	  * @param[in] scalar 
	  */
      PrimVar& operator=  (const double& scalar);
      
      /*!
	  * Operator adding a primitive variable vector into the current one
	  * @param[in] prim_var
	  */
      PrimVar& operator+= (const PrimVar& prim_var);
      
      /*!
	  * Function setting the current primitive variable as a linear sum of a two 
	  * entropy variables
	  * @param[in] prim1
	  * @param[in] prim2
	  * @param[in] c1  
	  * @param[in] c2 
	  */
      PrimVar& equ(const PrimVar& prim1, const PrimVar& prim2, 
                   const double& c1, const double& c2);
                   
      /*!
	  * Function setting the current primitive variable as the absolute difference of a two 
	  * entropy variables
	  * @param[in] prim1
	  * @param[in] prim2
	  */
      PrimVar& absdiff(const PrimVar& prim1, const PrimVar& prim2);             
      
      
      /*!
	  * Function adding the squared (component-wise) primitive variables into the given one 
	  * @param[in] prim
	  */ 
      PrimVar& sqadd(const PrimVar& prim);  
      
      /*!
	  * Function adding the scaled squared (component-wise) primitive variables into the given one 
	  * @param[in] prim
	  * @param[in] c1
	  */ 
      PrimVar& sqsadd(const PrimVar& prim, const double& c1);      
      
      /*!
	  * Function adding a scaled primitive variable to the current primitive variable
	  * @param[in] prim
	  * @param[in] c1  
	  */
      PrimVar& sadd(const PrimVar& prim, const double& c1);        
      
      /*!
	  * Function setting the current primitive variable as the minimum (component-wise) 
	  * of between itself and another primitive variable
	  * @param[in] p
	  */
      void min (const PrimVar& p);
      
      /*!
	  * Function setting the current primitive variable as the maximum (component-wise) 
	  * of between itself and another primitive variable
	  * @param[in] p
	  */
      void max (const PrimVar& p);
};


//------------------------------------------------------------------------------
// Multiply given primitive by scalar
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::operator*= (const double& scalar)
{
   temperature *= scalar;
   velocity    *= scalar; 
   pressure    *= scalar;

   return *this;
}

//------------------------------------------------------------------------------
// Set a scalar value
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::operator= (const double& scalar)
{
   temperature = scalar;
   velocity    = scalar; 
   pressure    = scalar;

   return *this;
}

//------------------------------------------------------------------------------
// Add primitive variable to given primitive variable
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::operator+= (const PrimVar& prim_var)
{
   temperature += prim_var.temperature;
   velocity    += prim_var.velocity;
   pressure    += prim_var.pressure;

   return *this;
}

//------------------------------------------------------------------------------
// Update PrimVar = min(PrimVar, p)
//------------------------------------------------------------------------------
inline
void PrimVar::min (const PrimVar& p)
{
   temperature = std::min(temperature, p.temperature);
   velocity.x  = std::min(velocity.x,  p.velocity.x);
   velocity.y  = std::min(velocity.y,  p.velocity.y);
   velocity.z  = std::min(velocity.z,  p.velocity.z);
   pressure    = std::min(pressure,    p.pressure);
}

//------------------------------------------------------------------------------
// Update PrimVar = max(PrimVar, p)
//------------------------------------------------------------------------------
inline
void PrimVar::max (const PrimVar& p)
{
   temperature = std::max(temperature, p.temperature);
   velocity.x  = std::max(velocity.x,  p.velocity.x);
   velocity.y  = std::max(velocity.y,  p.velocity.y);
   velocity.z  = std::max(velocity.z,  p.velocity.z);
   pressure    = std::max(pressure,    p.pressure);
}

//------------------------------------------------------------------------------
// Set current primitive variable vector as component-wise product of two primitive vectors
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::cdot(const PrimVar& prim_var1, const PrimVar& prim_var2)
{

   temperature = prim_var1.temperature * prim_var2.temperature;
   velocity.x  = prim_var1.velocity.x  * prim_var2.velocity.x;
   velocity.y  = prim_var1.velocity.y  * prim_var2.velocity.y;
   velocity.z  = prim_var1.velocity.z  * prim_var2.velocity.z;
   pressure    = prim_var1.pressure    * prim_var2.pressure;

   return *this;
}


//------------------------------------------------------------------------------
// Add linear sum of two primitive vectors into the current one
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::addcdot(const PrimVar& prim_var1, const PrimVar& prim_var2)
{

   temperature += prim_var1.temperature * prim_var2.temperature;
   velocity.x  += prim_var1.velocity.x  * prim_var2.velocity.x;
   velocity.y  += prim_var1.velocity.y  * prim_var2.velocity.y;
   velocity.z  += prim_var1.velocity.z  * prim_var2.velocity.z;
   pressure    += prim_var1.pressure    * prim_var2.pressure;

   return *this;
}


//------------------------------------------------------------------------------
// Set current primitive variable vector as linear sum of two primitive vectors
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::equ(const PrimVar& prim1, const PrimVar&  prim2, 
                     const double& c1, const double& c2)
{
   temperature  = prim1.temperature*c1  + prim2.temperature*c2;
   velocity.equ(prim1.velocity,prim2.velocity,c1,c2);
   pressure     = prim1.pressure*c1   + prim2.pressure*c2;
   
   return *this;
}

//------------------------------------------------------------------------------
// Absolute difference of two primitive variables
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::absdiff(const PrimVar& prim1, const PrimVar&  prim2)
{
   temperature  = fabs(prim1.temperature - prim2.temperature);
   velocity.x   = fabs(prim1.velocity.x - prim2.velocity.x);
   velocity.y   = fabs(prim1.velocity.y - prim2.velocity.y);
   velocity.z   = fabs(prim1.velocity.z - prim2.velocity.z);
   pressure     = fabs(prim1.pressure - prim2.pressure);
   
   return *this;
}

//------------------------------------------------------------------------------
// Adding the squared (component-wise) primitive variables into  given primitive variable, 
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::sqadd(const PrimVar& prim)
{
   temperature  += prim.temperature*prim.temperature;
   velocity.sqadd(prim.velocity);
   pressure     += prim.pressure*prim.pressure;
   
   return *this;
}

//------------------------------------------------------------------------------
// Adding the scaled squared (component-wise) primitive variables into  given primitive variable, 
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::sqsadd(const PrimVar& prim, const double& c1)
{
   temperature  += prim.temperature*prim.temperature*c1;
   velocity.sqsadd(prim.velocity,c1);
   pressure     += prim.pressure*prim.pressure*c1;
   
   return *this;
}

//------------------------------------------------------------------------------
// Add scaled primitive variable to given primitive variable, 
//------------------------------------------------------------------------------
inline
PrimVar& PrimVar::sadd(const PrimVar& prim, const double& c1)
{
   temperature  += prim.temperature*c1;
   velocity.sadd(prim.velocity,c1);
   pressure     += prim.pressure*c1;
   
   return *this;
}

#endif
