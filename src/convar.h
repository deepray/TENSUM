#ifndef __CONVAR_H__
#define __CONVAR_H__

#include "vec.h"
#include "flux.h"

//------------------------------------------------------------------------------
//! Conserved Variable class 
/*!
  This is describes conserved variable as an object, along with various operators
  and functions
*/
//------------------------------------------------------------------------------
class ConVar
{
   public:
      ConVar () {density =0; energy=0; momentum=0;};
      ~ConVar () {};
      
      /*!
	  * Operator copying a conserved variable vector into the current one
	  * @param[in] con_var 
	  */
      ConVar& operator=  (const ConVar& con_var);
      
      /*!
	  * Operator adding a conserved variable vector into the current one
	  * @param[in] con_var 
	  */
      ConVar& operator+= (const ConVar& con_var);
      
      /*!
	  * Operator multiplying the current conserved variable by a scalar
	  * @param[in] scalar 
	  */
      ConVar& operator*= (const double& scalar);
      
      /*!
	  * Function setting the current conserved variable as a linear sum of a conserved 
	  * variable and a Flux variable
	  * @param[in] con 
	  * @param[in] flux
	  * @param[in] c1  
	  * @param[in] c2 
	  */
      ConVar& equ(const ConVar& con, const Flux&   flux, 
                   const double& c1, const double& c2);
                   
      /*!
	  * Function adding a linear sum of a conserved variable and a Flux variable to the
	  * current conserved variable
	  * @param[in] con 
	  * @param[in] flux
	  * @param[in] c1  
	  * @param[in] c2 
	  */             
      ConVar& sadd(const ConVar& con, const Flux&   flux, 
                   const double& c1, const double& c2);             
      
      /*!
	  * Function adding a scaled conserved variable to the current conserved variable
	  * @param[in] con 
	  * @param[in] c1  
	  */
      ConVar& sadd(const ConVar& con, const double& c1);             

      double density, energy;
      Vector momentum;

};

//------------------------------------------------------------------------------
// Assign conserved variable
//------------------------------------------------------------------------------
inline
ConVar& ConVar::operator= (const ConVar& con_var)
{
   density  = con_var.density;
   momentum = con_var.momentum;
   energy   = con_var.energy;

   return *this;
}

//------------------------------------------------------------------------------
// Add conserved variable to given variable
//------------------------------------------------------------------------------
inline
ConVar& ConVar::operator+= (const ConVar& con_var)
{
   density  += con_var.density;
   momentum += con_var.momentum;
   energy   += con_var.energy;

   return *this;
}

//------------------------------------------------------------------------------
// Adding scaled conserved variable and flux into given conserved variable
//------------------------------------------------------------------------------
inline
ConVar& ConVar::equ(const ConVar& con, const Flux&   flux, 
                     const double& c1, const double& c2)
{
   density  = con.density*c1  + flux.mass_flux*c2;
   momentum.equ(con.momentum,flux.momentum_flux,c1,c2);
   energy   = con.energy*c1   + flux.energy_flux*c2;
   
   return *this;
}

//------------------------------------------------------------------------------
// Add scaled conserved variable and flux to  given conserved variable
//------------------------------------------------------------------------------
inline
ConVar& ConVar::sadd(const ConVar& con, const Flux&   flux, 
                     const double& c1, const double& c2)
{
   density  += con.density*c1  + flux.mass_flux*c2;
   momentum.sadd(con.momentum,flux.momentum_flux,c1,c2);
   energy   += con.energy*c1   + flux.energy_flux*c2;
   
   return *this;
}


//------------------------------------------------------------------------------
// Add scaled conserved variable to given conserved variable
//------------------------------------------------------------------------------
inline
ConVar& ConVar::sadd(const ConVar& con, const double& c1)
{
   density  += con.density*c1;
   momentum.sadd(con.momentum,c1);
   energy   += con.energy*c1;
   
   return *this;
}

//------------------------------------------------------------------------------
// Multiply given conserved variable by scalar
//------------------------------------------------------------------------------
inline
ConVar& ConVar::operator *= (const double& scalar)
{
   density  *= scalar;
   momentum *= scalar; 
   energy   *= scalar;

   return *this;
}



#endif
