#ifndef __FLUX_H__
#define __FLUX_H__

#include "vec.h"

//------------------------------------------------------------------------------
//! Flux class 
/*!
  This is describes flux variable as an object, along with various operators
  and functions
*/
//------------------------------------------------------------------------------
class Flux
{
   public:
      double mass_flux;
      Vector momentum_flux;
      double energy_flux;
      
      /*!
	  * Operator adding an flux variable vector into the current one
	  * @param[in] flux
	  */
      Flux& operator+= (const Flux& flux);
      
      /*!
	  * Operator subtracting an flux variable vector into the current one
	  * @param[in] flux
	  */
      Flux& operator-= (const Flux& flux);
      
      /*!
	  * Operator multiplying the current flux variable by a scalar
	  * @param[in] scalar 
	  */
      Flux& operator*= (const double& scalar);
      
      /*!
	  * Function adding a scaled flux variable to the current flux variable
	  * @param[in] flux
	  * @param[in] c1  
	  */
      Flux& sadd(const Flux& flux, const double& scalar);
      
      /*!
	  * Function adding a linear sum of a three flux variables to the current one.
	  * @param[in] fl1
	  * @param[in] fl2
	  * @param[in] fl3
	  * @param[in] c1  
	  * @param[in] c2 
	  * @param[in] c3
	  */
      Flux& sadd(const Flux& fl1, const Flux& fl2, const Flux& fl3,
                 const double& c1, const double& c2, const double& c3 ); // needed in LUSGS
      
      /*!
	  * Function setting the current flux variable to zero.
	  */
      void zero ();
};

//------------------------------------------------------------------------------
// Set all flux components to zero
//------------------------------------------------------------------------------
inline
void Flux::zero ()
{
   mass_flux     = 0.0;
   momentum_flux = 0.0;
   energy_flux   = 0.0;
}

//------------------------------------------------------------------------------
// Add flux to given flux
//------------------------------------------------------------------------------
inline
Flux& Flux::operator+= (const Flux& flux)
{
   mass_flux     += flux.mass_flux;
   momentum_flux += flux.momentum_flux;
   energy_flux   += flux.energy_flux;

   return *this;
}

//------------------------------------------------------------------------------
// Subtract flux from given flux
//------------------------------------------------------------------------------
inline
Flux& Flux::operator-= (const Flux& flux)
{
   mass_flux     -= flux.mass_flux;
   momentum_flux -= flux.momentum_flux;
   energy_flux   -= flux.energy_flux;

   return *this;
}

//------------------------------------------------------------------------------
// Multiply given flux by a scalar
//------------------------------------------------------------------------------
inline
Flux& Flux::operator*= (const double& scalar)
{
   mass_flux     *= scalar;
   momentum_flux *= scalar;
   energy_flux   *= scalar;

   return *this;
}

//------------------------------------------------------------------------------
// Adding a scaled flux to the given flux
//------------------------------------------------------------------------------
inline
Flux& Flux::sadd(const Flux& flux, const double& scalar)
{
   mass_flux     += flux.mass_flux*scalar;
   momentum_flux.sadd(flux.momentum_flux,scalar);
   energy_flux   += flux.energy_flux*scalar;

   return *this;
}

//------------------------------------------------------------------------------
// Adding scaled fluxes to the given flux
//------------------------------------------------------------------------------
inline
Flux& Flux::sadd(const Flux& fl1, const Flux& fl2, const Flux& fl3,
             const double& c1, const double& c2, const double& c3 )
{
   mass_flux     += fl1.mass_flux*c1 + fl2.mass_flux*c2 + fl3.mass_flux*c3;
   momentum_flux.sadd(fl1.momentum_flux,fl2.momentum_flux,fl3.momentum_flux,
                      c1,c2,c3);
   energy_flux   += fl1.energy_flux*c1 + fl2.energy_flux*c2 + fl3.energy_flux*c3;

   return *this;
}


#endif
