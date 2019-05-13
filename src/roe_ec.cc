#include <algorithm>
#include <cmath>
#include "material.h"

using namespace std;

//------------------------------------------------------------------------------
// ROE-EC flux
//------------------------------------------------------------------------------
void Material::roe_ec_flux (const PrimVar& left,
						   const PrimVar& right,
						   const Vector& normal,
						   //const FluxData& data,
						   Flux& flux) const
{  
   double area = normal.norm();
   Vector unit_normal;
   unit_normal.equ(normal, 1.0/area);
   double rhol = Density(left);
   double rhor = Density(right);
   double pl   = left.pressure;
   double pr   = right.pressure;
   double bl   = sqrt(rhol/pl);   
   double br   = sqrt(rhor/pr);
   
   // arithmetic average of parameter vector
   double z1b = 0.5*(bl+br);
   Vector zvb;
   zvb.equ(left.velocity,right.velocity,0.5*bl,0.5*br);
   double z5b = 0.5*(bl*pl + br*pr);
   double zvb_normal = zvb * unit_normal;
   
   //logarithmic avg of parameter vector
   double z1l = logavg( bl, br );
   double z5l = logavg( bl*pl, br*pr );

   // central flux
   flux.mass_flux = z5l * zvb_normal;
   flux.momentum_flux.equ(unit_normal,zvb,z5b/z1b,flux.mass_flux/z1b);
   flux.energy_flux = 0.5*( (gamma+1.0)/(gamma-1.0)* flux.mass_flux/z1l + 
                      flux.momentum_flux * zvb)/z1b;
                      
   //flux *= area;                   
}
