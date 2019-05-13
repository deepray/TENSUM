#include <algorithm>
#include <cmath>
#include "material.h"

using namespace std;

//------------------------------------------------------------------------------
// SIMPLE average flux
//------------------------------------------------------------------------------
void Material::simple_avg_flux (const PrimVar& left,
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
   double El = 0.5*rhol*left.velocity.square() + left.pressure/(gamma-1);
   double Er = 0.5*rhor*right.velocity.square() + right.pressure/(gamma-1);
   
   
   Flux flux_l, flux_r;
   double vel_normal;
   
   //Left flux
   vel_normal = left.velocity * unit_normal;
   flux_l.mass_flux = rhol*vel_normal;
   flux_l.momentum_flux.equ(unit_normal,left.velocity,left.pressure,flux_l.mass_flux);
   flux_l.energy_flux = (El+left.pressure)*vel_normal;
   
   //Right flux
   vel_normal = right.velocity * unit_normal;
   flux_r.mass_flux = rhor*vel_normal;
   flux_r.momentum_flux.equ(unit_normal,right.velocity,right.pressure,flux_r.mass_flux);
   flux_r.energy_flux = (Er+right.pressure)*vel_normal;

   // central flux
   flux.zero();
   flux.sadd(flux_l,0.5);
   flux.sadd(flux_r,0.5);
                      
   //flux *= area;                   
}
