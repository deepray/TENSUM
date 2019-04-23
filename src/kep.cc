#include <algorithm>
#include <cmath>
#include "material.h"

using namespace std;

//------------------------------------------------------------------------------
// KEP flux
//------------------------------------------------------------------------------
void Material::kep_flux (const PrimVar& left,
						   const PrimVar& right,
						   const Vector& normal,
						   const FluxData& data,
						   Flux& flux) const
{  
   double area = normal.norm();
   Vector unit_normal;
   unit_normal.equ(normal, 1.0/area);
   double rhol = Density(left);
   double rhor = Density(right);
   double Hl = pow(sound_speed(left),2.0)/(gamma-1.0) + 0.5*left.velocity.square();
   double Hr = pow(sound_speed(right),2.0)/(gamma-1.0) + 0.5*right.velocity.square();
   
   double p = 0.5*(left.pressure + right.pressure);
   double rho = 0.5*( rhol + rhor );
   double H = 0.5*(Hl+Hr);
   
   Vector vel;
   vel.equ(left.velocity,right.velocity,0.5,0.5);
   double vel_normal = vel * unit_normal;

   // central flux
   flux.mass_flux = rho * vel_normal;
   flux.momentum_flux.equ(unit_normal,vel,p,flux.mass_flux);
   flux.energy_flux = rho*H*vel_normal;
                      
   //flux *= area;                   
}
