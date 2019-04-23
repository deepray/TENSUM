#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include "material.h"
#include "constants.h"
#include "mpi_utils.h"

using namespace std;

extern Dimension dim;

//------------------------------------------------------------------------------
// Do some initializations
//------------------------------------------------------------------------------
void Material::initialize ()
{
   Cp = gamma * gas_const / (gamma - 1.0);

   if(model == euler && mu_ref > 0.0)
         MPI_ERR("read_material:: mu_ref must be 0 for the model=euler:");      
}

//------------------------------------------------------------------------------
// Numerical flux function
//------------------------------------------------------------------------------
void Material::num_flux (const Vector& cl,
                         const Vector& cr,
                         const PrimVar& left0, 
                         const PrimVar& right0, 
                         const PrimVar& left, 
                         const PrimVar& right, 
                         const std::vector<Vector>& dE_l,
                         const std::vector<Vector>& dE_r,
                         const Vector& normal, 
                         const FluxData& data,
                         Flux& flux) const
{
   switch (flux_scheme)
   {
      case kepes_tecno_roe:
         kepes_tecno_roe_flux (cl,cr,left0, right0, dE_l, dE_r, normal, data, flux);
         break;
      case kepes_tecno_rusanov:
         kepes_tecno_rusanov_flux (cl,cr,left0, right0, dE_l, dE_r, normal, data, flux);
         break;   
      case kepec:
         kepec_flux (left0, right0, normal, data, flux);
         break;
      case kep:
         kep_flux (left0, right0, normal, data, flux);
         break;
      case roe_ec:
         roe_ec_flux (left0, right0, normal, data, flux);
         break;  
      case avg:
         simple_avg_flux (left0, right0, normal, data, flux);
         break;       
      default:
         cout << "num_flux: unknown flux " << flux_scheme << endl;
         exit (0);
   }
   
   // if(order_correction)
//       central_flux_correction(cl,cr,left0, right0, dE_l, dE_r, normal, data, flux);
   double area = normal.norm();
   
   // cout<<flux.mass_flux<<" "
//        <<flux.momentum_flux.x<<" "
//        <<flux.momentum_flux.y<<" "
//        <<flux.momentum_flux.z<<" "
//        <<flux.energy_flux<<endl;
   flux *= area;
}

//------------------------------------------------------------------------------
// Flux on slip walls
//------------------------------------------------------------------------------
void Material::slip_flux (const PrimVar& state,
                          const Vector&  normal,
                          Flux&          flux) const
{
   flux.mass_flux     = 0.0;
   flux.momentum_flux.equ(normal,state.pressure);
   flux.energy_flux   = 0.0;
}

//------------------------------------------------------------------------------
// Euler Flux Calculation
//------------------------------------------------------------------------------
void Material::euler_flux (const PrimVar& prim, 
                           const Vector&  normal,
                           Flux&          flux) const
{
   // Enthalpy 
   double h  = gamma * gas_const * prim.temperature / (gamma-1.0) + 
               0.5 * prim.velocity.square();
   // Normal velocity
   double vn = prim.velocity * normal;
   double density = Density (prim);

   flux.mass_flux = density * vn;
   flux.momentum_flux.equ(normal,prim.velocity,prim.pressure,flux.mass_flux);
   flux.energy_flux = h * flux.mass_flux;
}

//------------------------------------------------------------------------------
// Inviscid boundary flux
//------------------------------------------------------------------------------
void Material::inviscid_bflux (const PrimVar& prim,
                               const PrimVar& primb,
                               const Vector&  normal,
                               Flux& flux) const
{
   // Normal velocity
   double vnb = primb.velocity * normal;
   double vn = prim.velocity*normal;
   double density = Density (prim);

   flux.mass_flux = density * vnb;
   // flux.momentum_flux = normal * primb.pressure + 
//                         prim.velocity * flux.mass_flux;
   flux.momentum_flux.equ(normal,prim.velocity,primb.pressure,flux.mass_flux);
   flux.energy_flux = (0.5*density*prim.velocity.square() + 
                       gamma * prim.pressure / (gamma-1.0))*vnb + 
                      (primb.pressure - prim.pressure)*vn;
}

// NEW METHOD: TO BE TESTED
//void Material::inviscid_bflux (const PrimVar& prim,
                              //const Vector& velb,
                              //const double& pb,
                              //const double& etab, 
                              //const Vector&  normal,
                              //Flux& flux) const
//{
   //// Normal velocity
   //double vnb = velb * normal;
   //double vn = prim.velocity*normal;
   //double density = Density (prim);
   //double eta = -Density(prim)*Entropy(prim)/(gamma - 1.0); 

   //flux.mass_flux = density * vnb;
   //flux.momentum_flux = normal * pb + 
                        //prim.velocity * flux.mass_flux;
   //flux.energy_flux = (0.5*density*prim.velocity.square() + 
                       //gamma * prim.pressure / (gamma-1.0))*vnb + 
                      //(pb - prim.pressure)*vn +
                      //(eta - etab)*vn*prim.pressure/Density(prim);
//}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// viscous flux: 
// This is for interior faces.
//------------------------------------------------------------------------------
void Material::viscous_flux (const EntVar& state, 
                             const std::vector<Vector>& dE,
                             const Vector&  normal0, Flux& flux0,
                             const Vector&  normal1, Flux& flux1,
                             const Vector&  normal2, Flux& flux2
                             ) const
{
   
   double mu = viscosity (-1/(gas_const*state.entl));
   double k = mu * Cp / prandtl;
   
   if(!visc_active)
      mu = 0.0;
   if(!heat_active)
      k = 0.0;
   
   // Terms appearing in matrices K_ij
   double e1 = state.entv.x/state.entl;
   double e2 = state.entv.y/state.entl;
   double e3 = state.entv.z/state.entl;
   double e4 = 1.0/state.entl;
   double c1 = 4.0*mu/3.0;
   double c2 = k/gas_const;
   double c3 = 2*mu/3.0;
   
   // Matrix gradient products (Matrices do not include extra e4 factor) 
   double K11dV1[5], K12dV2[5], K21dV1[5], K22dV2[5];
   K11dV1[0] = 0.0;
   K11dV1[1] = c1*(-dE[1].x + e1*dE[4].x); 
   K11dV1[2] = mu*(-dE[2].x + e2*dE[4].x);
   K11dV1[3] = mu*(-dE[3].x + e3*dE[4].x);
   K11dV1[4] = c1*e1*dE[1].x + mu*(e2*dE[2].x + e3*dE[3].x) + 
               (-c1*e1*e1 - mu*(e2*e2 + e3*e3) + c2*e4)*dE[4].x; 
   
   K12dV2[0] = 0.0;
   K12dV2[1] = c3*(dE[2].y - e2*dE[4].y); 
   K12dV2[2] = mu*(-dE[1].y + e1*dE[4].y);
   K12dV2[3] = 0.0; 
   K12dV2[4] = mu*e2*dE[1].y - c3*e1*dE[2].y -(mu-c3)*e1*e2*dE[4].y;
   
   K21dV1[0] = 0.0;
   K21dV1[1] = mu*(-dE[2].x + e2*dE[4].x); 
   K21dV1[2] = c3*(dE[1].x - e1*dE[4].x);
   K21dV1[3] = 0.0;
   K21dV1[4] = -c3*e2*dE[1].x +mu*e1*dE[2].x -(mu - c3)*e1*e2*dE[4].x;
   
   K22dV2[0] = 0.0;
   K22dV2[1] = mu*(-dE[1].y + e1*dE[4].y); 
   K22dV2[2] = c1*(-dE[2].y + e2*dE[4].y);
   K22dV2[3] = mu*(-dE[3].y + e3*dE[4].y);
   K22dV2[4] = mu*(e1*dE[1].y + e3*dE[3].y) + c1*e2*dE[2].y + 
               (-mu*(e1*e1 + e3*e3) - c1*e2*e2 + c2*e4)*dE[4].y; 

   // Evaluating viscous flux contributions
   flux0.mass_flux       = 0.0;
   flux0.momentum_flux.x = -((K11dV1[1] + K12dV2[1])* normal0.x + (K21dV1[1] + K22dV2[1])* normal0.y)*e4;
   flux0.momentum_flux.y = -((K11dV1[2] + K12dV2[2])* normal0.x + (K21dV1[2] + K22dV2[2])* normal0.y)*e4;
   flux0.momentum_flux.z = -((K11dV1[3] + K12dV2[3])* normal0.x + (K21dV1[3] + K22dV2[3])* normal0.y)*e4;
   flux0.energy_flux     = -((K11dV1[4] + K12dV2[4])* normal0.x + (K21dV1[4] + K22dV2[4])* normal0.y)*e4;

   flux1.mass_flux       = 0.0;
   flux1.momentum_flux.x = -((K11dV1[1] + K12dV2[1])* normal1.x + (K21dV1[1] + K22dV2[1])* normal1.y)*e4;
   flux1.momentum_flux.y = -((K11dV1[2] + K12dV2[2])* normal1.x + (K21dV1[2] + K22dV2[2])* normal1.y)*e4;
   flux1.momentum_flux.z = -((K11dV1[3] + K12dV2[3])* normal1.x + (K21dV1[3] + K22dV2[3])* normal1.y)*e4;
   flux1.energy_flux     = -((K11dV1[4] + K12dV2[4])* normal1.x + (K21dV1[4] + K22dV2[4])* normal1.y)*e4;
   
   flux2.mass_flux       = 0.0;
   flux2.momentum_flux.x = -((K11dV1[1] + K12dV2[1])* normal2.x + (K21dV1[1] + K22dV2[1])* normal2.y)*e4;
   flux2.momentum_flux.y = -((K11dV1[2] + K12dV2[2])* normal2.x + (K21dV1[2] + K22dV2[2])* normal2.y)*e4;
   flux2.momentum_flux.z = -((K11dV1[3] + K12dV2[3])* normal2.x + (K21dV1[3] + K22dV2[3])* normal2.y)*e4;
   flux2.energy_flux     = -((K11dV1[4] + K12dV2[4])* normal2.x + (K21dV1[4] + K22dV2[4])* normal2.y)*e4;

}

//------------------------------------------------------------------------------
// viscous flux:
// This is for boundary faces only.
//------------------------------------------------------------------------------
void Material::viscous_bflux (const EntVar& state, 
                              const bool&    adiabatic,
                              const std::vector<Vector>& dE,
                              const Vector&  normal, Flux& flux
                             ) const
{
   double mu = viscosity (-1/(gas_const*state.entl));
   double k = mu * Cp / prandtl;
   
   if(!visc_active)
      mu = 0.0;
   if(!heat_active)
      k = 0.0;

   // Terms appearing in matrices K_ij
   double e1 = state.entv.x/state.entl;
   double e2 = state.entv.y/state.entl;
   double e3 = state.entv.z/state.entl;
   double e4 = 1.0/state.entl;
   double c1 = 4.0*mu/3.0;
   double c2 = k/gas_const;
   double c3 = 2*mu/3.0;
   
   // Correction term for no-slip conditions if boundary heat-flux mentioned
   double Cor;
   if(adiabatic)
      Cor = -1.0*(0.0 + c2*e4*e4*(dE[4].x*normal.x + dE[4].y*normal.y));
   else
      Cor = 0.0;
   
   // Matrix gradient products
   double K11dV1[5], K12dV2[5], K21dV1[5], K22dV2[5];
   K11dV1[0] = 0.0;
   K11dV1[1] = c1*(-dE[1].x + e1*dE[4].x); 
   K11dV1[2] = mu*(-dE[2].x + e2*dE[4].x);
   K11dV1[3] = mu*(-dE[3].x + e3*dE[4].x);
   K11dV1[4] = c1*e1*dE[1].x + mu*(e2*dE[2].x + e3*dE[3].x) + 
               (-c1*e1*e1 - mu*(e2*e2 + e3*e3) + c2*e4)*dE[4].x; 
   
   K12dV2[0] = 0.0;
   K12dV2[1] = c3*(dE[2].y - e2*dE[4].y); 
   K12dV2[2] = mu*(-dE[1].y + e1*dE[4].y);
   K12dV2[3] = 0.0; 
   K12dV2[4] = mu*e2*dE[1].y - c3*e1*dE[2].y -(mu-c3)*e1*e2*dE[4].y;
   
   K21dV1[0] = 0.0;
   K21dV1[1] = mu*(-dE[2].x + e2*dE[4].x); 
   K21dV1[2] = c3*(dE[1].x - e1*dE[4].x);
   K21dV1[3] = 0.0;
   K21dV1[4] = -c3*e2*dE[1].x +mu*e1*dE[2].x -(mu - c3)*e1*e2*dE[4].x;
   
   K22dV2[0] = 0.0;
   K22dV2[1] = mu*(-dE[1].y + e1*dE[4].y); 
   K22dV2[2] = c1*(-dE[2].y + e2*dE[4].y);
   K22dV2[3] = mu*(-dE[3].y + e3*dE[4].y);
   K22dV2[4] = mu*(e1*dE[1].y + e3*dE[3].y) + c1*e2*dE[2].y + 
               (-mu*(e1*e1 + e3*e3) - c1*e2*e2 + c2*e4)*dE[4].y;
   

   // Evaluating viscous flux contributions
   flux.mass_flux       = 0.0;
   flux.momentum_flux.x = -((K11dV1[1] + K12dV2[1])* normal.x + (K21dV1[1] + K22dV2[1])* normal.y)*e4;
   flux.momentum_flux.y = -((K11dV1[2] + K12dV2[2])* normal.x + (K21dV1[2] + K22dV2[2])* normal.y)*e4;
   flux.momentum_flux.z = -((K11dV1[3] + K12dV2[3])* normal.x + (K21dV1[3] + K22dV2[3])* normal.y)*e4;
   flux.energy_flux     = -((K11dV1[4] + K12dV2[4])* normal.x + (K21dV1[4] + K22dV2[4])* normal.y)*e4 - Cor;

}
