#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include <iostream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
#include "vec.h"
#include "primvar.h"
#include "flux.h"
#include "convar.h"
#include "entvar.h"
#include "mpi_utils.h"

#define   NVAR   5

//------------------------------------------------------------------------------
//!  Logarithmic average: (a - b)/log(a/b) 
/*!
  Numerically stable alogorithm taken from Ismail and Roe
*/
//------------------------------------------------------------------------------
inline
double logavg(double a, double b)
{
   double xi = b/a;
   double f = (xi - 1.0) / (xi + 1.0);
   double u = f * f;

   double F;
   if (u < 1.0e-2)
   {
      double u2 = u * u;
      double u3 = u2 * u;
      F = 1.0 + u/3.0 + u2/5.0 + u3/7.0;
   }
   else
      F = log(xi)/2.0/f;

   return 0.5*(a+b)/F;
}

//------------------------------------------------------------------------------
//!  Used to evaluate Liou's fix 
//------------------------------------------------------------------------------
struct FluxData
{
   double ssw;
};

//------------------------------------------------------------------------------
//!  Material class 
/*!
  Deals with the material properties of the flow, including the flux expressions
*/
//------------------------------------------------------------------------------
class Material
{
   public:
      Material ()
      {
         mu_ref = 0.0;
         tvb_minmod_coef = 0.0;
         ev_mod_alpha = 0.0;
         ev_mod_beta = 0.0;
         nd_active = true;
         visc_active = true;
         heat_active = true;
      };
      
      enum FlowModel {euler,  /*!< Euler model */ 
                      ns      /*!< Navier-Stokes model */
                      };
      
      /*! Add new flux names here */  
      enum FluxScheme { kepes_tecno_roe,     /*!< Kinetic Energy and entropy stable flux with
                                                  Roe type dissipation */  
                        kepes_tecno_rusanov, /*!< Kinetic Energy and entropy stable flux with
                                                  Rusanov type dissipation */  
                        kepec,               /*!< Kinetic Energy and entropy conservative 
                                                  flux */ 
                        roe_ec,              /*!< Roe and Ismail's entropy conservative 
                                                  flux */
                        kep,                  /*!< Jameson's kinetic energy preserving 
                                                  flux */
                        avg                   /*!< Simple average flux */                                                                
                      };                     
      enum MuModel {mu_constant,    /*!< Constant viscosity model */  
                    mu_sutherland,  /*!< Sutherland's viscosity model */  
                    mu_power        /*!< Power law viscosity model */    
                    };     
      
      /*! Reconstruction of entropy variable is dissipation term of entropy stable fluxes */              
      enum ReconstructScheme 
      {                       
         first,        /*!< First order i.e., no reconstruction */ 
         second,       /*!< Second order reconstruction using central slope difference */
         second_muscl, /*!< Second order reconstruction using MUSCL scheme */
         van_albada,   /*!< Second order (limiter) reconstruction using MUSCL and Van Albada limiter */
         tvd_minmod,   /*!< Second order (limiter) reconstruction using minmod */
         tvb_minmod,   /*!< Second order (limiter) reconstruction using minmod (TVB version) */
         eno2          /*!< Second order (limiter) reconstruction using ENO2 */
      };                  
                        
      
      double gamma;               /*!< Ratio of specific heats */
      double gas_const;           /*!< Gas constant */
      double prandtl;             /*!< Prandtl number */
      double Cp;                  /*!< Specific heat at constant pressure */
      double T_0,                 /*!< Sutherland's temperature */
             T_ref,               /*!< Sutherland's law: reference temperature */
             mu_ref;              /*!< Sutherland's law: reference viscosity */
      double omega;               /*!< Exponent in poer law for viscosity */
      double ev_mod_alpha,        /*!< Used in Roe and Ismail's EC1 modification of eigen_values */
             ev_mod_beta;         /*!< Used in Roe and Ismail's EC2 modification of eigen_values */
      bool   nd_active,           /*!< If true, artificial numerical dissipation is active */
             visc_active,         /*!< If true, physical dissipation is active */
             heat_active;         /*!< If true, heat flux is present */
//       bool   order_correction;    /*!< If true, flux correction for the central flux is use */
      FlowModel model;      
      FluxScheme flux_scheme;
      MuModel mu_model;
      ReconstructScheme reconstruct_scheme;
      
      double tvb_minmod_coef;     /*!< Coefficient required for TVB minmod limiting */
      
      
	  
	  /*!
	  * Reconstruction function
	  * @param[in] recon_type type of reconstruction.
	  * @param[in] h nodal distance (needed for TVB minmod).
	  * @param[in] statel constant left state.
	  * @param[in] stater constant right state.
	  * @param[in] gradl_dot_dr left directional derivative along the line joining nodes
	  * @param[in] gradr_dot_dr right directional derivative along the line joining nodes
	  * @param[out] rec_statel  left reconstructed state.
	  * @param[out] rec_stater  right reconstructed state.
	  */
	  void reconstruct  (const unsigned int& recon_type,
						 const double& h,
						 const double (& statel) [NVAR],
						 const double (& stater) [NVAR],
						 const double (& gradl_dot_dr) [NVAR],
						 const double (& gradr_dot_dr)[NVAR],
						 double (& rec_statel)[NVAR],
						 double (& rec_stater) [NVAR]) const;
      
      /*!
	  * First order reconstruction function
	  * @param[in] statel constant left state.
	  * @param[in] stater constant right state.
	  * @param[out] rec_statel  left reconstructed state.
	  * @param[out] rec_stater  right reconstructed state.
	  */
      void reconstruct_first (const double (& statel) [NVAR],
                              const double (& stater) [NVAR],
                              double (& rec_statel)[NVAR],
                              double (& rec_stater) [NVAR]) const ;
                              
      /*!
	  * Second order (central slope) reconstruction function
	  * @param[in] statel constant left state.
	  * @param[in] stater constant right state.
	  * @param[in] gradl_dot_dr left directional derivative along the line joining nodes
	  * @param[in] gradr_dot_dr right directional derivative along the line joining nodes
	  * @param[out] rec_statel  left reconstructed state.
	  * @param[out] rec_stater  right reconstructed state.
	  */                        
      void reconstruct_second (const double (& statel) [NVAR],
                               const double (& stater) [NVAR],
                               const double (& gradl_dot_dr) [NVAR],
                               const double (& gradr_dot_dr)[NVAR],
                               double (& rec_statel)[NVAR],
                               double (& rec_stater) [NVAR]) const;
      
      /*!
	  * Second order (MUSCL) reconstruction function
	  * @param[in] statel constant left state.
	  * @param[in] stater constant right state.
	  * @param[in] gradl_dot_dr left directional derivative along the line joining nodes
	  * @param[in] gradr_dot_dr right directional derivative along the line joining nodes
	  * @param[out] rec_statel  left reconstructed state.
	  * @param[out] rec_stater  right reconstructed state.
	  */ 
      void reconstruct_second_muscl (const double (& statel) [NVAR],
                                     const double (& stater) [NVAR],
                                     const double (& gradl_dot_dr) [NVAR],
                                     const double (& gradr_dot_dr)[NVAR],
                                     double (& rec_statel)[NVAR],
                                     double (& rec_stater) [NVAR]) const;
      
      /*!
	  * Second order (MUSCL + limiter) reconstruction function
	  * @param[in] h nodal distance (needed for TVB minmod).
	  * @param[in] statel constant left state.
	  * @param[in] stater constant right state.
	  * @param[in] gradl_dot_dr left directional derivative along the line joining nodes
	  * @param[in] gradr_dot_dr right directional derivative along the line joining nodes
	  * @param[out] rec_statel  left reconstructed state.
	  * @param[out] rec_stater  right reconstructed state.
	  */ 
      void reconstruct_muscl_limited(const unsigned int& recon_type,
                                     const double& h,
                                     const double (& statel) [NVAR],
                                     const double (& stater) [NVAR],
                                     const double (& gradl_dot_dr) [NVAR],
                                     const double (& gradr_dot_dr)[NVAR],
                                     double (& rec_statel)[NVAR],
                                     double (& rec_stater) [NVAR]) const;                      
      
      /*!
	  * Returns limited slope
	  * @param[in] ul left slope.
	  * @param[in] ur right slope.
	  * @param[in] recon_type type of reconstruction.
	  * @param[in] h nodal distance (needed for TVB minmod).
	  */ 
      double limited_slope (const double& ul, 
                            const double& ur,
                            const unsigned int& recon_type,
                            const double& h) const;
      
      /*!
	  * Basic intializations
	  */
      void initialize ();
      
      /*!
	  * Converts primary variable to conserved variable
	  */
      ConVar  prim2con (const PrimVar& prim_var);
      
      /*!
	  * Converts conserved to primary variable
	  */
      PrimVar con2prim (const ConVar&  con_var);
      
      /*!
	  * Converts primary variable to entropy variable
	  */
      EntVar prim2ent(const PrimVar &prim_var);
     
      /*!
	  * Evaluates numerical flux at a face
	  * @param[in] cl coordinates of left node.
	  * @param[in] cr coordinates of right node.
	  * @param[in] left0 constant left (primary) state.
	  * @param[in] right0 constant right (primary) state.
	  * @param[in] left reconstructed left (primary) state.
	  * @param[in] right reconstructed right (primary) state.
	  * @param[in] dE_l gradient at left node.
	  * @param[in] dE_l gradient at right node.
	  * @param[in] normal dual face normal.
	  * @param[in] data information about shock indicator (needed for Liou fix).
	  * @param[out] flux  flux evaluated at face.
	  */
      void num_flux(const Vector& cl,
                    const Vector& cr,
                    const PrimVar& left0, 
                    const PrimVar& right0, 
                    const PrimVar& left, 
                    const PrimVar& right, 
                    const std::vector<Vector>& dE_l,
                    const std::vector<Vector>& dE_r,
                    const Vector& normal, 
                    const FluxData& data,
                    Flux& flux) const;
                    
      /*!
	  * Kinetic Energy and entropy stable flux with Roe type dissipation
	  * @param[in] cl coordinates of left node.
	  * @param[in] cr coordinates of right node.
	  * @param[in] left constant left (primary) state.
	  * @param[in] right constant right (primary) state.
	  * @param[in] ent_grads_l entropy gradient at left node.
	  * @param[in] ent_grads_r entropy gradient at right node.
	  * @param[in] normal dual face normal.
	  * @param[in] data information about shock indicator (needed for Liou fix).
	  * @param[out] flux  flux evaluated at face.
	  */              
      void   kepes_tecno_roe_flux (const Vector& cl,
                                   const Vector& cr,
                                   const PrimVar& left, 
                                   const PrimVar& right, 
                                   const std::vector<Vector>& ent_grads_l, 
                                   const std::vector<Vector>& ent_grads_r, 
                                   const Vector& normal, 
                                   const FluxData& data,
                                   Flux& flux) const;      
                                   
      /*!
	  * Kinetic Energy and entropy stable flux with Rusanov type dissipation
	  * @param[in] cl coordinates of left node.
	  * @param[in] cr coordinates of right node.
	  * @param[in] left constant left (primary) state.
	  * @param[in] right constant right (primary) state.
	  * @param[in] ent_grads_l entropy gradient at left node.
	  * @param[in] ent_grads_r entropy gradient at right node.
	  * @param[in] normal dual face normal.
	  * @param[in] data information about shock indicator (needed for Liou fix).
	  * @param[out] flux  flux evaluated at face.
	  */                                                      
      void   kepes_tecno_rusanov_flux (const Vector& cl,
                                   const Vector& cr,
                                   const PrimVar& left, 
                                   const PrimVar& right, 
                                   const std::vector<Vector>& ent_grads_l, 
                                   const std::vector<Vector>& ent_grads_r, 
                                   const Vector& normal, 
                                   //const FluxData& data,
                                   Flux& flux) const;  
                                   
      /*!
	  * Kinetic Energy and entropy conservative flux
	  * @param[in] left constant left (primary) state.
	  * @param[in] right constant right (primary) state.
	  * @param[in] normal dual face normal.
	  * @param[out] flux  flux evaluated at face.
	  */                                                         
      void   kepec_flux (const PrimVar& left,
                         const PrimVar& right,
                         const Vector& normal,
                         //const FluxData& data,
                         Flux& flux) const;                                              
      
      
      /*!
	  * Jameson's kinetic energy preserving flux
	  * @param[in] left constant left (primary) state.
	  * @param[in] right constant right (primary) state.
	  * @param[in] normal dual face normal.
	  * @param[out] flux  flux evaluated at face.
	  */                                                         
      void   kep_flux   (const PrimVar& left,
                         const PrimVar& right,
                         const Vector& normal,
                         //const FluxData& data,
                         Flux& flux) const;
      
      /*!
	  * Roe and Ismail's entropy conservative flux
	  * @param[in] left constant left (primary) state.
	  * @param[in] right constant right (primary) state.
	  * @param[in] normal dual face normal.
	  * @param[out] flux  flux evaluated at face.
	  */                                                         
      void   roe_ec_flux   (const PrimVar& left,
                         const PrimVar& right,
                         const Vector& normal,
                         //const FluxData& data,
                         Flux& flux) const;    
                         
      /*!
	  * Simple average flux flux
	  * @param[in] left constant left (primary) state.
	  * @param[in] right constant right (primary) state.
	  * @param[in] normal dual face normal.
	  * @param[out] flux  flux evaluated at face.
	  */                                                         
      void   simple_avg_flux   (const PrimVar& left,
                         const PrimVar& right,
                         const Vector& normal,
                         //const FluxData& data,
                         Flux& flux) const;                                  
      
      /*!
	  * Central flux correction as described by Katz and Sankaran
	  * @param[in] cl coordinates of left node.
	  * @param[in] cr coordinates of right node.
	  * @param[in] left constant left (primary) state.
	  * @param[in] right constant right (primary) state.
	  * @param[in] ent_grads_l entropy gradient at left node.
	  * @param[in] ent_grads_r entropy gradient at right node.
	  * @param[in] normal dual face normal.
	  * @param[in] data information about shock indicator (needed for Liou fix).
	  * @param[out] flux  flux evaluated at face.
	  */ 
      // void central_flux_correction(const Vector& cl,
//                                    const Vector& cr,
//                                    const PrimVar& left, 
//                                    const PrimVar& right, 
//                                    const std::vector<Vector>& ent_grads_l, 
//                                    const std::vector<Vector>& ent_grads_r, 
//                                    const Vector& normal, 
//                                    const FluxData& data,
//                                    Flux& flux) const;
      
      /*!
	  * Numerical flux at slip boundary
	  * @param[in] state primary state.
	  * @param[in] normal dual face normal.
	  * @param[out] flux  flux evaluated at face.
	  */
      void    slip_flux (const PrimVar& state, 
                         const Vector& normal, 
                         Flux& flux) const;
      
      /*!
	  * Evaluation of exact Euler flux
	  * @param[in] state primary state.
	  * @param[in] normal dual face normal.
	  * @param[out] flux  flux evaluated at face.
	  */
      void    euler_flux (const PrimVar& prim, 
                          const Vector&  normal,
                          Flux& flux) const;
     
      /*!
	  * Entropy stable numerical boundary flux
	  * @param[in] state primary state.
	  * @param[in] normal dual face normal.
	  * @param[out] flux  flux evaluated at face.
	  */
      void    inviscid_bflux (const PrimVar& prim,
                              const PrimVar& primb,
                              const Vector&  normal,
                              Flux& flux) const; 
      // NEW BOUNDARY FLUX TECHNIQUE: TO BE TESTED
      //void    inviscid_bflux (const PrimVar& prim,
      //                        const Vector& velb,
      //                        const double& pb,
      //                        const double& etab, 
      //                        const Vector&  normal,
      //                        Flux& flux) const; 
      //-----------------------------------------------                        
      
      /*!
	  * Evaluation of viscous flux for inner faces
	  * @param[in] state entropy variable at cell.
	  * @param[in] dE entropy variable gradient at cell.
	  * @param[in] normal0 primary face normal0.
	  * @param[in] normal1 primary face normal1.
	  * @param[in] normal2 primary face normal2.
	  * @param[out] flux0  flux contribution to vertex0.
	  * @param[out] flux1  flux contribution to vertex1.
	  * @param[out] flux2  flux contribution to vertex2.
	  */
      void viscous_flux (const EntVar& state, 
                         const std::vector<Vector>& dE,
                         const Vector&  normal0, Flux& flux0,
                         const Vector&  normal1, Flux& flux1,
                         const Vector&  normal2, Flux& flux2
                         ) const;
                         
      /*!
	  * Boundary viscous flux
	  * @param[in] state averaged entropy variable at boundary face.
	  * @param[in] adiabtic bool variable to check whether wall is adiabtic or not.
	  * @param[in] dE entropy gradient in adjoining boundary cell  .
	  * @param[in] normal face normal.
	  * @param[out] flux  boundary viscous flux .
	  */                   
      void viscous_bflux (const EntVar& state, 
                         const bool&    adiabatic,
                         const std::vector<Vector>& dE,
                         const Vector&  normal, Flux& flux
                         ) const;                                                              
      
      /*!
	  * Evaluating viscosity at temperature T
	  */
      double viscosity (const double T) const;
      
      /*!
	  * Evaluating total energy E
	  */
      double total_energy (const PrimVar& state) const;
      
      /*!
	  * Evaluating speed of sound
	  */
      double sound_speed (const PrimVar& state) const;
      
      /*!
	  * Evaluating Density
	  */
      double Density (const PrimVar& state) const;
      
      /*!
	  * Evaluating Mach number
	  */
      double Mach (const PrimVar& state) const;
      
      /*!
	  * Evaluating physical entropy
	  */
      double Entropy (const PrimVar& state) const;
      
      /*!
	  * Evaluating vorticity at node
	  */
      double Vorticity (const PrimVar& state,  const std::vector<Vector>& dE) const;

};

//------------------------------------------------------------------------------
// Convert primitive to conserved
//------------------------------------------------------------------------------
inline
ConVar Material::prim2con(const PrimVar& prim_var)
{
   ConVar con_var;
   double density = Density (prim_var);

   con_var.density  = density;
   con_var.momentum.equ(prim_var.velocity,density);
   con_var.energy   = prim_var.pressure/(gamma - 1.0) +
                        0.5 * prim_var.velocity.square() * density;

   return con_var;
}

//------------------------------------------------------------------------------
// Convert conserved to primitive
//------------------------------------------------------------------------------
inline
PrimVar Material::con2prim (const ConVar& con_var)
{
   PrimVar prim_var;

   prim_var.velocity.equ(con_var.momentum,1.0/ con_var.density);
   prim_var.pressure = (gamma - 1.0) * 
        ( con_var.energy - 0.5 * con_var.momentum.square() / con_var.density );
   prim_var.temperature = prim_var.pressure / (gas_const * con_var.density);

   return prim_var;
}

//------------------------------------------------------------------------------
// Converting primary to entropy variables
//------------------------------------------------------------------------------
inline
EntVar Material::prim2ent(const PrimVar &prim_var)
{
   EntVar ent_var;
   double density = Density (prim_var);
   double beta = density/(2*prim_var.pressure);
   double s = Entropy (prim_var);
   ent_var.entf = (gamma - s)/(gamma - 1) - (prim_var.velocity.square())*beta;
   ent_var.entv.equ(prim_var.velocity,(2*beta));
   ent_var.entl = -2*beta;
   
   return ent_var;   	
}

//------------------------------------------------------------------------------
// Viscosity coefficient according to sutherland law
//------------------------------------------------------------------------------
inline
double Material::viscosity (const double T) const
{
   switch (mu_model)
   {
      case mu_constant:
         return mu_ref;

      case mu_sutherland:
         return mu_ref * std::pow(T, 1.5)/ (T + T_0);

      case mu_power:
         return mu_ref * std::pow(T/T_ref, omega);

      default:
         MPI_ERR("viscosity: unknown model " << mu_model);
   }
}

//------------------------------------------------------------------------------
// Total energy per unit volume
//------------------------------------------------------------------------------
inline
double Material::total_energy (const PrimVar& state) const
{
   double density = Density (state);
   return state.pressure / (gamma - 1.0) + 
          0.5 * density * state.velocity.square();
}

//------------------------------------------------------------------------------
// sound speed
//------------------------------------------------------------------------------
inline
double Material::sound_speed (const PrimVar& state) const
{
   return sqrt(gamma * gas_const * state.temperature);
}

//------------------------------------------------------------------------------
// Density
//------------------------------------------------------------------------------
inline
double Material::Density (const PrimVar& state) const
{
   return state.pressure / (gas_const * state.temperature);
}

//------------------------------------------------------------------------------
// Mach number
//------------------------------------------------------------------------------
inline
double Material::Mach (const PrimVar& state) const
{
   double sonic = sound_speed (state);
   return state.velocity.norm() / sonic;
}

//------------------------------------------------------------------------------
// Physical entropy
//------------------------------------------------------------------------------
inline
double Material::Entropy (const PrimVar& state) const
{
   return log(state.pressure) - gamma*log(Density (state));
}


//------------------------------------------------------------------------------
// Vorticity
//------------------------------------------------------------------------------
inline
double Material::Vorticity (const PrimVar& state, const std::vector<Vector>& dE) const
{
   return gas_const*state.temperature*(dE[2].x + state.velocity.y*dE[4].x
								       - dE[1].y - state.velocity.x*dE[4].y);
}


#endif
