#include <algorithm>
#include <cmath>
#include "material.h"
#include "limiter.h"

using namespace std;

//------------------------------------------------------------------------------
// KEPEP_TECNO_ROE flux
//------------------------------------------------------------------------------
void Material::kepes_tecno_roe_flux (const Vector& cl,
									 const Vector& cr,
									 const PrimVar& left,
									 const PrimVar& right,
									 const vector<Vector>& dE_l, 
									 const vector<Vector>& dE_r, 
									 const Vector& normal,
									 const FluxData& data,
									 Flux& flux) const
{
   //static const double BETA = 0.0; //1.0/6.0;
   
   double area = normal.norm();
   Vector unit_normal;
   unit_normal.equ(normal, 1.0/area);

   double rhol = Density(left);
   double rhor = Density(right);
   double rho = logavg( rhol, rhor );
   //double vel2= 0.5 * (left.velocity.square() + right.velocity.square());
   double betal = 0.5 / (gas_const * left.temperature);
   double betar = 0.5 / (gas_const * right.temperature);
   double beta = logavg(betal, betar);
   double a   = sqrt(0.5 * gamma / beta);

   double p     = 0.5 * (rhol + rhor) / (betal + betar);

   //Vector vel = (left.velocity + right.velocity) / 2.0;
   Vector vel;
   vel.equ(left.velocity,right.velocity,0.5,0.5);
   double vel_normal = vel * unit_normal;

   // central flux
   // flux.mass_flux = rho * vel_normal;
//    flux.momentum_flux.equ(unit_normal,vel,p,flux.mass_flux);
//    flux.energy_flux = 0.5 * ( 1.0/((gamma-1.0)*beta) - vel2) * flux.mass_flux + 
//                       flux.momentum_flux *  vel;
    
   kepec_flux (left, right, normal, flux); 


   if(nd_active)
   {	
	   // entropy dissipation
	   // eigenvectors
	   double H = a*a/(gamma-1.0) + 0.5*vel.square();
	   double v1 = vel.x * unit_normal.y - vel.y * unit_normal.x;
	   double v2 = vel.z * unit_normal.x - vel.x * unit_normal.z;
   
	   double R[][5] = {
		  {                1         ,          1         ,           0         ,          0       ,             1           },
		  {vel.x - a*unit_normal.x   ,         vel.x      ,    unit_normal.y    , -unit_normal.z   , vel.x + a*unit_normal.x }, 
		  {vel.y - a*unit_normal.y   ,         vel.y      ,    -unit_normal.x   ,          0       , vel.y + a*unit_normal.y },
		  {vel.z - a*unit_normal.z   ,         vel.z      ,           0         ,    unit_normal.x , vel.z + a*unit_normal.z },
		  {H     - a*vel_normal      ,    0.5*vel.square(),          v1         ,          v2      ,  H     + a*vel_normal   } 
	   };
   
	   double const1 = sqrt(rho/gamma);
	   double psqrt = sqrt(p);
	   double S_sqrt[] = { sqrt(0.5)*const1, sqrt(gamma-1.0)*const1, psqrt, psqrt, sqrt(0.5)*const1 };
	   for(unsigned int i=0; i<5; ++i)
		  for(unsigned int j=0; j<5; ++j)
			 R[j][i] *= S_sqrt[i];
   
	   // eigenvalues
	   double vnl = left.velocity  * unit_normal;
	   double vnr = right.velocity * unit_normal;
	   double al  = sound_speed (left);
	   double ar  = sound_speed (right);
   
	   // ROE
	   double LambdaL[] = { vnl - al, vnl, vnl, vnl, vnl + al };
	   double LambdaR[] = { vnr - ar, vnr, vnr, vnr, vnr + ar };
	   double l2, l3;
	   l2 = l3 = fabs(vel_normal);
   
	   if(data.ssw > 0) l3 = max(l3, a); // Liou modification for carbuncle
	   double Lambda[]  = { (1+ev_mod_beta)*fabs(vel_normal - a) + ev_mod_alpha*fabs(LambdaL[0]-LambdaR[0]), 
							l2,
							l3,
							l3,
							(1+ev_mod_beta)*fabs(vel_normal + a) + ev_mod_alpha*fabs(LambdaL[4]-LambdaR[4])};

	   // entropy
	   double sl    = Entropy(left);
	   double sr    = Entropy(right); 
   
	   // entropy variables
	   betal = 0.5 / (gas_const * left.temperature);
	   betar = 0.5 / (gas_const * right.temperature);
	   double Vl[] = {(gamma-sl)/(gamma-1.0) - betal*left.velocity.square(),
					  2*left.velocity.x*betal,
					  2*left.velocity.y*betal,
					  2*left.velocity.z*betal,
					  -2*betal};
	   double Vr[] = {(gamma-sr)/(gamma-1.0) - betar*right.velocity.square(),
					  2*right.velocity.x*betar,
					  2*right.velocity.y*betar,
					  2*right.velocity.z*betar,
					  -2*betar};    
   
	   // Gradient dot dr
	   double rad_vec[] = {cr.x - cl.x, cr.y - cl.y, cr.z - cl.z};
	   double rad_vec_norm = sqrt(pow(rad_vec[0],2) + pow(rad_vec[1],2) + pow(rad_vec[2],2) );
	   double gradl_dot_dr[5], gradr_dot_dr[5];
	   for(unsigned int i=0; i<5; ++i)
	   { 
		  gradl_dot_dr[i] = dE_l[i].x*rad_vec[0] + dE_l[i].y*rad_vec[1] + dE_l[i].z*rad_vec[2];
		  gradr_dot_dr[i] = dE_r[i].x*rad_vec[0] + dE_r[i].y*rad_vec[1] + dE_r[i].z*rad_vec[2];
	   }
   
	   double Diff[5];             
                                   
	   // Scaled entropy variables
	   double Zl[] = {0.0,0.0,0.0,0.0,0.0};
	   double Zr[] = {0.0,0.0,0.0,0.0,0.0};
	   for(unsigned int i=0; i<5; ++i)
		 for(unsigned int j=0; j<5; ++j)
		 {
			Zl[i] += R[j][i]*Vl[j];
			Zr[i] += R[j][i]*Vr[j];
		 }

	   // Scaled grad_dot_dr
	   double sc_gradl_dot_dr[] =  {0.0, 0.0, 0.0, 0.0, 0.0}; 
	   double sc_gradr_dot_dr[] =  {0.0, 0.0, 0.0, 0.0, 0.0};
	   for(unsigned int i=0; i<5; ++i)
	   {
		 for(unsigned int j=0; j<5; ++j)
		 {
			sc_gradl_dot_dr[i] += R[j][i]*gradl_dot_dr[j];
			sc_gradr_dot_dr[i] += R[j][i]*gradr_dot_dr[j];
		 }
	   }      

	   // Reconstructed states
	   double Zij[5],Zji[5];     
	   reconstruct(reconstruct_scheme,rad_vec_norm,Zl,Zr,sc_gradl_dot_dr,sc_gradr_dot_dr,Zij,Zji); 

	   // Jump in reconstructed states
	   double dZij[5];
	   for(unsigned int i=0; i<5; ++i)
	   {   
		  dZij[i] = Zji[i] - Zij[i] ; 
		  Diff[i] = 0;
	   }    

	   // diffusive flux = R * Lambda * dZij
	   for(unsigned int i=0; i<5; ++i)
	   {
		 for(unsigned int j=0; j<5; ++j)
			 Diff[i] += R[i][j] * Lambda[j] * dZij[j];
	   }    
		
	   flux.mass_flux       -= 0.5 * Diff[0];
	   flux.momentum_flux.x -= 0.5 * Diff[1];
	   flux.momentum_flux.y -= 0.5 * Diff[2];
	   flux.momentum_flux.z -= 0.5 * Diff[3];
	   flux.energy_flux     -= 0.5 * Diff[4];
   }	   

   //flux *= area;
}
