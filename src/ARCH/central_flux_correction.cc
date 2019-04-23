#include <algorithm>
#include <cmath>
#include "material.h"
#include "limiter.h"

using namespace std;

//------------------------------------------------------------------------------
// Central flux correction as proposed by Katz and Sankaran
//------------------------------------------------------------------------------
void Material::central_flux_correction (const Vector& cl,
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
   Vector unitn,vl,vr;
   unitn.equ(normal, 1.0/area);
   
   double rhol = Density(left);
   double rhor = Density(right);
   double pl = left.pressure;
   double pr = right.pressure;
   vl = left.velocity;
   vr = right.velocity;
   double betal = 0.5 / (gas_const * left.temperature);
   double betar = 0.5 / (gas_const * right.temperature);
   double vnl = left.velocity * unitn;
   double vnr = right.velocity * unitn;
   double v2l = left.velocity.square();
   double v2r = right.velocity.square(); 
   double bvnl = betal*vnl;
   double bvnr = betar*vnr;

   // Gradient dot dr
   double rad_vec[] = {cr.x - cl.x, cr.y - cl.y, cr.z - cl.z};
   //.double rad_vec_norm = sqrt(pow(rad_vec[0],2) + pow(rad_vec[1],2) + pow(rad_vec[2],2) );
   double dir_gradl[5], dir_gradr[5];
   for(unsigned int i=0; i<5; ++i)
   { 
	  dir_gradl[i] = dE_l[i].x*rad_vec[0] + dE_l[i].y*rad_vec[1] + dE_l[i].z*rad_vec[2];
	  dir_gradr[i] = dE_r[i].x*rad_vec[0] + dE_r[i].y*rad_vec[1] + dE_r[i].z*rad_vec[2];
   }

   double Diff[5], Diff1[5], Diff2[5]; 
   double C1l = -betal*v2l - gamma/(gamma-1.0);
   double C1r = -betar*v2r - gamma/(gamma-1.0);   
   double C2l = C1l*C1l + 2*betal*v2l + gamma/(gamma-1);      
   double C2r = C1r*C1r + 2*betar*v2r + gamma/(gamma-1); 
   
   // Jacobians 
   double Anl[][5] = {
      { 2*bvnl*pl                  , (unitn.x + 2*bvnl*vl.x)*pl                           ,  (unitn.y + 2*bvnl*vl.y)*pl                          , (unitn.z + 2*bvnl*vl.z)*pl                          , -vnl*C1l*pl                                   },
      { (unitn.x + 2*bvnl*vl.x)*pl , (2*vl.x*unitn.x + (1 + 2*betal*vl.x*vl.x)*vnl)*pl    ,  (unitn.x*vl.y + unitn.y*vl.x + 2*bvnl*vl.x*vl.y)*pl , (unitn.x*vl.z + unitn.z*vl.x + 2*bvnl*vl.x*vl.z)*pl ,(C1l*unitn.x/(-2*betal) - vl.x*vnl*(C1l-1))*pl },
	  { (unitn.y + 2*bvnl*vl.y)*pl , (unitn.x*vl.y + unitn.y*vl.x + 2*bvnl*vl.x*vl.y)*pl  ,  (2*vl.y*unitn.y + (1 + 2*betal*vl.y*vl.y)*vnl)*pl   , (unitn.y*vl.z + unitn.z*vl.y + 2*bvnl*vl.y*vl.z)*pl ,(C1l*unitn.y/(-2*betal) - vl.y*vnl*(C1l-1))*pl },		
	  {	(unitn.z + 2*bvnl*vl.z)*pl , (unitn.x*vl.z + unitn.z*vl.x + 2*bvnl*vl.x*vl.z)*pl  ,  (unitn.y*vl.z + unitn.z*vl.y + 2*bvnl*vl.y*vl.z)*pl , (2*vl.z*unitn.z + (1 + 2*betal*vl.z*vl.z)*vnl)*pl   ,(C1l*unitn.z/(-2*betal) - vl.z*vnl*(C1l-1))*pl },
	  { -vnl*C1l*pl                , (C1l*unitn.x/(-2*betal) - vl.x*vnl*(C1l-1))*pl       ,  (C1l*unitn.y/(-2*betal) - vl.y*vnl*(C1l-1))*pl      , (C1l*unitn.z/(-2*betal) - vl.z*vnl*(C1l-1))*pl      ,(-C2l*vnl/(-2*betal))*pl                       }
	  };
	  
   double Anr[][5] = {
      { 2*bvnr*pr                  , (unitn.x + 2*bvnr*vr.x)*pr                           ,  (unitn.y + 2*bvnr*vr.y)*pr                          , (unitn.z + 2*bvnr*vr.z)*pr                          , -vnr*C1r*pr                                   },
      { (unitn.x + 2*bvnr*vr.x)*pr , (2*vr.x*unitn.x + (1 + 2*betar*vr.x*vr.x)*vnr)*pr    ,  (unitn.x*vr.y + unitn.y*vr.x + 2*bvnr*vr.x*vr.y)*pr , (unitn.x*vr.z + unitn.z*vr.x + 2*bvnr*vr.x*vr.z)*pr ,(C1r*unitn.x/(-2*betar) - vr.x*vnr*(C1r-1))*pr },
	  { (unitn.y + 2*bvnr*vr.y)*pr , (unitn.x*vr.y + unitn.y*vr.x + 2*bvnr*vr.x*vr.y)*pr  ,  (2*vr.y*unitn.y + (1 + 2*betar*vr.y*vr.y)*vnr)*pr   , (unitn.y*vr.z + unitn.z*vr.y + 2*bvnr*vr.y*vr.z)*pr ,(C1r*unitn.y/(-2*betar) - vr.y*vnr*(C1r-1))*pr },		
	  {	(unitn.z + 2*bvnr*vr.z)*pr , (unitn.x*vr.z + unitn.z*vr.x + 2*bvnr*vr.x*vr.z)*pr  ,  (unitn.y*vr.z + unitn.z*vr.y + 2*bvnr*vr.y*vr.z)*pr , (2*vr.z*unitn.z + (1 + 2*betar*vr.z*vr.z)*vnr)*pr   ,(C1r*unitn.z/(-2*betar) - vr.z*vnr*(C1r-1))*pr },
	  { -vnr*C1r*pr                , (C1r*unitn.x/(-2*betar) - vr.x*vnr*(C1r-1))*pr       ,  (C1r*unitn.y/(-2*betar) - vr.y*vnr*(C1r-1))*pr      , (C1r*unitn.z/(-2*betar) - vr.z*vnr*(C1r-1))*pr      ,(-C2r*vnr/(-2*betar))*pr                       }
	  };
	  
   
   for(unsigned int i=0; i<5; ++i)
   {
	  Diff[i] = 0;
	  Diff1[i] = 0;
	  Diff2[i] = 0;
	  for(unsigned int j=0; j<5; ++j)
	  {
		  Diff[i] += Anr[i][j]*dir_gradr[j] - Anl[i][j]*dir_gradl[j];
		  // Diff1[i] += dir_gradr[j];
// 		  Diff2[i] += dir_gradl[j];
      }		  
   }	  	  
   // Jump in directed gradients
   // Diff[0] = 2*bvnr*pr*dir_gradr[0] - 2*bvnl*pl*dir_gradl[0]
//              +(unitn.x + 2*bvnr*right.velocity.x)*pr*dir_gradr[1] - (unitn.x + 2*bvnl*left.velocity.x)*pl*dir_gradl[1]
//              +(unitn.y + 2*bvnr*right.velocity.y)*pr*dir_gradr[2] - (unitn.y + 2*bvnl*left.velocity.y)*pl*dir_gradl[2]
//              +(unitn.z + 2*bvnr*right.velocity.z)*pr*dir_gradr[3] - (unitn.z + 2*bvnl*left.velocity.z)*pl*dir_gradl[3]
//              +(-vnr*C1r)*pr*dir_gradr[4] - (-vnl*C1l)*pl*dir_gradl[4];
//              
//    Diff[1] = (unitn.x + 2*bvnr*right.velocity.x)*pr*dir_gradr[0] - (unitn.x + 2*bvnl*left.velocity.x)*pl*dir_gradl[0]
//              +(2*right.velocity.x*unitn.x + (1 + 2*betar*pow(right.velocity.x,2.0))*vnr)*pr*dir_gradr[1] 
//                - (2*left.velocity.x*unitn.x + (1 + 2*betal*pow(left.velocity.x,2.0))*vnl)*pl*dir_gradl[1]
//              +(unitn.x*right.velocity.y + unitn.y*right.velocity.x + 2*bvnr*right.velocity.x*right.velocity.y)*pr*dir_gradr[2] 
//                - (unitn.x*left.velocity.y + unitn.y*left.velocity.x + 2*bvnl*left.velocity.x*left.velocity.y)*pl*dir_gradl[2]
//              +(unitn.x*right.velocity.z + unitn.z*right.velocity.x + 2*bvnr*right.velocity.x*right.velocity.z)*pr*dir_gradr[3] 
//                - (unitn.x*left.velocity.z + unitn.z*left.velocity.x + 2*bvnl*left.velocity.x*left.velocity.z)*pl*dir_gradl[3]  
//              +(C1r*unitn.x/(-2*betar) - right.velocity.x*vnr*(C1r-1))*pr*dir_gradr[4] 
//                -(C1l*unitn.x/(-2*betal) - left.velocity.x*vnr*(C1l-1))*pl*dir_gradl[4];
//    
//    Diff[2] = (unitn.y + 2*bvnr*right.velocity.y)*pr*dir_gradr[0] - (unitn.y + 2*bvnl*left.velocity.y)*pl*dir_gradl[0]
//              +(unitn.x*right.velocity.y + unitn.y*right.velocity.x + 2*bvnr*right.velocity.x*right.velocity.y)*pr*dir_gradr[1] 
//                - (unitn.x*left.velocity.y + unitn.y*left.velocity.x + 2*bvnl*left.velocity.x*left.velocity.y)*pl*dir_gradl[1]
//              +(2*right.velocity.y*unitn.y + (1 + 2*betar*pow(right.velocity.y,2.0))*vnr)*pr*dir_gradr[2] 
//                -(2*left.velocity.y*unitn.y + (1 + 2*betal*pow(left.velocity.y,2.0))*vnl)*pl*dir_gradl[2]
//              +(unitn.y*right.velocity.z + unitn.z*right.velocity.y + 2*bvnr*right.velocity.y*right.velocity.z)*pr*dir_gradr[3] 
//                - (unitn.y*left.velocity.z + unitn.z*left.velocity.y + 2*bvnl*left.velocity.y*left.velocity.z)*pl*dir_gradl[3]  
//              +(C1r*unitn.y/(-2*betar) - right.velocity.y*vnr*(C1r-1))*pr*dir_gradr[4] 
//                -(C1l*unitn.y/(-2*betal) - left.velocity.y*vnr*(C1l-1))*pl*dir_gradl[4];   
//    
//    Diff[3] = (unitn.z + 2*bvnr*right.velocity.z)*pr*dir_gradr[0] - (unitn.z + 2*bvnl*left.velocity.z)*pl*dir_gradl[0]
//              +(unitn.x*right.velocity.z + unitn.z*right.velocity.x + 2*bvnr*right.velocity.x*right.velocity.z)*pr*dir_gradr[1] 
//                - (unitn.x*left.velocity.z + unitn.z*left.velocity.x + 2*bvnl*left.velocity.x*left.velocity.z)*pl*dir_gradl[1]
//              +(unitn.y*right.velocity.z + unitn.z*right.velocity.y + 2*bvnr*right.velocity.y*right.velocity.z)*pr*dir_gradr[2] 
//                - (unitn.y*left.velocity.z + unitn.z*left.velocity.y + 2*bvnl*left.velocity.y*left.velocity.z)*pl*dir_gradl[2]  
//              +(2*right.velocity.z*unitn.z + (1 + 2*betar*pow(right.velocity.z,2.0))*vnr)*pr*dir_gradr[3] 
//                -(2*left.velocity.z*unitn.z + (1 + 2*betal*pow(left.velocity.z,2.0))*vnl)*pl*dir_gradl[3]
//              +(C1r*unitn.z/(-2*betar) - right.velocity.z*vnr*(C1r-1))*pr*dir_gradr[4] 
//                -(C1l*unitn.z/(-2*betal) - left.velocity.z*vnr*(C1l-1))*pl*dir_gradl[4]; 
//                
//    Diff[4] = (-vnr*C1r)*pr*dir_gradr[0] - (-vnl*C1l)*pl*dir_gradl[0]   
//              +(C1r*unitn.x/(-2*betar) - right.velocity.x*vnr*(C1r-1))*pr*dir_gradr[1] 
//                -(C1l*unitn.x/(-2*betal) - left.velocity.x*vnr*(C1l-1))*pl*dir_gradl[1]  
//              +(C1r*unitn.y/(-2*betar) - right.velocity.y*vnr*(C1r-1))*pr*dir_gradr[2] 
//                -(C1l*unitn.y/(-2*betal) - left.velocity.y*vnr*(C1l-1))*pl*dir_gradl[2]  
//              +(C1r*unitn.z/(-2*betar) - right.velocity.z*vnr*(C1r-1))*pr*dir_gradr[3] 
//                -(C1l*unitn.z/(-2*betal) - left.velocity.z*vnr*(C1l-1))*pl*dir_gradl[3]   
//              +(-C2r*vnr/(-2*betar))*pr*dir_gradr[4] - (-C2l*vnl/(-2*betal))*pl*dir_gradl[4];                              
             
	
   flux.mass_flux       -= 0.25 * Diff[0];
   flux.momentum_flux.x -= 0.25 * Diff[1];
   flux.momentum_flux.y -= 0.25 * Diff[2];
   flux.momentum_flux.z -= 0.25 * Diff[3];
   flux.energy_flux     -= 0.25 * Diff[4];   
   
   // cout<<cr.x<<" "<<cr.y<<" "<<cr.z<<" <-- "<<cl.x<<" "<<cl.y<<" "<<cl.z<<endl;
//    cout<<dir_gradl[0]<<" "<<dir_gradl[1]<<" "<<dir_gradl[2]<<" "<<dir_gradl[3]<<" "<<dir_gradl[4]<<endl;
//    cout<<dir_gradr[0]<<" "<<dir_gradr[1]<<" "<<dir_gradr[2]<<" "<<dir_gradr[3]<<" "<<dir_gradr[4]<<endl;

   //flux *= area;
}
