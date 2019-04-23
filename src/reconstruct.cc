#include <cmath>
#include "fv.h"
#include "limiter.h"

#define EPSILON  1.0e-14
#define limit_albada(a,b)  max(0.0, (2*a*b + EPSILON)/(a*a + b*b + EPSILON))

using namespace std;

//------------------------------------------------------------------------------
// First order Reconstruct left and right states
//------------------------------------------------------------------------------
void Material::reconstruct_first (const double (& statel) [NVAR],
								  const double (& stater) [NVAR],
								  double (& rec_statel)[NVAR],
								  double (& rec_stater) [NVAR]) const
{
   for(unsigned int i=0; i<NVAR; ++i)
   {
      rec_statel[i] = statel[i];
      rec_stater[i] = stater[i];
   }   
}


//------------------------------------------------------------------------------
// Second order Reconstruct left and right states
//------------------------------------------------------------------------------
void Material::reconstruct_second (const double (& statel) [NVAR],
								   const double (& stater) [NVAR],
								   const double (& gradl_dot_dr) [NVAR],
								   const double (& gradr_dot_dr)[NVAR],
								   double (& rec_statel)[NVAR],
								   double (& rec_stater) [NVAR]) const
{
   for(unsigned int i=0; i<NVAR; i++)
   {
       rec_statel[i] = statel[i] + 0.5*gradl_dot_dr[i]; 
       rec_stater[i] = stater[i] - 0.5*gradr_dot_dr[i];
   }
}

//------------------------------------------------------------------------------
// Second order Reconstruct left and right states (unlimited MUSCL)
//------------------------------------------------------------------------------
void Material::reconstruct_second_muscl ( const double (& statel) [NVAR],
										  const double (& stater) [NVAR],
										  const double (& gradl_dot_dr) [NVAR],
										  const double (& gradr_dot_dr)[NVAR],
										  double (& rec_statel)[NVAR],
										  double (& rec_stater) [NVAR]) const
{
   double dstate[NVAR], dstatel[NVAR], dstater[NVAR];
   for(unsigned int i=0; i<NVAR; ++i)
   {
      dstate[i]    = stater[i] - statel[i] ;
      
      // left state
      dstatel[i]   = 2.0 *gradl_dot_dr[i] - dstate[i];
      rec_statel[i] = statel[i] + (dstatel[i]*(1.0-KKK) + dstate[i]*(1.0+KKK))*0.25;
      
      // right state
      dstater[i]   = 2.0 *gradr_dot_dr[i] - dstate[i];
      rec_stater[i] = stater[i] - (dstater[i]*(1.0-KKK) + dstate[i]*(1.0+KKK))*0.25;
   }
}

//------------------------------------------------------------------------------
// Reconstruct left and right states (limited MUSCL)
//------------------------------------------------------------------------------
void Material::reconstruct_muscl_limited(const unsigned int& recon_type,
										 const double& h,
										 const double (& statel) [NVAR],
										 const double (& stater) [NVAR],
										 const double (& gradl_dot_dr) [NVAR],
										 const double (& gradr_dot_dr)[NVAR],
										 double (& rec_statel)[NVAR],
										 double (& rec_stater) [NVAR]) const
{
   double dstate[NVAR], dstatel[NVAR], dstater[NVAR];
   double si,sj;
   for(unsigned int i=0; i<NVAR; ++i)
   {
      dstate[i]    = stater[i] - statel[i] ;
      
      // left state
      dstatel[i]   = 2.0 *gradl_dot_dr[i] - dstate[i];
      si = limited_slope(dstatel[i],dstate[i],recon_type,h);
      rec_statel[i] = statel[i] + si;
      
      // right state
      dstater[i]   = 2.0 *gradr_dot_dr[i] - dstate[i];
      sj = limited_slope(dstater[i],dstate[i],recon_type,h);
      rec_stater[i] = stater[i] - sj;
   }
}

//------------------------------------------------------------------------------
// Computed limited slope 
//------------------------------------------------------------------------------
double Material::limited_slope (const double& ul, 
                                const double& ur,
                                const unsigned int& recon_type,
                                const double& h) const
{
   double s,result;
   switch(recon_type)
   {
      case van_albada:
         s = limit_albada(ul,ur);
         result = 0.25*s*( (1.0-KKK*s)*ul + (1.0+KKK*s)*ur );
         break;
         
      case tvd_minmod:
         result = 0.5*minmod_tvd(ul,ur);
         break;

      case eno2:
         result = 0.5*minabs(ul,ur);
         break;
         
      case tvb_minmod:
         result = 0.5*minmod_tvb(ul,ur,h,tvb_minmod_coef);
         break;         
         
      default:
         cout << "reconstruct: unknown reconstruction scheme = " 
              << reconstruct_scheme << endl;
         exit (0);
   }
   

   return result;
}

//------------------------------------------------------------------------------
// Reconstruct left and right states
//------------------------------------------------------------------------------
void Material::reconstruct (const unsigned int& recon_type,
                            const double& h,
                            const double (& statel) [NVAR],
                            const double (& stater) [NVAR],
                            const double (& gradl_dot_dr) [NVAR],
                            const double (& gradr_dot_dr)[NVAR],
                            double (& rec_statel)[NVAR],
                            double (& rec_stater) [NVAR]) const
{
   switch(recon_type)
   {
      // First order
      case first:
         reconstruct_first (statel,stater,rec_statel,rec_stater);
         break;

      // Second order
      case second:
         reconstruct_second (statel,stater,gradl_dot_dr,gradr_dot_dr,rec_statel,rec_stater);
         break;
         
      // Second order, MUSCL
      case second_muscl:
         reconstruct_second_muscl (statel,stater,gradl_dot_dr,gradr_dot_dr,rec_statel,rec_stater);
         break;   

      // MUSCL with van Albada limiter
      case van_albada:
         reconstruct_muscl_limited (recon_type,h,statel,stater,gradl_dot_dr,gradr_dot_dr,rec_statel,rec_stater);
         break;
    
      // MUSCL with minmod limiter
      case tvd_minmod:
         reconstruct_muscl_limited (recon_type,h,statel,stater,gradl_dot_dr,gradr_dot_dr,rec_statel,rec_stater);
         break;
      
      // MUSCL with eno2/minabs limiter
      case eno2:
         reconstruct_muscl_limited (recon_type,h,statel,stater,gradl_dot_dr,gradr_dot_dr,rec_statel,rec_stater);
         break;
      
      // MUSCL with TVB minmod limiter
      case tvb_minmod:
         reconstruct_muscl_limited (recon_type,h,statel,stater,gradl_dot_dr,gradr_dot_dr,rec_statel,rec_stater);
         break;      

      default:
         MPI_LOC_ERR("reconstruct: unknown reconstruction scheme = "<<reconstruct_scheme);
   }
}
