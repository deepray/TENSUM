#ifndef __LIMITER_H__
#define __LIMITER_H__

#include <cmath>

//!  Limiters. 
/*!
  Defines various limiter functions used for reconstruction
*/

//------------------------------------------------------------------------------
// minmod function
//------------------------------------------------------------------------------
inline
double minmod_tvd (const double &a, const double &b)
{
   double result;

   if( a < 0.0 && b < 0.0 )
      result = -std::min(-a, -b);
   else if( a > 0.0 && b > 0.0 )
      result = std::min(a, b);
   else
      result = 0.0;

   return result;
}

//------------------------------------------------------------------------------
// minabs function == eno2
//------------------------------------------------------------------------------
inline
double minabs(const double &a, const double &b)
{
   double result;

   if( fabs(a) < fabs(b) )
      result = a;
   else 
      result = b;

   return result;
}

//------------------------------------------------------------------------------
// minmod function for three arguments
//------------------------------------------------------------------------------
inline
double minmod3 (const double &a, const double &b, const double &c)
{
   double result;

   if( a < 0.0 && b < 0.0  && c < 0.0)
      result = -std::min(std::min(-a, -b),-c);
   else if( a >= 0.0 && b >= 0.0 && c >= 0.0)
      result = std::min(std::min(a, b),c);
   else
      result = 0.0;

   return result;
}

//------------------------------------------------------------------------------
// TVB minmod function
//------------------------------------------------------------------------------
inline
double minmod_tvb (const double &a, const double &b, const double &h, const double &M)
{
   double result,c;
   c = 0.5*(a+b);
   if( fabs(c) < M*pow(h,2))
      result = c;   
   else
      result = minmod_tvd(a,b);

   return result;
}

//------------------------------------------------------------------------------
// Update phi using Barth-Jespersen scheme
//------------------------------------------------------------------------------
inline
void minmax (const double& T, 
             const double& u, 
             const double& v, 
             const double& w, 
             const double& p, 
             const PrimVar& pmin, 
             const PrimVar& pmax, 
             const PrimVar& prim, 
             PrimVar& phi)
{
   // temperature
   if(T > pmax.temperature)
   {
      double fact = (pmax.temperature - prim.temperature) / 
                    (T - prim.temperature);
      phi.temperature = std::min(phi.temperature, fact);
   }
   else if(T < pmin.temperature)
   {
      double fact = (pmin.temperature - prim.temperature) / 
                    (T - prim.temperature);
      phi.temperature = std::min(phi.temperature, fact);
   }

   // x velocity
   if(u > pmax.velocity.x)
   {
      double fact = (pmax.velocity.x - prim.velocity.x) / 
                    (u - prim.velocity.x);
      phi.velocity.x = std::min(phi.velocity.x, fact);
   }
   else if(u < pmin.velocity.x)
   {
      double fact = (pmin.velocity.x - prim.velocity.x) / 
                    (u - prim.velocity.x);
      phi.velocity.x = std::min(phi.velocity.x, fact);
   }

   // y velocity
   if(v > pmax.velocity.y)
   {
      double fact = (pmax.velocity.y - prim.velocity.y) / 
                    (v - prim.velocity.y);
      phi.velocity.y = std::min(phi.velocity.y, fact);
   }
   else if(v < pmin.velocity.y)
   {
      double fact = (pmin.velocity.y - prim.velocity.y) / 
                    (v - prim.velocity.y);
      phi.velocity.y = std::min(phi.velocity.y, fact);
   }

   // z velocity
   if(w > pmax.velocity.z)
   {
      double fact = (pmax.velocity.z - prim.velocity.z) / 
                    (w - prim.velocity.z);
      phi.velocity.z = std::min(phi.velocity.z, fact);
   }
   else if(w < pmin.velocity.z)
   {
      double fact = (pmin.velocity.z - prim.velocity.z) / 
                    (w - prim.velocity.z);
      phi.velocity.z = std::min(phi.velocity.z, fact);
   }

   // pressure
   if(p > pmax.pressure)
   {
      double fact = (pmax.pressure - prim.pressure) / 
                    (p - prim.pressure);
      phi.pressure = std::min(phi.pressure, fact);
   }
   else if(p < pmin.pressure)
   {
      double fact = (pmin.pressure - prim.pressure) / 
                    (p - prim.pressure);
      phi.pressure = std::min(phi.pressure, fact);
   }

}

#endif
