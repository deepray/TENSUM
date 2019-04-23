#ifndef __EXACT_H__
#define __EXACT_H__

#include <iostream>
#include <string>
#include "vec.h"
#include "material.h"
#include "fparser.h"

//------------------------------------------------------------------------------
//!  Exact Solution class 
/*!
  This is the class used to specify exact solution
*/
//------------------------------------------------------------------------------
class ExactSolution
{
   public:
      ExactSolution (){};
      void    add (std::string, std::string);
      PrimVar value (const Vector& p, const double& t);

   private:
      FParser temperature;
      FParser xvelocity;
      FParser yvelocity;
      FParser zvelocity;
      FParser pressure;
};

//------------------------------------------------------------------------------
// Add function as defined in "fun"
//------------------------------------------------------------------------------
inline
void ExactSolution::add (std::string variable, std::string fun)
{
   if(variable == "temperature")
      temperature.FParse (fun);
   else if(variable == "xvelocity")
      xvelocity.FParse (fun);
   else if(variable == "yvelocity")
      yvelocity.FParse (fun);
   else if(variable == "zvelocity")
      zvelocity.FParse (fun);
   else if(variable == "pressure")
      pressure.FParse (fun);
   else
      MPI_LOC_ERR("ExactSolution::add: Unknown variable " << variable);
}

//------------------------------------------------------------------------------
// Evaluate primitive variables for given point
//------------------------------------------------------------------------------
inline
PrimVar ExactSolution::value (const Vector& p, const double& t)
{
   PrimVar result;

   double vals[3] = {p.x, p.y, t};
   result.temperature= temperature.Eval (vals);
   result.velocity.x = xvelocity.Eval (vals);
   result.velocity.y = yvelocity.Eval (vals);
   result.velocity.z = zvelocity.Eval (vals);
   result.pressure   = pressure.Eval (vals);

   return result;
}

#endif
