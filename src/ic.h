#ifndef __IC_H__
#define __IC_H__

#include <iostream>
#include <string>
#include <dlfcn.h>
#include "vec.h"
#include "material.h"
#include "fparser.h"

//------------------------------------------------------------------------------
//!  Initial Condition class 
/*!
  This is the class used to specify and apply initial conditions
*/
//------------------------------------------------------------------------------
class InitialCondition
{
   public:
      InitialCondition ()
         : has_lib (false) {};
      void    add (std::string, std::string);
      void    add (std::string);
      PrimVar value (const Vector& p);

   private:
      bool has_lib;
      FParser temperature;
      FParser xvelocity;
      FParser yvelocity;
      FParser zvelocity;
      FParser pressure;

      void INITIAL_CONDITION(double,double,double&,double&,double&,double&,double&);
      void (*initial_condition)(double,double,double&,double&,double&,double&,double&);
};

//------------------------------------------------------------------------------
// Load function from external shared library
//------------------------------------------------------------------------------
inline
void InitialCondition::add (std::string lib_file)
{
   void *lib_handle = dlopen(lib_file.c_str(), RTLD_LAZY);
   if(lib_handle == NULL)
      MPI_LOC_ERR(dlerror());

   initial_condition = 
      (void (*)(double,double,double&,double&,double&,double&,double&))
      dlsym (lib_handle, "INITIAL_CONDITION");
   has_lib = true;
   
}

//------------------------------------------------------------------------------
// Add function as defined in "fun"
//------------------------------------------------------------------------------
inline
void InitialCondition::add (std::string variable, std::string fun)
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
      MPI_LOC_ERR("InitialCondition::add: Unknown variable " << variable);
}

//------------------------------------------------------------------------------
// Evaluate primitive variables for given point
//------------------------------------------------------------------------------
inline
PrimVar InitialCondition::value (const Vector& p)
{
   PrimVar result;

   if (has_lib)
   {
      (*initial_condition)(p.x, p.y,
                           result.temperature, 
                           result.velocity.x,
                           result.velocity.y,
                           result.velocity.z,
                           result.pressure);
   }
   else
   {
      double t = 0.0;
      double vals[3] = {p.x, p.y, t};
      result.temperature= temperature.Eval (vals);
      result.velocity.x = xvelocity.Eval (vals);
      result.velocity.y = yvelocity.Eval (vals);
      result.velocity.z = zvelocity.Eval (vals);
      result.pressure   = pressure.Eval (vals);
   }

   return result;
}

#endif
