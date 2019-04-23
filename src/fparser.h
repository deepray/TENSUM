#ifndef __FPARSER_H__
#define __FPARSER_H__


#include <cmath>
#include <map>
#include "fparser.hh"
#include "perturb.h"

extern std::map<std::string, double> constants;


//------------------------------------------------------------------------------
//!  FParser class 
/*!
  Parses functions and constants for evaluating variables. This is used in setting
  initial and boundary conditions
*/
//------------------------------------------------------------------------------
class FParser: public FunctionParser
{
   public:
      FParser()
      {
         AddConstant("pi", M_PI);
         AddFunction("pert",PERT,4);
      }
      void FParse (const std::string& fun)
      {
         std::map<std::string,double>::iterator it;
         for ( it=constants.begin() ; it != constants.end(); it++ )
            AddConstant ((*it).first, (*it).second);
         Parse (fun, "x,y,t");
      }
};

#endif
