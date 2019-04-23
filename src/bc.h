#ifndef __BC_H__
#define __BC_H__

#include <iostream>
#include <string>
#include "face.h"
#include "fparser.h"
#include "primvar.h"
#include "material.h"
#include "mpi_utils.h"

extern bool PERIODIC;

namespace BC
{
   enum BCType { none, 
                 slip, 
                 noslip, 
                 //maxwell, 
                 farfield, 
                 inlet, 
                 outlet, 
                 pressure, 
                 periodic };
}

//------------------------------------------------------------------------------
//! Boundary Condition class 
/*!
  This is the class used to specify and apply boundary conditions
*/
//------------------------------------------------------------------------------
class BoundaryCondition
{
   public:
      BoundaryCondition () {};
      
      /*!
	  * Main class constructor intializing various elements
	  * @param[in] material Material object
	  * @param[in] bc_type type of boundary condition
	  * @param[in] variable variables whose values are specified for the boundary condition
	  * @param[in] function the functional expressions for the boundary conditions
	  */
      BoundaryCondition (Material                 &material,
                         std::string              &bc_type,
                         std::vector<std::string> &variable,
                         std::vector<std::string> &function);
      
      /*!
	  * Applying boundary condition at the vertex corresponding to a boundary face
	  * @param[in] vertex face vertex
	  * @param[in] time current simulation time
	  * @param[in] face boundary face in question
	  * @param[in,out] state (vector of) state which are modified based on the specified 
	  				boundary condition
	  */
      void apply (const Vector         &vertex,
                  const double         &time,
                  const Face           &face,
                  std::vector<PrimVar> &state);
      
      /*!
	  * Find variable/parameter values for evaluation of boundary inviscid flux
	  * @param[in] face face in question
	  * @param[in] state vector of states used to assign the parameter values
	  * @param[out] vel velocity paramater
	  * @param[out] p pressure parameter
	  * @param[out] eta mathematical entropy parameter
	  */
      void inviscid_states (const Face           &face,
                            const std::vector<PrimVar> &state,
                            Vector &vel,
                            double &p,
                            double &eta);            
      
      /*!
	  * Applying boundary condition at the vertex corresponding to a boundary face
	  * @param[in] vertex face vertex
	  * @param[in] time current simulation time
	  * @param[in] face boundary face in question
	  * @param[in,out] state state which is modified based on the specified 
	  				boundary condition
	  */
      void apply (const Vector &vertex,
                  const double         &time,
                  const Face   &face,
                  PrimVar      &state);
                  
      /*!
	  * Applying slip boundary condition
	  * @param[in] face boundary face in question
	  * @param[in,out] state (vector of) state which are modified based on the specified 
	  				boundary condition
	  */            
      void apply_slip (const Face           &face,
                       std::vector<PrimVar> &state);
                       
      /*!
	  * Applying slip boundary condition
	  * @param[in] face boundary face in question
	  * @param[in,out] state state which is modified based on the specified 
	  				boundary condition
	  */                     
      void apply_slip (const Face &face,
                       PrimVar    &state);
                       
      /*!
	  * Applying no-slip boundary condition
	  * @param[in] vertex face vertex
	  * @param[in] time current simulation time
	  * @param[in,out] state (vector of) state which are modified based on the specified 
	  				boundary condition
	  */                     
      void apply_noslip (const Vector         &vertex,
                         const double         &time,
                         std::vector<PrimVar> &state);
      
      /*!
	  * Applying no-slip boundary condition
	  * @param[in] vertex face vertex
	  * @param[in] time current simulation time
	  * @param[in,out] state state which is modified based on the specified 
	  				boundary condition
	  */                     
      void apply_noslip (const Vector &vertex,
                         const double         &time,
                         PrimVar      &state);
      
      /*!
	  * Applying Maxwell boundary condition
	  * @param[in] face boundary face in question
	  * @param[in] time current simulation time
	  * @param[in,out] state (vector of) state which are modified based on the specified 
	  				boundary condition
	  */  
      // void apply_maxwell (const Face           &face,
//                           const double         &time,
//                           std::vector<PrimVar> &state);
      
      /*!
	  * Applying pressure boundary condition
	  * @param[in] vertex face vertex
	  * @param[in] time current simulation time
	  * @param[in,out] state (vector of) state which are modified based on the specified 
	  				boundary condition
	  */
      void apply_pressure (const Vector         &vertex,
                           const double         &time,
                           std::vector<PrimVar> &state);
      
      /*!
	  * Applying pressure boundary condition
	  * @param[in] vertex face vertex
	  * @param[in] time current simulation time
	  * @param[in,out] state state which is modified based on the specified 
	  				boundary condition
	  */
      void apply_pressure (const Vector         &vertex,
                           const double         &time,
                           PrimVar &state);                     
      
      /*!
	  * Applying inlet boundary condition
	  * @param[in] vertex face vertex
	  * @param[in] time current simulation time
	  * @param[in,out] state (vector of) state which are modified based on the specified 
	  				boundary condition
	  */
      void apply_inlet (const Vector         &vertex,
                        const double         &time,
                        std::vector<PrimVar> &state);
      
      /*!
	  * Applying inlet boundary condition
	  * @param[in] vertex face vertex
	  * @param[in] time current simulation time
	  * @param[in,out] state state which is modified based on the specified 
	  				boundary condition
	  */
      void apply_inlet (const Vector &vertex,
                        const double         &time,
                        PrimVar      &state);
      
      /*!
	  * Applying outlet boundary condition
	  * @param[in,out] state (vector of) state which are modified based on the specified 
	  				boundary condition
	  */
      void apply_outlet (std::vector<PrimVar> &state);
      
      /*!
	  * Applying farfield boundary condition
	  * @param[in] vertex face vertex
	  * @param[in] time current simulation time
	  * @param[in,out] state (vector of) state which are modified based on the specified 
	  				boundary condition
	  */
      void apply_farfield (const Vector         &vertex,
                           const double         &time,
                           std::vector<PrimVar> &state);
      
      std::string    name;		/*!< stores bc_type */ 
      BC::BCType     type;      /*!< finds corresponding BC tag */ 
      bool           adiabatic; /*!< true if wall is adiabatic */ 

   private:
      Material*      material;
      FParser        xvelocity;  /*!< parses functional expression for x-velocity */ 
      FParser        yvelocity;  /*!< parses functional expression for y-velocity */ 
      FParser        zvelocity;  /*!< parses functional expression for z-velocity */ 
      FParser        pressure;   /*!< parses functional expression for pressure */ 
      FParser        temperature;/*!< parses functional expression for temperature */ 
};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
inline
BoundaryCondition::BoundaryCondition (Material                 &material,
                                      std::string              &bc_type,
                                      std::vector<std::string> &variable,
                                      std::vector<std::string> &function)
:
   name (bc_type),
   adiabatic (false),
   material(&material)
{
   // Set to none for safety purpose
   type = BC::none;
   // Slip bc, no state is required
   if(bc_type == "slip")
   {
      MPI_ASSERT(variable.size() == 0);
      type = BC::slip;
      adiabatic = true;
   }
   // Slip bc, no state is required
   else if(bc_type == "periodic")
   {
      MPI_ASSERT(variable.size() == 0);
      type = BC::periodic;
      PERIODIC = true;
   }
   // noslip bc: velocity is specified. If temperature is also specified
   // then it is also used. In this case, adiabatic bc is not used
   else if(bc_type == "noslip")
   {
      MPI_ASSERT(variable.size() == 3 || variable.size() == 4);
      type = BC::noslip;
      adiabatic = true;
      bool has_xvelocity   = false;
      bool has_yvelocity   = false;
      bool has_zvelocity   = false;
      for(unsigned int i=0; i<variable.size(); ++i)
      {
         if(variable[i] == "xvelocity")
         {
            has_xvelocity = true;
            xvelocity.FParse (function[i]);
         }
         else if(variable[i] == "yvelocity")
         {
            has_yvelocity = true;
            yvelocity.FParse (function[i]);
         }
         else if(variable[i] == "zvelocity")
         {
            has_zvelocity = true;
            zvelocity.FParse (function[i]);
         }
         else if(variable[i] == "temperature")
         {
            temperature.FParse (function[i]);
            adiabatic = false;
         }
      }
      MPI_ASSERT(has_xvelocity && has_yvelocity && has_zvelocity);
   }
   // else if(bc_type == "maxwell")
//    {
//       MPI_ASSERT(variable.size() == 4);
//       type = BC::maxwell;
//       adiabatic = false;
//       bool has_xvelocity   = false;
//       bool has_yvelocity   = false;
//       bool has_zvelocity   = false;
//       bool has_temperature = false;
//       for(unsigned int i=0; i<variable.size(); ++i)
//       {
//          if(variable[i] == "xvelocity")
//          {
//             has_xvelocity = true;
//             xvelocity.FParse (function[i]);
//          }
//          else if(variable[i] == "yvelocity")
//          {
//             has_yvelocity = true;
//             yvelocity.FParse (function[i]);
//          }
//          else if(variable[i] == "zvelocity")
//          {
//             has_zvelocity = true;
//             zvelocity.FParse (function[i]);
//          }
//          else if(variable[i] == "temperature")
//          {
//             has_temperature = true;
//             temperature.FParse (function[i]);
//          }
//       }
//       MPI_ASSERT(has_xvelocity && has_yvelocity && has_zvelocity && has_temperature);
//    }
   // In this case only pressure is specified
   else if(bc_type == "pressure")
   {
      MPI_ASSERT (variable.size() == 1);
      type = BC::pressure;
      MPI_ASSERT (variable[0] == "pressure");
      pressure.FParse (function[0]);
   }
   // All values are specified
   else if(bc_type == "inlet" || bc_type == "farfield")
   {
      MPI_ASSERT (variable.size() == 5);
      if(bc_type == "inlet")
         type = BC::inlet;
      else
         type = BC::farfield;
      bool has_temperature = false;
      bool has_xvelocity   = false;
      bool has_yvelocity   = false;
      bool has_zvelocity   = false;
      bool has_pressure    = false;
      for(unsigned int i=0; i<variable.size(); ++i)
      {
         if(variable[i] == "temperature")
         {
            has_temperature = true;
            temperature.FParse (function[i]);
         }
         else if(variable[i] == "xvelocity")
         {
            has_xvelocity = true;
            xvelocity.FParse (function[i]);
         }
         else if(variable[i] == "yvelocity")
         {
            has_yvelocity = true;
            yvelocity.FParse (function[i]);
         }
         else if(variable[i] == "zvelocity")
         {
            has_zvelocity = true;
            zvelocity.FParse (function[i]);
         }
         else if(variable[i] == "pressure")
         {
            has_pressure = true;
            pressure.FParse (function[i]);
         }
      }
      MPI_ASSERT (has_temperature && has_xvelocity && has_yvelocity && has_zvelocity &&
              has_pressure);
   }
   // At outflow nothing is specified
   else if(bc_type == "outlet")
   {
      MPI_ASSERT (variable.size() == 0);
      type = BC::outlet;
   }
   else
      MPI_ERR("BoundaryCondition: Unknown boundary condition " << bc_type);

   if(type == BC::none)
      MPI_ERR("BoundaryCondition: unknown bc for " << bc_type);
}

//------------------------------------------------------------------------------
// Normal velocity is zero
// u[1] = u[0] - (u[0] * n) n
// Then u[1] * n = u[0] * n = 0
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_slip(const Face           &face,
                                   std::vector<PrimVar> &state)
{
   Vector unit_normal;
   unit_normal.equ(face.normal,1.0 / face.measure);
   double fact = state[0].velocity * unit_normal;
   state[0].velocity.sadd(unit_normal,-fact);

   state[1] = state[0];
}

inline
void BoundaryCondition::apply_slip(const Face &face,
                                   PrimVar    &state)
{
   Vector unit_normal;
   unit_normal.equ(face.normal,1.0 / face.measure);
   double fact = state.velocity * unit_normal;
   state.velocity.sadd(unit_normal,-fact);
}

//------------------------------------------------------------------------------
// Velocity is specified
// If temperature is specified, then pressure is computed using the temperature
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_noslip(const Vector         &vertex,
                                     const double         &time,
                                     std::vector<PrimVar> &state)
{
   double point[3]  = {vertex.x, vertex.y,time};
   state[0].velocity.x = xvelocity.Eval(point);
   state[0].velocity.y = yvelocity.Eval(point);
   state[0].velocity.z = zvelocity.Eval(point);

   state[1].velocity = state[0].velocity;
   state[1].pressure = state[0].pressure;

   if(adiabatic)
      state[1].temperature  = state[0].temperature;
   else
   {
      double T = temperature.Eval(point);
      state[0].temperature = T;
      state[1].temperature = T;
   }
}

//------------------------------------------------------------------------------
// Velocity is modified, others kept same.
// Used to enforce noslip after update.
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_noslip(const Vector &vertex,
                                     const double &time,
                                     PrimVar      &state)
{
   double point[3]  = {vertex.x, vertex.y,time};
   state.velocity.x = xvelocity.Eval(point);
   state.velocity.y = yvelocity.Eval(point);
   state.velocity.z = zvelocity.Eval(point);
   if(!adiabatic)
   {
      state.temperature = temperature.Eval(point);
   }
}

//------------------------------------------------------------------------------
// state[1] = wall state
//------------------------------------------------------------------------------
// inline
// void BoundaryCondition::apply_maxwell(const Face           &face,
//                                       const double         &time,
//                                       std::vector<PrimVar> &state)
// {
//    // wall values
//    double point[3]  = {face.centroid.x, face.centroid.y, time};
//    state[1].velocity.x = xvelocity.Eval(point);
//    state[1].velocity.y = yvelocity.Eval(point);
//    state[1].velocity.z = zvelocity.Eval(point);
//    state[1].temperature = temperature.Eval(point);
// 
//    // normal velocity should be zero
//    // state[1] must already have zero normal velocity
//    Vector unit_normal;
//    unit_normal.equ(face.normal,1.0 / face.measure);
//    //state[0].velocity -= unit_normal * (state[0].velocity * unit_normal);
//    double fact = state[0].velocity * unit_normal;
//    state[0].velocity.sadd(unit_normal,-fact);
// 
//    double density = material->Density(state[0]);
// 
//    double taunn= 0.0;
//    double rhow = density * 
//                  sqrt(state[0].temperature / state[1].temperature) * 
//                  (1.0 - 0.5 * taunn/state[0].pressure);
//    state[1].pressure = rhow * material->gas_const * state[1].temperature;
// }

//------------------------------------------------------------------------------
// Reset pressure value. Other states remain same
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_pressure (const Vector         &vertex,
                                        const double         &time,
                                        std::vector<PrimVar> &state)
{
   double point[3]  = {vertex.x, vertex.y, time};

   state[1].pressure  = pressure.Eval(point);
}

inline
void BoundaryCondition::apply_pressure (const Vector         &vertex,
                                        const double         &time,
                                        PrimVar &state)
{
   double point[3]  = {vertex.x, vertex.y,time};

   state.pressure  = pressure.Eval(point);
}

//------------------------------------------------------------------------------
// At inlet all values are specified
// Both states are set to inlet values
// e.g., supersonic inlet
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_inlet (const Vector         &vertex,
                                     const double         &time,
                                     std::vector<PrimVar> &state)
{
   double point[3]  = {vertex.x, vertex.y,time};
   state[0].temperature= temperature.Eval(point);
   state[0].velocity.x = xvelocity.Eval(point);
   state[0].velocity.y = yvelocity.Eval(point);
   state[0].velocity.z = zvelocity.Eval(point);
   state[0].pressure   = pressure.Eval(point);

   state[1] = state[0];
}

inline
void BoundaryCondition::apply_inlet (const Vector &vertex,
                                     const double &time,
                                     PrimVar      &state)
{
   double point[3]  = {vertex.x, vertex.y,time};
   state.temperature= temperature.Eval(point);
   state.velocity.x = xvelocity.Eval(point);
   state.velocity.y = yvelocity.Eval(point);
   state.velocity.z = zvelocity.Eval(point);
   state.pressure   = pressure.Eval(point);
}

//------------------------------------------------------------------------------
// At outlet all all values are from inside values
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_outlet (std::vector<PrimVar> &state)
{
   state[1] = state[0];
}

//------------------------------------------------------------------------------
// At farfield all values are specified
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply_farfield (const Vector         &vertex,
                                        const double         &time,
                                        std::vector<PrimVar> &state)
{
   double point[3]     = {vertex.x, vertex.y, time};
   state[1].temperature= temperature.Eval(point);
   state[1].velocity.x = xvelocity.Eval(point);
   state[1].velocity.y = yvelocity.Eval(point);
   state[1].velocity.z = zvelocity.Eval(point);
   state[1].pressure   = pressure.Eval(point);
}

//------------------------------------------------------------------------------
// Apply boundary condition based on type
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply(const Vector         &vertex,
                              const double         &time,
                              const Face           &face,
                              std::vector<PrimVar> &state)
{
   switch(type)
   {
      case BC::slip:
         apply_slip (face, state);
         break;

      case BC::noslip:
         apply_noslip (vertex, time, state);
         break;

      // case BC::maxwell:
//          apply_maxwell (face, time, state);
//          break;

      case BC::pressure:
         apply_pressure (vertex, time, state);
         break;

      case BC::inlet:
         apply_inlet (vertex, time, state);
         break;

      case BC::outlet:
         apply_outlet (state);
         break;

      case BC::farfield:
         apply_farfield (vertex, time,  state);
         break;
      case BC::periodic:
         // Nothing to be done
         break;   

      default:
         MPI_LOC_ERR("Unknown boundary condition: " << name);
   }
}

//------------------------------------------------------------------------------
// Apply boundary condition based on type
//------------------------------------------------------------------------------
inline
void BoundaryCondition::apply(const Vector  &vertex,
                              const double  &time,
                              const Face    &face,
                              PrimVar       &state)
{
   switch(type)
   {
      case BC::slip:
      // case BC::maxwell:
         apply_slip (face, state);
         break;

      case BC::noslip:
         apply_noslip (vertex, time, state);
         break;

      case BC::pressure:
         apply_pressure (vertex, time, state);
         break;

      case BC::inlet:
         apply_inlet (vertex, time, state);
         break;

      case BC::outlet:
         // Nothing to be done
         break;

      case BC::farfield:
         //apply_farfield (vertex, state);
         break;
      case BC::periodic:
         // Nothing to be done
         break;     

      default:
         MPI_LOC_ERR("Unknown boundary condition: " << name);
   }
}

//------------------------------------------------------------------------------
// Find parameter values for inviscid flux
//------------------------------------------------------------------------------
inline
void BoundaryCondition::inviscid_states(const Face &face,
                            const std::vector<PrimVar> &state,
                            Vector &vel,
                            double &p,
                            double &eta)
{
   switch(type)
   {
      case BC::slip:
         vel = state[0].velocity;
         break;
      
      // case BC::maxwell:
//          // TO DO
//          break;
         
      case BC::noslip:
         vel = state[0].velocity;
         break;
      
      case BC::inlet:
      case BC::outlet:
      case BC::farfield: 
      case BC::pressure:
      {
         double vn = state[1].velocity * face.normal;
         double a = material->sound_speed(state[1]);
         if(vn < 0.0) //INFLOW AS face.normal is outward normal 
         {
			 vel = state[1].velocity;
		     eta = -material->Density(state[1])*material->Entropy(state[1])/(material->gamma - 1.0);
			 if(std::abs(vn)>a) //Supersonic
			 {
				 p = state[1].pressure;
			 }
		 }
		 else //OUTFLOW
		 {
             if(std::abs(vn)<a) //Subsonic
			 {
				 p = state[1].pressure;
			 }
		 }
         break;
         
      }
      case BC::periodic:
         // Nothing to be done
         break;    
      default:
         MPI_LOC_ERR("Unknown boundary condition: " << name);
   }
}

#endif
