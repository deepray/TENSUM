#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <sstream>
#include "fv.h"
#include "parameter.h"
#include "mpi_utils.h"

using namespace std;
//------------------------------------------------------------------------------
// Reduction of error norm from neighbouring processors
//------------------------------------------------------------------------------
void FiniteVolume::reduce_err_norm()
{
    double SBUF[NVAR+3], RBUF[NVAR+3];
    
    // L1 norm
	SBUF[0] = prim_err_L1.temperature;
	SBUF[1] = prim_err_L1.velocity.x;
	SBUF[2] = prim_err_L1.velocity.y;
	SBUF[3] = prim_err_L1.velocity.z;
	SBUF[4] = prim_err_L1.pressure;
	SBUF[5] = density_err_L1;
	SBUF[6] = entropy_err_L1;
	SBUF[7] = mach_err_L1;  
    
    MPI_Allreduce(&SBUF, &RBUF, NVAR+3, MPI_DOUBLE, MPI_SUM, grid.run_comm);
    
    prim_err_L1.temperature = RBUF[0];
	prim_err_L1.velocity.x  = RBUF[1];
	prim_err_L1.velocity.y  = RBUF[2];
	prim_err_L1.velocity.z  = RBUF[3];
	prim_err_L1.pressure    = RBUF[4];
	density_err_L1          = RBUF[5];
	entropy_err_L1          = RBUF[6];
	mach_err_L1             = RBUF[7]; 
	
	// L2 norm
	SBUF[0] = prim_err_L2.temperature;
	SBUF[1] = prim_err_L2.velocity.x;
	SBUF[2] = prim_err_L2.velocity.y;
	SBUF[3] = prim_err_L2.velocity.z;
	SBUF[4] = prim_err_L2.pressure;
	SBUF[5] = density_err_L2;
	SBUF[6] = entropy_err_L2;
	SBUF[7] = mach_err_L2;  
    
    MPI_Allreduce(&SBUF, &RBUF, NVAR+3, MPI_DOUBLE, MPI_SUM, grid.run_comm);
    
    prim_err_L2.temperature = RBUF[0];
	prim_err_L2.velocity.x  = RBUF[1];
	prim_err_L2.velocity.y  = RBUF[2];
	prim_err_L2.velocity.z  = RBUF[3];
	prim_err_L2.pressure    = RBUF[4];
	density_err_L2          = RBUF[5];
	entropy_err_L2          = RBUF[6];
	mach_err_L2             = RBUF[7];  
	
	// Linf norm
	SBUF[0] = prim_err_Linf.temperature;
	SBUF[1] = prim_err_Linf.velocity.x;
	SBUF[2] = prim_err_Linf.velocity.y;
	SBUF[3] = prim_err_Linf.velocity.z;
	SBUF[4] = prim_err_Linf.pressure;
	SBUF[5] = density_err_Linf;
	SBUF[6] = entropy_err_Linf;
	SBUF[7] = mach_err_Linf;  
    
    MPI_Allreduce(&SBUF, &RBUF, NVAR+3, MPI_DOUBLE, MPI_SUM, grid.run_comm);
    
    prim_err_Linf.temperature = RBUF[0];
	prim_err_Linf.velocity.x  = RBUF[1];
	prim_err_Linf.velocity.y  = RBUF[2];
	prim_err_Linf.velocity.z  = RBUF[3];
	prim_err_Linf.pressure    = RBUF[4];
	density_err_Linf          = RBUF[5];
	entropy_err_Linf          = RBUF[6];
	mach_err_Linf             = RBUF[7];
  
}