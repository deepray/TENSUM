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
// Reduction of square of residual norm from neighbouring processors
//------------------------------------------------------------------------------
void FiniteVolume::reduce_res_norm()
{
    t_res_norm_red.start_time();
    
    double SBUF[NVAR], RBUF[NVAR];
	SBUF[0] = residual_norm.mass_flux;
	SBUF[1] = residual_norm.momentum_flux.x;
	SBUF[2] = residual_norm.momentum_flux.y;
	SBUF[3] = residual_norm.momentum_flux.z;
	SBUF[4] = residual_norm.energy_flux;  
    
    MPI_Allreduce(&SBUF, &RBUF, NVAR, MPI_DOUBLE, MPI_SUM, grid.run_comm);
    
    residual_norm.mass_flux = RBUF[0];
	residual_norm.momentum_flux.x = RBUF[1];
	residual_norm.momentum_flux.y = RBUF[2];
	residual_norm.momentum_flux.z = RBUF[3];
	residual_norm.energy_flux = RBUF[4];  
	
	t_res_norm_red.add_time();   
}