#include "mpi_utils.h"
#include <iostream>

using namespace std;

/*!
 * Initializes MPI
 */
void mpi_init(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
}

/*!
 * Returns processor ID
 */
unsigned int get_proc_id()
{
   int MPI_PROC_ID;
   MPI_Comm_rank ( MPI_COMM_WORLD, &MPI_PROC_ID );
   return MPI_PROC_ID;
}

/*!
 * Returns total number of processors in use
 */
unsigned int get_comm_size()
{
   int MPI_PROC_TOTAL_NUM;
   MPI_Comm_size ( MPI_COMM_WORLD, &MPI_PROC_TOTAL_NUM );
   return MPI_PROC_TOTAL_NUM;
}

/*!
 * Create group communicators for processors involved in a single sample run
 * @param[in,out] group_comm group communicator.
 */
void create_group_comm(MPI_Comm& group_comm)
{
    int rank  = get_proc_id();
    int color = rank / NPART;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &group_comm);  
}


/*!
 * Leads to an MPI Barrier for the particular group
 */
void mpi_barrier(MPI_Comm group_comm)
{
   MPI_Barrier(group_comm); 
}

/*!
 * Finalizes MPI
 */
void mpi_finalize()
{
   MPI_Barrier(MPI_COMM_WORLD); 
   MPI_Finalize(); 
}

/*!
 * Finalizes MPI when an error is encountered and exits code
 */
void mpi_err_finalize()
{
   MPI_Barrier(MPI_COMM_WORLD); 
   MPI_Finalize(); 
   exit(0);
}

/*!
 * Returns true if the current processor is a base processor for the group it belongs to
 */
bool check_group_base()
{
   unsigned int rank = get_proc_id();
   if(rank % NPART == 0 && rank < NPART*N_MC_GROUPS)
       return true;
   else 
       return false;
       
}

/*!
 * Check if the current processor is active
 */
bool active_proc()
{
   unsigned int rank = get_proc_id();
   if(rank < NPART*N_MC_GROUPS)
       return true;
   else 
       return false;
       
}

/*!
 * Returns local processor ID in its group
 */
unsigned int get_proc_loc_id()
{
   unsigned int rank = get_proc_id();
   if(rank >= NPART*N_MC_GROUPS)
   {    MPI_ERR("Processor "<<rank<<" does not belong to a mc group, thus does "
               <<"not have a valid local id"); 
   }            
   else 
       return rank % NPART;
}

/*!
 * Returns global ID of the group base for the current processor belongs to
 */
unsigned int get_gbase_glob_id()
{
   unsigned int rank = get_proc_id();
   return (rank/NPART)*NPART; // remember all these are integers
}

/*!
 * Returns global id of the processor
 */
unsigned int get_proc_glob_id(int id)
{
   unsigned int b_rank = get_gbase_glob_id();
   return b_rank + id;
}
