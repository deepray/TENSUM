#ifndef __MPI_UTILS_H__
#define __MPI_UTILS_H__

#include <mpi.h>
#include <map>
#include <vector>
#include <cstdlib>
#include <set>
#include <string>
#include <sstream>
extern unsigned int NPART;
extern unsigned int N_MC_GROUPS;

#define MPI_ASSERT(__condition__) \
if(!(__condition__)){\
   if(get_proc_id() == 0){std::cout<<"In "<<__FILE__<<" line "<<__LINE__\
                                   <<"; assertion failed: ("<<#__condition__<<")"<<std::endl;}\
   MPI_Barrier(MPI_COMM_WORLD); \
   MPI_Finalize(); \
   exit(0);};

#define MPI_LOC_ASSERT(__condition__) \
if(!(__condition__)){\
   std::cout<<"In "<<__FILE__<<" line "<<__LINE__\
            <<"; assertion failed: ("<<#__condition__<<")"<<std::endl;\
   MPI_Abort(MPI_COMM_WORLD,1);};   
   
#define MPI_ERR(__message__) \
{std::stringstream ss; \
ss <<__message__; \
if(get_proc_id() == 0){std::cout<<ss.str().c_str()<<std::endl;}\
MPI_Barrier(MPI_COMM_WORLD); \
MPI_Finalize(); \
exit(0);};

#define MPI_LOC_ERR(__message__) \
{std::stringstream ss; \
ss <<__message__; \
std::cout<<ss.str().c_str()<<std::endl;\
MPI_Abort(MPI_COMM_WORLD,1);};

#define MPI_DISP(__message__,__status__) \
if(get_proc_id() == 0 && (__status__)){\
      cout<<__message__<<std::endl;};

/*!
  Class used for creating variable sized buffer array for MPI communication.
*/
class VARARRAY
{
  public:
    double* buf;
    VARARRAY(){ buf = NULL;}
  ~VARARRAY(){ delete[] buf;}
};

/*!
 * Structure for group communicator tags
*/
struct GROUP_COMM
{
   MPI_Group new_group;
   MPI_Comm  new_comm;
};

void mpi_init(int argc, char *argv[]);
void send_all_groupings(int sender, int receiver, std::set< std::set<int> >& list,int tag);
void get_all_groupings(std::set< std::set<int> >& list);
void create_group_comm(MPI_Comm& group_comm);
void mpi_barrier(MPI_Comm group_comm);
void mpi_err_finalize();
void mpi_finalize();

unsigned int get_proc_id();
unsigned int get_proc_loc_id();
unsigned int get_proc_glob_id(int id);
unsigned int get_gbase_glob_id();
unsigned int get_comm_size();

bool check_group_base();
bool active_proc();



#endif