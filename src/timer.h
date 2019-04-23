#ifndef __TIMER_H__
#define __TIMER_H__

#include <iostream>
#include <stdint.h>
#include "mpi_utils.h"

class Timer
{
    public:
       Timer(){t_total=t_min=t_max=t_avg=0; count=0;};
       ~Timer(){};
       double t_start, t_end;
       double t_total;
       double t_min, t_max, t_avg;
       uint64_t count;
       void start_time();
       void end_time();
       void add_time();
       double hours();
       void reduce(MPI_Comm group_comm);
       void reset();
       
};

inline
void Timer::start_time()
{
   t_start = MPI_Wtime();
   ++count;
}

inline
void Timer::end_time()
{
   t_end = MPI_Wtime();
}

inline
void Timer::add_time()
{
   end_time();
   t_total  += (t_end-t_start)/60.0;
}

inline
void Timer::reduce(MPI_Comm group_comm)
{
   MPI_Allreduce(&t_total,&t_min,1,MPI_DOUBLE,MPI_MIN,group_comm);
   MPI_Allreduce(&t_total,&t_max,1,MPI_DOUBLE,MPI_MAX,group_comm);
   MPI_Allreduce(&t_total,&t_avg,1,MPI_DOUBLE,MPI_SUM,group_comm);
   int size;
   MPI_Comm_size(group_comm, &size); 
   t_avg /= double(size);
   
}

inline
double Timer::hours()
{
   return t_total/60.0;
}


inline
void Timer::reset()
{
   t_total = t_min = t_max = t_avg = 0;
   count = 0;
}

#endif
