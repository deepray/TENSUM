#include <map>
#include <stdlib.h>
#include <iostream>
#include "mpi_utils.h"

using namespace std;

void send_all_groupings(int sender, int receiver, std::set< std::set<int> >& list,int tag)
{
    int lsize=list.size();
    int gsize[lsize+1];
    gsize[0] = lsize;
    int sum_of_elems = 0;
    int ind = 1;
    for(set<set<int> >::iterator j=list.begin();j!=list.end();++j)
    {
       sum_of_elems += (*j).size();
       gsize[ind++] = (*j).size();
    }   
    int total_size = 1+lsize+sum_of_elems;
    int pack_vec[total_size];
    
    //PACKING
    pack_vec[0] = lsize;
    ind = 1;
    while(ind < lsize+1)
    {
       pack_vec[ind] = gsize[ind];
       ind++;
    }   
    for(set<set<int> >::iterator j=list.begin();j!=list.end();++j) 
       for(set<int>::iterator k=(*j).begin();k!=(*j).end();++k)
          pack_vec[ind++] = *k;      
    
    //SENDING
    //cout<<get_proc_id()<<" -> " <<receiver<<endl;
    MPI_Ssend(pack_vec,total_size,MPI_INT,receiver,tag,MPI_COMM_WORLD);
    //cout<<get_proc_id()<<" -> " <<receiver<<"...done"<<endl;
}

void get_all_groupings(set< set<int> >& list)
{
    MPI_Status status;
    
    // Probe for an incoming message from process zero
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    int number_amount;
    
    // When probe returns, the status object has the size and other
    // attributes of the incoming message. Get the message size
    MPI_Get_count(&status, MPI_INT, &number_amount);
    int tag = status.MPI_TAG;
    int source = status.MPI_SOURCE;

    // Allocate a buffer to hold the incoming numbers
    int* unpack_list = (int*)malloc(sizeof(int) * number_amount);

    // Now receive the message with the allocated buffer
    //cout<<get_proc_id()<<" <- " <<source<<endl;
    MPI_Recv(unpack_list, number_amount, MPI_INT, source, tag,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);   
    //cout<<get_proc_id()<<" <- " <<source<<"...done"<<endl;                     
     
    //UNPACKING
    int lsize = unpack_list[0];
    int ind = lsize+1;
    for(int j=1;j<=lsize;j++)
    {
       set<int> dummy;
       int indc = ind;
       while(indc < ind+unpack_list[j])
          dummy.insert(unpack_list[indc++]);
       ind = indc;
       list.insert(dummy);    
    }      
    free(unpack_list);
}
