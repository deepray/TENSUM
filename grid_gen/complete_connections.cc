#include <iostream>
#include <cstdlib>
#include <string>
#include <set>
#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <cassert>
#include "grid.h"

using namespace std;

//------------------------------------------------------------------------------
// Completes node periodic connections
//------------------------------------------------------------------------------
template <int DIM>
void Grid<DIM>::complete_connections()
{
   vector<vector<int> >::iterator itvv1, itvv2,itvv3;
   vector<int>::iterator itv1,itv2,itv3;

   // for(itvv1 = periodic_node_sets.begin(); itvv1!=periodic_node_sets.end();++itvv1)
//    {
//       for(itv1 = (*itvv1).begin(); itv1 != (*itvv1).end(); ++itv1)
//          cout<<*itv1<<" ";
//       cout<<endl;   
//    }
   
   int next_top =0;
   bool connected = false;
   itvv1 = periodic_node_sets.begin()+next_top;
   if(periodic_node_sets.size()>1)
   {
      itvv2 = itvv1+1;   
	  while(!connected)
	  {
		 int union_size = (*itvv1).size() + (*itvv2).size();
		 vector<int> dvec(union_size);
		 sort((*itvv1).begin(),(*itvv1).end());
		 sort((*itvv2).begin(),(*itvv2).end());
		 itv1 = set_difference ((*itvv1).begin(), (*itvv1).end(), (*itvv2).begin(), (*itvv2).end(), dvec.begin());
		 dvec.resize(itv1-dvec.begin());
	  
		 if(dvec.size() == (*itvv1).size())
		 {
			itvv2++;
		 }
		 else
		 {
			dvec.resize(union_size,0);
			itv1 = set_union ((*itvv1).begin(), (*itvv1).end(), (*itvv2).begin(), (*itvv2).end(), dvec.begin());
			dvec.resize(itv1-dvec.begin());
		 
			// cout<<"New graph connection vector:\n";
// 			for(itv3 = dvec.begin(); itv3 != dvec.end(); ++itv3)
// 			   cout<<*itv3<<" ";
// 			cout<<endl;   
		 
			(*itvv2).resize(0);
			(*itvv2) = dvec;
			periodic_node_sets.erase(itvv1);
			itvv1 = periodic_node_sets.begin()+next_top;
			itvv2 = itvv1+1;
		 
			// cout<<"\nGraph connection list changed to :" <<endl;
// 			for(itvv3 = periodic_node_sets.begin(); itvv3!=periodic_node_sets.end();++itvv3)
// 			{
// 			   for(itv1 = (*itvv3).begin(); itv1 != (*itvv3).end(); ++itv1)
// 				  cout<<*itv1<<" ";
// 			   cout<<endl;   
// 			}
		   
		 } 
		 if(itvv2==periodic_node_sets.end())
		 {
		    next_top++;
		    if(periodic_node_sets.size()>1)
		    {
			   itvv1 = periodic_node_sets.begin()+next_top;
			   if(itvv1 != periodic_node_sets.end() - 1)
			   {
			      itvv2 = itvv1+1;
			   } 
			   else // reached last position of listing
			   {
				  connected = true;
			   }
			}      
			else // listing only has one connection now
			{
			   connected = true;
			}
		 }
	  }
   }
   
   //Sorting lists
   for(itvv1 = periodic_node_sets.begin(); itvv1!=periodic_node_sets.end();++itvv1)
   {
      sort((*itvv1).begin(),(*itvv1).end());
   }
   
   // for(itvv1 = periodic_node_sets.begin(); itvv1!=periodic_node_sets.end();++itvv1)
//    {
//       for(itv1 = (*itvv1).begin(); itv1 != (*itvv1).end(); ++itv1)
//          cout<<*itv1<<" ";
//       cout<<endl;   
//    }
   
}

template class Grid<2>;
template class Grid<3>;

