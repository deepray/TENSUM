// #include <iostream>
// #include <iomanip>
// #include <string>
// #include <cstdlib>
// #include <sstream>
// #include <algorithm>
// #include "grid.h"
// 
// using namespace std;
// 
// // ------------------------------------------------------------------------------
// // Compare faces
// // ------------------------------------------------------------------------------
// template <int DIM>
// int Grid<DIM>::compare_faces (Face<DIM>& master, Face<DIM>& slave)
// {
//     int val = 0;
//     map<int,vector<int> >::iterator it0, it1;
//     if(DIM == 2)
//     {
//        it0 = periodic_nodes.find(slave.vertex[0]);
//        it1 = periodic_nodes.find(slave.vertex[1]);
//        
//        if(find((it0->second).begin(),(it0->second).end(),master.vertex[0])!=(it0->second).end()
//        && find((it1->second).begin(),(it1->second).end(),master.vertex[1])!=(it1->second).end())
//           val = 1;
//        
//        if(find((it0->second).begin(),(it0->second).end(),master.vertex[1])!=(it0->second).end()
//        && find((it1->second).begin(),(it1->second).end(),master.vertex[0])!=(it1->second).end())
//           val = -1;
//        
//        // A fix to avoid matching of the (*) faces when all four corners are periodic
//        // wrt to each other. NOTE THIS FIX WILL NOT WORK FOR GENERAL DOMAINS
//        //     
//        //     X*************X-------------X
//        //     |             |             |
//        //     |             |             |
//        //     |             |             |  
//        //     X-------------X-------------X       INCORRECT MATCH
//        //     |             |             |
//        //     |             |             |
//        //     |             |             | 
//        //     X-------------X*************X
//        //
//        //     X*************X-------------X
//        //     |             |             |
//        //     |             |             |
//        //     |             |             |  
//        //     X-------------X-------------X       CORRECT MATCH
//        //     |             |             |
//        //     |             |             |
//        //     |             |             | 
//        //     X*************X-------------X
//        //
//        // Also there could be issues with the orientation of the face
//        //
//        //     X-------------X-------------X
//        //     ^             |             v
//        //     ^             |             v       INCORRECT ORIENTATION
//        //     ^             |             v  
//        //     X-------------X-------------X       
//        //
//        //     X-------------X-------------X
//        //     ^             |             ^
//        //     ^             |             ^       CORRECT ORIENTATION
//        //     ^             |             ^  
//        //     X-------------X-------------X 
//        //
//        // These issue would occur for a very coarse mesh, which is never considered in 
//        // practice
//        double eps = 1.0e-13;
//        if(val !=0 && (it0->second).size()>1)
//        {
//           int sv0_in = (1 - val)/2;
//           int sv1_in = (1 + val)/2;
//           int mv0 = master.vertex[0];
//           int mv1 = master.vertex[1];
//           int sv0 = slave.vertex[sv0_in];
//           int sv1 = slave.vertex[sv1_in];
//           
//           Vector dr1, dr2,dr3,dr4;
//           dr1.equ(vertex[mv0].coord,vertex[sv0].coord,1.0,-1.0);
//           double d1 = dr1.norm();
//           
//           dr2.equ(vertex[mv1].coord,vertex[sv1].coord,1.0,-1.0);
//           double d2 = dr2.norm();
//           
//           if(abs(d1 - d2) > eps)
//              val = 0;
//           
//           sv0 = slave.vertex[0];
//           sv1 = slave.vertex[1];
//           dr3.equ(vertex[mv0].coord,vertex[mv1].coord,1.0,-1.0);
//           dr4.equ(vertex[sv0].coord,vertex[sv1].coord,1.0,-1.0);
//           if(dr3*dr4 > 0)
//              val = 1;     
//              
//        }
//     }
//     else if(DIM == 3)
//     {
// 	   cout<<"WARNING: 3D face matching not implemented yet. Need to create edges as well"<<endl;
// 	   exit(0);
// 	   // vector<int> master_nodes, slave_nodes;
// // 	   for(int i=0; i< DIM; ++i)
// // 	   {
// // 		  master_nodes.push_back(master.vertex[i]);
// // 		  slave_nodes.push_back(periodic_nodes.find(slave.vertex[i])->second);
// // 	   }
// // 	   sort(master_nodes.begin(), master_nodes.end());
// // 	   sort(slave_nodes.begin(),  slave_nodes.end());
// // 	   if(master_nodes == slave_nodes)
// // 	   {
// // 		  val = 1;
// // 		 
// // 		  // Check orientation
// // 		  int m_v1 = master.vertex[0];
// // 		  int m_v2 = master.vertex[1];
// // 		  int i=0;
// // 		  bool found = false;
// // 		  while(!found)
// // 		  {
// // 			 if(periodic_nodes.find(slave.vertex[i++])->second == m_v1)
// // 				found = true;
// // 			 else
// // 				i++;   
// // 		  } 
// // 		  if(periodic_nodes.find(slave.vertex[i%DIM])->second != m_v2)
// // 			 val = -1;
// // 	   }     
//     }
//     else
//     {
//        cout<<"Error: Only 2D and 3D problems can be solved"<<endl;
//        exit(0);  
//     }
//     return val;
// }
// 
// template class Grid<2>;
// template class Grid<3>;