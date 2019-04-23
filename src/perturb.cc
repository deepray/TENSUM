#include "perturb.h"
#include "mpi_utils.h"
#include <iostream>
#include <vector>

using namespace std;

extern vector <double> rnd_nos;  // This corresponds to the vector of random numbers
                                      // available for the given sample
                                      
double PERT1(const double* input); 
double PERT2(const double* input); 
double PERT3(const double* input);     
//double PERT4(const double* input);                                  

double PERT(const double* input)
{
   int pert_type = (int) input[0];
   double pert_val = 0.0;
   double s_input[] = {input[1],input[2],input[3]};
   switch(pert_type)
   {
      case 1:
      {
         pert_val = PERT1(s_input);
         break;
      }   
      case 2:
      {
         pert_val = PERT2(s_input);
         break;
      }   
      case 3:
      {
         pert_val = PERT3(s_input);
         break;
      } 
      // case 4:
//       {
//          pert_val = PERT4(s_input);
//          break;
//       }   
      default:
         MPI_LOC_ERR("Unknown perturbation type: " << pert_type);
   }
        
   return pert_val;   
}
  
//------------------------------------------------------------------------------
// 	Multi-modal perturbation based on random inputs and number of modes 
//  M==rnd_nos.size()/2
//  The perturbation parameter corresponds to the fraction of the boundary
//  (centered) where the perturbation must be active.
//------------------------------------------------------------------------------    
double PERT1(const double* input)
{
   int M = (int)(floor(rnd_nos.size()/2.0));
   double y = input[1];
   double frac = max(min(input[2],1.0),0.0); // Making sure the parameter is between 0 and 1
   double llim = (1.0-frac)*0.5;
   double ulim = (frac + 1.0)*0.5;
   
   double a[M], b[M];
   double sum = 0.0;
   double pert_val = 0.0;
   int i=0;
   while(i<M)
   {
      a[i] = rnd_nos[i];
      sum += a[i];
      //cout<<a[i]<<" ";
      i++;
   }
   //cout<<endl;
   i=0;
   while(i<M)
   {   
      b[i] =rnd_nos[M+i];
      //cout<<b[i]<<" ";
      i++;
   }
   //cout<<endl<<endl;
   
   if(y > llim && y < ulim)
   {
	  for(int i=0; i<M; ++i)
		 pert_val += (a[i]/sum)*sin(2*M_PI*y + b[i]);
   }   
   
   //cout<<pert_val<<endl;     
   return pert_val;   
}

//------------------------------------------------------------------------------
// M partitions == rnd_nos.size() in y direction of total length H.
//  The perturbation parameter corresponds to the fraction of the boundary
//  (centered) where the perturbation must be active.
//------------------------------------------------------------------------------
// double PERT2(const double* input)
// {
//    int M = rnd_nos.size();
//    double y = input[1]; // scaled to 1
//    double div_size = (1.0/double(M));
//    int div = (int)floor(y/div_size);
//    if(div == M)
//       div--;
//    int sgn = pow(-1,div);
//    double pert_val =rnd_nos[div]*sgn;
//    //cout <<pert_val<<endl;
//    return pert_val;   
// }

//------------------------------------------------------------------------------
// M partitions == rnd_nos.size() in y direction of total length H.
//  The perturbation parameter corresponds to the fraction of the boundary
//  (centered) where the perturbation must be active.
//------------------------------------------------------------------------------
double PERT2(const double* input)
{
   int M = rnd_nos.size();
   double y = input[1]; // scaled to 1
   double frac = max(min(input[2],1.0),0.0); // Making sure the parameter is between 0 and 1
   double llim = (1.0-frac)*0.5;
   double ulim = (frac + 1.0)*0.5;
   
   double div_size = (frac/double(M));
   double pert_val;
   if(y>llim && y<ulim)
   {
	  y = y - llim;
	  int div = (int)floor(y/div_size);
	  if(div == M)
		 div--;
	  int sgn = pow(-1,div);
	  pert_val =rnd_nos[div]*sgn;
   }
   else
      pert_val = 0.0;
   //cout <<pert_val<<endl;
   return pert_val;   
}

//------------------------------------------------------------------------------
// 	Multi-modal perturbation based on random input numbers and number of modes 
//  M==rnd_nos.size()/2. In this case different wavenumber are also used
//  The perturbation parameter corresponds to the fraction of the boundary
//  (centered) where the perturbation must be active.
//------------------------------------------------------------------------------    
double PERT3(const double* input)
{
   int M = (int)(floor(rnd_nos.size()/2.0));
   double y = input[1];
   double frac = max(min(input[2],1.0),0.0); // Making sure the parameter is between 0 and 1
   double llim = (1.0-frac)*0.5;
   double ulim = (frac + 1.0)*0.5;
   
   double a[M], b[M];
   double sum = 0.0;
   double pert_val = 0.0;
   int i=0;
   while(i<M)
   {
      a[i] = rnd_nos[i];
      sum += a[i];
      //cout<<a[i]<<" ";
      i++;
   }
   //cout<<endl;
   i=0;
   while(i<M)
   {   
      b[i] =rnd_nos[M+i];
      //cout<<b[i]<<" ";
      i++;
   }
   //cout<<endl<<endl;
   
   if(y > llim && y < ulim)
   {
	  for(int i=0; i<M; ++i)
		 pert_val += (a[i]/sum)*sin(2*i*M_PI*y + b[i]);
   }   
   //cout<<pert_val<<endl;     
   return pert_val;   
}

