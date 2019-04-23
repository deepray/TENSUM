#ifndef __PERTURB_H__
#define __PERTURB_H__

#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>

// input should have two input parameters:
// input[0] -> type of perturbation
// input[1] -> x coordinate scaled to lie in [0 1]
// input[2] -> y coordinate scaled to lie in [0 1]
// input[3] -> additional perturbation parameter
double PERT(const double* input);
#endif
