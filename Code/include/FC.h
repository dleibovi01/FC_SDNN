/* FC related functions */

#ifndef FC_H
#define FC_H

#include <iostream>
#include "mkl.h"
#include "mkl_dfti.h"
#include <complex.h>
#include <fstream>
#include <string>

void read_FC_Data(double *A, double *Q, int d, int C, std::string filename_A, std::string filename_Q);
//void build_Cont_Mat(const double *A, const double *Q, double *AQ, double *FAQF);


#endif 