/* 2D FC spatial differentiation schemes */

#ifndef FC_2D_H
#define FC_2D_H

#include "FC.h"
#include "FC_1D.h"
#include "Patch2D.h"
#include "SpatDiffScheme.h"
#include <string.h>
#include "printing.h"


class FC_2D : public SpatDiffScheme<Patch2D>{

FC_1D * FC_x;
FC_1D * FC_y;
int Nx;
int Ny;

public:

   FC_2D(const char* xx, const char*  yy, int _Nx, int dx, int Cx, double hx, 
        int _Ny, int dy, int Cy, double hy, int alpha, int p);

   FC_2D(const char*  xx, const char*  yy, int _Nx, int dx, int Cx, double hx, 
        int _Ny, int dy, int Cy, double hy, double delta);

   ~FC_2D();

   FC_1D * getFCx() const {return FC_x;}
   FC_1D * getFCy() const {return FC_y;}
   int getNx() const {return Nx;}
   int getNy() const {return Ny;}

   void diff_y(const std::complex<double> *y_hat, double * y_der,
      double* y_der_2) const;

   void diff_y(const double* y, double *y_der, double h, const double* bc_d,
      const double* bc_u) const;

   void diff_y(const double* y, double * y_der, double* y_der_2, double h, 
      const double* bc_d, const double* bc_u) const;

   void filter_y(VectorField *v, const std::vector<int> &unknowns, 
      std::vector<std::complex<double> *> *ffts, 
      const std::vector<int> &fft_loc, double h, 
      const std::vector<double* > &bc_d, const std::vector<double* > &bc_u) const;       

   void diff_x(const std::complex<double> *y_hat, double * y_der,
      double* y_der_2) const;

   void diff_x(const double* y, double *y_der, double h, const double* bc_l,
      const double* bc_r) const;

   void diff_x(const double* y, double * y_der, double* y_der_2, double h, 
      const double* bc_l, const double* bc_r) const;

   void filter_x(VectorField *v, const std::vector<int> &unknowns, 
      std::vector<std::complex<double> *> *ffts, 
      const std::vector<int> &fft_loc, double h, 
      const std::vector<double* > &bc_l, const std::vector<double* > &bc_r) const; 

};



#endif