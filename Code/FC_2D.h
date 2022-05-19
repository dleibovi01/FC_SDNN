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


   FC_2D(const std::string xx, const std::string yy, int _Nx, int dx, int Cx,
      double hx, int _Ny, int dy, int Cy, double hy, int alpha, int p);

   FC_2D(const std::string xx, const std::string yy, int _Nx, int dx, int Cx,
      double hx, int _Ny, int dy, int Cy, double hy, double delta);

   ~FC_2D();

   FC_1D * getFCx() const {return FC_x;}

   FC_1D * getFCy() const {return FC_y;}

   int getNx() const {return Nx;}

   int getNy() const {return Ny;}

   int getCx() const {return FC_x->getC();}

   int getCy() const {return FC_y->getC();}

   int getDx() const {return FC_x->getD();}

   int getDy() const {return FC_y->getD();}

   double getFourPtsDblx() const {return FC_x->getFourPts_dbl();}

   double getFourPtsDbly() const {return FC_y->getFourPts_dbl();}

   double * getAQx() const {return FC_x->getAQ();}

   double * getFAQFx() const {return FC_x->getFAQF();} 

   double * getAQy() const {return FC_y->getAQ();}

   double * getFAQFy() const {return FC_y->getFAQF();} 

   std::complex<double>* getShiftCoeffsx() const{return FC_x->getShiftCoeffs();}

   std::complex<double>* getShiftCoeffsy() const{return FC_y->getShiftCoeffs();}

   DFTI_DESCRIPTOR_HANDLE getDescHandlex() const {return FC_x->getDescHandle();}

   DFTI_DESCRIPTOR_HANDLE getDescHandley() const {return FC_y->getDescHandle();}

   void diff_y(const std::complex<double> *y_hat, double * y_der,
      double* y_der_2) const;

   void diff_y(const double* y, double *y_der, double h, const double* bc_d,
      const double* bc_u) const;

   void diff_y(const double* y, double * y_der, double* y_der_2, double h, 
      const double* bc_d, const double* bc_u) const;

   void filter_y(VectorField *v, const std::vector<int> &unknowns, 
      std::vector<std::complex<double> *> *ffts, 
      const std::vector<int> &fft_loc, double h, 
      const std::vector<double* > &bc_d, const std::vector<double* > &bc_u) 
      const;      

   void filter_y(VectorField *v, const std::vector<int> &unknowns, double h, 
      const std::vector<double* > &bc_d, const std::vector<double* > &bc_u) 
      const;   

   void diff_x(const std::complex<double> *y_hat, double * y_der,
      double* y_der_2) const;

   void diff_x(const double* y, double *y_der, double h, const double* bc_l,
      const double* bc_r) const;

   void diff_x(const double* y, double * y_der, double* y_der_2, double h, 
      const double* bc_l, const double* bc_r) const;

   void filter_x(VectorField *v, const std::vector<int> &unknowns, 
      std::vector<std::complex<double> *> *ffts, 
      const std::vector<int> &fft_loc, double h, 
      const std::vector<double* > &bc_l, const std::vector<double* > &bc_r)
      const;

   void filter_x(VectorField *v, const std::vector<int> &unknowns, double h, 
      const std::vector<double* > &bc_l, const std::vector<double* > &bc_r)
      const;

   void shift1D_x(const double* y, double * y_shift, double h, double bc_l,
      double bc_r) const {FC_x->shift(y, y_shift, h, bc_l,bc_r);}

   void shift1D_y(const double* y, double * y_shift, double h, double bc_d,
      double bc_u) const {FC_y->shift(y, y_shift, h, bc_d, bc_u);}

   void filter(VectorField *v, const std::vector<int> &unknowns, double hx,
      double hy, const std::vector<double* > &bc_l,
      const std::vector<double* > &bc_r, const std::vector<double* > &bc_d, 
      const std::vector<double* > &bc_u) const;

};

#endif