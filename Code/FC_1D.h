/* 1D FC spatial differentiation schemes */

#ifndef FC_1D_H
#define FC_1D_H

#include <complex>
#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#include "mkl.h"
#include "VectorField.h"
#include "Mesh.h"
#include "Patch1D.h"
#include <algorithm> 
#include <vector>
#include "FC.h"
#include "FC_data.h"
#include "SpatDiffScheme.h"

class FC_1D : public SpatDiffScheme<Patch1D>{

protected:

   Patch1D* patch;
   double * der_coeffs;
   double * der_coeffs_2;
   double * filter_coeffs;
   std::complex<double> * shift_coeffs;
   int N;
   int d; 
   int C; 
   double h;
   int fourPts;
   double fourPts_dbl;
   double prd;
   double * AQ;
   double * FAQF;
   DFTI_DESCRIPTOR_HANDLE desc_handle;

public:

   FC_1D(int _N, int _d, int _C, Patch1D *_patch) : patch{_patch}
   {
      init(_N, _d, _C);
   }

   FC_1D(int _N, int _d, int _C, Patch1D *_patch, int alpha, int p) : patch{_patch}
   {
      init(_N, _d, _C);
      getFiltCoeffs(filter_coeffs, fourPts, alpha, p); 
   }

   FC_1D(int _N, int _d, int _C, Patch1D *_patch, double delta);

   FC_1D(const FC_1D &slv);

   ~FC_1D();

   int getN() const {return N;}
   int getD() const {return d;}
   int getC() const {return C;}
   double getH() const {return h;}
   int getFourPts() const {return fourPts;}
   double getFourPts_dbl() const {return fourPts_dbl;}
   double getPrd() const {return prd;}
   double * getFilterCoeffs() const {return filter_coeffs;}
   double * getDerCoeffs() const {return der_coeffs;}
   double * getDerCoeffs2() const {return der_coeffs_2;}
   std::complex<double> * getShiftCoeffs() const {return shift_coeffs;}
   double * getAQ() const {return AQ;}
   double * getFAQF() const {return FAQF;}
   DFTI_DESCRIPTOR_HANDLE getDescHandle() const {return desc_handle;}
   Patch1D* getPatch() const {return patch;}

   void init(int _N, int _d, int _C);

   virtual void filter(VectorField *v, const std::vector<int> &unknowns, 
      std::vector<std::complex<double> *> *ffts, 
      const std::vector<int> &fft_loc, double h, 
      const std::vector<double> &bc_l, const std::vector<double> &bc_r) const {}

   virtual void diff(const std::complex<double> *y_hat, double * y_der,
      double* y_der_2) const {}

   virtual void diff(const double* y, double *y_der, double h, double bc_l,
      double bc_r) const {}

   virtual void diff(const double* y, double * y_der, double* y_der_2, double h, 
      double bc_l, double bc_r) const {}

};


// FC scheme with left and right Dirichlet boundary conditions
class FC_1D_DD : public FC_1D{

public:

   FC_1D_DD(int _N, int _d, int _C, Patch1D *_patch);

   FC_1D_DD(int _N, int _d, int _C, Patch1D *_patch, int alpha, int p);

   FC_1D_DD(int _N, int _d, int _C, Patch1D *_patch, double delta);

   FC_1D_DD(const FC_1D_DD &slv) : FC_1D(slv) {}

   void diff(const std::complex<double> *y_hat, double * y_der,
      double* y_der_2) const;

   void diff(const double* y, double *y_der, double h, double bc_l,
      double bc_r) const;

   void diff(const double* y, double * y_der, double* y_der_2, double h, 
      double bc_l, double bc_r) const;

   void filter(VectorField *v, const std::vector<int> &unknowns, 
      std::vector<std::complex<double> *> *ffts, 
      const std::vector<int> &fft_loc, double h, 
      const std::vector<double> &bc_l, const std::vector<double> &bc_r) const;

private :

   void set_FC_Data(double* A, double* Q, int d, int C);
};




// FC scheme with left Dirichlet right Neumann boundary conditions
class FC_1D_DN : public FC_1D{

public:

   FC_1D_DN(int _N, int _d, int _C, Patch1D *_patch);

   FC_1D_DN(int _N, int _d, int _C, Patch1D *_patch, int alpha, int p);

   FC_1D_DN(int _N, int _d, int _C, Patch1D *_patch, double delta);

   FC_1D_DN(const FC_1D_DN &slv) : FC_1D(slv) {}

   void diff(const std::complex<double> *y_hat, double * y_der,
      double* y_der_2) const;

   void diff(const double* y, double *y_der, double h, double bc_l,
      double bc_r) const;

   void diff(const double* y, double * y_der, double* y_der_2, double h, 
      double bc_l, double bc_r) const;

   void filter(VectorField *v, const std::vector<int> &unknowns, 
      std::vector<std::complex<double> *> *ffts, 
      const std::vector<int> &fft_loc, double h, 
      const std::vector<double> &bc_l, const std::vector<double> &bc_r) const;

private :

   void set_FC_Data(double* A, double* Q, double* Q_tilde, int d, int C);

};

// FC scheme with left Neumann right Dirichlet boundary conditions
class FC_1D_ND : public FC_1D{

public:

   FC_1D_ND(int _N, int _d, int _C, Patch1D *_patch);

   FC_1D_ND(int _N, int _d, int _C, Patch1D *_patch, int alpha, int p);

   FC_1D_ND(int _N, int _d, int _C, Patch1D *_patch, double delta) ;

   FC_1D_ND(const FC_1D_ND &slv) : FC_1D(slv) {}

   void diff(const std::complex<double> *y_hat, double * y_der,
      double* y_der_2) const;

   void diff(const double* y, double *y_der, double h, double bc_l,
      double bc_r) const;

   void diff(const double* y, double * y_der, double* y_der_2, double h, 
      double bc_l, double bc_r) const;

   void filter(VectorField *v, const std::vector<int> &unknowns, 
      std::vector<std::complex<double> *> *ffts, 
      const std::vector<int> &fft_loc, double h, 
      const std::vector<double> &bc_l, const std::vector<double> &bc_r) const;

private :

   void set_FC_Data(double* A, double* Q, double* Q_tilde, int d, int C);

};


// FC scheme with left and right Neumann boundary conditions
class FC_1D_NN : public FC_1D{

public:

   FC_1D_NN(int _N, int _d, int _C, Patch1D *_patch);

   FC_1D_NN(int _N, int _d, int _C, Patch1D *_patch, int alpha, int p);

   FC_1D_NN(int _N, int _d, int _C, Patch1D *_patch, double delta);

   FC_1D_NN(const FC_1D_NN &slv) : FC_1D(slv) {}

   void diff(const std::complex<double> *y_hat, double * y_der,
      double* y_der_2) const;

   void diff(const double* y, double *y_der, double h, double bc_l,
      double bc_r) const;

   void diff(const double* y, double * y_der, double* y_der_2, double h, 
      double bc_l, double bc_r) const;

   void filter(VectorField *v, const std::vector<int> &unknowns, 
      std::vector<std::complex<double> *> *ffts, 
      const std::vector<int> &fft_loc, double h, 
      const std::vector<double> &bc_l, const std::vector<double> &bc_r) const;

private :

   void set_FC_Data(double* A, double* Q, int d, int C);

};

#endif