/* Spatial Differentiation scheme */

#ifndef SPADIFFSCHEME_H
#define SPADIFFSCHEME_H

#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#include "mkl.h"
#include "VectorField1D.h"
#include "Mesh.h"
#include <algorithm> 
#include <vector>
#include "FC.h"
#include "FC_data.h"
// #define MKL_Complex16 std::complex<double>

template<typename Patch>
class SpatDiffScheme{

};

template<typename Patch>
class FD1_1D : public SpatDiffScheme<Patch>{
 
Patch* patch;
double h;

public:

FD1_1D(Patch *_patch) : patch{_patch}, h{_patch->getH()} {};
VectorField1D  diff(const VectorField1D & v)
{
   VectorField1D u = circshift(v, -1);
   u = (v - u) / h;
   return u;
}

};


template<typename Patch, typename VectorField>
class WENO_5Z_1D : public SpatDiffScheme<Patch>{

Patch* patch;
double epsilon;
double h;
// std::vector<VectorField> data;

static constexpr double d0p = 0.3;
static constexpr double d1p = 0.6;
static constexpr double d2p = 0.1;

static constexpr double d0n = 0.1;
static constexpr double d1n = 0.6;
static constexpr double d2n = 0.3;


public :

WENO_5Z_1D(Patch *_patch, double _epsilon) : patch{_patch}, epsilon{_epsilon}, 
   h{_patch->getH()} 
{
   // int ndata = 11;
   // for(int i = 0; i < ndata; i++)
   // {
   //    data.push_back(VectorField1D(patch->getUnknowns(), patch->getNnodes()));
   // }
}

void diff(VectorField1D *fluxL, VectorField1D* fluxR, 
   std::vector<VectorField1D>* data) const
{
   const int N = fluxR->getLength();
   const int unknowns = fluxR->getUnknowns();

   // Setting some vectors to hold the data
   // int ndata = 11;
   // std::vector<VectorField1D> data;
   // for(int i = 0; i < ndata; i++)
   // {
   //    data.push_back(VectorField1D(unknowns, N));
   // }

   ///////////////////////////////////////////////////////////
   // Left flux treatment
   //////////////////////////////////////////////////////////

   data->at(0) = circshift(*fluxL, -2);
   data->at(1) = circshift(*fluxL, -1);
   data->at(2) = circshift(*fluxL, 1);
   data->at(3) = circshift(*fluxL, 2);


   // Polynomials
   // data->at(4) = (1.0/6.0)*(5.0*data->at(1) - data->at(0) + 2.0*(*fluxL));
   // data->at(5) = (1.0/6.0)*(2.0*data->at(1) + 5.0*(*fluxL) - data->at(2));
   // data->at(6) = (1.0/6.0)*(11.0*(*fluxL) - 7.0*data->at(2) + 2.0*data->at(3));
   linComb3(data->at(1), data->at(0), *fluxL, &(data->at(4)), 5.0/6.0, -1.0/6.0,
      2.0/6.0);
   linComb3(data->at(1), *fluxL, data->at(2), &(data->at(5)), 2.0/6.0, 5.0/6.0,
      - 1.0/6.0);
   linComb3(*fluxL, data->at(2), data->at(3), &(data->at(6)), 11.0/6.0, 
      - 7.0/6.0, 2.0/6.0);


   // Smooth indicators
   // data->at(7) = (13.0/12.0)*sqr(data->at(0) - 2.0*data->at(1) + (*fluxL))
   //    + 0.25*sqr(data->at(0) - 4.0*data->at(1) + 3.0*(*fluxL));
   // data->at(8) = (13.0/12.0)*sqr(data->at(1) - 2.0*(*fluxL) + data->at(2))
   //    + 0.25*sqr(data->at(1) - data->at(2));
   // data->at(9) = (13.0/12.0)*sqr((*fluxL) - 2.0*data->at(2) + data->at(3))
   //    + 0.25*sqr(3.0*(*fluxL) - 4.0*data->at(2) + data->at(3));
   linComb3(data->at(0), data->at(1), *fluxL, &(data->at(7)), 1.0, -2.0, 1.0);
   linComb3(data->at(0), data->at(1), *fluxL, &(data->at(0)), 1.0, -4.0, 3.0);
   linComb2(data->at(7).sqr(), data->at(0).sqr(), &(data->at(7)), 13.0/12.0,
      0.25);
   linComb3(data->at(1), *fluxL, data->at(2), &(data->at(8)), 1.0, -2.0, 1.0);
   linComb2(data->at(1), data->at(2), &(data->at(1)), 1.0, -1.0);
   linComb2(data->at(8).sqr(), data->at(1).sqr(), &(data->at(8)), 13.0/12.0,
      0.25);
   linComb3(*fluxL, data->at(2), data->at(3), &(data->at(9)), 1.0, -2.0, 1.0);
   linComb3(*fluxL, data->at(2), data->at(3), &(data->at(2)), 3.0, -4.0, 1.0);
   linComb2(data->at(9).sqr(), data->at(2).sqr(), &(data->at(9)), 13.0/12.0,
      0.25);

   // tau5
   // data->at(10) = abs(data->at(7) - data->at(9));
   linComb2(data->at(7), data->at(9), &(data->at(10)), 1.0, -1.0);
   (data->at(10)).abs();


   // alpha weights
   // data->at(7) = d0p*(1.0 + data->at(10)/(data->at(7) + epsilon));
   // data->at(8) = d1p*(1.0 + data->at(10)/(data->at(8) + epsilon));
   // data->at(9) = d2p*(1.0 + data->at(10)/(data->at(9) + epsilon));
   // data->at(10) = data->at(7) + data->at(8) + data->at(9);
   data->at(7) += epsilon;
   data->at(7) = data->at(10)/data->at(7);
   data->at(7) += 1.0;
   data->at(7) *= d0p;
   data->at(8) += epsilon;
   data->at(8) = data->at(10)/data->at(8);
   data->at(8) += 1.0;
   data->at(8) *= d1p;
   data->at(9) += epsilon;
   data->at(9) = data->at(10)/data->at(9);
   data->at(9) += 1.0;
   data->at(9) *= d2p;
   linComb3(data->at(7), data->at(8), data->at(9), &(data->at(10)), 1.0, 1.0,
      1.0);


   // ENO stencils weights
   data->at(7) = data->at(7)/data->at(10);
   data->at(8) = data->at(8)/data->at(10);
   data->at(9) = data->at(9)/data->at(10);


   // Numerical flux at cell boundary
   // fluxL->setFlow((data->at(7)*data->at(4) + data->at(8)*data->at(5) + data->at(9)*data->at(6)).getFlow());
   data->at(7) *= data->at(4);
   data->at(8) *= data->at(5);
   data->at(9) *= data->at(6);
   linComb3(data->at(7), data->at(8), data->at(9), fluxL, 1.0, 1.0, 1.0);




   ///////////////////////////////////////////////////////////
   // Right flux treatment
   //////////////////////////////////////////////////////////

   // Rotating vectors
   data->at(0) = circshift(*fluxR, -2);
   data->at(1) = circshift(*fluxR, -1);
   data->at(2) = circshift(*fluxR, 1);
   data->at(3) = circshift(*fluxR, 2);

   // Polynomials
   // data->at(4) = (1.0/6.0)*(2.0*data->at(0) - 7.0*data->at(1) + 11.0*(*fluxR));
   // data->at(5) = (1.0/6.0)*(5.0*(*fluxR) - data->at(1) + 2.0*data->at(2));
   // data->at(6) = (1.0/6.0)*(2.0*(*fluxR) + 5.0*data->at(2) - data->at(3));
   linComb3(data->at(1), data->at(0), *fluxR, &(data->at(4)), -7.0/6.0, 2.0/6.0,
      11.0/6.0);
   linComb3(data->at(1), *fluxR, data->at(2), &(data->at(5)), -1.0/6.0, 5.0/6.0,
      2.0/6.0);
   linComb3(*fluxR, data->at(2), data->at(3), &(data->at(6)), 2.0/6.0, 5.0/6.0,
      -1.0/6.0);

   // Smooth indicators
   // data->at(7) = (13.0/12.0)*sqr(data->at(0) - 2.0*data->at(1) + (*fluxR))
   //    + 0.25*sqr(data->at(0) - 4.0*data->at(1) + 3.0*(*fluxR));
   // data->at(8) = (13.0/12.0)*sqr(data->at(1) - 2.0*(*fluxR) + data->at(2))
   //    + 0.25*sqr(data->at(1) - data->at(2));
   // data->at(9) = (13.0/12.0)*sqr((*fluxR)- 2.0*data->at(2) + data->at(3))
   //    + 0.25*sqr(3*(*fluxR) - 4.0*data->at(2) + data->at(3));
   linComb3(data->at(0), data->at(1), *fluxR, &(data->at(7)), 1.0, -2.0, 1.0);
   linComb3(data->at(0), data->at(1), *fluxR, &(data->at(0)), 1.0, -4.0, 3.0);
   linComb2(data->at(7).sqr(), data->at(0).sqr(), &(data->at(7)), 13.0/12.0,
      0.25);
   linComb3(data->at(1), *fluxR, data->at(2), &(data->at(8)), 1.0, -2.0, 1.0);
   linComb2(data->at(1), data->at(2), &(data->at(1)), 1.0, -1.0);
   linComb2(data->at(8).sqr(), data->at(1).sqr(), &(data->at(8)), 13.0/12.0,
      0.25);
   linComb3(*fluxR, data->at(2), data->at(3), &(data->at(9)), 1.0, -2.0, 1.0);
   linComb3(*fluxR, data->at(2), data->at(3), &(data->at(2)), 3.0, -4.0, 1.0);
   linComb2(data->at(9).sqr(), data->at(2).sqr(), &(data->at(9)), 13.0/12.0,
      0.25);

   // tau5
   // data->at(10) = abs(data->at(7) - data->at(9));
   linComb2(data->at(7), data->at(9), &(data->at(10)), 1.0, -1.0);
   (data->at(10)).abs();

   // alpha weights
   // data->at(7) = d0n*(1.0 + data->at(10)/(data->at(7) + epsilon));
   // data->at(8) = d1n*(1.0 + data->at(10)/(data->at(8) + epsilon));
   // data->at(9) = d2n*(1.0 + data->at(10)/(data->at(9) + epsilon));
   // data->at(10) = data->at(7) + data->at(8) + data->at(9);
   data->at(7) += epsilon;
   data->at(7) = data->at(10)/data->at(7);
   data->at(7) += 1.0;
   data->at(7) *= d0n;
   data->at(8) += epsilon;
   data->at(8) = data->at(10)/data->at(8);
   data->at(8) += 1.0;
   data->at(8) *= d1n;
   data->at(9) += epsilon;
   data->at(9) = data->at(10)/data->at(9);
   data->at(9) += 1.0;
   data->at(9) *= d2n;
   linComb3(data->at(7), data->at(8), data->at(9), &(data->at(10)), 1.0, 1.0,
      1.0);

   // ENO stencils weights
   data->at(7) = data->at(7)/data->at(10);
   data->at(8) = data->at(8)/data->at(10);
   data->at(9) = data->at(9)/data->at(10);

   // Numerical flux at cell boundary
   // fluxR->setFlow((data->at(7)*data->at(4) + data->at(8)*data->at(5)
   //    + data->at(9)*data->at(6)).getFlow());
   data->at(7) *= data->at(4);
   data->at(8) *= data->at(5);
   data->at(9) *= data->at(6);
   linComb3(data->at(7), data->at(8), data->at(9), fluxR, 1.0, 1.0, 1.0);

   // Differentiation
   // fluxR->setFlow(((1.0/h) * ((*fluxR) - circshift(*fluxR, -1) +
   //    (*fluxL) - circshift(*fluxL, -1))).getFlow());
   linComb4(*fluxR, circshift(*fluxR, -1), *fluxL, circshift(*fluxL, -1), fluxR,
      1.0/h, -1.0/h, 1.0/h, -1.0/h);

}


};

template<typename Patch>
class FC_1D : public SpatDiffScheme<Patch>{

Patch* patch;
// std::complex<double> * der_coeffs;
// std::complex<double> * der_coeffs_2;
// std::complex<double> * filter_coeffs;
double * der_coeffs;
double * der_coeffs_2;
double * filter_coeffs;
std::complex<double> * shift_coeffs;
int N;
int d; 
int C; 
int fourPts;
double fourPts_dbl;
double prd;
double * AQ;
double * FAQF;
DFTI_DESCRIPTOR_HANDLE desc_handle;

public:

FC_1D(int _N, int _d, int _C, Patch *_patch, std::string filename_A, 
   std::string filename_Q) : patch{_patch}
{
   N = _N;
   d = _d;
   C = _C;
   fourPts = N + C;
   fourPts_dbl = double(fourPts);
   double h = patch->getH();
   prd = fourPts_dbl*h;
   // der_coeffs = new std::complex<double> [N+C];
   // der_coeffs_2 = new std::complex<double> [N+C];
   // filter_coeffs = new std::complex<double> [N+C];
   der_coeffs = double[N+C];
   der_coeffs_2 = double[N+C];
   filter_coeffs = new double [N+C];
   double A[C*d];
   double Q[d*d];   
   AQ = new double [C*d];
   FAQF = new double [C*d];
   MKL_LONG status;
   const double pi = std::acos(-1);
   const std::complex<double> I(0, 1); 
   int k[N + C];
   getK(k, N + C);
   // for(int j = 0; j < N + C; j++)
   // {
   //    der_coeffs[j] = 2.0*pi/prd*I*double(k[j]);
   // } 
   for(int j = 0; j < N + C; j++)
   {
      der_coeffs[j] = 2.0*pi/prd*double(k[j]);
   }         
   for(int j = 0; j < N + C; j++)
   {
      filter_coeffs[j] = 1.0;
   }    
   status = DftiCreateDescriptor(&desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1,
      fourPts); 
   status = DftiCommitDescriptor(desc_handle);
   read_FC_Data(A, Q, d, C, filename_A, filename_Q);
   build_Cont_Mat(A, Q, d, C, AQ, FAQF);
}

FC_1D(int _N, int _d, int _C, Patch *_patch) : patch{_patch}
{
   init(_N, _d, _C);
}

FC_1D(int _N, int _d, int _C, Patch *_patch, int alpha, int p) : patch{_patch}
{
   init(_N, _d, _C);
   getFiltCoeffs(filter_coeffs, fourPts, alpha, p); 
}

FC_1D(int _N, int _d, int _C, Patch *_patch, double delta) : patch{_patch}
{
   init(_N, _d, _C);
   int k[fourPts];
   getK(k, fourPts);
   const double pi = std::acos(-1);
   const std::complex<double> I(0, 1);     
   for(int i = 0; i < fourPts; i++)
   {
      shift_coeffs[i] = std::exp(2.0*pi*I* double(k[i]) * delta / 
         double(fourPts));
   }

}


FC_1D(const FC_1D &slv)
{
   patch = slv.getPatch();
   N = slv.getN();
   d = slv.getD();
   C = slv.getC();
   fourPts = slv.getFourPts();
   fourPts_dbl = slv.fourPts_dbl;
   prd = slv.getPrd();
   desc_handle = slv.getDescHandle();
   // der_coeffs = new std::complex<double> [N+C];
   // der_coeffs_2 = new std::complex<double> [N+C];
   // filter_coeffs = new std::complex<double> [N+C]; 
   der_coeffs = new double [N+C];
   der_coeffs_2 = new double [N+C];
   filter_coeffs = new double [N+C]; 
   shift_coeffs = new std::complex<double> [N+C];
   AQ = new double [C*d];
   FAQF = new double [C*d];
   // std::complex<double> * d_c = slv.getDerCoeffs();
   double * d_c = slv.getDerCoeffs();
   std::copy(d_c, d_c + N + C, der_coeffs);
   // std::complex<double> * d_c_2 = slv.getDerCoeffs2();
   double * d_c_2 = slv.getDerCoeffs2();
   std::copy(d_c_2, d_c_2 + N + C, der_coeffs_2);
   // std::complex<double> * f_c = slv.getFilterCoeffs();
   double * f_c = slv.getFilterCoeffs();
   std::copy(f_c, f_c + N + C, filter_coeffs);  
   std::complex<double> * s_c = slv.getShiftCoeffs();
   std::copy(s_c, s_c + N + C, shift_coeffs);  
   double * aq = slv.getAQ();
   std::copy(aq, aq + C*d, AQ);    
   double * faqf = slv.getFAQF();
   std::copy(faqf, faqf + C*d, FAQF);      
}


~FC_1D()
{
   delete [] AQ;
   delete [] FAQF;
   delete [] filter_coeffs;
   delete [] der_coeffs;
   delete [] der_coeffs_2;
   delete [] shift_coeffs;
}


int getN() const {return N;}
int getD() const {return d;}
int getC() const {return C;}
int getFourPts() const {return fourPts;}
double getFourPts_dbl() const {return fourPts_dbl;}
double getPrd() const {return prd;}
// std::complex<double> * getFilterCoeffs() const {return filter_coeffs;}
double * getFilterCoeffs() const {return filter_coeffs;}
// std::complex<double> * getDerCoeffs() const {return der_coeffs;}
// std::complex<double> * getDerCoeffs2() const {return der_coeffs_2;}
double * getDerCoeffs() const {return der_coeffs;}
double * getDerCoeffs2() const {return der_coeffs_2;}
std::complex<double> * getShiftCoeffs() const {return shift_coeffs;}
double * getAQ() const {return AQ;}
double * getFAQF() const {return FAQF;}
DFTI_DESCRIPTOR_HANDLE getDescHandle() const {return desc_handle;}
Patch* getPatch() const {return patch;}

// VectorField1D  diff(const VectorField1D & v) const
// {
//    VectorField1D u{v.getUnknowns(), v.getLength()};
//    int unknowns = v.getUnknowns();
//    double data[v.getLength()];
//    for(int i = 0; i < unknowns; i++)
//    {
//       FC_Der(v.getField(i), data, der_coeffs, filter_coeffs, N, d, C,
//          fourPts_dbl, AQ, FAQF, desc_handle);
//       u.setField(i, N, data);
      
//    }
//    return u;
// }

void diff(const double* y, double *y_der) const
{
   // FC_Der(y, y_der, der_coeffs, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle);
   std::complex<double> y_hat[N + C];
   Fcont_Gram_Blend(y, y_hat, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle);
   FC_Der(y_der, y_hat, der_coeffs, N, C, desc_handle);
}

void diff(const double* y, double * y_der, double* y_der_2) const
{
   std::complex<double> y_hat[N + C];
   Fcont_Gram_Blend(y, y_hat, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle);
   diff(y_hat, y_der, y_der_2);
}

void diff(const std::complex<double> *y_hat, double * y_der,
   double* y_der_2) const
{
   FC_Der(y_der, y_hat, der_coeffs, N, C, desc_handle);
   FC_Der(y_der_2, y_hat, der_coeffs_2, N, C, desc_handle, true);
}


void set_FC_Data(double* A, double* Q, int d, int C)
{
   if(C == 27)
   {
      if(d == 5)
      {
         std::copy(Ad5C27_data.data(), Ad5C27_data.data() + d*C, A);
         std::copy(Qd5C27_data.data(), Qd5C27_data.data() + d*d, Q);
      }
   }
}

// VectorField1D filter(const VectorField1D & v) const
// {
//    VectorField1D u{v};
//    int unknowns = v.getUnknowns();
//    int length = v.getLength();   
//    double data[length];  
//    for(int i = 0; i < unknowns; i++)
//    {
//       FC_Der(v.getField(i), data, filter_coeffs, filter_coeffs, length, d, C,
//          fourPts_dbl, AQ, FAQF, desc_handle);
//       u.setField(i, length, data);
//    } 
//    return u;
// }

// template<typename VectorField>
// void filter(VectorField *v, const std::vector<int> &unknowns) const
// {
//    int length = v->getLength();    
//    for(int i = 0; i < unknowns.size(); i++)
//    {
//       FC_Der(v->getField(unknowns[i]), v->getField(unknowns[i]), filter_coeffs,
//          filter_coeffs, length, d, C, fourPts_dbl, AQ, FAQF, desc_handle);
//    } 
// }

template<typename VectorField>
void filter(VectorField *v, const std::vector<int> &unknowns, 
   std::vector<std::complex<double> *> *ffts, 
   const std::vector<int> &fft_loc) const
{
   int length = v->getLength();  
   std::complex<double> f_ext[length + C]; 
   for(int i = 0; i < unknowns.size(); i++)
   {
      Fcont_Gram_Blend(v->getField(unknowns[i]), ffts->at(fft_loc[i]), N, d, C,
         fourPts_dbl, AQ, FAQF, desc_handle);
      // VectorMul(N + C, ffts->at(fft_loc[i]), filter_coeffs, ffts->at(fft_loc[i]));
      VectorMulReCmp(N + C, filter_coeffs, ffts->at(fft_loc[i]), ffts->at(fft_loc[i]));
      std::copy(ffts->at(fft_loc[i]), ffts->at(fft_loc[i]) + N + C, f_ext);
      int status = DftiComputeBackward(desc_handle, f_ext);
      for (int j = 0; j < length; j++)
      {
         v->getField(unknowns[i])[j] = f_ext[j].real();
      }     

   } 
}


private:

   void init(int _N, int _d, int _C)
   {
      N = _N;
      d = _d;
      C = _C;
      fourPts = N + C;
      fourPts_dbl = double(fourPts);
      double h = patch->getH();
      prd = fourPts_dbl*h;
      // der_coeffs = new std::complex<double> [N+C];
      // der_coeffs_2 = new std::complex<double> [N+C];
      // filter_coeffs = new std::complex<double> [N+C];
      der_coeffs = new double [N+C];
      der_coeffs_2 = new double [N+C];
      filter_coeffs = new double [N+C];
      shift_coeffs = new std::complex<double> [N+C];
      double A[C*d];
      double Q[d*d];   
      AQ = new double [C*d];
      FAQF = new double [C*d];
      MKL_LONG status;
      const double pi = std::acos(-1);
      const std::complex<double> I(0, 1); 
      int k[N + C];
      getK(k, N + C);
      for(int j = 0; j < N + C; j++)
      {
         // der_coeffs[j] = 2.0*pi/prd*I*double(k[j]);
         der_coeffs[j] = 2.0*pi/prd*double(k[j]);
         der_coeffs_2[j] = -4.0*pi*pi/prd/prd*double(k[j])*double(k[j]);
      }      
      for(int j = 0; j < N + C; j++)
      {
         filter_coeffs[j] = 1.0;
      }    
      status = DftiCreateDescriptor(&desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1,
         fourPts); 
      status = DftiCommitDescriptor(desc_handle);
      // read_FC_Data(A, Q, d, C, filename_A, filename_Q);
      set_FC_Data(A, Q, d, C);
      build_Cont_Mat(A, Q, d, C, AQ, FAQF);      
   }

};


#endif 