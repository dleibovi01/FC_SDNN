/* Spatial Differentiation scheme */

#ifndef SPADIFFSCHEME_H
#define SPADIFFSCHEME_H

#include "VectorField1D.h"
#include "Mesh.h"
#include <algorithm> 
#include "FC.h"
#include "FC_data.h"

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

template<typename Patch>
class FC_1D : public SpatDiffScheme<Patch>{

Patch* patch;
std::complex<double> * der_coeffs;
std::complex<double> * filter_coeffs;
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
   der_coeffs = new std::complex<double> [N+C];
   filter_coeffs = new std::complex<double> [N+C];
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
      der_coeffs[j] = 2.0*pi/prd*I*double(k[j]);
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
   // std::complex<double> shift_coeffs[fourPts];
   // std::complex<double> filter_coeffs[fourPts];
   const double pi = std::acos(-1);
   const std::complex<double> I(0, 1);     
   for(int i = 0; i < fourPts; i++)
   {
      shift_coeffs[i] = std::exp(2.0*pi*I* double(k[i]) * delta / 
         double(fourPts));
   }
   // Print_Mat(shift_coeffs, fourPts, 1);
   // std::cout << std::endl;
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
   der_coeffs = new std::complex<double> [N+C];
   filter_coeffs = new std::complex<double> [N+C]; 
   shift_coeffs = new std::complex<double> [N+C];
   AQ = new double [C*d];
   FAQF = new double [C*d];
   std::complex<double> * d_c = slv.getDerCoeffs();
   std::copy(d_c, d_c + N + C, der_coeffs);
   std::complex<double> * f_c = slv.getFilterCoeffs();
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
   delete [] shift_coeffs;
}


int getN() const {return N;}
int getD() const {return d;}
int getC() const {return C;}
int getFourPts() const {return fourPts;}
double getPrd() const {return prd;}
std::complex<double> * getFilterCoeffs() const {return filter_coeffs;}
std::complex<double> * getDerCoeffs() const {return der_coeffs;}
std::complex<double> * getShiftCoeffs() const {return shift_coeffs;}
double * getAQ() const {return AQ;}
double * getFAQF() const {return FAQF;}
DFTI_DESCRIPTOR_HANDLE getDescHandle() const {return desc_handle;}
Patch* getPatch() const {return patch;}

VectorField1D  diff(const VectorField1D & v) const
{
   VectorField1D u{v.getUnknowns(), v.getLength()};
   int unknowns = v.getUnknowns();
   double data[v.getLength()];
   for(int i = 0; i < unknowns; i++)
   {
      FC_Der(v.getField(i), data, der_coeffs, filter_coeffs, N, d, C,
         fourPts_dbl, AQ, FAQF, desc_handle);
      u.setField(i, N, data);
      
   }
   return u;
}

void diff(const double* y, double *y_der) const
{
   // std::cout << "fourPts = " << fourPts << std::endl;
   // std::cout << "fourPts_dbl = " << fourPts_dbl << std::endl;
   FC_Der(y, y_der, der_coeffs, filter_coeffs, N, d, C, fourPts_dbl, AQ,
      FAQF, desc_handle);
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

VectorField1D filter(const VectorField1D & v) const
{
   VectorField1D u{v};
   int unknowns = v.getUnknowns();
   int length = v.getLength();   
   double data[length];  
   for(int i = 0; i < unknowns; i++)
   {
      FC_Der(v.getField(i), data, filter_coeffs, filter_coeffs, length, d, C,
         fourPts_dbl, AQ, FAQF, desc_handle);
      u.setField(i, length, data);
   } 
   return u;
}

template<typename VectorField>
void filter(VectorField *v, const std::vector<int> &unknowns) const
{
   int length = v->getLength();   
   double data[length];  
   for(int i = 0; i < unknowns.size(); i++)
   {
      FC_Der(v->getField(unknowns[i]), data, filter_coeffs, filter_coeffs,
         length, d, C, fourPts_dbl, AQ, FAQF, desc_handle);
      v->setField(unknowns[i], length, data);
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
      der_coeffs = new std::complex<double> [N+C];
      filter_coeffs = new std::complex<double> [N+C];
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
         der_coeffs[j] = 2.0*pi/prd*I*double(k[j]);
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