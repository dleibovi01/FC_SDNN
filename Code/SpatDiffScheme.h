/* Spatial Differentiation scheme */

#ifndef SPADIFFSCHEME_H
#define SPADIFFSCHEME_H

// #define MKL_Complex16 std::complex<double>
// #define MKL_Complex8 std::complex<float>
// #include "mkl.h"
#include "VectorField.h"
#include "Mesh.h"
#include "Patch1D.h"
#include <algorithm> 
#include <vector>
#include "FC.h"

template<typename Patch>
class SpatDiffScheme{

};

// template<typename Patch>
class FD1_1D : public SpatDiffScheme<Patch1D>{
 
Patch1D* patch;
double h;

public:

FD1_1D(Patch1D *_patch) : patch{_patch}, h{_patch->getH()} {};
VectorField  diff(const VectorField & v)
{
   VectorField u = circshift(v, -1);
   u = (v - u) / h;
   return u;
}

};


// template<typename Patch>
class WENO_5Z_1D : public SpatDiffScheme<Patch1D>{

Patch1D* patch;
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

   WENO_5Z_1D(Patch1D *_patch, double _epsilon) : patch{_patch}, epsilon{_epsilon}, 
      h{_patch->getH()} {}

   void diff(VectorField *fluxL, VectorField* fluxR, 
      std::vector<VectorField>* data) const
   {
      const int N = fluxR->getLength();
      const int unknowns = fluxR->getUnknowns();

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


#endif 