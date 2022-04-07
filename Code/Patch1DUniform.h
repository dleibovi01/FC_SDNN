/* Uniform 1d patch */

#ifndef PATCH1DUNIFORM_H
#define PATCH1DUNIFORM_H

#include "Patch1D.h"

class Patch1DUniform : public Patch1D {

   double a;
   double b;
public:

   Patch1DUniform(int N, int unknowns, double _a, double _b, int lb, int rb,
      int intrbl, int intrbr);
   double getH() const {return (b - a) / (double(Nnodes - 1));}

};




#endif