/* 2D patch */

#ifndef PATCH2D_H
#define PATCH2D_H

#include "Node.h"
#include <iostream>
#include "Patch.h"
#include <algorithm>

// #include "printing.h"

class Patch2D : public Patch{

protected:

int Nx;
int Ny;
double hx;
double hy;
int fringe;

public:

    Patch2D(){};
    Patch2D(const Patch2D &patch) : Patch(patch), Nx{patch.Nx}, Ny{patch.Ny}, 
        fringe{patch.fringe} {}
    Patch2D(int _Nx, int _Ny, int _unknowns, double a, double b, double c,
        double d, bool lb, bool rb, bool db, bool ub, int fringe);
    Patch2D & operator=(const Patch2D &patch);
    int getNx() const {return Nx;}
    int getNy() const{return Ny;}
    double getHx() const {return hx;}
    double getHy() const {return hy;}
    double getH() const {return std::min(hx, hy);}
    int getFringe() const {return fringe;}
    void setInnerBdry(Patch2D* patch, int position, int fringe, 
        int overlap, int unknowns, int stage);
 
};

#endif 