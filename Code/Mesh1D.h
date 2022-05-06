/* 1D Mesh */

#ifndef MESH_1D_H
#define MESH_1D_H

#include "Mesh.h"


class Mesh1DUniform : public Mesh<Patch1D>
{
public:

    int overlap;

    /* Constructor */   
    Mesh1DUniform(double a, double b, int n_patches, int patchsize, 
        int _overlap, int intrb, int unknowns, int l_b, int r_b);  

    int getOverlap() const {return overlap;}

    double getH() const {return (patches[0]->getH());}

    void setIntraPatchBC(int unknowns, int stage);
};


#endif