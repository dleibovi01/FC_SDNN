/* 2D Mesh */

#ifndef MESH_2D_H
#define MESH_2D_H

#include "Mesh.h"


class Mesh2DUniform : public Mesh<Patch2D>
{

protected: 

    /* number of patches in the x-direction. */
    int Npatches_x;

    /* number of patches in the y-direction. */
    int Npatches_y;

    /* overlap between the patches. */
    int overlap;

public:

    /* Constructor */   
    Mesh2DUniform(double a, double b, double c, double d, int npatches_x,
        int npatches_y, int patchsize_x, int patchsize_y, int _overlap, 
        int fringe, int unknowns, bool lb, bool rb, bool db, bool ub);  

    int getOverlap() const {return overlap;}

    int getNpatches_x() const {return Npatches_x;}

    int getNpatches_y() const {return Npatches_y;}

    void setIntraPatchBC(int unknowns, int stage);
};


#endif