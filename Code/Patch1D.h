/* 1D patch */

#ifndef PATCH1D_H
#define PATCH1D_H

#include "Node.h"
#include "VectorField.h"
#include <iostream>
#include "Patch.h"


class Patch1D : public Patch{

protected:

int intra_patch_nodes_l;
int intra_patch_nodes_r;
double h;



public:

    Patch1D(){};
    Patch1D(const VectorField &v1, std::vector<Node*> n, std::vector<int> p,
            int i_l, int i_r);
    Patch1D(const Patch1D &patch) : Patch(patch), h{patch.h},
        intra_patch_nodes_l{patch.intra_patch_nodes_l}, 
        intra_patch_nodes_r{patch.intra_patch_nodes_r} {} 
    Patch1D(int N, int _unknowns, double _a, double _b, int lb, int rb,
        int intrbl, int intrbr);
    Patch1D & operator=(const Patch1D &patch);
    double getH() const {return h;}
    // ~Patch1D();
    int getIntraPatchNodesL() const {return intra_patch_nodes_l;}
    int getIntraPatchNodesR() const {return intra_patch_nodes_r;}
    void setInnerBdry(Patch1D* patch_l, Patch1D* patch_r, int overlap_l,
        int overlap_r, int unknowns, int stage);
    void setInnerBdry(Patch1D* patch, int overlap, int unknowns, int stage,
        bool direction);

};

#endif 