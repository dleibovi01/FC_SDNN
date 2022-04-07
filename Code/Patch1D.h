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
    Patch1D(const Patch1D &patch);
    Patch1D & operator=(const Patch1D &patch);
    ~Patch1D();
    int getIntraPatchNodesL() const {return intra_patch_nodes_l;}
    int getIntraPatchNodesR() const {return intra_patch_nodes_r;}

    
    template<typename Patch>
    void setInnerBdry(const Patch &patch_l, const Patch &patch_r, int overlap_l, 
        int overlap_r, int unknowns, int stage)
    {
        int Nnodes_l = patch_l->getNnodes();
        int Nnodes_r = patch_r->getNnodes();
        for(int i = 0; i < intra_patch_nodes_l; i++)
        {
            for(int j = 0; j < unknowns; j++)
            {
                v.setFieldValue(stage*unknowns + j, i, 
                    patch_l->getFlow().getFieldValue(stage*unknowns + j,
                    Nnodes_l - overlap_l + i));
            }
        }
        for(int i = 0; i < intra_patch_nodes_r; i++)
        {
            for(int j = 0; j < unknowns; j++)
            {
                v.setFieldValue(stage*unknowns + j, Nnodes_r - 1 - i,
                    patch_r->getFlow().getFieldValue(stage*unknowns + j,
                    overlap_r - 1 - i));
            }
        }        
    }

    template<typename Patch>
    void setInnerBdry(const Patch &patch, int overlap, int unknowns, int stage,
        bool direction)
    {
        int Nnodes_n = patch->getNnodes();
        if(direction)
        {
            for(int i = 0; i < intra_patch_nodes_l; i++)
            {
                for(int j = 0; j < unknowns; j++)
                {
                    v.setFieldValue(stage*unknowns + j, i, patch->getFlow().
                        getFieldValue(stage*unknowns + j, Nnodes_n - overlap +
                        i));
                }
            }
        }
        else
        {
            for(int i = 0; i < intra_patch_nodes_r; i++)
            {
                for(int j = 0; j < unknowns; j++)
                {
                    v.setFieldValue(stage*unknowns + j, Nnodes - 1 - i,
                        patch->getFlow().getFieldValue(stage*unknowns + j,
                        overlap- 1 - i));
                }
            }   
        }
    }
};

#endif 