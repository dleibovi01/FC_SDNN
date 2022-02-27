/* 1D patch */

#ifndef PATCH1D_H
#define PATCH1D_H

#include "Node1D.h"
#include "VectorField1D.h"
#include <iostream>


class Patch1D{

protected:

VectorField1D v;
std::vector<Node1D*> nodes;
int Nnodes;
std::vector<int> phys_bdry_nodes;
int intra_patch_nodes_l;
int intra_patch_nodes_r;


public:

    Patch1D(){};
    Patch1D(const VectorField1D &v1, std::vector<Node1D*> n, std::vector<int> p,
            int i_l, int i_r);
    Patch1D(const Patch1D &patch);
    Patch1D & operator=(const Patch1D &patch);
    ~Patch1D();
    const VectorField1D & getFlow() const {return v;}
    VectorField1D & getFlowRef() {return v;}
    // VectorField1D* getVectorField() const {return &v;}
    std::vector<Node1D*> getNodes() const {return nodes;}
    Node1D* getNode(int i) const {return nodes[i];}
    int getNnodes() const {return Nnodes;}
    const std::vector<int> & getPhysBdryNodes() const {return phys_bdry_nodes;}
    int getIntraPatchNodesL() const {return intra_patch_nodes_l;}
    int getIntraPatchNodesR() const {return intra_patch_nodes_r;}

    void setField(const VectorField1D &_flow){v = _flow;}
    void setFlowValue(int i, int j, double d){v.setFieldValue(i, j, d);};
    void NodesToVectorField();
    void VectorFieldToNodes();
    
    template<typename Patch>
    void setInnerBdry(const Patch &patch_l, const Patch &patch_r, int overlap_l, 
        int overlap_r, int unknowns)
    {
        int Nnodes_l = patch_l->getNnodes();
        int Nnodes_r = patch_r->getNnodes();
        for(int i = 0; i < intra_patch_nodes_l; i++)
        {
            for(int j = 0; j < unknowns; j++)
            {
                v.setFieldValue(j, i, patch_l->getFlow().getFieldValue(
                    j, Nnodes_l - overlap_l + i));
            }
        }
        for(int i = 0; i < intra_patch_nodes_r; i++)
        {
            for(int j = 0; j < unknowns; j++)
            {
                v.setFieldValue(j, Nnodes_r - 1 - i, patch_r->getFlow().
                    getFieldValue(j, overlap_r - 1 - i));
            }
        }        
    }

    template<typename Patch>
    void setInnerBdry(const Patch &patch, int overlap, int unknowns, 
        bool direction)
    {
        int Nnodes_n = patch->getNnodes();
        if(direction)
        {
            for(int i = 0; i < intra_patch_nodes_l; i++)
            {
                for(int j = 0; j < unknowns; j++)
                {
                    v.setFieldValue(j, i, patch->getFlow().
                        getFieldValue(j, Nnodes_n - overlap + i));
                }
            }
        }
        else
        {
            for(int i = 0; i < intra_patch_nodes_r; i++)
            {
                for(int j = 0; j < unknowns; j++)
                {
                    v.setFieldValue(j, Nnodes - 1 - i, patch->getFlow().
                        getFieldValue(j, overlap- 1 - i));
                }
            }   
        }
    }
};

#endif 