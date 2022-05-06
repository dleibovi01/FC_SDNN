/* Patch 1D */

#include "Patch1D.h"

Patch1D::Patch1D(const VectorField &v1, std::vector<Node*> n, std::vector<int> p,
            int i_l, int i_r)
{
    v = v1;
    Nnodes = n.size();
    for(int i = 0; i < Nnodes; i++)
    {
        nodes.push_back(n[i]);
    }
    phys_bdry_nodes = p;
    intra_patch_nodes_l = i_l;
    intra_patch_nodes_r = i_r;
}

Patch1D::Patch1D(int N, int _unknowns, double _a, double _b, int lb, int rb,
    int intrbl, int intrbr)
{
    double a = _a;
    double b = _b;
    Nnodes = N;
    h = (b - a)/double(N - 1);
    std::vector<double> position;
    position.push_back(0.0);
    for(int i = 0; i < N; i++)
    {
        nodes.push_back(new Node(_unknowns));
        nodes[i]->setIndex(i);
        position[0] = a + i*h;
        nodes[i]->setPosition(position);
    }
    for(int i = 0; i < lb; i++)
    {
        phys_bdry_nodes.push_back(i);
    }
    for(int i = 0; i < rb; i++)
    {
        phys_bdry_nodes.push_back(N - rb + i);
    }
    intra_patch_nodes_l = intrbl;
    intra_patch_nodes_r = intrbr;
    v = VectorField{_unknowns, N};
}



Patch1D & Patch1D::operator=(const Patch1D &patch)
{
    v = patch.getFlow();
    std::vector<Node*> n = patch.getNodes();
    Nnodes = n.size();    
    nodes.clear();
    for(int i = 0; i < Nnodes; i++)
    {
        nodes.push_back(n[i]);
    }
    phys_bdry_nodes = patch.getPhysBdryNodes();
    intra_patch_nodes_l = patch.getIntraPatchNodesL();
    intra_patch_nodes_r = patch.getIntraPatchNodesR();   
    return *this; 
}


void Patch1D::setInnerBdry(Patch1D* patch_l, Patch1D* patch_r, int overlap_l,
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

void Patch1D::setInnerBdry(Patch1D* patch, int overlap, int unknowns, int stage,
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

