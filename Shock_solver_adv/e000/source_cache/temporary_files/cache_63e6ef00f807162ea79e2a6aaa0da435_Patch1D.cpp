/* Patch 1D */

#include "Patch1D.h"

Patch1D::Patch1D(const VectorField1D &v1, std::vector<Node1D*> n, std::vector<int> p,
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

Patch1D::Patch1D(const Patch1D &patch)
{
    v = patch.getFlow();
    std::vector<Node1D*> n = patch.getNodes();
    Nnodes = n.size();
    for(int i = 0; i < Nnodes; i++)
    {
        nodes.push_back(n[i]);
    }
    phys_bdry_nodes = patch.getPhysBdryNodes();
    intra_patch_nodes_l = patch.getIntraPatchNodesL();
    intra_patch_nodes_r = patch.getIntraPatchNodesR();
}

Patch1D & Patch1D::operator=(const Patch1D &patch)
{
    v = patch.getFlow();
    std::vector<Node1D*> n = patch.getNodes();
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

Patch1D::~Patch1D()
{
    while(!nodes.empty())
    {
        delete nodes.back();
        Nnodes--;
        nodes.pop_back();
    }
}


void Patch1D::NodesToVectorField()
{
    int unknowns;
    for(int i = 0; i < Nnodes; i++)
    {
        unknowns = nodes[i]->getUnknowns();
        for(int j = 0; j < unknowns; j++)
        {
            v.setFieldValue(j, i, (nodes[i]->getValues())[j]);
        }
    }
}


void Patch1D::VectorFieldToNodes()
{
    int unknowns = v.getUnknowns();
    for(int i = 0; i < Nnodes; i++)
    {
        for(int j = 0; j < unknowns; j++)
        {
            (nodes[i])->setValue(j, (v.getField(j))[i]);
        }
    }
}

