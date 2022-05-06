/* Patch */

#include "Patch.h"



Patch::Patch(const Patch &patch) : v{patch.v}, Nnodes{patch.Nnodes}, 
    phys_bdry_nodes{patch.phys_bdry_nodes}, inner_nodes{patch.inner_nodes}
{
    std::vector<Node*> n = patch.getNodes();
    for(int i = 0; i < Nnodes; i++)
    {
        nodes.push_back(n[i]);
    }
}


Patch::~Patch()
{
    while(!nodes.empty())
    {
        delete nodes.back();
        Nnodes--;
        nodes.pop_back();
    }
}


void Patch::NodesToVectorField()
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


void Patch::VectorFieldToNodes()
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