/* Patch */

#include "Patch.h"

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