
#ifndef TESTINGMESH_H
#define TESTINGMESH_H

#include "Mesh2D.h"
#include "printing.h"

void Test2DMesh()
{
    std::cout << "Testing 2D Mesh" << std::endl;
    int npatches_x = 3;
    int npatches_y = 3;
    int patchsize_x = 16;
    int patchsize_y = 16;
    int overlap = 4;
    int fringe = 2;
    int unknowns = 1;
    Mesh2DUniform mesh{0.0, 1.0, 0.0, 1.0, npatches_x, npatches_y, patchsize_x, 
        patchsize_y, overlap, fringe, unknowns, true, true, true, true}; 
    for(int j = 0; j < npatches_x; j++)
    {
        for(int i = 0; i < npatches_y; i++)
        {
            for(int k = 0; k < patchsize_x * patchsize_y; k++)
            {
                mesh.getPatches()[j*npatches_y + i]->getNode(k)->
                    setValue(0, double(j*npatches_y + i));
            }
        }
    }
    for(int i = 0; i < npatches_x * npatches_y; i++)
    {
        mesh.getPatches()[i]->NodesToVectorField();
    }
    mesh.setIntraPatchBC(1, 0);
     for(int i = 0; i < npatches_x * npatches_y; i++)
    {
        mesh.getPatches()[i]->VectorFieldToNodes();
    }   
    Print_Mesh2DUniform(mesh);


    
}


#endif