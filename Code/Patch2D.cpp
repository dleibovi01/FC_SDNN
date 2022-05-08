/* 2D Patch */

#include "Patch2D.h"
#include "printing.h"
#include "VectorField.h"

Patch2D & Patch2D::operator=(const Patch2D &patch)
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
    hx = patch.getHx();
    hy = patch.getHy();
    Nx = patch.getNx();
    Ny = patch.getNy();
    inner_nodes = patch.inner_nodes;
    return *this; 
}


Patch2D::Patch2D(int _Nx, int _Ny, int _unknowns, double a, double b, 
    double c, double d, bool lb, bool rb, bool db, bool ub, int _fringe) : 
    Nx{_Nx}, Ny{_Ny}, fringe{_fringe}
{
    Nnodes = Nx * Ny;
    hx = (b - a) / double(Nx - 1);
    hy = (d - c) / double(Ny - 1);
    std::vector<double> position(2);
    for(int i = 0; i < Nx; i++)
    {
        for(int j = 0; j < Ny; j++)
        {
            nodes.push_back(new Node(_unknowns));
            nodes[Ny * i + j]->setIndex(Ny * i + j);
            position[0] = a + i*hx;
            position[1] = c + j*hy;
            nodes[Ny * i + j]->setPosition(position);
        }
    }
    if(lb)
    {
        for(int j = 0; j < Ny; j++)
        {
            phys_bdry_nodes.push_back(j);
        }
    }
    if(rb)
    {
        for(int j = 0; j < Ny; j++)
        {
            phys_bdry_nodes.push_back(Ny * (Nx - 1) + j);
        }       
    }
    if(db)
    {
         for(int i = 0; i < Nx; i++)
        {
            phys_bdry_nodes.push_back(Ny * i);
        }          
    }
    if(ub)
    {
         for(int i = 0; i < Nx; i++)
        {
            phys_bdry_nodes.push_back(Ny * (i + 1) - 1);
        }          
    }
    v = VectorField{_unknowns, Nnodes};    
}


void Patch2D::setInnerBdry(Patch2D* patch, int position, int fringe, 
        int overlap, int unknowns, int stage)
{
    //  Position: position of the patch from which the boundary information is
    //  is extracted relative to the current patch
    //                      1     2     3
    //                      4   Patch   6
    //                      7     8     9
    //
    int Nx2 = patch->getNx();
    int Nnodes2 = patch->getNnodes();

    switch (position)
    {
        case 1:
            for(int i = 0; i < fringe; i++)
            {
                for(int j = 0; j < fringe; j++)
                {
                    for(int k = 0; k < unknowns; k++)
                    {
                        v.setFieldValue(stage*unknowns + k, j*Ny + i, 
                            patch->getFlow().getFieldValue(stage*unknowns + k,
                            (Nx - overlap + j)*Ny + Ny - overlap + i));
                    }
                }
            }
            break;

        case 2:
            for(int i = 0; i < fringe; i++)
            {
                for(int j = 0; j < Nx; j++)
                {
                    for(int k = 0; k < unknowns; k++)
                    {
                        v.setFieldValue(stage*unknowns + k, j*Ny + i, 
                            patch->getFlow().getFieldValue(stage*unknowns + k,
                            j*Ny + Ny - overlap + i));
                    }
                }
            }
            break;

        case 3:
            for(int i = 0; i < fringe; i++)
            {
                for(int j = 0; j < fringe; j++)
                {
                    for(int k = 0; k < unknowns; k++)
                    { 
                        v.setFieldValue(stage*unknowns + k, (Nx - fringe + j)*Ny
                            + i, patch->getFlow().getFieldValue(stage*unknowns
                            + k, (overlap - fringe + j)*Ny + Ny - overlap + i));
                    }
                }
            }
            break;

        case 4:
            for(int i = 0; i < Ny; i++)
            {
                for(int j = 0; j < fringe; j++)
                {
                    for(int k = 0; k < unknowns; k++)
                    {
                        v.setFieldValue(stage*unknowns + k, j*Ny + i, 
                            patch->getFlow().getFieldValue(stage*unknowns + k,
                            (Nx - overlap + j)*Ny + i));
                    }
                }
            }
            break;

        case 6:
            for(int i = 0; i < Ny; i++)
            {
                for(int j = 0; j < fringe; j++)
                {
                    for(int k = 0; k < unknowns; k++)
                    {                      
                        v.setFieldValue(stage*unknowns + k, (Nx - fringe + j)*Ny
                            + i, 
                            patch->getFlow().getFieldValue(stage*unknowns + k,
                            (overlap - fringe + j)*Ny + i));
                    }
                }
            }
            break;

        case 7:
            for(int i = 0; i < fringe; i++)
            {
                for(int j = 0; j < fringe; j++)
                {
                    for(int k = 0; k < unknowns; k++)
                    { 
                        v.setFieldValue(stage*unknowns + k, j*Ny + Ny - fringe
                            + i, patch->getFlow().getFieldValue(stage*unknowns
                            + k, (Nx - overlap + j)*Ny + overlap - fringe + i));
                    }
                }
            }
            break;

        case 8:
            for(int i = 0; i < fringe; i++)
            {
                for(int j = 0; j < Nx; j++)
                {
                    for(int k = 0; k < unknowns; k++)
                    {
                        v.setFieldValue(stage*unknowns + k, j*Ny + Ny - fringe
                            + i, patch->getFlow().getFieldValue(stage*unknowns
                            + k, j*Ny + overlap - fringe + i));
                    }
                }
            }
            break;

        case 9:
            for(int i = 0; i < fringe; i++)
            {
                for(int j = 0; j < fringe; j++)
                {
                    for(int k = 0; k < unknowns; k++)
                    {                        
                        v.setFieldValue(stage*unknowns + k, (Nx - fringe + j)*Ny
                            + Ny - fringe + i, 
                            patch->getFlow().getFieldValue(stage*unknowns
                            + k, (overlap - fringe + j)*Ny 
                            + overlap - fringe + i));
                    }
                }
            }
            break;
    }
}