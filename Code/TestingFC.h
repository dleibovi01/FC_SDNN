/* Tests for FC and FC spatial differentiation schemes routines. */

#ifndef TESTINGFC_H
#define TESTINGFC_H

#include "Mesh2D.h"
#include "FC_2D.h"
#include "IC.h"
#include <iostream>

void Testing2DFC()
{
    std::cout << "Testing 2D FC differentiation" << std::endl;

    // Setting up a 2D mesh
    int npatches_x = 1;
    int npatches_y = 1;
    int patchsize_x = 20;
    int patchsize_y = 20;
    int overlap = 4;
    int fringe = 2;
    int unknowns = 1;
    Mesh2DUniform mesh{0.0, 1.0, 0.0, 1.0, npatches_x, npatches_y, patchsize_x, 
        patchsize_y, overlap, fringe, unknowns, true, true, true, true}; 

    // Setting up an initial condition
    std::string problem = "2D_test";
    IC ic{problem, 1};
    ic(&mesh);
    std::cout << "Initial condition" << std::endl;
    Print_Mesh2DUniform(mesh);

    // Setting up an FC_2D spatial differentiation scheme
    int dx = 5;
    int dy = 5;
    int Cx = 27;
    int Cy = 27;
    double delta = 0.1;
    double hx = 1.0/(double(patchsize_x) - 1.0);
    double hy = 1.0/(double(patchsize_y) - 1.0);
    double bc_d[patchsize_x];
    double bc_u[patchsize_x];
    double bc_l[patchsize_y];
    double bc_r[patchsize_y];
    memset(bc_d, 0.0, patchsize_x);
    memset(bc_u, 0.0, patchsize_x);
    memset(bc_l, 0.0, patchsize_y);
    memset(bc_r, 0.0, patchsize_y);

    FC_2D fc_2d{"DD", "DD", patchsize_x, dx, Cx, hx, patchsize_y, dy, Cy, hy,
        delta};

    // Differentiating the data on the patch
    auto patches = mesh.getPatches();
    // fc_2d.diff_y(patches[0]->getFlowPtr()->getField(0), 
    //     patches[0]->getFlowPtr()->getField(0), hy, bc_d, bc_u);
    fc_2d.diff_x(patches[0]->getFlowPtr()->getField(0), 
        patches[0]->getFlowPtr()->getField(0), hx, bc_l, bc_r);

    // Setting the values of the nodes and printing
    patches[0]->VectorFieldToNodes();
    mesh.setPatches(patches);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "After y-differentiation" << std::endl;
    Print_Mesh2DUniform(mesh);
        

}

#endif