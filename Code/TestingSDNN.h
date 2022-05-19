/* Tests for SDNN routines. */

#ifndef TESTINGSDNN_H
#define TESTINGSDNN_H

#include "Mesh2D.h"
#include "FC_2D.h"
#include "IC.h"
#include <iostream>
#include "BC_Euler_2D.h"
#include "SDNN.h"

void Testing2DSDNN()
{
    std::cout << "Testing 2D SDNN" << std::endl;

    // Setting up a 2D mesh
    int npatches_x = 1;
    int npatches_y = 1;
    int patchsize_x = 10;
    int patchsize_y = 10;
    int overlap = 4;
    int fringe = 2;
    int unknowns = 4;
    int stages = 5;
    double x_a = 0.0;
    double x_b = 1.2;
    double y_a = 0.0;
    double y_b = 1.2;
    int mesh_unknowns = unknowns*stages + 3;
    Mesh2DUniform mesh{x_a, x_b, y_a, y_b, npatches_x, npatches_y, patchsize_x, 
        patchsize_y, overlap, fringe, mesh_unknowns, true, true, true, true}; 

    // Setting up an initial condition
    std::string problem = "Euler2D_Riemann4";
    IC ic{problem, unknowns};
    ic(&mesh);
    std::cout << "Initial condition" << std::endl;
    Print_Mesh2DUniform(mesh, unknowns);


    // Setting up an FC_2D spatial differentiation scheme
    int dx = 5;
    int dy = 5;
    int Cx = 27;
    int Cy = 27;

    double alpha0 = 10.0;
    double p_0 = 2.0;
    double p = 14.0;    
    double delta = 0.1;
    double hx = 1.0/(double(patchsize_x) - 1.0);
    double hy = 1.0/(double(patchsize_y) - 1.0);


    FC_2D fc_2d{"NN", "NN", patchsize_x, dx, Cx, hx, patchsize_y, dy, Cy, hy,
        delta};

    // Creating a PDE
    BC_Euler2D_Riemann4_NNNN bc;
    double alpha = 1.0;
    double gamma = 7.0/5.0;
    double w1 = 2.0;
    double w2 = 1.0;
    double w3 = 0.0;
    double w4 = 0.0;
    double discard_noise = 0.01;
    double T = 0.25;
    SDNN sdnn{discard_noise, w1, w2, w3, w4, alpha};
    Euler2D_SDNN<BC_Euler2D_Riemann4_NNNN> pde {ic, bc, T, gamma, sdnn};

    

    // Update tau
    sdnn.updateTau(mesh.getPatches()[0], &fc_2d, pde, unknowns, stages);

    std::cout << "Printing tau" << std::endl;
    Print_Mat(mesh.getPatches()[0]->getFlowPtr()->getField(unknowns*stages), 
        patchsize_y, patchsize_x);

    // Set viscosity normalization coefficients
    SVW_mesh svw_m{mesh};

    // Update artificial viscosity
    sdnn.updateVisc(&mesh, svw_m, pde, unknowns, stages, 0);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "Printing viscosity" << std::endl;
    Print_Mat(mesh.getPatches()[0]->getFlowPtr()->getField(unknowns*stages + 1), 
        patchsize_y, patchsize_x);

    std::cout << std::endl;
}

#endif