/* Tests for solver routines. */

#ifndef TESTINGSOLVER_H
#define TESTINGSOLVER_H

#include "Mesh2D.h"
#include "FC_2D.h"
#include "IC.h"
#include <iostream>
#include "BC_Euler_2D.h"
#include "SDNN.h"
#include "PDE.h"


void Testing1DEulerSDNN()
{
   // Create a PDE
    double T = 0.2;
    double dt = 0.000005;
    // std::string problem = "smooth_LA";
    std::string problem = "Euler1D_Sod";
    // BC bc{problem, 3};
    BC_Euler_Sod_NN bc;
    IC ic{problem, 3};
    double alpha = 1.0;
    double gamma = 7.0/5.0;

    double w1 = 2.0;
    double w2 = 1.0;
    double w3 = 0.0;
    double w4 = 0.0;
    double discard_noise = 0.01;
    
    // ANN ann{alpha};  
    SDNN sdnn{discard_noise, w1, w2, w3, w4, alpha};
    Euler1D_SDNN<BC_Euler_Sod_NN> pde {ic, bc, T, gamma, sdnn};      

    int intrb = 6;
    int npatches = 4;
    int overlap = 12;

    int N = 101;


   // Create a time-stepping scheme
    SSPRK_4 TS{3, N};
    int stages = TS.getStages();
    
    // Create a mesh
    int unknowns = pde.getPhysUnknowns();
    Mesh1DUniform mesh{0, 1.0, npatches, N, overlap, intrb, unknowns*stages + 3, 
        1, 1};
    ic(&mesh);
    //Print_Mesh1D(mesh);   

    int C = 27;
    int d = 5;

    int stage = 0;
    double alpha0 = 10.0;
    double p_0 = 2.0;
    double p = 14.0;
    double delta = 0.1;

    std::vector<FC_1D * > diff_schemes;
    std::vector<FC_1D * > filters;
    initialize1DFCSDNN(&diff_schemes, &filters, mesh, npatches, "NN", N, d, C,
        delta, alpha0, p_0, p);


    // Create a solver
    Solver<Mesh1DUniform, Euler1D_SDNN<BC_Euler_Sod_NN>, SSPRK_4, FC_1D*,
        FC_1D* > slv{mesh, pde, TS, diff_schemes, filters};


    // Run the solver
    double CFL = 2.0;
    bool visc = true;
    bool adaptive = true;  
    slv.solve_sdnn(dt, CFL, adaptive, visc);  

    // Print the solution
    Mesh1DUniform mesh1 = slv.getMesh();

    std::cout << "Solution" << std::endl;
    Print_Mesh1D(mesh1);

    std::string result_file = "result.txt";
    Print_Mesh1D(mesh1, unknowns, intrb, result_file);

    freeFCSDNN(&diff_schemes, &filters);

    std::cout << "Memory was freed" << std::endl;
}



void Testing2DFCSDNN()
{

    std::cout << "Testing 2D FC-SDNN" << std::endl;

    // Setting up a 2D mesh
    int npatches_x = 1;
    int npatches_y = 1;
    int patchsize_x = 101;
    int patchsize_y = 101;
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

    // Setting up an FC_2D spatial differentiation scheme
    int dx = 2;
    int dy = 2;
    int Cx = 27;
    int Cy = 27;

    double alpha0 = 10.0;
    double p_0 = 2.0;
    double p = 14.0;    
    double delta = 0.1;
    double hx = 1.0/(double(patchsize_x) - 1.0);
    double hy = 1.0/(double(patchsize_y) - 1.0);

    // Setting up a time differentiation scheme
    SSPRK_4 TS{unknowns, patchsize_x * patchsize_y};

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

    
    // Creating spatial diff schemes and filters
    std::vector<FC_2D* > diff_schemes;
    std::vector<FC_2D* > filters;

    initialize2DFCSDNN(&diff_schemes, &filters, mesh, npatches_x, npatches_y,
        "N", "N", "N", "N", patchsize_x, patchsize_y, dx, dy, Cx, Cy, delta,
        alpha0, p_0, p);

    // Create a solver
    Solver<Mesh2DUniform, Euler2D_SDNN<BC_Euler2D_Riemann4_NNNN>, SSPRK_4,
        FC_2D*, FC_2D* > slv{mesh, pde, TS, diff_schemes, filters};

    // Run the solver
    double CFL = 2.0;
    bool visc = true;
    bool adaptive = true;  
    double dt = 0.001;
    slv.solve2DSDNN(dt, CFL, adaptive, visc);  

    // Printing solution
    std::cout << "Solution" << std::endl;
    // Print_Mesh2DUniform(mesh, unknowns);

    // Free the memory of the spat-diff and filters
    freeFCSDNN(&diff_schemes, &filters);
    std::cout << "Memory was freed" << std::endl;
}


#endif