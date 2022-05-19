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
    std::string problem = "Burgers_Square";
    IC ic{problem, 1};
    ic(&mesh);
    std::cout << "Initial condition" << std::endl;
    Print_Mesh2DUniform(mesh);

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
    std::vector<double *> bc_d;
    std::vector<double *> bc_u;
    std::vector<double *> bc_r;
    std::vector<double *> bc_l;
    bc_d.push_back(new double[patchsize_x]);
    bc_u.push_back(new double[patchsize_x]);
    bc_l.push_back(new double[patchsize_y]);
    bc_r.push_back(new double[patchsize_y]);
    // [patchsize_x];
    // double bc_u[patchsize_x];
    // double bc_l[patchsize_y];
    // double bc_r[patchsize_y];
    memset(bc_d[0], 0.0, patchsize_x);
    memset(bc_u[0], 0.0, patchsize_x);
    memset(bc_l[0], 0.0, patchsize_y);
    memset(bc_r[0], 0.0, patchsize_y);

    FC_2D fc_2d{"DD", "DD", patchsize_x, dx, Cx, hx, patchsize_y, dy, Cy, hy,
        alpha0, p};

    // Differentiating the data on the patch
    auto patches = mesh.getPatches();
    // // fc_2d.diff_y(patches[0]->getFlowPtr()->getField(0), 
    // //     patches[0]->getFlowPtr()->getField(0), hy, bc_d[0], bc_u[0]);
    // fc_2d.diff_x(patches[0]->getFlowPtr()->getField(0), 
    //     patches[0]->getFlowPtr()->getField(0), hx, bc_l[0], bc_r[0]);

    // testing filtering, der and 2nd der

    std::vector<std::complex<double> * > fft_data_x;
    std::vector<std::complex<double> * > fft_data_y;
    std::vector<std::vector<int> > fft_locs;
    std::vector<int> loc;
    for(int j = 0; j < unknowns; j++)
    {
        fft_data_x.push_back(
            new std::complex<double> [mesh.getPatches()[0]->getNnodes() + Cx]);
        fft_data_y.push_back(
            new std::complex<double> [mesh.getPatches()[0]->getNnodes() + Cy]);
        loc.push_back(j);
    }
    fft_locs.push_back(loc);
    loc.clear();
    std::vector<int> filt_unknowns;
    for(int i = 0; i < unknowns; i++)
    {
        filt_unknowns.push_back(i);
    }

    fc_2d.filter_y(patches[0]->getFlowPtr(), filt_unknowns, &fft_data_y,
        fft_locs[0], hy, bc_d, bc_u);


    // Setting the values of the nodes and printing
    patches[0]->VectorFieldToNodes();
    mesh.setPatches(patches);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "After y-filtering" << std::endl;
    Print_Mesh2DUniform(mesh);

    fc_2d.filter_x(patches[0]->getFlowPtr(), filt_unknowns, &fft_data_x,
        fft_locs[0], hx, bc_l, bc_r);

    // Setting the values of the nodes and printing
    patches[0]->VectorFieldToNodes();
    mesh.setPatches(patches);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "After y and x filtering" << std::endl;
    Print_Mesh2DUniform(mesh);
        
        

}

#endif