/* Boundary conditions for 2D Euler equations */


#ifndef BC_Euler_2D_H
#define BC_Euler_2D_H

#include "BC.h"
#include <vector>


class BC_Euler2D_Riemann4_NNNN : public BC
{

    std::vector<double> bc_l;
    std::vector<double> bc_r;
    std::vector<double> bc_d;
    std::vector<double> bc_u;
    std::vector<int> enforceable_bdries;

public:

    BC_Euler2D_Riemann4_NNNN() : BC{4}
    { 

        // bc_d.push_back(new double[Nx]);
        // bc_u.push_back(new double[Nx]);
        // bc_l.push_back(new double[Ny]);
        // bc_r.push_back(new double[Ny]);
        // memset(bc_d[0], 0.0, Nx);
        // memset(bc_d[1], 0.0, Nx);
        // memset(bc_d[2], 0.0, Nx);
        // memset(bc_d[3], 0.0, Nx);
        // memset(bc_u[0], 0.0, Nx);
        // memset(bc_u[1], 0.0, Nx);
        // memset(bc_u[2], 0.0, Nx);
        // memset(bc_u[3], 0.0, Nx);
        // memset(bc_l[0], 0.0, Ny);
        // memset(bc_l[1], 0.0, Ny);
        // memset(bc_l[2], 0.0, Ny);
        // memset(bc_l[3], 0.0, Ny);
        // memset(bc_r[0], 0.0, Ny);
        // memset(bc_r[1], 0.0, Ny);
        // memset(bc_r[2], 0.0, Ny);
        // memset(bc_r[3], 0.0, Ny);

        bc_d.push_back(0.0);
        bc_d.push_back(0.0);
        bc_d.push_back(0.0);
        bc_d.push_back(0.0);
        bc_u.push_back(0.0);
        bc_u.push_back(0.0);
        bc_u.push_back(0.0);
        bc_u.push_back(0.0);
        bc_l.push_back(0.0);
        bc_l.push_back(0.0);
        bc_l.push_back(0.0);
        bc_l.push_back(0.0);
        bc_r.push_back(0.0);
        bc_r.push_back(0.0);
        bc_r.push_back(0.0);
        bc_r.push_back(0.0);
    }


    std::vector<double> enforceBC(const Node* node, const double t) 
    {
        std::vector<double> bdry_values;
        std::vector<double> x = node->getPosition();
        if(x[0] == 0.0)
            return bc_l;
        else if(x[0] == 1.2)
            return bc_r;
        else if(x[1] == 0.0)
            return bc_d;
        else if(x[1] == 1.2)
            return bc_u;
        return bdry_values;
    }

    void getBC_L(std::vector<double* > * BC_l, const Patch2D * patch,
        const double t) const 
    {
        int Ny = patch->getNy();
        memset((*BC_l)[0], 0.0, Ny*sizeof(double));
        memset((*BC_l)[1], 0.0, Ny*sizeof(double));
        memset((*BC_l)[2], 0.0, Ny*sizeof(double));
        memset((*BC_l)[3], 0.0, Ny*sizeof(double));
    }

    void getBC_R(std::vector<double* > * BC_r, const Patch2D * patch,
        const double t) const 
    {
        int Ny = patch->getNy();
        memset((*BC_r)[0], 0.0, Ny*sizeof(double));
        memset((*BC_r)[1], 0.0, Ny*sizeof(double));
        memset((*BC_r)[2], 0.0, Ny*sizeof(double));
        memset((*BC_r)[3], 0.0, Ny*sizeof(double));
    }

    void getBC_D(std::vector<double* > * BC_d, const Patch2D * patch,
        const double t) const 
    {
        int Nx = patch->getNx();
        memset((*BC_d)[0], 0.0, Nx*sizeof(double));
        memset((*BC_d)[1], 0.0, Nx*sizeof(double));
        memset((*BC_d)[2], 0.0, Nx*sizeof(double));
        memset((*BC_d)[3], 0.0, Nx*sizeof(double));
    }

    void getBC_U(std::vector<double* > * BC_u, const Patch2D * patch,
        const double t) const 
    {
        int Nx = patch->getNx();
        memset((*BC_u)[0], 0.0, Nx*sizeof(double));
        memset((*BC_u)[1], 0.0, Nx*sizeof(double));
        memset((*BC_u)[2], 0.0, Nx*sizeof(double));
        memset((*BC_u)[3], 0.0, Nx*sizeof(double));
    }

    std::vector<int> getEnforceableBdries(const Node* node,
        const double t) const {return enforceable_bdries;}
};


#endif