/* SDNN functions */

#ifndef SDNN_H
#define SDNN_H


#include "FC.h"
#include "ANN.h"
#include <algorithm>
#include <cmath>
#include <string>
#include "VectorField.h"
#include "Mesh.h"
#include "SVW.h"
#include "MVOperations.h"
#include "VectorOperations.h"

class SDNN{

ANN* ann;
double discard_noise;
double w1;
double w2;
double w3;
double w4;

public :

    SDNN(const ANN &_ann, double _discard_noise, double _w1, double _w2,
        double _w3, double _w4);

    SDNN(double _discard_noise, double _w1, double _w2, double _w3, double _w4,
        double alpha);

    SDNN(const SDNN &sdnn);

    SDNN & operator= (const SDNN &sdnn);

    ~SDNN() {delete ann;}

    template<typename Mesh, typename PDE>
    void updateVisc(Mesh *mesh, const SVW_mesh &svw_mesh, const PDE &pde,
        int phys_unknowns, int stages, int stage) const
    {
        for(int i = 0; i < mesh->getPatches().size(); i++)
            updateViscPatch(mesh->getPatches()[i], mesh, svw_mesh.getSVWs()[i],
                pde, phys_unknowns, stages, stage);
    }

    template<typename SpatDiffScheme, typename PDE>
    void updateTau (Patch1D * patch, const SpatDiffScheme &sp, const PDE &pde,
        int unknowns, int stages) const
    {
        MKL_LONG status;
        int fourPts;
        auto v = patch->getFlowPtr();
        double h = patch->getH();
        
        int N = v->getLength();
        double proxy[N];
        pde.getProxy(*v, proxy);
        int tau[N];
        bool discard[N];

        int d = sp->getD();
        int C = sp->getC();
        double* AQ = sp->getAQ();
        double* FAQF = sp->getFAQF();
        DFTI_DESCRIPTOR_HANDLE desc_handle = sp->getDescHandle();

        // Form the continuation
        fourPts = N + C;
        double proxy_ext[N+C];   
        // Fcont_shift(proxy, proxy_ext, sp->getShiftCoeffs(), N, d, C,
        //     double(fourPts), AQ, FAQF, desc_handle);
        sp->shift(proxy, proxy_ext, h, 0.0, 0.0);

        // Form the stencils 
        double stencils[s*N];
        form_stencils(N, C, s, proxy_ext, stencils);

        // Pre-process the stencils
        preprocess_stencils(N, stencils, discard, tau);

        // get classification (forward propagation)
        getRegularity(tau, discard, stencils, N, s);

        double tau_dbl[N];  

        for(int i = 0; i < N; i++)
            tau_dbl[i] = double(tau[i]);

        v->setField(unknowns*stages, N, tau_dbl);
    }


    template<typename SpatDiffScheme, typename PDE>
    void updateTau (Patch2D * patch, const SpatDiffScheme &sp, const PDE &pde,
        int unknowns, int stages) const
    {
        MKL_LONG status;
        int fourPts;
        auto v = patch->getFlowPtr();
        
        int Nx = patch->getNx();
        int Ny = patch->getNy();
        int Cx = sp->getCx();
        int Cy = sp->getCy();
        double hx = patch->getHx();
        double hy = patch->getHy();

        int dx = sp->getDx();
        int dy = sp->getDy();

        double fourPtsx = sp->getFourPtsDblx();
        double fourPtsy = sp->getFourPtsDbly();

        double * AQx = sp->getAQx();
        double * FAQFx = sp->getFAQFx();
        double * AQy = sp->getAQy();
        double * FAQFy = sp->getFAQFy();
        DFTI_DESCRIPTOR_HANDLE desc_handlex = sp->getDescHandlex();
        DFTI_DESCRIPTOR_HANDLE desc_handley = sp->getDescHandley();

        double proxy[Nx*Ny];
        int tau[2*Nx*Ny];
        bool discard[2*Nx*Ny];

        pde.getProxy(*v, proxy);

        double proxy_ext_x[Nx + Cx];
        double proxy_ext_y[Ny + Cy];

        double stencils[2*s*Nx*Ny];

        // Building pre-processed stencils in the y-direction
        for(int i = 0; i < Nx; i++)
        {
            // Forming continuation
            sp->shift1D_y(proxy + i*Ny, proxy_ext_y, hy, 0.0, 0.0);

            // Form the stencils 
            form_stencils(Ny, Cy, s, proxy_ext_y, stencils + i*s*Ny);
        }


        // Transpose proxy
        mkl_dimatcopy('C', 'T', Ny, Nx, 1.0, proxy, Ny, Nx);
        // 
        for(int i = 0; i < Ny; i++)
        {
            // Forming continuation
            sp->shift1D_x(proxy + i*Nx, proxy_ext_x, hx, 0.0, 0.0);

            // Form the stencils 
            form_stencils(Nx, Cx, s, proxy_ext_x, stencils + s*Ny*Nx + i*s*Nx);
        }

        // Pre-process the stencils
        preprocess_stencils(2*Nx*Ny, stencils, discard, tau);

        // get classification (forward propagation)
        getRegularity(tau, discard, stencils, 2*Nx*Ny, s);

        double tau_dbl_x[Nx*Ny];  
        double tau_dbl_y[Nx*Ny];  

        for(int i = 0; i < Nx*Ny; i++)
        {
            tau_dbl_y[i] = double(tau[i]);
            tau_dbl_x[i] = double(tau[Nx*Ny + i]);
        }
        
        // get the maximum between the two partial classifications
        mkl_dimatcopy('C', 'T', Nx, Ny, 1.0, tau_dbl_x, Nx, Ny);
    
        for(int i = 0; i < Nx*Ny; i++)
        {
            tau_dbl_y[i] = std::min(tau_dbl_y[i], tau_dbl_x[i]);
        }

        v->setField(unknowns*stages, Nx*Ny, tau_dbl_y);

    }

private :

    void formRegStencils(int N, int s, const double * stencils, 
        const bool * discard, double * regStencils, int * indices) const;

    double Q(double tau) const;

    void getRegularity(int * tau, bool * discard, double * stencils, int N,
        int s) const;

    std::vector<double> Visc_weights(int N, const double* tau) const;

    void form_stencils(int N, int C, int s, const double* proxy_ext,
        double* stencils) const;

    void preprocess_stencils(int N, double *stencils, bool* discard, int* tau)
        const;

    void getMaxedMWSB(Patch1D * patch, int s, const double * MWSB,
        double * MWSB_maxed) const;

    void getMaxedMWSB(Patch2D * patch, int s, const double * MWSB,
        double * MWSB_maxed) const;

    void form2DLocStencil(const double * y, double * stencil, int Nx, int Ny,
        int s, int i, int j) const;

    template<typename Patch, typename Mesh, typename SVW, typename PDE>
    void updateViscPatch(Patch *patch, Mesh *mesh, const SVW &svw,
        const PDE &pde, int phys_unknowns, int stages, int stage) const
    {
        int N = patch->getNnodes();
        auto W = svw->getPatchSVWS();
        int M = W.size();

        auto patchIds = svw->getPatchIds();
        auto patches = mesh->getPatches();
        
    
        double MWSB[N];
        std::vector<double> weighted_tau;
        std::vector<double> zeros(N, 0.0);
        double* mu = zeros.data();
        double h = patch->getH();
        
        double alpha = 1.0;
        double beta = 1.0;  
        int status;
        int N0;
        // Compute the contribution of each smooth wave function
        for(int i = 0; i < M; i++)
        {
            weighted_tau = 
                Visc_weights(patches[patchIds[i]]->getFlowRef().getLength(),
                patches[patchIds[i]]->getFlowRef().getField(phys_unknowns*stages));

            W[i]->SMatVect(weighted_tau.data(), mu);
        }
    
        // Need to form stencils and get maximums of MWSB
        int s = 7;

        pde.getMWSB(patch->getFlow(), MWSB);
        vdAbs(N, MWSB, MWSB);
        double MWSB_maxed[N];
        getMaxedMWSB(patch, s, MWSB, MWSB_maxed);

        double y[N];
        for(int i = 0; i < N; i++)
        {
            y[i] = 0.0;
        }
        // Muiltiply by the wave speed
        VectorMul(N, MWSB_maxed, mu, mu);
        // Muiltiply by h
        cblas_daxpy(N, h, mu, 1, y, 1); 
        patch->getFlowPtr()->setField(phys_unknowns*stages + 1, N, y);
    }
};


#endif