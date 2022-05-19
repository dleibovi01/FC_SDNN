/* Partial differential equation object */

#ifndef PDE_H
#define PDE_H

#include "Mesh.h"
#include "IC.h"
#include "BC.h"
#include "VectorField.h"
#include <string.h>
#include <fstream>
#include <cmath>
#include "mkl_operations.h"
#include "SDNN.h"
#include "VectorOperations.h"
#include "FC.h"
#include "printing.h"
#include <array>
#include <vector>


template <typename BC>
class PDE{

protected:

/* Initial condition */
IC ic;
/* Boundary condition */
BC bc;
/* Final computation time */
double T;
/* Number of physical unknowns */
int phys_unknowns;

public:

    PDE(const IC &_ic, const BC &_bc, double &_T, int _pu) : ic{_ic}, bc{_bc}, 
        T{_T}, phys_unknowns{_pu} {}
    
    PDE(const PDE<BC> &pde) : ic{pde.ic}, bc{pde.bc}, T{pde.T}, 
        phys_unknowns{pde.phys_unknowns} {}
 
    virtual ~PDE() {};

    IC getIC() const {return ic;}
    BC getBC() const {return bc;}
    double getT() const {return T;}
    int getPhysUnknowns() const {return phys_unknowns;}


};


class LA1D : public PDE<BC>{

double a;

public:

    LA1D(const IC &_ic, const BC &_bc, double &_T, double _a) : 
        PDE<BC>{_ic, _bc, _T, 1}, a{_a} {};

    VectorField Prim_to_cons(const VectorField &v) {return v;}
    VectorField Cons_to_prim(const VectorField &v) {return v;}
    VectorField Prim_to_flux(const VectorField &v) {return v;}
    VectorField Flux_to_prim(const VectorField &v) {return v;}
    VectorField Cons_to_flux(const VectorField &v) {return v;}

};

template<typename BC>
class SDNN_flux : public PDE<BC>{

protected: 

SDNN sdnn;

public:

    SDNN_flux(const IC &_ic, const BC &_bc, double &_T, int _pu, 
        const SDNN &_sdnn) : PDE<BC>{_ic, _bc, _T, _pu}, sdnn{_sdnn} {}

    // Need a copy constructor and copy-assignment
    SDNN_flux(const SDNN_flux<BC> & flux) : PDE<BC>(flux), 
        sdnn{flux.sdnn} {}    

    const SDNN& getSDNN() const {return sdnn;}

    template<typename Mesh>
    double getAdaptiveTimeStep(const Mesh &mesh, int unknowns, int stages, 
        double CFL) const
    {
        auto patches = mesh.getPatches();
        auto patch = patches[0];
        auto v = patch->getFlow();
        int npatches = patches.size();
        double timestep_min;
        double timestep;
        int N_origin;
        double h;
        double MWSB_max;
        double mu_max;
        const double pi = std::acos(-1);
        for(int i = 0; i < npatches; i++)
        {
            patch = patches[i];
            v = patch->getFlow();
            MWSB_max = getMWSBMax(v);
            mu_max = getMuMax(v, unknowns, stages);
            h = patch->getH();
            timestep = CFL / (pi*(MWSB_max/h + mu_max/h/h));
            if(i == 0)
            {
                timestep_min = timestep;
            }
            else
            {
                if(timestep < timestep_min)
                {
                    timestep_min = timestep;
                }
            }
        }
        return timestep_min;
    }

virtual void getMWSB(const VectorField &v, double* MWSB) const {}



double getMWSBMax(const VectorField &v) const
{
    int N = v.getLength();
    double MWSB[N];
    getMWSB(v, MWSB);
    double M = *(std::max_element(MWSB, MWSB + N));
    return M;
}

double getMuMax(const VectorField &v, int unknowns, int stages) const
{
    int N = v.getLength();
    double* mu = v.getField(stages*unknowns + 1);
    double M =  *(std::max_element(mu, mu + N));   
    return M;
}


};

template <typename BC> 
class LA1D_SDNN : public SDNN_flux<BC>{

public:


    LA1D_SDNN(const IC &_ic, const BC &_bc, double &_T, double _a, 
        const SDNN &sdnn) : SDNN_flux<BC>{_ic, _bc, _T, 1, sdnn}, 
        a{_a} {}

    LA1D_SDNN(const LA1D_SDNN<BC> & flux) : SDNN_flux<BC>(flux), 
        a{flux.a} {}


    VectorField Prim_to_cons(const VectorField &v, int stages, int stage) {return v;}
    VectorField Cons_to_prim(const VectorField &v, int stages, int stage) {return v;}

    template<typename Sp_diff>
    void Cons_to_flux(const VectorField &v, VectorField* flux, 
        const Sp_diff &sp, int stages, int stage) const
    {
        int N = v.getLength();
        double data[N];
        sp.diff(v.getField(stage), data);
        vdMul(N, v.getField(stages + 1), data, data);
        cblas_dscal(N, -1.0, data, 1);
        cblas_daxpy(N, a, v.getField(stage), 1, data, 1);
        flux->setField(0, N, data);
    }


    template<typename Sp_diff>
    void Cons_to_der_flux(const VectorField &v, VectorField* flux, 
        const Sp_diff &sp, int stages, int stage, const double* mux,
        const double t) const
    {
        const int N = v.getLength();
        double k1x[N];
        double k1xx[N];

        std::vector<double> bc_l = bc.getBC_L(t);
        std::vector<double> bc_r = bc.getBC_R(t);

        // Differentiating
        sp.diff(v.getField(stage), k1x, k1xx, bc_l[0], bc_r[0]);

        computeDerFlux(v, flux, stages, stage, mux, k1x, k1xx);
    }

    template<typename Sp_diff>
    void Cons_to_der_flux(const VectorField &v, VectorField* flux, 
        const Sp_diff &sp, int stages, int stage, const double * mux, 
        const std::vector<std::complex<double> *>  & u_hat,
        const std::vector<int> & fft_loc) const
    {
        const int N = v.getLength();
        double k1x[N];
        double k1xx[N];

        // Differentiating
        sp.diff(u_hat[fft_loc[0]], k1x, k1xx);

        computeDerFlux(v, flux, stages, stage, mux, k1x, k1xx);
    }

    void computeDerFlux(const VectorField &v, VectorField* flux, 
        int stages, int stage, const double * mux, const double *k1x, 
        const double *k1xx) const
    {
        const int N = v.getLength();
        double data[N];
        double data2[N];
        VectorMul(N, k1x, mux, data);
        VectorMul(N, v.getField(stages + 1), k1xx, data2);
        VectorAdd(N, data2, data, data);
        cblas_dscal(N, -1.0, data, 1);
        cblas_daxpy(N, a, k1x, 1, data, 1);
        flux->setField(0, N, data);

    }


    void getMWSB(const VectorField &v, double * MWSB) const
    {
        int N = v.getLength();
        for(int i = 0; i < N; i++)
        {
            MWSB[i] = a;
        }
    }

    void getProxy(const VectorField &v, double* proxy) const
    {
        int N = v.getLength();
        std::copy(v.getField(0), v.getField(0) + N, proxy);    
    }


private:

    double a;

};

template<typename BC>
class Euler1D : public PDE<BC>{

    using PDE<BC>::phys_unknowns;

public:


    Euler1D(const IC &_ic, const BC &_bc, double &_T, double _gamma) : 
        PDE<BC>{_ic, _bc, _T, 3}, gamma{_gamma} {}

    Euler1D(const Euler1D<BC> & flux) : PDE<BC>(flux), 
        gamma{flux.gamma} {}

    void Cons_to_flux(const VectorField &v, VectorField* flux, int stages,
        int stage) const
    {
        const int N = v.getLength();

        double data1[N];
        double data2[N];
        double data3[N];
        double vel[N];
        double kin[N];
        double data_temp[N];      

        vdDiv(N, v.getField(stage*phys_unknowns + 1), 
            v.getField(stage*phys_unknowns), vel);
        vdMul(N, vel, v.getField(stage*phys_unknowns + 1), kin);

        // 1st element
        std::copy(v.getField(stage*phys_unknowns + 1), 
            v.getField(stage*phys_unknowns + 1) + N, data1);
        // cblas_dscal(N, -1.0, data1, 1);

        // 2nd element
        std::copy(kin, kin + N, data2);
        cblas_dscal(N, 1.5 - 0.5*gamma, data2, 1);
        cblas_daxpy(N, gamma - 1.0, v.getField(stage*phys_unknowns + 2), 1, 
            data2, 1);
        // cblas_dscal(N, -1.0, data2, 1);

        // 3rd element
        std::copy(kin, kin + N, data3);
        cblas_dscal(N, -0.5*(gamma - 1.0), data3, 1);
        cblas_daxpy(N, gamma, v.getField(stage*phys_unknowns + 2), 1, 
            data3, 1);
        vdMul(N, vel, data3, data3);



        // Differentiating
        flux->setField(0, N, data1);
        flux->setField(1, N, data2);
        flux->setField(2, N, data3);
    }


protected :

    double gamma;

};


template<typename BC>
class Euler1D_LF : public Euler1D<BC>{

    // using Euler1D<BC>::phys_unknowns;
    using Euler1D<BC>::gamma;
    using Euler1D<BC>::Cons_to_flux;

    public:

    Euler1D_LF(const IC &_ic, const BC &_bc, double &_T, double _gamma) : 
        Euler1D<BC>{_ic, _bc, _T, _gamma} {}

    Euler1D_LF(const Euler1D_LF<BC> & flux) : Euler1D<BC>(flux) {}

    template<typename Sp_diff>
    void Cons_to_der_flux(const VectorField &v, VectorField* flux, 
        const Sp_diff &sp, int stages, int stage, 
        std::vector<VectorField>* data) const
    {

        int phys_unknowns = 3;
        VectorField flux2 = *flux;
        const int N = v.getLength();
        double vel[N];
        double vel_sq[N];
        double a[N];
        double lambda;
        
        // // get the velocity
        vdDiv(N, v.getField(1), v.getField(0), vel);
        // Forming the speed of sound
        vdMul(N, vel, vel, vel_sq);
        vdDiv(N, v.getField(2),v.getField(0), a);
        cblas_daxpy(N, -0.5, vel_sq, 1, a, 1);
        vdSqrt(N, a, a);
        cblas_dscal(N, std::sqrt(gamma * (gamma - 1)), a, 1);
        // get lambda
        vdAbs(N, vel, vel);
        cblas_daxpy(N, 1.0, vel, 1, a, 1);
        lambda= *(std::max_element(a, a + N));

        // build the two fluxes
        std::vector<int> extractions;
        for(int i = 0; i < phys_unknowns; i++)
        {
            extractions.push_back(stage*phys_unknowns + i);
        }
        Cons_to_flux(v, flux, stages, stage);
                
        // flux2 = 0.5*(*flux - lambda*v.extract(extractions));
        linComb2(*flux, v.extract(extractions), &flux2, 0.5, -0.5*lambda);
        flux2.circshift(1);
        
        // flux->setFlow((0.5*(*flux + lambda*v.extract(extractions))).getFlow());
        linComb2(*flux, v.extract(extractions), flux, 0.5, 0.5*lambda);

        // obtain the flux derivative
        sp.diff(&flux2, flux, data);
    }
    

};


template<typename BC>
class Euler1D_SDNN : public SDNN_flux<BC>{

    using SDNN_flux<BC>::bc;
    using SDNN_flux<BC>::phys_unknowns;

public:


    Euler1D_SDNN(const IC &_ic, const BC &_bc, double &_T, double _gamma, 
        const SDNN &sdnn) : SDNN_flux<BC>{_ic, _bc, _T, 3, sdnn}, 
        gamma{_gamma} {}

    Euler1D_SDNN(const Euler1D_SDNN<BC> & flux) : SDNN_flux<BC>(flux), 
        gamma{flux.gamma} {}

    template<typename Sp_diff>
    void Cons_to_flux(const VectorField &v, VectorField* flux, 
        const Sp_diff &sp, int stages, int stage) const
    {

        // int mode = vmlSetMode(VML_EP);
        int N = v.getLength();
        double data1[N];
        double data2[N];
        double vel[N];
        double kin[N];
        double data_temp[N];
        double data3[N];

        vdDiv(N, v.getField(stage*phys_unknowns + 1), 
            v.getField(stage*phys_unknowns), vel);
        vdMul(N, vel, v.getField(stage*phys_unknowns + 1), kin);

        // Differentiating
        sp.diff(v.getField(stage*phys_unknowns), data1);
        sp.diff(v.getField(stage*phys_unknowns + 1), data2);
        sp.diff(v.getField(stage*phys_unknowns + 2), data3);   

        // 1st element
        vdMul(N, v.getField(stages*phys_unknowns + 1), data1, data1);
        cblas_daxpy(N, -1.0, v.getField(stage*phys_unknowns + 1), 1, data1, 1);
        cblas_dscal(N, -1.0, data1, 1);

        // 2nd element
        std::copy(kin, kin + N, data_temp);
        cblas_dscal(N, 1.0 - 0.5*(gamma - 1.0), data_temp, 1);
        cblas_daxpy(N, gamma - 1.0, v.getField(stage*phys_unknowns + 2), 1, 
            data_temp, 1);
        vdMul(N, v.getField(stages*phys_unknowns + 1), data2, data2);
        cblas_daxpy(N, -1.0, data_temp, 1, data2, 1);
        cblas_dscal(N, -1.0, data2, 1);

        // 3rd element
        std::copy(kin, kin + N, data_temp);
        cblas_dscal(N, -0.5*(gamma - 1.0), data_temp, 1);
        cblas_daxpy(N, gamma, v.getField(stage*phys_unknowns + 2), 1, 
            data_temp, 1);
        vdMul(N, vel, data_temp, data_temp);
        vdMul(N, v.getField(stages*phys_unknowns + 1), data3, data3);
        cblas_daxpy(N, -1.0, data_temp, 1, data3, 1);
        cblas_dscal(N, -1.0, data3, 1);

        // Differentiating
        flux->setField(0, N, data1);
        flux->setField(1, N, data2);
        flux->setField(2, N, data3);
    }


    template<typename Sp_diff>
    void Cons_to_der_flux(const VectorField &v, VectorField* flux, 
        const Sp_diff &sp, int stages, int stage, const double * mux, 
        const double h, const double t) const
    {
        const int N = v.getLength();

        double k1x[N];
        double k2x[N];
        double k3x[N];
        double k1xx[N];
        double k2xx[N];
        double k3xx[N];   

        std::vector<double> bc_l = bc.getBC_L(t);
        std::vector<double> bc_r = bc.getBC_R(t);

        // Differentiating
        sp->diff(v.getField(stage*phys_unknowns), k1x, k1xx, h, bc_l[0], bc_r[0]);
        sp->diff(v.getField(stage*phys_unknowns + 1), k2x, k2xx, h, bc_l[1],
            bc_r[1]);
        sp->diff(v.getField(stage*phys_unknowns + 2), k3x, k3xx, h, bc_l[2],
            bc_r[2]);


        computeDerFlux(v, flux, stages, stage, mux, k1x, k2x, k3x, k1xx, 
            k2xx, k3xx);
    }


    template<typename Sp_diff>
    void Cons_to_der_flux(const VectorField &v, VectorField* flux, 
        const Sp_diff &sp, int stages, int stage, const double * mux, 
        const std::vector<std::complex<double> *>  & u_hat,
        const std::vector<int> & fft_loc) const
    {
        const int N = v.getLength();

        double k1x[N];
        double k2x[N];
        double k3x[N];
        double k1xx[N];
        double k2xx[N];
        double k3xx[N];   

        // Differentiating
        sp->diff(u_hat[fft_loc[0]], k1x, k1xx);
        sp->diff(u_hat[fft_loc[1]], k2x, k2xx);
        sp->diff(u_hat[fft_loc[2]], k3x, k3xx);

        computeDerFlux(v, flux, stages, stage, mux, k1x, k2x, k3x, k1xx, 
            k2xx, k3xx);
    }


    void computeDerFlux(const VectorField &v, VectorField* flux, 
        int stages, int stage, const double * mux, const double *k1x, 
        const double *k2x, const double *k3x, const double *k1xx, 
        const double *k2xx, const double *k3xx) const
    {
        const int N = v.getLength();
        double data1[N];
        double data1_temp[N];

        double data2[N];
        double data2_temp1[N];

        double vel[N];
        double vel2[N];
        double e[N];

        double data3[N];
        double data3_temp1[N];
        double data3_temp2[N]; 

        // Pre-computing some useful fields
        // vel = k2/k1
        vdDiv(N, v.getField(stage*phys_unknowns + 1), 
            v.getField(stage*phys_unknowns), vel);
        // vel2 = (k2/k1)^2
        VectorMul(N, vel, vel, vel2);
        // e = k3/k1
        vdDiv(N, v.getField(stage*phys_unknowns + 2), 
            v.getField(stage*phys_unknowns), e);

        // 1st element
        VectorMul(N, mux, k1x, data1);
        VectorMul(N, v.getField(stages*phys_unknowns + 1), k1xx, data1_temp);
        VectorAdd(N, data1_temp, data1, data1);
        cblas_dscal(N, -1.0, data1, 1);
        VectorAdd(N, k2x, data1, data1);

        // 2nd element
        // get the viscous term
        VectorMul(N, mux, k2x, data2);
        VectorMul(N, v.getField(stages*phys_unknowns + 1), k2xx, data2_temp1);
        VectorAdd(N, data2_temp1, data2, data2);
        cblas_dscal(N, -1.0, data2, 1);
        // get the flux term
        VectorMul(N, vel, k2x, data2_temp1);
        cblas_daxpy(N, 3.0 - gamma, data2_temp1, 1, data2, 1);
        VectorMul(N, vel2, k1x, data2_temp1);
        cblas_daxpy(N, -0.5*(3.0 - gamma), data2_temp1, 1, data2, 1);
        cblas_daxpy(N, gamma - 1.0, k3x, 1, data2, 1);

        // 3rd element
         // get the viscous term
        // vdMul(N, mux, k3x, data3);
        VectorMul(N, mux, k3x, data3);
        VectorMul(N, v.getField(stages*phys_unknowns + 1), k3xx, data3_temp1);
        VectorAdd(N, data3_temp1, data3, data3);
        cblas_dscal(N, -1.0, data3, 1);
        // get the flux term
        // 1st part (gamma*k3*k2/k1)_x
        VectorMul(N, vel, k3x, data3_temp1);
        VectorMul(N, e, k2x, data3_temp2);
        VectorAdd(N, data3_temp1, data3_temp2, data3_temp1);
        VectorMul(N, e, vel, data3_temp2);
        VectorMul(N, data3_temp2, k1x, data3_temp2);
        cblas_daxpy(N, -1.0, data3_temp2, 1, data3_temp1, 1);
        cblas_daxpy(N, gamma, data3_temp1, 1, data3, 1);
        // 2nd part -0.5*(gamma - 1)*k2^3/k1^2
        VectorMul(N, vel2, k2x, data3_temp1);
        cblas_dscal(N, 3.0, data3_temp1, 1);
        VectorMul(N, vel2, vel, data3_temp2);
        VectorMul(N, data3_temp2, k1x, data3_temp2);
        cblas_daxpy(N, -2.0, data3_temp2, 1, data3_temp1, 1);
        cblas_daxpy(N, -0.5*(gamma - 1.0), data3_temp1, 1, data3, 1);


        // Setting the full flux
        flux->setField(0, N, data1);
        flux->setField(1, N, data2);
        flux->setField(2, N, data3);
    }



    void getMWSB(const VectorField &v, double * MWSB) const
    {
        int N = v.getLength();
        double vel[N];
        double vel_sq[N];

        // Obtaining v
        vdDiv(N, v.getField(1), v.getField(0), vel);
        // Obtaining v^2
        vdMul(N, vel, vel, vel_sq);
        // Forming the speed of sound
        vdDiv(N, v.getField(2), v.getField(0), MWSB);
        cblas_daxpy(N, -0.5, vel_sq, 1, MWSB, 1);
        vdSqrt(N, MWSB, MWSB);
        cblas_dscal(N, std::sqrt(gamma * (gamma - 1)), MWSB, 1);
        // Getting the MWSB
        vdAbs(N, vel, vel);
        cblas_daxpy(N, 1.0, vel, 1, MWSB, 1);
    }

    void getProxy(const VectorField &v, double* proxy) const
    {
        int N = v.getLength();
        double vel[N];
        double vel_sq[N];
        // Obtaining v
        vdDiv(N, v.getField(1), v.getField(0), vel);
        // Obtaining v^2
        vdMul(N, vel, vel, vel_sq);
        // Forming the speed of sound
        vdDiv(N, v.getField(2), v.getField(0), proxy);
        cblas_daxpy(N, -0.5, vel_sq, 1, proxy, 1);
        vdSqrt(N, proxy, proxy);
        cblas_dscal(N, std::sqrt(gamma * (gamma - 1)), proxy, 1); 
        // Obtaining Ma
        vdAbs(N, vel, vel);
        vdDiv(N, vel, proxy, proxy);
    }


private:

    double gamma;

};


template<typename BC>
class Euler2D_SDNN : public SDNN_flux<BC>{

    using SDNN_flux<BC>::bc;
    using SDNN_flux<BC>::phys_unknowns;

public:


    Euler2D_SDNN(const IC &_ic, const BC &_bc, double &_T, double _gamma, 
        const SDNN &sdnn) : SDNN_flux<BC>{_ic, _bc, _T, 4, sdnn}, 
        gamma{_gamma} {}

    Euler2D_SDNN(const Euler2D_SDNN<BC> & flux) : SDNN_flux<BC>(flux), 
        gamma{flux.gamma} {}


    template<typename Sp_diff>
    void Cons_to_der_flux(const VectorField &v, VectorField* flux, 
        const Sp_diff &sp, int stages, int stage, const double * mux, 
        const double *muy, const double hx, const double hy, const double t,
        const std::vector<double *> &bc_l, const std::vector<double *> &bc_r,
        const std::vector<double *> &bc_d, const std::vector<double *> &bc_u)
        const
    {
        const int N = v.getLength();

        double k1x[N];
        double k2x[N];
        double k3x[N];
        double k4x[N];
        double k1xx[N];
        double k2xx[N];
        double k3xx[N]; 
        double k4xx[N];  

        double k1y[N];
        double k2y[N];
        double k3y[N];
        double k4y[N];
        double k1yy[N];
        double k2yy[N];
        double k3yy[N]; 
        double k4yy[N];  
        
        // std::vector<double* > bc_l = bc.getBC_L(t);
        // std::vector<double* > bc_r = bc.getBC_R(t);
        // std::vector<double* > bc_d = bc.getBC_D(t);
        // std::vector<double* > bc_u = bc.getBC_U(t);

        // Differentiating
        sp->diff_x(v.getField(stage*phys_unknowns), k1x, k1xx, hx, bc_l[0],
            bc_r[0]);
        sp->diff_x(v.getField(stage*phys_unknowns + 1), k2x, k2xx, hx, bc_l[1],
            bc_r[1]);
        sp->diff_x(v.getField(stage*phys_unknowns + 2), k3x, k3xx, hx, bc_l[2],
            bc_r[2]);
        sp->diff_x(v.getField(stage*phys_unknowns + 3), k4x, k4xx, hx, bc_l[3],
            bc_r[3]);

        sp->diff_y(v.getField(stage*phys_unknowns), k1y, k1yy, hy, bc_d[0],
            bc_u[0]);
        sp->diff_y(v.getField(stage*phys_unknowns + 1), k2y, k2yy, hy, bc_d[1],
            bc_u[1]);
        sp->diff_y(v.getField(stage*phys_unknowns + 2), k3y, k3yy, hy, bc_d[2],
            bc_u[2]);
        sp->diff_y(v.getField(stage*phys_unknowns + 3), k4y, k4yy, hy, bc_d[3],
            bc_u[3]);


        // Building the divergence of the flux

        double data1[N];
        double data1_temp[N];

        double data2[N];
        double data2_temp[N];
        double data2_temp2[N];
        double data2_temp3[N];

        double data3[N];
        double data3_temp[N];
        double data3_temp2[N]; 
        double data3_temp3[N];

        double data4[N];
        double data4_temp[N];
        double data4_temp2[N]; 
        double data4_temp3[N];
        double data4_temp4[N];


        // Pre-computing some useful fields

        double rho2[N]; // rho^2
        VectorMul(N, v.getField(stage*phys_unknowns),
            v.getField(stage*phys_unknowns), rho2);

        double rho3[N]; // rho^3
        VectorMul(N, v.getField(stage*phys_unknowns), rho2, rho3);

        //////////////////////////////////////////////////////////////////////
        // 1st element : (rho*u)_x + (rho*v)_y - mux*(rho)_x - mu*rho_xx
        //               - muy*rho_y - mu*rho_yy

        // get the viscous term
        VectorMul(N, mux, k1x, data1);
        VectorMul(N, v.getField(stages*phys_unknowns + 1), k1xx, data1_temp);
        VectorAdd(N, data1, data1_temp, data1);
        VectorMul(N, muy, k1y, data1_temp);
        VectorAdd(N, data1, data1_temp, data1);
        VectorMul(N, v.getField(stages*phys_unknowns + 1), k1yy, data1_temp);
        VectorAdd(N, data1, data1_temp, data1);
        cblas_dscal(N, -1.0, data1, 1);

        // the rest of the flux
        VectorAdd(N, k2x, data1, data1);
        VectorAdd(N, k3y, data1, data1);


        //////////////////////////////////////////////////////////////////////
        // 2nd element : (rho*u^2 + p)_x + (rho*u*v)_y - mux*(rho*u)_x
        //              - mu*(rho*u)_xx - muy*(rho*u)_y - mu*(rho*u)_yy

         // get the viscous term
        VectorMul(N, mux, k2x, data2);
        VectorMul(N, v.getField(stages*phys_unknowns + 1), k2xx, data2_temp);
        VectorAdd(N, data2, data2_temp, data2);
        VectorMul(N, muy, k2y, data2_temp);
        VectorAdd(N, data2, data2_temp, data2);
        VectorMul(N, v.getField(stages*phys_unknowns + 1), k2yy, data2_temp);
        VectorAdd(N, data2, data2_temp, data2);
        cblas_dscal(N, -1.0, data2, 1);

        // the 1st term of the flux
        VectorMul(N, k2x, v.getField(stage*phys_unknowns + 1), data2_temp);
        VectorMul(N, data2_temp, v.getField(stage*phys_unknowns), data2_temp);
        cblas_dscal(N, 2.0, data2_temp, 1);
        VectorMul(N, v.getField(stage*phys_unknowns + 1),
            v.getField(stage*phys_unknowns + 1), data2_temp2);
        VectorMul(N, data2_temp2, k1x, data2_temp2);
        cblas_daxpy(N, -1.0, data2_temp2, 1, data2_temp, 1);
        vdDiv(N, data2_temp, rho2, data2_temp);
        cblas_dscal(N, 0.5*(3.0 - gamma), data2_temp, 1);

        VectorMul(N, k3x, v.getField(stage*phys_unknowns + 2), data2_temp2);
        VectorMul(N, data2_temp2, v.getField(stage*phys_unknowns), data2_temp2);
        cblas_dscal(N, 2.0, data2_temp2, 1);
        VectorMul(N, v.getField(stage*phys_unknowns + 2),
            v.getField(stage*phys_unknowns + 2), data2_temp3);
        VectorMul(N, data2_temp3, k1x, data2_temp3);
        cblas_daxpy(N, -1.0, data2_temp3, 1, data2_temp2, 1);
        vdDiv(N, data2_temp2, rho2, data2_temp2);
        cblas_daxpy(N, -0.5*(gamma - 1.0), data2_temp2, 1, data2_temp, 1);

        cblas_daxpy(N, gamma -1.0, k4x, 1, data2_temp, 1);

        // the 2nd term of the flux
        VectorMul(N, k2y, v.getField(stage*phys_unknowns + 2), data2_temp2);
        VectorMul(N, k3y, v.getField(stage*phys_unknowns + 1), data2_temp3);
        VectorAdd(N, data2_temp2, data2_temp3, data2_temp2);
        VectorMul(N, data2_temp2, v.getField(stage*phys_unknowns), data2_temp2);

        VectorMul(N, k1y, v.getField(stage*phys_unknowns + 2), data2_temp3);
        VectorMul(N, data2_temp3, v.getField(stage*phys_unknowns + 1),
            data2_temp3);
        
        cblas_daxpy(N, -1.0, data2_temp3, 1, data2_temp2, 1);

        vdDiv(N, data2_temp2, rho2, data2_temp2);

        // summing up all terms
        VectorAdd(N, data2, data2_temp, data2);
        VectorAdd(N, data2, data2_temp2, data2);
        

        //////////////////////////////////////////////////////////////////////
        // 3rd element : (rho*v^2 + p)_y + (rho*u*v)_x - mux*(rho*v)_x
        //              - mu*(rho*v)_xx - muy*(rho*v)_y - mu*(rho*v)_yy

         // get the viscous term
        VectorMul(N, mux, k3x, data3);
        VectorMul(N, v.getField(stages*phys_unknowns + 1), k3xx, data3_temp);
        VectorAdd(N, data3, data3_temp, data3);
        VectorMul(N, muy, k3y, data3_temp);
        VectorAdd(N, data3, data3_temp, data3);
        VectorMul(N, v.getField(stages*phys_unknowns + 1), k3yy, data3_temp);
        VectorAdd(N, data3, data3_temp, data3);
        cblas_dscal(N, -1.0, data3, 1);

        // the 1st term of the flux
        VectorMul(N, k3y, v.getField(stage*phys_unknowns + 2), data3_temp);
        VectorMul(N, data3_temp, v.getField(stage*phys_unknowns), data3_temp);
        cblas_dscal(N, 2.0, data3_temp, 1);
        VectorMul(N, v.getField(stage*phys_unknowns + 2),
            v.getField(stage*phys_unknowns + 2), data3_temp2);
        VectorMul(N, data3_temp2, k1y, data3_temp2);
        cblas_daxpy(N, -1.0, data3_temp2, 1, data3_temp, 1);
        vdDiv(N, data3_temp, rho2, data3_temp);
        cblas_dscal(N, 0.5*(3.0 - gamma), data3_temp, 1);

        VectorMul(N, k2y, v.getField(stage*phys_unknowns + 1), data3_temp2);
        VectorMul(N, data3_temp2, v.getField(stage*phys_unknowns), data3_temp2);
        cblas_dscal(N, 2.0, data3_temp2, 1);
        VectorMul(N, v.getField(stage*phys_unknowns + 1),
            v.getField(stage*phys_unknowns + 1), data3_temp3);
        VectorMul(N, data3_temp3, k1y, data3_temp3);
        cblas_daxpy(N, -1.0, data3_temp3, 1, data3_temp2, 1);
        vdDiv(N, data3_temp2, rho2, data3_temp2);
        cblas_daxpy(N, -0.5*(gamma - 1.0), data3_temp2, 1, data3_temp, 1);

        cblas_daxpy(N, gamma -1.0, k4y, 1, data3_temp, 1);

        // the 2nd term of the flux
        VectorMul(N, k3x, v.getField(stage*phys_unknowns + 1), data3_temp2);
        VectorMul(N, k2x, v.getField(stage*phys_unknowns + 2), data3_temp3);
        VectorAdd(N, data3_temp2, data3_temp3, data3_temp2);
        VectorMul(N, data3_temp2, v.getField(stage*phys_unknowns), data3_temp2);

        VectorMul(N, k1x, v.getField(stage*phys_unknowns + 1), data3_temp3);
        VectorMul(N, data3_temp3, v.getField(stage*phys_unknowns + 2),
            data3_temp3);
        
        cblas_daxpy(N, -1.0, data3_temp3, 1, data3_temp2, 1);

        vdDiv(N, data3_temp2, rho2, data3_temp2);

        // summing up all terms
        VectorAdd(N, data3, data3_temp, data3);
        VectorAdd(N, data3, data3_temp2, data3);
        

        //////////////////////////////////////////////////////////////////////
        // 4th element : (u*(E + p))_x + (v*(E + p))_y - mux*(E)_x
        //              - mu*(E)_xx - muy*(E)_y - mu*(E)_yy

         // get the viscous term
        VectorMul(N, mux, k4x, data4);
        VectorMul(N, v.getField(stages*phys_unknowns + 1), k4xx, data4_temp);
        VectorAdd(N, data4, data4_temp, data4);
        VectorMul(N, muy, k4y, data4_temp);
        VectorAdd(N, data4, data4_temp, data4);
        VectorMul(N, v.getField(stages*phys_unknowns + 1), k4yy, data4_temp);
        VectorAdd(N, data4, data4_temp, data4);
        cblas_dscal(N, -1.0, data4, 1);

        // 1st term of the flux
        VectorMul(N, k2x, v.getField(stage*phys_unknowns + 3), data4_temp);
        VectorMul(N, k4x, v.getField(stage*phys_unknowns + 1), data4_temp2);
        VectorAdd(N, data4_temp, data4_temp2, data4_temp);
        VectorMul(N, data4_temp, v.getField(stage*phys_unknowns), data4_temp); // just added this
        VectorMul(N, k1x, v.getField(stage*phys_unknowns + 3), data4_temp2);
        VectorMul(N, data4_temp2, v.getField(stage*phys_unknowns + 1),
            data4_temp2);
        cblas_daxpy(N, -1.0, data4_temp2, 1, data4_temp, 1);
        vdDiv(N, data4_temp, rho2, data4_temp);
        cblas_dscal(N, gamma, data4_temp, 1);

        VectorMul(N, k2x, v.getField(stage*phys_unknowns), data4_temp2);
        cblas_dscal(N, 3.0, data4_temp2, 1);
        VectorMul(N, k1x, v.getField(stage*phys_unknowns + 1), data4_temp3);
        cblas_daxpy(N, -2.0, data4_temp3, 1, data4_temp2, 1);
        VectorMul(N, data4_temp2, v.getField(stage*phys_unknowns + 1),
            data4_temp2);
        VectorMul(N, data4_temp2, v.getField(stage*phys_unknowns + 1),
            data4_temp2);
        vdDiv(N, data4_temp2, rho3, data4_temp2);
        cblas_dscal(N, - 0.5*(gamma - 1.0), data4_temp2, 1);

        VectorMul(N, k2x, v.getField(stage*phys_unknowns + 2), data4_temp3);
        VectorMul(N, k3x, v.getField(stage*phys_unknowns + 1), data4_temp4);
        cblas_daxpy(N, 2.0, data4_temp4, 1, data4_temp3, 1);
        VectorMul(N, data4_temp3, v.getField(stage*phys_unknowns), data4_temp3);
        VectorMul(N, k1x, v.getField(stage*phys_unknowns + 2), data4_temp4);
        VectorMul(N, data4_temp4, v.getField(stage*phys_unknowns + 1),
            data4_temp4);
        cblas_daxpy(N, - 2.0, data4_temp4, 1, data4_temp3, 1);
        VectorMul(N, data4_temp3, v.getField(stage*phys_unknowns + 2),
            data4_temp3);
        vdDiv(N, data4_temp3, rho3, data4_temp3);
        cblas_dscal(N, - 0.5*(gamma - 1.0), data4_temp3, 1);

        VectorAdd(N, data4_temp, data4_temp2, data4_temp);
        VectorAdd(N, data4_temp, data4_temp3, data4_temp);

        // std::cout << "printing first term of the 4th coord of the flux" << std::endl;
        // Print_Mat(data4_temp, 10, 10);
        // std::cout << std::endl;

        // sum up with the viscosity already to free data_temp
        VectorAdd(N, data4_temp, data4, data4); 


        // 2nd term of the flux
        VectorMul(N, k3y, v.getField(stage*phys_unknowns + 3), data4_temp);
        VectorMul(N, k4y, v.getField(stage*phys_unknowns + 2), data4_temp2);
        VectorAdd(N, data4_temp, data4_temp2, data4_temp);
        VectorMul(N, data4_temp, v.getField(stage*phys_unknowns), data4_temp); // just added this
        // std::cout << "printing first bit second term of the 4th coord of the flux" << std::endl;
        // Print_Mat(data4_temp, 10, 10);
        // std::cout << std::endl;
        
        
        VectorMul(N, k1y, v.getField(stage*phys_unknowns + 3), data4_temp2);
        VectorMul(N, data4_temp2, v.getField(stage*phys_unknowns + 2),
            data4_temp2);
        // std::cout << "printing first bit second term of the 4th coord of the flux" << std::endl;
        // Print_Mat(data4_temp2, 10, 10);
        // std::cout << std::endl;
        cblas_daxpy(N, -1.0, data4_temp2, 1, data4_temp, 1);
        vdDiv(N, data4_temp, rho2, data4_temp);
        cblas_dscal(N, gamma, data4_temp, 1);
        // std::cout << "printing first bit second term of the 4th coord of the flux" << std::endl;
        // Print_Mat(data4_temp, 10, 10);
        // std::cout << std::endl;

        VectorMul(N, k3y, v.getField(stage*phys_unknowns), data4_temp2);
        cblas_dscal(N, 3.0, data4_temp2, 1);
        VectorMul(N, k1y, v.getField(stage*phys_unknowns + 2), data4_temp3);
        cblas_daxpy(N, -2.0, data4_temp3, 1, data4_temp2, 1);
        VectorMul(N, data4_temp2, v.getField(stage*phys_unknowns + 2),
            data4_temp2);
        VectorMul(N, data4_temp2, v.getField(stage*phys_unknowns + 2),
            data4_temp2);
        vdDiv(N, data4_temp2, rho3, data4_temp2);
        cblas_dscal(N, - 0.5*(gamma - 1.0), data4_temp2, 1);

        VectorMul(N, k3y, v.getField(stage*phys_unknowns + 1), data4_temp3);
        VectorMul(N, k2y, v.getField(stage*phys_unknowns + 2), data4_temp4);
        cblas_daxpy(N, 2.0, data4_temp4, 1, data4_temp3, 1);
        VectorMul(N, data4_temp3, v.getField(stage*phys_unknowns), data4_temp3);
        VectorMul(N, k1y, v.getField(stage*phys_unknowns + 1), data4_temp4);
        VectorMul(N, data4_temp4, v.getField(stage*phys_unknowns + 2),
            data4_temp4);
        cblas_daxpy(N, - 2.0, data4_temp4, 1, data4_temp3, 1);
        VectorMul(N, data4_temp3, v.getField(stage*phys_unknowns + 1),
            data4_temp3);
        vdDiv(N, data4_temp3, rho3, data4_temp3);
        cblas_dscal(N, - 0.5*(gamma - 1.0), data4_temp3, 1);

        VectorAdd(N, data4_temp, data4_temp2, data4_temp);
        VectorAdd(N, data4_temp, data4_temp3, data4_temp);

        // std::cout << "printing second term of the 4th coord of the flux" << std::endl;
        // Print_Mat(data4_temp, 10, 10);
        // std::cout << std::endl;

        // sum up with the rest
        VectorAdd(N, data4_temp, data4, data4); 


        //////////////////////////////////////////////////////////////////////
        // Setting the full flux
        flux->setField(0, N, data1);
        flux->setField(1, N, data2);
        flux->setField(2, N, data3);
        flux->setField(3, N, data4);
    }


    void getMWSB(const VectorField &v, double * MWSB) const
    {
        // std::cout << "Getting MWSB" << std::endl;
        int N = v.getLength();

        // MWSB = c + sqrt(u^2 + v^2) 
        //      = sqrt(gamma*(gamma - 1)*(E/rho - 0.5*(u^2 + v^2)))
        //          + sqrt(u^2 + v^2) 
        
        double u2[N];
        double u_norm[N];

        // Obtaining u
        vdDiv(N, v.getField(1), v.getField(0), u2);
        // Obtaining u^2
        VectorMul(N, u2, u2, u2);
        // std::cout << "printing u2" << std::endl;
        // Print_Mat(u2, 10, 10);
        // std::cout << std::endl;
        // Obtain v
        vdDiv(N, v.getField(2), v.getField(0), u_norm);
        // Obtaining v^2
        VectorMul(N, u_norm, u_norm, u_norm);
        // Obtaining ||u||
        cblas_daxpy(N, 1.0, u_norm, 1, u2, 1);
        vdSqrt(N, u2, u_norm);


        // Obtaining E/rho, store it in MWSB
        vdDiv(N, v.getField(3), v.getField(0), MWSB);
        // update MWSB by - 0.5*(u^2 + v^2)
        cblas_daxpy(N, -0.5, u2, 1, MWSB, 1);
        // update MWSB by scaling by gamma*(gamma - 1), and taking sqrt
        cblas_dscal(N, gamma * (gamma - 1), MWSB, 1);
        vdSqrt(N, MWSB, MWSB);

        cblas_daxpy(N, 1.0, u_norm, 1, MWSB, 1);

        // std::cout << "printing MWSB" << std::endl;
        // Print_Mat(MWSB, 10, 10);
        // std::cout << std::endl;
    }

    void getProxy(const VectorField &v, double* proxy) const
    {
        int N = v.getLength();
        double data1[N];
        double kin[N];
        memset(proxy, 0.0, N*sizeof(double));
        // Ma = sqrt(rho*(u^2 + v^2)/(gamma*(gamma - 1)*(E - 0.5*rho(u^2 + v^2))))

        // Obtaining (rho*u)^2
        VectorMul(N, v.getField(1), v.getField(1), data1);
        // Obtaining (rho*v)^2
        VectorMul(N, v.getField(2), v.getField(2), kin);
        // Obtaining rho*(u^2 + v^2)
        VectorAdd(N, data1, kin, kin);
        vdDiv(N, kin, v.getField(0), kin);
        // Obtaining E - 0.5*rho*(u^2 + v^2) and storing it in proxy
        cblas_daxpy(N, 1.0, v.getField(3), 1, proxy, 1);
        cblas_daxpy(N, -0.5, kin, 1, proxy, 1);
        // Multiplying proxy by gamma*(gamma - 1)
        cblas_dscal(N, gamma*(gamma - 1), proxy, 1);

        // std::cout << "printing intermediary proxy" << std::endl;
        // Print_Mat(proxy, 20, 20);
        // Obtaining Ma and storing it in proxy
        vdDiv(N, kin, proxy, proxy);
        vdSqrt(N, proxy, proxy);
    }


private:

    double gamma;

};

#endif
