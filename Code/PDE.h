/* Partial differential equation object */

#ifndef PDE_H
#define PDE_H

#include "Mesh.h"
#include "IC.h"
#include "BC.h"
#include "VectorField1D.h"
#include <string>
#include <fstream>
#include <cmath>
#include "mkl_operations.h"
#include "SDNN.h"
#include "VectorOperations.h"
#include "FC.h"
#include "printing.h"
#include <array>


template <typename VectorField>
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
    
    PDE(const PDE<VectorField> &pde) : ic{pde.ic}, bc{pde.bc}, T{pde.T}, 
        phys_unknowns{pde.phys_unknowns} {}
 
    virtual ~PDE() {};

    IC getIC() const {return ic;}
    BC getBC() const {return bc;}
    double getT() const {return T;}
    int getPhysUnknowns() const {return phys_unknowns;}


};


class LA1D : public PDE<VectorField1D>{

double a;

public:

    LA1D(const IC &_ic, const BC &_bc, double &_T, double _a) : 
        PDE<VectorField1D>{_ic, _bc, _T, 1}, a{_a} {};

    VectorField1D Prim_to_cons(const VectorField1D &v) {return v;}
    VectorField1D Cons_to_prim(const VectorField1D &v) {return v;}
    VectorField1D Prim_to_flux(const VectorField1D &v) {return v;}
    VectorField1D Flux_to_prim(const VectorField1D &v) {return v;}
    VectorField1D Cons_to_flux(const VectorField1D &v) {return v;}

};

template<typename VectorField>
class SDNN_flux : public PDE<VectorField>{

protected: 

ANN ann;

public:


    SDNN_flux(const IC &_ic, const BC &_bc, double &_T, int _pu, 
        const ANN &_ann) : PDE<VectorField>{_ic, _bc, _T, _pu}, ann{_ann} {}

    // Need a copy constructor and copy-assignment
    SDNN_flux(const SDNN_flux<VectorField> & flux) : PDE<VectorField>(flux), 
        ann{flux.ann} {}    

    const ANN & getANN() const {return ann;}

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
            // std::cout <<"mu_max = " << mu_max << std::endl;
            // std::cout <<"MWSB_max = " << MWSB_max << std::endl;
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


class LA1D_SDNN : public SDNN_flux<VectorField1D>{

public:


    LA1D_SDNN(const IC &_ic, const BC &_bc, double &_T, double _a, 
        const ANN &ann) : SDNN_flux<VectorField1D>{_ic, _bc, _T, 1, ann}, 
        a{_a} {}

    LA1D_SDNN(const LA1D_SDNN & flux) : SDNN_flux<VectorField1D>(flux), 
        a{flux.a} {}


    VectorField1D Prim_to_cons(const VectorField1D &v, int stages, int stage) {return v;}
    VectorField1D Cons_to_prim(const VectorField1D &v, int stages, int stage) {return v;}

    template<typename Sp_diff>
    void Cons_to_flux(const VectorField1D &v, VectorField1D* flux, 
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

    void getMWSB(const VectorField1D &v, double * MWSB) const
    {
        int N = v.getLength();
        for(int i = 0; i < N; i++)
        {
            MWSB[i] = a;
        }
    }

    void getProxy(const VectorField1D &v, double* proxy) const
    {
        int N = v.getLength();
        std::copy(v.getField(0), v.getField(0) + N, proxy);    
    }


private:

    double a;

};



class Euler1D_SDNN : public SDNN_flux<VectorField1D>{

public:


    Euler1D_SDNN(const IC &_ic, const BC &_bc, double &_T, double _gamma, 
        const ANN &ann) : SDNN_flux<VectorField1D>{_ic, _bc, _T, 3, ann}, 
        gamma{_gamma} {}

    Euler1D_SDNN(const Euler1D_SDNN & flux) : SDNN_flux<VectorField1D>(flux), 
        gamma{flux.gamma} {}


    VectorField1D Prim_to_cons(const VectorField1D &v, int stages, int stage) {return v;}
    VectorField1D Cons_to_prim(const VectorField1D &v, int stages, int stage) {return v;}

    template<typename Sp_diff>
    void Cons_to_flux(const VectorField1D &v, VectorField1D* flux, 
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


    template<typename Sp_diff, std::size_t fourPts>
    void Cons_to_der_flux(const VectorField1D &v, VectorField1D* flux, 
        const Sp_diff &sp, int stages, int stage, 
        const std::vector<std::array<std::complex<double>, fourPts> > & ffts_old, 
        std::vector<std::array<std::complex<double> , fourPts> > *ffts_flux) const
    {
        const int N = v.getLength();
        const int C = sp.getC();
        int d = sp.getD();
        std::complex<double> * der_coeffs = sp.getDerCoeffs();
        double fourPts_dbl = sp.getFourPts_dbl();
        DFTI_DESCRIPTOR_HANDLE desc_handle = sp.getDescHandle();

        double data1[N];
        double data2[N];
        double data2_temp[N];
        double data3[N];
        double vel[N];
        double kin[N];
        double data_temp[N];
 

        double k1x[N];
        double k2x[N];
        double k3x[N];       

        std::complex<double> fft_temp[N+C]; 

        vdDiv(N, v.getField(stage*phys_unknowns + 1), 
            v.getField(stage*phys_unknowns), vel);
        vdMul(N, vel, v.getField(stage*phys_unknowns + 1), kin);

        // Differentiating
        FC_Der(k1x, ffts_old[0].data(), der_coeffs, N, C, desc_handle);
        FC_Der(k2x, ffts_old[1].data(), der_coeffs, N, C, desc_handle);
        FC_Der(k3x, ffts_old[2].data(), der_coeffs, N, C, desc_handle);

        // 1st element
        vdMul(N, v.getField(stages*phys_unknowns + 1), k1x, data1);
        cblas_daxpy(N, -1.0, v.getField(stage*phys_unknowns + 1), 1, data1, 1);
        cblas_dscal(N, -1.0, data1, 1);        
        // Fcont_Gram_Blend(data1, ffts_flux->at(0), N, d, C, fourPts_dbl, sp.getAQ(),
        //     sp.getFAQF(), desc_handle);
        Fcont_Gram_Blend(data1, fft_temp, N, d, C, fourPts_dbl, sp.getAQ(),
            sp.getFAQF(), desc_handle);
        std::copy(fft_temp, fft_temp + N + C, ffts_flux->at(0).begin());

        // 2nd element
        vdMul(N, v.getField(stages*phys_unknowns + 1), k2x, data2);
        std::copy(kin, kin + N, data_temp);
        cblas_dscal(N, 1.5 - 0.5*gamma, data_temp, 1);
        cblas_daxpy(N, gamma - 1.0, v.getField(stage*phys_unknowns + 2), 1, 
            data_temp, 1);
        cblas_daxpy(N, -1.0, data_temp, 1, data2, 1);
        cblas_dscal(N, -1.0, data2, 1);
        // Fcont_Gram_Blend(data2, ffts_flux->at(1), N, d, C, fourPts_dbl, sp.getAQ(),
        //     sp.getFAQF(), desc_handle);
        Fcont_Gram_Blend(data1, fft_temp, N, d, C, fourPts_dbl, sp.getAQ(),
            sp.getFAQF(), desc_handle);
        std::copy(fft_temp, fft_temp + N + C, ffts_flux->at(1).begin());

        // 3rd element
        vdMul(N, v.getField(stages*phys_unknowns + 1), k3x, data3);
        std::copy(kin, kin + N, data_temp);
        cblas_dscal(N, -0.5*(gamma - 1.0), data_temp, 1);
        cblas_daxpy(N, gamma, v.getField(stage*phys_unknowns + 2), 1, 
            data_temp, 1);
        vdMul(N, vel, data_temp, data_temp);
        cblas_daxpy(N, -1.0, data_temp, 1, data3, 1);
        cblas_dscal(N, -1.0, data3, 1);
        // Fcont_Gram_Blend(data3, ffts_flux->at(2), N, d, C, fourPts_dbl, sp.getAQ(),
        //     sp.getFAQF(), desc_handle);
        Fcont_Gram_Blend(data1, fft_temp, N, d, C, fourPts_dbl, sp.getAQ(),
            sp.getFAQF(), desc_handle);
        std::copy(fft_temp, fft_temp + N + C, ffts_flux->at(2).begin());

        // Differentiating
        flux->setField(0, N, data1);
        flux->setField(1, N, data2);
        flux->setField(2, N, data3);
    }


    template<typename Sp_diff>
    void Cons_to_der_flux(const VectorField1D &v, VectorField1D* flux, 
        const Sp_diff &sp, int stages, int stage, const double * mux) const
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

        double k1x[N];
        double k2x[N];
        double k3x[N];

        double k1xx[N];
        double k2xx[N];
        double k3xx[N];   

        // Differentiating
        sp.diff(v.getField(stage*phys_unknowns), k1x, k1xx);
        sp.diff(v.getField(stage*phys_unknowns + 1), k2x, k2xx);
        sp.diff(v.getField(stage*phys_unknowns + 2), k3x, k3xx);

        // Pre-computing some useful fields
        // vel = k2/k1
        vdDiv(N, v.getField(stage*phys_unknowns + 1), 
            v.getField(stage*phys_unknowns), vel);
        // vel2 = (k2/k1)^2
        // vdMul(N, vel, vel, vel2);
        VectorMul(N, vel, vel, vel2);
        // e = k3/k1
        vdDiv(N, v.getField(stage*phys_unknowns + 2), 
            v.getField(stage*phys_unknowns), e);

        // 1st element
        // vdMul(N, mux, k1x, data1);
        VectorMul(N, mux, k1x, data1);
        // vdMul(N, v.getField(stages*phys_unknowns + 1), k1xx, data1_temp);
        VectorMul(N, v.getField(stages*phys_unknowns + 1), k1xx, data1_temp);
        // vdAdd(N, data1_temp, data1, data1);
        VectorAdd(N, data1_temp, data1, data1);
        cblas_dscal(N, -1.0, data1, 1);
        vdAdd(N, k2x, data1, data1);

        // 2nd element
        // get the viscous term
        // vdMul(N, mux, k2x, data2);
        VectorMul(N, mux, k2x, data2);
        // vdMul(N, v.getField(stages*phys_unknowns + 1), k2xx, data2_temp1);
        VectorMul(N, v.getField(stages*phys_unknowns + 1), k2xx, data2_temp1);
        // vdAdd(N, data2_temp1, data2, data2);
        VectorAdd(N, data2_temp1, data2, data2);
        cblas_dscal(N, -1.0, data2, 1);
        // get the flux term
        // vdMul(N, vel, k2x, data2_temp1);
        VectorMul(N, vel, k2x, data2_temp1);
        cblas_daxpy(N, 3.0 - gamma, data2_temp1, 1, data2, 1);
        // vdMul(N, vel2, k1x, data2_temp1);
        VectorMul(N, vel2, k1x, data2_temp1);
        cblas_daxpy(N, -0.5*(3.0 - gamma), data2_temp1, 1, data2, 1);
        cblas_daxpy(N, gamma - 1.0, k3x, 1, data2, 1);

        // 3rd element
         // get the viscous term
        // vdMul(N, mux, k3x, data3);
        VectorMul(N, mux, k3x, data3);
        // vdMul(N, v.getField(stages*phys_unknowns + 1), k3xx, data3_temp1);
        VectorMul(N, v.getField(stages*phys_unknowns + 1), k3xx, data3_temp1);
        // vdAdd(N, data3_temp1, data3, data3);
        VectorAdd(N, data3_temp1, data3, data3);
        cblas_dscal(N, -1.0, data3, 1);
        // get the flux term
        // 1st part (gamma*k3*k2/k1)_x
        // vdMul(N, vel, k3x, data3_temp1);
        VectorMul(N, vel, k3x, data3_temp1);
        // vdMul(N, e, k2x, data3_temp2);
        VectorMul(N, e, k2x, data3_temp2);
        // vdAdd(N, data3_temp1, data3_temp2, data3_temp1);
        VectorAdd(N, data3_temp1, data3_temp2, data3_temp1);
        // vdMul(N, e, vel, data3_temp2);
        VectorMul(N, e, vel, data3_temp2);
        // vdMul(N, data3_temp2, k1x, data3_temp2);
        VectorMul(N, data3_temp2, k1x, data3_temp2);
        cblas_daxpy(N, -1.0, data3_temp2, 1, data3_temp1, 1);
        cblas_daxpy(N, gamma, data3_temp1, 1, data3, 1);
        // 2nd part -0.5*(gamma - 1)*k2^3/k1^2
        // vdMul(N, vel2, k2x, data3_temp1);
        VectorMul(N, vel2, k2x, data3_temp1);
        cblas_dscal(N, 3.0, data3_temp1, 1);
        // vdMul(N, vel2, vel, data3_temp2);
        VectorMul(N, vel2, vel, data3_temp2);
        // vdMul(N, data3_temp2, k1x, data3_temp2);
        VectorMul(N, data3_temp2, k1x, data3_temp2);
        cblas_daxpy(N, -2.0, data3_temp2, 1, data3_temp1, 1);
        cblas_daxpy(N, -0.5*(gamma - 1.0), data3_temp1, 1, data3, 1);


        // Setting the full flux
        flux->setField(0, N, data1);
        flux->setField(1, N, data2);
        flux->setField(2, N, data3);
    }

    void getMWSB(const VectorField1D &v, double * MWSB) const
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

    void getProxy(const VectorField1D &v, double* proxy) const
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

#endif
