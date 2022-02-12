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
#include "FC.h"
#include "printing.h"

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

    double * getProxy(const VectorField1D &v, int stage) const
    {
        return v.getField(stage);    
    }


private:

    double a;

};



class Euler1D_SDNN : public SDNN_flux<VectorField1D>{

public:


    Euler1D_SDNN(const IC &_ic, const BC &_bc, double &_T, double _gamma, 
    const ANN &ann) : SDNN_flux<VectorField1D>{_ic, _bc, _T, 1, ann}, 
    gamma{_gamma} {}

    Euler1D_SDNN(const Euler1D_SDNN & flux) : SDNN_flux<VectorField1D>(flux), 
    gamma{flux.gamma} {}


    VectorField1D Prim_to_cons(const VectorField1D &v, int stages, int stage) {return v;}
    VectorField1D Cons_to_prim(const VectorField1D &v, int stages, int stage) {return v;}

    template<typename Sp_diff>
    void Cons_to_flux(const VectorField1D &v, VectorField1D* flux, 
        const Sp_diff &sp, int stages, int stage) const
    {
        int N = v.getLength();
        double data1[N];
        double data2[N];
        double vel[N];
        double kin[N];
        double data3[N];

        vdDiv(N, v.getField(stages + 1), v.getField(stages), vel);
        vdMul(N, vel, v.getField(stages + 1), kin);


        // 1st element
        std::copy(v.getField(stage + 1), v,getField(stage + 1) + N, data1);

        // 2nd element
        std::copy(kin, kin + N, data2);
        cblas_dscal(N, 1.0 - 0.5*(gamma - 1.0), data2, 1);
        cblas_daxpy(N, gamma - 1.0, v.getField(stage + 2), 1, data2, 1);

        // 3rd element
        std::copy(kin, kin + N, data3);
        cblas_dscal(N, 0.5*(gamma - 1.0), data3, 1);
        cblas_daxpy(N, gamma, v.getField(stage + 2), 1, data3, 1);
        vdMul(N, vel, data3, data3);

        // Differentiating

        sp.diff(data1, data1);
        sp.diff(data2, data2);
        sp.diff(data3, data3);

        flux->setField(0, N, data1);
        flux->setField(1, N, data2);
        flux->setField(2, N, data3);
    }

    void getMWSB(const VectorField1D &v, double * MWSB) const
    {
        int N = v.getLength();
        for(int i = 0; i < N; i++)
        {
            MWSB[i] = a;
        }
    }

    double * getProxy(const VectorField1D &v, int stage) const
    {
        return v.getField(stage);    
    }


private:

    double gamma;

};

#endif