/* Time stepping scheme */

#ifndef TIMESTEPSCHEME_H
#define TIMESTEPSCHEME_H

#include "VectorField1D.h"
#include "Mesh.h"
#include "SpatDiffScheme.h"
#include <algorithm> 
#include "printing.h"
#include "PDE.h"

/*
template<typename Sp_Diff>
class TimeStepScheme{


};
*/

/*
template<typename Sp_Diff>
class Fwd_Euler : public TimeStepScheme<Sp_Diff>{

Sp_Diff sp_diff;
double dt;

public:

Fwd_Euler(const Sp_Diff & sd, double _dt, double _h) : sp_diff{sd}, dt{_dt}, 
h{_h} {}
*/


template<typename VectorField, typename Sp_Diff>
void diffLinComb(VectorField* v, 
    const std::vector<std::vector<int> > &extractions, 
    const std::vector<int> &target, const std::vector<double> & coeffs,
    const VectorField &flux, const Sp_Diff &sp)
{
    // extractions[ncoeffs, unknowns]
    // target [unknowns]
    int ncoeffs = coeffs.size();
    int N = v->getLength();
    double data[v->getLength()];
    int unknowns = v->getUnknowns();
    for(int i = 0; i < N; i++)
    {
        data[i] = 0.0;
    }
    for(int i = 0; i < extractions[0].size(); i++) 
    {
        sp.diff(flux.getField(i), data);
        cblas_dscal(N, coeffs[ncoeffs - 1], data, 1);
        for(int j = 0; j < ncoeffs - 1; j++)
        {
            cblas_daxpy(N, coeffs[j], v->getField(extractions[j][i]), 1, data, 1);
        }
        v->setField(target[i], N, data);
    }
};


template<typename VectorField, typename Sp_Diff>
void diffLinComb(VectorField* v, 
    const std::vector<std::vector<int> > &extractions, 
    const std::vector<int> &target, const std::vector<double> & coeffs,
    const std::vector<VectorField*> &fluxes, 
    const std::vector<double> &coeffs_diff, const Sp_Diff &sp)
{
    // extractions[ncoeffs, unknowns]
    // target [unknowns]
    int ncoeffs = coeffs.size();
    int N = v->getLength();
    double data[v->getLength()];
    double data_temp[v->getLength()];
    int unknowns = v->getUnknowns();
    for(int i = 0; i < N; i++)
    {
        data[i] = 0.0;
    }
    for(int i = 0; i < extractions[0].size(); i++) 
    {
        sp.diff(fluxes[0]->getField(i), data);
        cblas_dscal(N, coeffs_diff[0], data, 1);
        for(int j = 1; j < fluxes.size(); j++)
        {
            sp.diff(fluxes[j]->getField(i), data_temp);
            cblas_daxpy(N, coeffs_diff[j], data_temp, 1, data, 1);
        }    
        for(int j = 0; j < ncoeffs; j++)
        {
            cblas_daxpy(N, coeffs[j], v->getField(extractions[j][i]), 1, data, 1);
        }
        v->setField(target[i], N, data);
    }
};



class Fwd_Euler{

double dt;
int stages = 1;

public:

double getStages() {return stages;}
Fwd_Euler(double _dt) : dt{_dt} {}

template<typename Sp_Diff>
VectorField1D advance(const VectorField1D & v, const Sp_Diff &sp_diff)
{
    VectorField1D u = v - dt*sp_diff.diff(v);
    return u;
}

};


class SSPRK_4{


static constexpr double a11 = 0.391752226571890;

static constexpr double a21 = 0.444370493651235;
static constexpr double a22 = 0.555629506348765;
static constexpr double a23 = 0.368410593050371;

static constexpr double a31 = 0.620101851488403;
static constexpr double a32 = 0.379898148511597;
static constexpr double a33 = 0.251891774271694;

static constexpr double a41 = 0.178079954393132;
static constexpr double a42 = 0.821920045606868;
static constexpr double a43 = 0.544974750228521;

static constexpr double a51 = 0.517231671970585;
static constexpr double a52 = 0.096059710526147;
static constexpr double a53 = 0.063692468666290;
static constexpr double a54 = 0.386708617503269;
static constexpr double a55 = 0.226007483236906;


int stages = 5;

public: 

    int getStages() const {return stages;}


    template<typename VectorField, typename Sp_Diff, typename PDE>
    // void advance(VectorField1D* v, const Sp_Diff &sp_diff, const PDE &pde,
    void advance(VectorField* v, const Sp_Diff &sp_diff, const PDE &pde,
        const double dt)
    {
        const int unknowns = pde.getPhysUnknowns();
        const int N = v->getLength();

        std::vector<int> stage0;
        std::vector<int> stage1;
        std::vector<int> stage2;
        std::vector<int> stage3;
        std::vector<int> stage4;
        for(int i = 0; i < unknowns; i++)
        {
            stage0.push_back(i);
            stage1.push_back(i + unknowns);
            stage2.push_back(i + 2*unknowns);
            stage3.push_back(i + 3*unknowns);
            stage4.push_back(i + 4*unknowns);
        }

        //VectorField1D v = u;       
        // VectorField1D flux{unknowns, N};
        VectorField flux{unknowns, N};


        // 1st stage
        // VectorField1D v0 = v->extract(stage0);
        // pde.Cons_to_flux(*v, &flux, sp_diff, stages, 0);
        // v->setFields(v0 - a11*dt*sp_diff.diff(flux), stage1);

        // 2nd stage
        // VectorField1D v1 = v->extract(stage1);
        // pde.Cons_to_flux(*v, &flux, sp_diff, stages, 1);
        // v->setFields(a21*v0 + a22*v1 - a23*dt*sp_diff.diff(flux), stage2); 
        // Print_VectorField1D(*v);

        // 3rd stage
        // VectorField1D v2 = v->extract(stage2);
        // pde.Cons_to_flux(*v, &flux, sp_diff, stages, 2);
        // v->setFields(a31*v0 + a32*v2 - a33*dt*sp_diff.diff(flux), stage3); 

        // 4th stage
        // VectorField1D v3 = v->extract(stage3);
        // pde.Cons_to_flux(*v, &flux, sp_diff, stages, 3);
        // v->setFields(a41*v0 + a42*v3 - a43*dt*sp_diff.diff(flux), stage4);

        // Stepping   
        // VectorField1D flux3 = flux;
        // VectorField1D v4 = v->extract(stage4);
        // pde.Cons_to_flux(*v, &flux, sp_diff, stages, 4);
        // v->setFields(a51*v2 + a52*v3 - a53*dt*sp_diff.diff(flux3) + a54*v4 - 
        //     a55*dt*sp_diff.diff(flux), stage0); 

        // 1st stage
        std::vector<std::vector<int> > extractions1;
        extractions1.push_back(stage0);
        std::vector<double> coeffs1 = {1.0};
        std::vector<double> coeffs1_diff = {-a11*dt};
        // VectorField1D v0 = v->extract(stage0);
        pde.Cons_to_flux(*v, &flux, sp_diff, stages, 0);
        std::vector<VectorField*> fluxes = {&flux};
        diffLinComb(v, extractions1, stage1, coeffs1, fluxes, coeffs1_diff, sp_diff);

        // 2nd stage
        // VectorField1D v1 = v->extract(stage1);
        pde.Cons_to_flux(*v, &flux, sp_diff, stages, 1);
        std::vector<std::vector<int> > extractions2;
        extractions2.push_back(stage0);
        extractions2.push_back(stage1);
        std::vector<double> coeffs2 = {a21, a22, -a23*dt};
        diffLinComb(v, extractions2, stage2, coeffs2, flux, sp_diff);


        // 3rd stage
        // VectorField1D v2 = v->extract(stage2);
        pde.Cons_to_flux(*v, &flux, sp_diff, stages, 2);
        std::vector<std::vector<int> > extractions3;
        extractions3.push_back(stage0);
        extractions3.push_back(stage2);
        std::vector<double> coeffs3 = {a31, a32, -a33*dt};
        diffLinComb(v, extractions3, stage3, coeffs3, flux, sp_diff);

        // 4th stage
        // VectorField1D v3 = v->extract(stage3);
        pde.Cons_to_flux(*v, &flux, sp_diff, stages, 3);
        std::vector<std::vector<int> > extractions4;
        extractions4.push_back(stage0);
        extractions4.push_back(stage3);
        std::vector<double> coeffs4 = {a41, a42, -a43*dt};
        diffLinComb(v, extractions4, stage4, coeffs4, flux, sp_diff);
           
        // Stepping   
        VectorField1D flux3 = flux;
        // VectorField1D v4 = v->extract(stage4);
        pde.Cons_to_flux(*v, &flux, sp_diff, stages, 4);
        std::vector<std::vector<int> > extractions5;
        extractions5.push_back(stage2);
        extractions5.push_back(stage3);
        extractions5.push_back(stage4);
        std::vector<double> coeffs5 = {a51, a52, a54};
        std::vector<double> coeffs5_diff = {-a53*dt, -a55*dt};
        std::vector<VectorField1D*> fluxes5 = {&flux3, &flux};
        diffLinComb(v, extractions5, stage0, coeffs5, fluxes5, coeffs5_diff, sp_diff);
    }

};




#endif 