/* Time stepping scheme */

#ifndef TIMESTEPSCHEME_H
#define TIMESTEPSCHEME_H

#include "VectorField.h"
#include "Mesh.h"
#include "SpatDiffScheme.h"
#include "FC_1D.h"
#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#include "mkl.h"
#include <algorithm> 
#include "printing.h"
#include "PDE.h"
#include <array>

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


// template<typename VectorField>
// void linComb(VectorField* v, 
//     const std::vector<std::vector<int> > &extractions, 
//     const std::vector<int> &target, const std::vector<double> & coeffs,
//     const VectorField &flux)
// {
//     int ncoeffs = coeffs.size();
//     int N = v->getLength();
//     double data[v->getLength()];
//     int unknowns = v->getUnknowns();
//     for(int i = 0; i < N; i++)
//     {
//         data[i] = 0.0;
//     }
//     for(int i = 0; i < extractions[0].size(); i++) 
//     {
//         std::copy(flux.getField(i), flux.getField(i) + N, data);
//         cblas_dscal(N, coeffs[ncoeffs - 1], data, 1);
//         for(int j = 0; j < ncoeffs - 1; j++)
//         {
//             cblas_daxpy(N, coeffs[j], v->getField(extractions[j][i]), 1, data, 1);
//         }
//         v->setField(target[i], N, data);
//     }
// };



// void linComb(int N, const std::vector<std::complex<double> *> & inputs,
//     std::complex<double> * output, 
//     const std::vector<std::complex<double> > & coeffs)

// {
//     int ncoeffs = coeffs.size();
//     for(int i = 0; i < ncoeffs; i++) 
//     {
//         cblas_zaxpy(N, &coeffs[i], inputs[i], 1, output, 1);
//     }
// };




template<typename Sp_Diff>
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


template<typename Sp_Diff>
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
VectorField advance(const VectorField & v, const Sp_Diff &sp_diff)
{
    VectorField u = v - dt*sp_diff.diff(v);
    return u;
}

};


class SSPRK_4{

std::vector<int> stage0;
std::vector<int> stage1;
std::vector<int> stage2;
std::vector<int> stage3;
std::vector<int> stage4;
std::vector<std::vector<int> > extractions1;
std::vector<std::vector<int> > extractions2;
std::vector<std::vector<int> > extractions3;
std::vector<std::vector<int> > extractions4;
std::vector<std::vector<int> > extractions5;
std::vector<double> coeffs1;
std::vector<double> coeffs2;
std::vector<double> coeffs3;
std::vector<double> coeffs4;
std::vector<double> coeffs5;
std::vector<double> coeffs5_bis;
int stages;
int unknowns;
VectorField flux;
VectorField flux3;

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

public:

    SSPRK_4(int _unknowns, int _N) : unknowns{_unknowns}, flux{_unknowns, _N}, 
        flux3{_unknowns, _N}
    {
        for(int i = 0; i < unknowns; i++)
        {
            stage0.push_back(i);
            stage1.push_back(i + unknowns);
            stage2.push_back(i + 2*unknowns);
            stage3.push_back(i + 3*unknowns);
            stage4.push_back(i + 4*unknowns);
        }   
        stages = 5;   
        extractions1.push_back(stage0); 
        extractions2.push_back(stage0);
        extractions2.push_back(stage1);
        extractions3.push_back(stage0);
        extractions3.push_back(stage2);
        extractions4.push_back(stage0);
        extractions4.push_back(stage3);
        extractions5.push_back(stage2);
        extractions5.push_back(stage3);
        extractions5.push_back(stage4);

        coeffs1 = {1.0, 0.0};
        coeffs2 = {a21, a22, 0.0};
        coeffs3 = {a31, a32, 0.0};
        coeffs4 = {a41, a42, 0.0};
        coeffs5 = {a51, a52, a54, 0.0};
        coeffs5_bis = {1, 0.0};
    }


public: 

    int getStages() const {return stages;}

    template<typename Patch, typename Sp_Diff, typename PDE, typename D>
    void advance(Patch* patch, const Sp_Diff &sp_diff, const PDE &pde,
        const double dt, double t, D data)
    {
        auto v = patch->getFlowPtr();
        const int N = v->getLength();
        
        // auto flux = v->extract(stage0);

        // 1st stage
        coeffs1[1] = -a11*dt;
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 0, data);
        linComb(v, extractions1, stage1, coeffs1, flux);
        // Print_VectorField(flux);
        pde.getBC().setBC(patch, t, 1);


        // 2nd stage

        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 1, data);
        coeffs2[2] = -a23*dt;
        linComb(v, extractions2, stage2, coeffs2, flux);
        pde.getBC().setBC(patch, t, 2);


        // 3rd stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 2, data);
        coeffs3[2] = -a33*dt;
        linComb(v, extractions3, stage3, coeffs3, flux);
        pde.getBC().setBC(patch, t, 3);

        // 4th stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 3, data);
        coeffs4[2] = -a43*dt;
        linComb(v, extractions4, stage4, coeffs4, flux);
        pde.getBC().setBC(patch, t, 4);


        // Stepping   
        // VectorField flux3 = flux;
        // auto flux3 = flux;
        flux3 = flux;
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 4, data);
        coeffs5[3] = -a53*dt;
        coeffs5_bis[1] = -a55*dt;
        std::vector<VectorField*> fluxes5 = {&flux3, &flux};
        linComb(v, extractions5, stage0, coeffs5, flux3);   
        linComb(v, extractions1, stage0, coeffs5_bis, flux); 
        pde.getBC().setBC(patch, t);

    }       


    template<typename Sp_Diff, typename PDE>
    void advance(VectorField* v, const Sp_Diff &sp_diff, const PDE &pde,
        const double dt, bool flag)
    {
        const int N = v->getLength();
        const int C = sp_diff.getC();

        // VectorField flux{unknowns, N};

        // Pre-computing the derivative of the viscosity mu
        double mux[N];
        sp_diff->diff(v->getField(stages*unknowns + 1), mux);


        // 1st stage
        coeffs1[1] = -a11*dt;
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 0, mux);
        linComb(v, extractions1, stage1, coeffs1, flux);
        

        // 2nd stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 1, mux);
        coeffs2[2] = -a23*dt;
        linComb(v, extractions2, stage2, coeffs2, flux);

        // 3rd stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 2, mux);
        coeffs3[2] = -a33*dt;
        linComb(v, extractions3, stage3, coeffs3, flux);

        // 4th stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 3, mux);
        coeffs4[2] = -a43*dt;
        linComb(v, extractions4, stage4, coeffs4, flux);

        // Stepping   
        // VectorField flux3 = flux;
        flux3 = flux;
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 4, mux);
        coeffs5[3] = -a53*dt;
        coeffs5_bis[1] = -a55*dt;
        std::vector<VectorField*> fluxes5 = {&flux3, &flux};
        linComb(v, extractions5, stage0, coeffs5, flux3);   
        linComb(v, extractions1, stage0, coeffs5_bis, flux); 
    }  





    template<typename Patch, typename Sp_Diff, typename PDE>
    void advance_sdnn(Patch* patch, const Sp_Diff &sp_diff, const PDE &pde,
        const double dt, const double t, bool flag, 
        const std::vector<std::complex<double> *> & ffts, 
        const std::vector<int> & ffts_loc)
    {

        auto v = patch->getFlowPtr();
        const int N = v->getLength();
        const int C = sp_diff->getC();
        const double h = patch->getH();
        
        // auto flux = v->extract(stage0);

        // Pre-computing the derivative of the viscosity mu
        double mux[N];
        sp_diff->diff(v->getField(stages*unknowns + 1), mux, h, 0, 0);

        // 1st stage
        coeffs1[1] = -a11*dt;
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 0, mux, ffts, ffts_loc);
        linComb(v, extractions1, stage1, coeffs1, flux);


        // 2nd stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 1, mux, h, t);
        coeffs2[2] = -a23*dt;
        linComb(v, extractions2, stage2, coeffs2, flux);


        // 3rd stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 2, mux, h, t);
        coeffs3[2] = -a33*dt;
        linComb(v, extractions3, stage3, coeffs3, flux);


        // 4th stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 3, mux, h, t);
        coeffs4[2] = -a43*dt;
        linComb(v, extractions4, stage4, coeffs4, flux);

        // Stepping   
        // VectorField flux3 = flux;
        flux3 = flux;
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 4, mux, h, t);
        coeffs5[3] = -a53*dt;
        coeffs5_bis[1] = -a55*dt;
        std::vector<VectorField*> fluxes5 = {&flux3, &flux};
        linComb(v, extractions5, stage0, coeffs5, flux3);   
        linComb(v, extractions1, stage0, coeffs5_bis, flux); 
    }        

    template<typename Mesh, typename Sp_Diff, typename PDE>
    void advance_sdnn(Mesh* mesh, const std::vector<Sp_Diff> & sp_diff,
        const PDE &pde, const double dt, std::vector<double *> *mux, 
        const std::vector<std::complex<double> *> & ffts, 
        const std::vector<std::vector<int> > & ffts_loc)
    {
        auto patches = mesh->getPatchesPtr();
        int npatches = patches->size();
        

        auto v = patches->at(0)->getFlowPtr();
        // auto flux = v->extract(stage0);
        // auto flux3 = v->extract(stage0);
        int N = v->getLength();
        int C = sp_diff[0].getC();

        // 1st stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();

            // Pre-computing the derivative of the viscosity mu
            sp_diff[i].diff(v->getField(stages*unknowns + 1), mux->at(i));    

            coeffs1[1] = -a11*dt;
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 0, mux->at(i),
                ffts, ffts_loc[i]);
            linComb(v, extractions1, stage1, coeffs1, flux);
        }


        // 2nd stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();

            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 1, mux->at(i), t);
            coeffs2[2] = -a23*dt;
            linComb(v, extractions2, stage2, coeffs2, flux);
        }
        mesh->setIntraPatchBC(unknowns, 2);


        // 3rd stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();
            
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 2, mux->at(i), t);
            coeffs3[2] = -a33*dt;
            linComb(v, extractions3, stage3, coeffs3, flux);
        }
        mesh->setIntraPatchBC(unknowns, 3);


        // 4th stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();
    
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 3, mux->at(i), t);
            coeffs4[2] = -a43*dt;
            linComb(v, extractions4, stage4, coeffs4, flux);
        }
        mesh->setIntraPatchBC(unknowns, 4);


        // Stepping 
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();
      
            pde.Cons_to_der_flux(*v, &flux3, sp_diff[i], stages, 3, mux->at(i), t);
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 4, mux->at(i), t);
            coeffs5[3] = -a53*dt;
            coeffs5_bis[1] = -a55*dt;
            std::vector<VectorField*> fluxes5 = {&flux3, &flux};
            linComb(v, extractions5, stage0, coeffs5, flux3);   
            linComb(v, extractions1, stage0, coeffs5_bis, flux); 
        }
        mesh->setIntraPatchBC(unknowns, 0);

    }  

    template<typename Mesh, typename Sp_Diff, typename PDE, typename Filter>
    void advance_sdnn(Mesh* mesh, const std::vector<Sp_Diff> & sp_diff,
        const PDE &pde, const double dt, std::vector<double *> *mux, 
        const std::vector<std::complex<double> *> & ffts, 
        const std::vector<std::vector<int> > & ffts_loc,
        const std::vector<Filter> filters)
    {
        auto patches = mesh->getPatchesPtr();
        int npatches = patches->size();
        

        auto v = patches->at(0)->getFlowPtr();
        // auto flux = v->extract(stage0);
        // auto flux3 = v->extract(stage0);
        int N = v->getLength();
        int C = sp_diff[0].getC();

        // 1st stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();

            // Pre-computing the derivative of the viscosity mu
            sp_diff[i].diff(v->getField(stages*unknowns + 1), mux->at(i));    
            
            coeffs1[1] = -a11*dt;
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 0, mux->at(i),
                ffts, ffts_loc[i]);
            filters[npatches + i].filter(&flux, stage0);
            linComb(v, extractions1, stage1, coeffs1, flux);
        }


        // 2nd stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();

            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 1, mux->at(i), t);
            coeffs2[2] = -a23*dt;
            filters[npatches + i].filter(&flux, stage0);
            linComb(v, extractions2, stage2, coeffs2, flux);
        }
        mesh->setIntraPatchBC(unknowns, 2);


        // 3rd stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();
            
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 2, mux->at(i), t);
            coeffs3[2] = -a33*dt;
            filters[npatches + i].filter(&flux, stage0);
            linComb(v, extractions3, stage3, coeffs3, flux);
        }
        mesh->setIntraPatchBC(unknowns, 3);


        // 4th stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();
    
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 3, mux->at(i), t);
            coeffs4[2] = -a43*dt;
            filters[npatches + i].filter(&flux, stage0);
            linComb(v, extractions4, stage4, coeffs4, flux);
        }
        mesh->setIntraPatchBC(unknowns, 4);


        // Stepping 
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();
      
            pde.Cons_to_der_flux(*v, &flux3, sp_diff[i], stages, 3, mux->at(i), t);
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 4, mux->at(i), t);
            coeffs5[3] = -a53*dt;
            coeffs5_bis[1] = -a55*dt;
            std::vector<VectorField*> fluxes5 = {&flux3, &flux};
            filters[npatches + i].filter(&flux, stage0);
            filters[npatches + i].filter(&flux3, stage0);
            linComb(v, extractions5, stage0, coeffs5, flux3);   
            linComb(v, extractions1, stage0, coeffs5_bis, flux); 
        }
        mesh->setIntraPatchBC(unknowns, 0);

    }   


    template<typename Mesh, typename Sp_Diff, typename PDE>
    void advance_sdnn(Mesh* mesh, const std::vector<Sp_Diff> & sp_diff,
        const PDE &pde, const double dt, std::vector<double *> *mux)
    {
        auto patches = mesh->getPatchesPtr();
        int npatches = patches->size();

        auto v = patches->at(0)->getFlowPtr();
        // auto flux = v->extract(stage0);
        // auto flux3 = v->extract(stage0);
        int N = v->getLength();
        int C = sp_diff[0].getC();

        // 1st stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();

            // Pre-computing the derivative of the viscosity mu
            sp_diff[i].diff(v->getField(stages*unknowns + 1), mux->at(i), 0, 0);    
            
            coeffs1[1] = -a11*dt;
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 0, mux->at(i), t);
            linComb(v, extractions1, stage1, coeffs1, flux);
        }
        mesh->setIntraPatchBC(unknowns, 1);


        // 2nd stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();

            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 1, mux->at(i), t);
            coeffs2[2] = -a23*dt;
            linComb(v, extractions2, stage2, coeffs2, flux);
        }
        mesh->setIntraPatchBC(unknowns, 2);


        // 3rd stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();
            
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 2, mux->at(i), t);
            coeffs3[2] = -a33*dt;
            linComb(v, extractions3, stage3, coeffs3, flux);
        }
        mesh->setIntraPatchBC(unknowns, 3);


        // 4th stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();
    
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 3, mux->at(i), t);
            coeffs4[2] = -a43*dt;
            linComb(v, extractions4, stage4, coeffs4, flux);
        }
        mesh->setIntraPatchBC(unknowns, 4);


        // Stepping 
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();
      
            pde.Cons_to_der_flux(*v, &flux3, sp_diff[i], stages, 3, mux->at(i), t);
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 4, mux->at(i), t);
            coeffs5[3] = -a53*dt;
            coeffs5_bis[1] = -a55*dt;
            std::vector<VectorField*> fluxes5 = {&flux3, &flux};
            linComb(v, extractions5, stage0, coeffs5, flux3);   
            linComb(v, extractions1, stage0, coeffs5_bis, flux); 
        }
        mesh->setIntraPatchBC(unknowns, 0);

    } 








    // template<typename Mesh, typename Sp_Diff, typename PDE>
    // void advance_sdnnv2(Mesh* mesh, const std::vector<Sp_Diff> & sp_diff,
    //     const PDE &pde, const double dt, std::vector<double *> *mux, 
    //     const std::vector<std::complex<double> *> & ffts, 
    //     const std::vector<std::vector<int> > & ffts_loc)
    // {
    //     auto patches = mesh->getPatchesPtr();
    //     int npatches = patches->size();
        
    //     const int unknowns = pde.getPhysUnknowns();

    //     std::vector<int> stage0;
    //     std::vector<int> stage1;
    //     std::vector<int> stage2;
    //     std::vector<int> stage3;
    //     std::vector<int> stage4;
    //     for(int i = 0; i < unknowns; i++)
    //     {
    //         stage0.push_back(i);
    //         stage1.push_back(i + unknowns);
    //         stage2.push_back(i + 2*unknowns);
    //         stage3.push_back(i + 3*unknowns);
    //         stage4.push_back(i + 4*unknowns);
    //     }

    //     auto v = patches->at(0)->getFlowPtr();
    //     auto flux = v->extract(stage0);
    //     auto flux3 = v->extract(stage0);
    //     int N = v->getLength();
    //     int C = sp_diff[0].getC();

    //     // int i = 0;

    //     for(int i = 0; i < npatches; i++)
    //     {
    //         v = patches->at(i)->getFlowPtr();
    //         flux = v->extract(stage0);
    //         N = v->getLength();
    //         C = sp_diff[i].getC();

    //         // Pre-computing the derivative of the viscosity mu

    //         sp_diff[i].diff(v->getField(stages*unknowns + 1), mux->at(i));


    //         // 1st stage
    //         std::vector<std::vector<int> > extractions1;
    //         extractions1.push_back(stage0);
    //         std::vector<double> coeffs1 = {1.0, -a11*dt};
    //         // pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 0, mux->at(i));
    //         pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 0, mux->at(i),
    //             ffts, ffts_loc[i]);
    //         linComb(v, extractions1, stage1, coeffs1, flux);

    //         // Print_VectorField(v->extract(stage1), true);



    //         // 2nd stage
    //         pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 1, mux->at(i));
    //         std::vector<std::vector<int> > extractions2;
    //         extractions2.push_back(stage0);
    //         extractions2.push_back(stage1);
    //         std::vector<double> coeffs2 = {a21, a22, -a23*dt};
    //         linComb(v, extractions2, stage2, coeffs2, flux);

    //         Print_VectorField(v->extract(stage2), true);


    //         // 3rd stage
    //         pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 2, mux->at(i));
    //         std::vector<std::vector<int> > extractions3;
    //         extractions3.push_back(stage0);
    //         extractions3.push_back(stage2);
    //         std::vector<double> coeffs3 = {a31, a32, -a33*dt};
    //         linComb(v, extractions3, stage3, coeffs3, flux);
 

    //         // 4th stage
    //         pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 3, mux->at(i));
    //         std::vector<std::vector<int> > extractions4;
    //         extractions4.push_back(stage0);
    //         extractions4.push_back(stage3);
    //         std::vector<double> coeffs4 = {a41, a42, -a43*dt};
    //         linComb(v, extractions4, stage4, coeffs4, flux);

    //         // Stepping   
    //         flux3 = flux;
    //         pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 4, mux->at(i));
    //         std::vector<std::vector<int> > extractions5;
    //         extractions5.push_back(stage2);
    //         extractions5.push_back(stage3);
    //         extractions5.push_back(stage4); 
    //         std::vector<double> coeffs5 = {a51, a52, a54, -a53*dt};
    //         std::vector<double> coeffs5_bis = {1, -a55*dt};
    //         std::vector<VectorField*> fluxes5 = {&flux3, &flux};
    //         linComb(v, extractions5, stage0, coeffs5, flux3);   
    //         linComb(v, extractions1, stage0, coeffs5_bis, flux); 
    //     }




    // }    
    

};




#endif 