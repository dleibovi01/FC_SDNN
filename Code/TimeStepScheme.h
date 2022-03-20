/* Time stepping scheme */

#ifndef TIMESTEPSCHEME_H
#define TIMESTEPSCHEME_H

#include "VectorField1D.h"
#include "Mesh.h"
#include "SpatDiffScheme.h"
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
        pde.Cons_to_flux(*v, &flux, sp_diff, stages, 0);
        std::vector<VectorField*> fluxes = {&flux};
        diffLinComb(v, extractions1, stage1, coeffs1, fluxes, coeffs1_diff, sp_diff);

        // 2nd stage
        pde.Cons_to_flux(*v, &flux, sp_diff, stages, 1);
        std::vector<std::vector<int> > extractions2;
        extractions2.push_back(stage0);
        extractions2.push_back(stage1);
        std::vector<double> coeffs2 = {a21, a22, -a23*dt};
        diffLinComb(v, extractions2, stage2, coeffs2, flux, sp_diff);


        // 3rd stage
        pde.Cons_to_flux(*v, &flux, sp_diff, stages, 2);
        std::vector<std::vector<int> > extractions3;
        extractions3.push_back(stage0);
        extractions3.push_back(stage2);
        std::vector<double> coeffs3 = {a31, a32, -a33*dt};
        diffLinComb(v, extractions3, stage3, coeffs3, flux, sp_diff);

        // 4th stage
        pde.Cons_to_flux(*v, &flux, sp_diff, stages, 3);
        std::vector<std::vector<int> > extractions4;
        extractions4.push_back(stage0);
        extractions4.push_back(stage3);
        std::vector<double> coeffs4 = {a41, a42, -a43*dt};
        diffLinComb(v, extractions4, stage4, coeffs4, flux, sp_diff);

           
        // Stepping   
        VectorField1D flux3 = flux;
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


    template<typename VectorField, typename Sp_Diff, typename PDE>
    void advance_fast(VectorField* v, 
        std::vector<std::complex<double> *> * cons_fft,
        std::vector<std::complex<double> *> * flux_fft,
        const Sp_Diff &sp_diff, const PDE &pde, const double dt)
    {
        const int unknowns = pde.getPhysUnknowns();
        const int N = v->getLength();
        const int C = sp_diff.getC();
        int d = sp_diff.getD();
        double fourPts_dbl = sp_diff.getFourPts_dbl();

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


        VectorField flux{unknowns, N};

        std::vector<std::complex<double> * > inputs;
        std::complex<double> fft_temp[N + C];
        std::complex<double> temp_cons[N+C];
        std::complex<double> temp_cons2[N+C];
        std::complex<double> temp_flux[N+C];
        


        // for(int i = 0; i < unknowns; i++)
        // {
        //     Fcont_Gram_Blend(v->getField(0), fft_temp, N, d, C, fourPts_dbl, 
        //         sp_diff.getAQ(), sp_diff.getFAQF(), sp_diff.getDescHandle());
        //     cons_fft_1.push_back(new std::complex<double>[N+C])
        // }
        // 1st stage
        // std::vector<std::vector<int> > extractions1;
        // extractions1.push_back(stage0);
        // std::vector<double> coeffs1 = {1.0};
        // std::vector<double> coeffs1_diff = {-a11*dt};
        // pde.Cons_to_flux(*v, &flux, sp_diff, stages, 0);
        // std::vector<VectorField*> fluxes = {&flux};
        // diffLinComb(v, extractions1, stage1, coeffs1, fluxes, coeffs1_diff, sp_diff);
        std::vector<std::vector<int> > extractions1;
        extractions1.push_back(stage0);
        std::vector<double> coeffs1 = {1.0, -a11*dt};
        std::vector<std::complex<double> > zcoeffs1 = {1.0, -a11*dt};
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 0, cons_fft, & flux_fft);
        linComb(v, extractions1, stage1, coeffs1, flux);
        for(int i = 0; i < unknowns; i++)
        {
            std::copy(cons_fft[i].begin(), cons_fft[i].end(), temp_cons);
            inputs.push_back(temp_cons);
            std::copy(flux_fft[i].begin(), flux_fft[i].end(), temp_flux);
            inputs.push_back(temp_flux);
            linComb(N + C, inputs, fft_temp, zcoeffs1);
            std::copy(fft_temp, fft_temp + N + C, cons_fft[i + unknowns].begin());
        }     



        // // 2nd stage
        // pde.Cons_to_flux(*v, &flux, sp_diff, stages, 1);
        // std::vector<std::vector<int> > extractions2;
        // extractions2.push_back(stage0);
        // extractions2.push_back(stage1);
        // std::vector<double> coeffs2 = {a21, a22, -a23*dt};
        // diffLinComb(v, extractions2, stage2, coeffs2, flux, sp_diff);
        std::vector<std::vector<int> > extractions2;
        extractions2.push_back(stage0);
        extractions2.push_back(stage1);
        std::vector<double> coeffs2 = {a21, a22, -a23*dt};
        std::vector<std::complex<double> > zcoeffs2 = {a21, a22, -a23*dt};
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 1, cons_fft, & flux_fft);
        linComb(v, extractions1, stage1, coeffs1, flux);
        for(int i = 0; i < unknowns; i++)
        {
            // std::copy(cons_fft[i].begin(), cons_fft[i].end(), temp_cons);
            // inputs.push_back(temp_cons);
            std::copy(cons_fft[i + unknowns].begin(), 
                cons_fft[i + unknowns].end(), temp_cons2);
            inputs.push_back(temp_cons2);
            std::copy(flux_fft[i].begin(), flux_fft[i].end(), temp_flux);
            inputs.push_back(temp_flux);
            linComb(N + C, inputs, fft_temp, zcoeffs1);
            std::copy(fft_temp, fft_temp + N + C, cons_fft[i + 2*unknowns].begin());
        }   


        // // 3rd stage
        // pde.Cons_to_flux(*v, &flux, sp_diff, stages, 2);
        // std::vector<std::vector<int> > extractions3;
        // extractions3.push_back(stage0);
        // extractions3.push_back(stage2);
        // std::vector<double> coeffs3 = {a31, a32, -a33*dt};
        // diffLinComb(v, extractions3, stage3, coeffs3, flux, sp_diff);
        std::vector<std::vector<int> > extractions3;
        extractions3.push_back(stage0);
        extractions3.push_back(stage2);
        std::vector<double> coeffs3 = {a31, a32, -a33*dt};
        std::vector<std::complex<double> > zcoeffs3 = {a31, a32, -a33*dt};
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 2, cons_fft, & flux_fft);
        linComb(v, extractions2, stage2, coeffs2, flux);
        for(int i = 0; i < unknowns; i++)
        {
            // std::copy(cons_fft[i].begin(), cons_fft[i].end(), temp_cons);
            // inputs.push_back(temp_cons);
            inputs.pop_back();
            std::copy(cons_fft[i + 2*unknowns].begin(), 
                cons_fft[i + 2*unknowns].end(), temp_cons2);
            inputs.push_back(temp_cons2);
            std::copy(flux_fft[i].begin(), flux_fft[i].end(), temp_flux);
            inputs.push_back(temp_flux);
            linComb(N + C, inputs, fft_temp, zcoeffs2);
            std::copy(fft_temp, fft_temp + N + C, cons_fft[i + 3*unknowns].begin());
        }   

        // // 4th stage
        // pde.Cons_to_flux(*v, &flux, sp_diff, stages, 3);
        // std::vector<std::vector<int> > extractions4;
        // extractions4.push_back(stage0);
        // extractions4.push_back(stage3);
        // std::vector<double> coeffs4 = {a41, a42, -a43*dt};
        // diffLinComb(v, extractions4, stage4, coeffs4, flux, sp_diff);
        std::vector<std::vector<int> > extractions3;
        extractions3.push_back(stage0);
        extractions3.push_back(stage3);
        std::vector<double> coeffs4 = {a41, a42, -a43*dt};
        std::vector<std::complex<double> > zcoeffs4 = {a41, a42, -a43*dt};
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 3, cons_fft, & flux_fft);
        linComb(v, extractions3, stage3, coeffs3, flux);
        for(int i = 0; i < unknowns; i++)
        {
            inputs.pop_back();
            // std::copy(cons_fft[i].begin(), cons_fft[i].end(), temp_cons);
            // inputs.push_back(temp_cons);
            std::copy(cons_fft[i + 3*unknowns].begin(), 
                cons_fft[i + 3*unknowns].end(), temp_cons2);
            inputs.push_back(temp_cons2);
            std::copy(flux_fft[i].begin(), flux_fft[i].end(), temp_flux);
            inputs.push_back(temp_flux);
            linComb(N + C, inputs, fft_temp, zcoeffs2);
            std::copy(fft_temp, fft_temp + N + C, cons_fft[i + 4*unknowns].begin());
        }  

           
        // // Stepping   
        // VectorField1D flux3 = flux;
        // pde.Cons_to_flux(*v, &flux, sp_diff, stages, 4);
        // std::vector<std::vector<int> > extractions5;
        // extractions5.push_back(stage2);
        // extractions5.push_back(stage3);
        // extractions5.push_back(stage4);
        // std::vector<double> coeffs5 = {a51, a52, a54};
        // std::vector<double> coeffs5_diff = {-a53*dt, -a55*dt};
        // std::vector<VectorField1D*> fluxes5 = {&flux3, &flux};
        // diffLinComb(v, extractions5, stage0, coeffs5, fluxes5, coeffs5_diff, sp_diff);
    }

    // template<typename VectorField, typename Sp_Diff, typename PDE>
    // void advance(VectorField* v, const Sp_Diff &sp_diff, const PDE &pde,
    //     const double dt, double t)
    template<typename Patch, typename Sp_Diff, typename PDE, typename D>
    void advance(Patch* patch, const Sp_Diff &sp_diff, const PDE &pde,
        const double dt, double t, D data)
    {
        auto v = patch->getFlowPtr();
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

        // VectorField flux{unknowns, N};
        auto flux = v->extract(stage0);



        // 1st stage
        std::vector<std::vector<int> > extractions1;
        extractions1.push_back(stage0);
        std::vector<double> coeffs1 = {1.0, -a11*dt};
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 0, data);
        linComb(v, extractions1, stage1, coeffs1, flux);
        // Print_VectorField1D(flux);
        pde.getBC().setBC(patch, t, 1);

        // std::cout << std::endl;
        // std::cout << std::endl;
        // std::cout << "1st stage" << std::endl;
        // Print_VectorField1D(v->extract(stage1), true, 17);

        // 2nd stage

        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 1, data);
        std::vector<std::vector<int> > extractions2;
        extractions2.push_back(stage0);
        extractions2.push_back(stage1);
        std::vector<double> coeffs2 = {a21, a22, -a23*dt};
        linComb(v, extractions2, stage2, coeffs2, flux);
        pde.getBC().setBC(patch, t, 2);

        // std::cout << std::endl;
        // std::cout << std::endl;
        // std::cout << "2nd stage" << std::endl;
        // Print_VectorField1D(v->extract(stage2), true, 17);

        // 3rd stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 2, data);
        std::vector<std::vector<int> > extractions3;
        extractions3.push_back(stage0);
        extractions3.push_back(stage2);
        std::vector<double> coeffs3 = {a31, a32, -a33*dt};
        linComb(v, extractions3, stage3, coeffs3, flux);
        pde.getBC().setBC(patch, t, 3);

        // std::cout << std::endl;
        // std::cout << std::endl;
        // std::cout << "3rd stage" << std::endl;
        // Print_VectorField1D(v->extract(stage3), true, 17);

        // Print_VectorField1D(*v);

        // 4th stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 3, data);
        std::vector<std::vector<int> > extractions4;
        extractions4.push_back(stage0);
        extractions4.push_back(stage3);
        std::vector<double> coeffs4 = {a41, a42, -a43*dt};
        linComb(v, extractions4, stage4, coeffs4, flux);
        pde.getBC().setBC(patch, t, 4);

        // std::cout << std::endl;
        // std::cout << std::endl;
        // std::cout << "4th stage" << std::endl;
        // Print_VectorField1D(*v);

        // Stepping   
        // VectorField1D flux3 = flux;
        auto flux3 = flux;
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 4, data);
        std::vector<std::vector<int> > extractions5;
        extractions5.push_back(stage2);
        extractions5.push_back(stage3);
        extractions5.push_back(stage4); 
        std::vector<double> coeffs5 = {a51, a52, a54, -a53*dt};
        std::vector<double> coeffs5_bis = {1, -a55*dt};
        std::vector<VectorField1D*> fluxes5 = {&flux3, &flux};
        linComb(v, extractions5, stage0, coeffs5, flux3);   
        linComb(v, extractions1, stage0, coeffs5_bis, flux); 
        pde.getBC().setBC(patch, t);

        // std::cout << " after first step " << std::endl;
        // Print_VectorField1D(*v, true);
    }       


    template<typename VectorField, typename Sp_Diff, typename PDE>
    void advance(VectorField* v, const Sp_Diff &sp_diff, const PDE &pde,
        const double dt, bool flag)
    {
        const int unknowns = pde.getPhysUnknowns();
        const int N = v->getLength();
        const int C = sp_diff.getC();
  

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

        VectorField flux{unknowns, N};



        // Pre-computing the derivative of the viscosity mu
        double mux[N];
        sp_diff.diff(v->getField(stages*unknowns + 1), mux);


        // 1st stage
        std::vector<std::vector<int> > extractions1;
        extractions1.push_back(stage0);
        std::vector<double> coeffs1 = {1.0, -a11*dt};
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 0, mux);
        linComb(v, extractions1, stage1, coeffs1, flux);
        

        // 2nd stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 1, mux);
        std::vector<std::vector<int> > extractions2;
        extractions2.push_back(stage0);
        extractions2.push_back(stage1);
        std::vector<double> coeffs2 = {a21, a22, -a23*dt};
        linComb(v, extractions2, stage2, coeffs2, flux);

        // 3rd stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 2, mux);
        std::vector<std::vector<int> > extractions3;
        extractions3.push_back(stage0);
        extractions3.push_back(stage2);
        std::vector<double> coeffs3 = {a31, a32, -a33*dt};
        linComb(v, extractions3, stage3, coeffs3, flux);

        // 4th stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 3, mux);
        std::vector<std::vector<int> > extractions4;
        extractions4.push_back(stage0);
        extractions4.push_back(stage3);
        std::vector<double> coeffs4 = {a41, a42, -a43*dt};
        linComb(v, extractions4, stage4, coeffs4, flux);

        // Stepping   
        VectorField1D flux3 = flux;
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 4, mux);
        std::vector<std::vector<int> > extractions5;
        extractions5.push_back(stage2);
        extractions5.push_back(stage3);
        extractions5.push_back(stage4); 
        std::vector<double> coeffs5 = {a51, a52, a54, -a53*dt};
        std::vector<double> coeffs5_bis = {1, -a55*dt};
        std::vector<VectorField1D*> fluxes5 = {&flux3, &flux};
        linComb(v, extractions5, stage0, coeffs5, flux3);   
        linComb(v, extractions1, stage0, coeffs5_bis, flux); 
    }  





    template<typename VectorField, typename Sp_Diff, typename PDE>
    void advance_sdnn(VectorField* v, const Sp_Diff &sp_diff, const PDE &pde,
        const double dt, bool flag, 
        const std::vector<std::complex<double> *> & ffts, 
        const std::vector<int> & ffts_loc)
    {
        const int unknowns = pde.getPhysUnknowns();
        const int N = v->getLength();
        const int C = sp_diff.getC();
  

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

        VectorField flux{unknowns, N};

        

        // Pre-computing the derivative of the viscosity mu
        double mux[N];
        sp_diff.diff(v->getField(stages*unknowns + 1), mux);


        // 1st stage
        std::vector<std::vector<int> > extractions1;
        extractions1.push_back(stage0);
        std::vector<double> coeffs1 = {1.0, -a11*dt};
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 0, mux, ffts, ffts_loc);
        linComb(v, extractions1, stage1, coeffs1, flux);
        // Print_VectorField1D(v->extract(stage1), true);
        // Print_VectorField1D(flux);
        // Print_VectorField1D(*v);

        // 2nd stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 1, mux);
        std::vector<std::vector<int> > extractions2;
        extractions2.push_back(stage0);
        extractions2.push_back(stage1);
        std::vector<double> coeffs2 = {a21, a22, -a23*dt};
        linComb(v, extractions2, stage2, coeffs2, flux);
        // Print_VectorField1D(v->extract(stage2), true);

        // 3rd stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 2, mux);
        std::vector<std::vector<int> > extractions3;
        extractions3.push_back(stage0);
        extractions3.push_back(stage2);
        std::vector<double> coeffs3 = {a31, a32, -a33*dt};
        linComb(v, extractions3, stage3, coeffs3, flux);

        // 4th stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 3, mux);
        std::vector<std::vector<int> > extractions4;
        extractions4.push_back(stage0);
        extractions4.push_back(stage3);
        std::vector<double> coeffs4 = {a41, a42, -a43*dt};
        linComb(v, extractions4, stage4, coeffs4, flux);

        // Stepping   
        VectorField1D flux3 = flux;
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 4, mux);
        std::vector<std::vector<int> > extractions5;
        extractions5.push_back(stage2);
        extractions5.push_back(stage3);
        extractions5.push_back(stage4); 
        std::vector<double> coeffs5 = {a51, a52, a54, -a53*dt};
        std::vector<double> coeffs5_bis = {1, -a55*dt};
        std::vector<VectorField1D*> fluxes5 = {&flux3, &flux};
        linComb(v, extractions5, stage0, coeffs5, flux3);   
        linComb(v, extractions1, stage0, coeffs5_bis, flux); 
        // Print_VectorField1D(v->extract(stage0), true);
        // std::cout << std::endl;
    }        

    template<typename Mesh, typename Sp_Diff, typename PDE>
    void advance_sdnn(Mesh* mesh, const std::vector<Sp_Diff> & sp_diff,
        const PDE &pde, const double dt, std::vector<double *> *mux, 
        const std::vector<std::complex<double> *> & ffts, 
        const std::vector<std::vector<int> > & ffts_loc)
    {
        auto patches = mesh->getPatchesPtr();
        int npatches = patches->size();
        
        const int unknowns = pde.getPhysUnknowns();

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

        auto v = patches->at(0)->getFlowPtr();
        auto flux = v->extract(stage0);
        auto flux3 = v->extract(stage0);
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

            std::vector<std::vector<int> > extractions1;
            extractions1.push_back(stage0);
            std::vector<double> coeffs1 = {1.0, -a11*dt};
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 0, mux->at(i),
                ffts, ffts_loc[i]);
            linComb(v, extractions1, stage1, coeffs1, flux);
        }

        // std::cout << "Patch 0" << std::endl;
        // Print_VectorField1D(patches->at(0)->getFlowPtr()->extract(stage1), true);
        // std::cout << "Patch 1" << std::endl;
        // Print_VectorField1D(patches->at(1)->getFlowPtr()->extract(stage1), true);
        // mesh->setIntraPatchBC(unknowns, 1);
        // std::cout << std::endl;
        // std::cout << "Patch 0" << std::endl;
        // Print_VectorField1D(patches->at(0)->getFlowPtr()->extract(stage1), true);
        // std::cout << "Patch 1" << std::endl;
        // Print_VectorField1D(patches->at(1)->getFlowPtr()->extract(stage1), true);


        // 2nd stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();

            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 1, mux->at(i));
            std::vector<std::vector<int> > extractions2;
            extractions2.push_back(stage0);
            extractions2.push_back(stage1);
            std::vector<double> coeffs2 = {a21, a22, -a23*dt};
            linComb(v, extractions2, stage2, coeffs2, flux);
        }
        mesh->setIntraPatchBC(unknowns, 2);

        // Print_VectorField1D(patches->at(0)->getFlowPtr()->extract(stage2), true);
            // Print_VectorField1D(v->extract(stage2), true);


        // 3rd stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();
            
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 2, mux->at(i));
            std::vector<std::vector<int> > extractions3;
            extractions3.push_back(stage0);
            extractions3.push_back(stage2);
            std::vector<double> coeffs3 = {a31, a32, -a33*dt};
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
    
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 3, mux->at(i));
            std::vector<std::vector<int> > extractions4;
            extractions4.push_back(stage0);
            extractions4.push_back(stage3);
            std::vector<double> coeffs4 = {a41, a42, -a43*dt};
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
      
            pde.Cons_to_der_flux(*v, &flux3, sp_diff[i], stages, 3, mux->at(i));
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 4, mux->at(i));
            std::vector<std::vector<int> > extractions5;
            std::vector<std::vector<int> > extractions1;
            extractions1.push_back(stage0);
            extractions5.push_back(stage2);
            extractions5.push_back(stage3);
            extractions5.push_back(stage4); 
            std::vector<double> coeffs5 = {a51, a52, a54, -a53*dt};
            std::vector<double> coeffs5_bis = {1, -a55*dt};
            std::vector<VectorField1D*> fluxes5 = {&flux3, &flux};
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
        
        const int unknowns = pde.getPhysUnknowns();

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

        auto v = patches->at(0)->getFlowPtr();
        auto flux = v->extract(stage0);
        auto flux3 = v->extract(stage0);
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
            
            std::vector<std::vector<int> > extractions1;
            extractions1.push_back(stage0);
            std::vector<double> coeffs1 = {1.0, -a11*dt};
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 0, mux->at(i),
                ffts, ffts_loc[i]);
            filters[npatches + i].filter(&flux, stage0);
            linComb(v, extractions1, stage1, coeffs1, flux);
        }

        // std::cout << "Patch 0" << std::endl;
        // Print_VectorField1D(patches->at(0)->getFlowPtr()->extract(stage1), true);
        // std::cout << "Patch 1" << std::endl;
        // Print_VectorField1D(patches->at(1)->getFlowPtr()->extract(stage1), true);
        // mesh->setIntraPatchBC(unknowns, 1);
        // std::cout << std::endl;
        // std::cout << "Patch 0" << std::endl;
        // Print_VectorField1D(patches->at(0)->getFlowPtr()->extract(stage1), true);
        // std::cout << "Patch 1" << std::endl;
        // Print_VectorField1D(patches->at(1)->getFlowPtr()->extract(stage1), true);


        // 2nd stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();

            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 1, mux->at(i));
            std::vector<std::vector<int> > extractions2;
            extractions2.push_back(stage0);
            extractions2.push_back(stage1);
            std::vector<double> coeffs2 = {a21, a22, -a23*dt};
            filters[npatches + i].filter(&flux, stage0);
            linComb(v, extractions2, stage2, coeffs2, flux);
        }
        mesh->setIntraPatchBC(unknowns, 2);

        // Print_VectorField1D(patches->at(0)->getFlowPtr()->extract(stage2), true);
            // Print_VectorField1D(v->extract(stage2), true);


        // 3rd stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();
            
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 2, mux->at(i));
            std::vector<std::vector<int> > extractions3;
            extractions3.push_back(stage0);
            extractions3.push_back(stage2);
            std::vector<double> coeffs3 = {a31, a32, -a33*dt};
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
    
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 3, mux->at(i));
            std::vector<std::vector<int> > extractions4;
            extractions4.push_back(stage0);
            extractions4.push_back(stage3);
            std::vector<double> coeffs4 = {a41, a42, -a43*dt};
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
      
            pde.Cons_to_der_flux(*v, &flux3, sp_diff[i], stages, 3, mux->at(i));
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 4, mux->at(i));
            std::vector<std::vector<int> > extractions5;
            std::vector<std::vector<int> > extractions1;
            extractions1.push_back(stage0);
            extractions5.push_back(stage2);
            extractions5.push_back(stage3);
            extractions5.push_back(stage4); 
            std::vector<double> coeffs5 = {a51, a52, a54, -a53*dt};
            std::vector<double> coeffs5_bis = {1, -a55*dt};
            std::vector<VectorField1D*> fluxes5 = {&flux3, &flux};
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
        
        const int unknowns = pde.getPhysUnknowns();

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

        auto v = patches->at(0)->getFlowPtr();
        auto flux = v->extract(stage0);
        auto flux3 = v->extract(stage0);
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
            
            std::vector<std::vector<int> > extractions1;
            extractions1.push_back(stage0);
            std::vector<double> coeffs1 = {1.0, -a11*dt};
            // pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 0, mux->at(i),
            //     ffts, ffts_loc[i]);
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 0, mux->at(i));
            linComb(v, extractions1, stage1, coeffs1, flux);
        }
        mesh->setIntraPatchBC(unknowns, 1);

        // std::cout << "Patch 0" << std::endl;
        // Print_VectorField1D(patches->at(0)->getFlowPtr()->extract(stage1), true);
        // std::cout << "Patch 1" << std::endl;
        // Print_VectorField1D(patches->at(1)->getFlowPtr()->extract(stage1), true);
        // mesh->setIntraPatchBC(unknowns, 1);
        // std::cout << std::endl;
        // std::cout << "Patch 0" << std::endl;
        // Print_VectorField1D(patches->at(0)->getFlowPtr()->extract(stage1), true);
        // std::cout << "Patch 1" << std::endl;
        // Print_VectorField1D(patches->at(1)->getFlowPtr()->extract(stage1), true);


        // 2nd stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();

            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 1, mux->at(i));
            std::vector<std::vector<int> > extractions2;
            extractions2.push_back(stage0);
            extractions2.push_back(stage1);
            std::vector<double> coeffs2 = {a21, a22, -a23*dt};
            linComb(v, extractions2, stage2, coeffs2, flux);
        }
        mesh->setIntraPatchBC(unknowns, 2);

        // Print_VectorField1D(patches->at(0)->getFlowPtr()->extract(stage2), true);
            // Print_VectorField1D(v->extract(stage2), true);


        // 3rd stage
        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();
            
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 2, mux->at(i));
            std::vector<std::vector<int> > extractions3;
            extractions3.push_back(stage0);
            extractions3.push_back(stage2);
            std::vector<double> coeffs3 = {a31, a32, -a33*dt};
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
    
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 3, mux->at(i));
            std::vector<std::vector<int> > extractions4;
            extractions4.push_back(stage0);
            extractions4.push_back(stage3);
            std::vector<double> coeffs4 = {a41, a42, -a43*dt};
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
      
            pde.Cons_to_der_flux(*v, &flux3, sp_diff[i], stages, 3, mux->at(i));
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 4, mux->at(i));
            std::vector<std::vector<int> > extractions5;
            std::vector<std::vector<int> > extractions1;
            extractions1.push_back(stage0);
            extractions5.push_back(stage2);
            extractions5.push_back(stage3);
            extractions5.push_back(stage4); 
            std::vector<double> coeffs5 = {a51, a52, a54, -a53*dt};
            std::vector<double> coeffs5_bis = {1, -a55*dt};
            std::vector<VectorField1D*> fluxes5 = {&flux3, &flux};
            linComb(v, extractions5, stage0, coeffs5, flux3);   
            linComb(v, extractions1, stage0, coeffs5_bis, flux); 
        }
        mesh->setIntraPatchBC(unknowns, 0);

    } 








    template<typename Mesh, typename Sp_Diff, typename PDE>
    void advance_sdnnv2(Mesh* mesh, const std::vector<Sp_Diff> & sp_diff,
        const PDE &pde, const double dt, std::vector<double *> *mux, 
        const std::vector<std::complex<double> *> & ffts, 
        const std::vector<std::vector<int> > & ffts_loc)
    {
        auto patches = mesh->getPatchesPtr();
        int npatches = patches->size();
        
        const int unknowns = pde.getPhysUnknowns();

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

        auto v = patches->at(0)->getFlowPtr();
        auto flux = v->extract(stage0);
        auto flux3 = v->extract(stage0);
        int N = v->getLength();
        int C = sp_diff[0].getC();

        // int i = 0;

        for(int i = 0; i < npatches; i++)
        {
            v = patches->at(i)->getFlowPtr();
            flux = v->extract(stage0);
            N = v->getLength();
            C = sp_diff[i].getC();

            // Pre-computing the derivative of the viscosity mu

            sp_diff[i].diff(v->getField(stages*unknowns + 1), mux->at(i));


            // 1st stage
            std::vector<std::vector<int> > extractions1;
            extractions1.push_back(stage0);
            std::vector<double> coeffs1 = {1.0, -a11*dt};
            // pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 0, mux->at(i));
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 0, mux->at(i),
                ffts, ffts_loc[i]);
            linComb(v, extractions1, stage1, coeffs1, flux);

            // Print_VectorField1D(v->extract(stage1), true);



            // 2nd stage
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 1, mux->at(i));
            std::vector<std::vector<int> > extractions2;
            extractions2.push_back(stage0);
            extractions2.push_back(stage1);
            std::vector<double> coeffs2 = {a21, a22, -a23*dt};
            linComb(v, extractions2, stage2, coeffs2, flux);

            Print_VectorField1D(v->extract(stage2), true);


            // 3rd stage
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 2, mux->at(i));
            std::vector<std::vector<int> > extractions3;
            extractions3.push_back(stage0);
            extractions3.push_back(stage2);
            std::vector<double> coeffs3 = {a31, a32, -a33*dt};
            linComb(v, extractions3, stage3, coeffs3, flux);
 

            // 4th stage
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 3, mux->at(i));
            std::vector<std::vector<int> > extractions4;
            extractions4.push_back(stage0);
            extractions4.push_back(stage3);
            std::vector<double> coeffs4 = {a41, a42, -a43*dt};
            linComb(v, extractions4, stage4, coeffs4, flux);

            // Stepping   
            flux3 = flux;
            pde.Cons_to_der_flux(*v, &flux, sp_diff[i], stages, 4, mux->at(i));
            std::vector<std::vector<int> > extractions5;
            extractions5.push_back(stage2);
            extractions5.push_back(stage3);
            extractions5.push_back(stage4); 
            std::vector<double> coeffs5 = {a51, a52, a54, -a53*dt};
            std::vector<double> coeffs5_bis = {1, -a55*dt};
            std::vector<VectorField1D*> fluxes5 = {&flux3, &flux};
            linComb(v, extractions5, stage0, coeffs5, flux3);   
            linComb(v, extractions1, stage0, coeffs5_bis, flux); 
        }




    }    
    

};




#endif 