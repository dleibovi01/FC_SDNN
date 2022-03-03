/* PDE solver */

#ifndef SOLVER_H
#define SOLVER_H

#include "SpatDiffScheme.h"
#include "TimeStepScheme.h"
#include "Mesh.h"
#include "PDE.h"
#include "printing.h"
#include "SDNN.h"
#include <chrono>
#include <array>

template<typename Meshtype, typename PDE, typename TS, typename Sp_diff, 
    typename Filter>
class Solver{
    Meshtype mesh;
    PDE pde;
    TS ts;
    const std::vector<Sp_diff> sp_diffs;
    const std::vector<Filter> filters;

    // ArtVisc artvisc;
    public:

    Solver(const Meshtype &_mesh, const PDE &_pde, TS &_ts, 
        const std::vector<Sp_diff> &_sp_diff, 
        const std::vector<Filter> &_filter) : mesh{_mesh}, pde{_pde}, ts{_ts},
         sp_diffs{_sp_diff}, filters{_filter} {}

    // const Meshtype& getMesh() const {return mesh;}
    const Meshtype& getMesh() const {return mesh;}


    /*
    void solve(double dt)
    {
        double t = 0.0;
        double T = pde.getT();
        // std::vector<Patch*> patches = mesh.getPatches();
        auto patches = mesh.getPatches();
        int npatches = patches.size();
        pde.getIC()(&mesh);
        while (t < T)
        {
            patches = mesh.getPatches();
            for(int i = 0; i < npatches; i++)
            {
                // VectorField v = patches[i]->getFlow();
                auto v = patches[i]->getFlow();
                v = pde.Prim_to_cons(v);
                v = ts.advance(v, sp_diff);
                v = pde.Cons_to_prim(v);
                patches[i]->setField(v);
            }
            mesh.setPatches(patches);
            if(t + dt > T)
            {
                t = T; 
            }
            else
            {
                t += dt;
            }
            mesh.setIntraPatchBC();
            pde.getBC()(&mesh, t);
        }
    }
    */

    void solve_sdnn(double dt, double CFL, bool adaptive, bool visc)
    {
        std::cout << "Solver" << std::endl;
        double t = 0.0;
        double T = pde.getT();
        auto patches = mesh.getPatches();
        int npatches = sp_diffs.size();
        int phys_unknowns = pde.getPhysUnknowns();
        int stages = ts.getStages();
        pde.getIC()(&mesh);

        // Filter data      
        std::vector<int> filt_unknowns;
        for(int i = 0; i < phys_unknowns; i++)
        {
            filt_unknowns.push_back(i);
        }

        // Set viscosity normalization coefficients
        SVW_mesh svw_m{mesh};

        // Pre-allocate a container to save the FFTS
        std::vector<std::complex<double> * > fft_data;
        std::vector<std::vector<int> > fft_locs;
        std::vector<int> loc;
        for(int i = 0; i < npatches; i++)
        {
            for(int j = 0; j < phys_unknowns; j++)
            {
                fft_data.push_back(
                    new std::complex<double> [mesh.getPatches()[i]->getNnodes() + 
                    sp_diffs[i].getC()]);
                for(int k = 0; k < mesh.getPatches()[i]->getNnodes() + 
                    sp_diffs[i].getC(); k++)
                {
                    fft_data[i*npatches + j][k] = 0.0;
                }
                loc.push_back(i*phys_unknowns + j);
            }
            fft_locs.push_back(loc);
            loc.clear();
        }



        auto t1 = std::chrono::high_resolution_clock::now();
      
        
        // Time loop
        while (t < T)
        {
            patches = mesh.getPatches();
            // Smooth viscosity assignment
            if(visc)
            {
                for(int i = 0; i < npatches; i++)
                {
                    updateTau(patches[i], sp_diffs[i], pde, 
                        sp_diffs[i].getShiftCoeffs(), phys_unknowns, stages);
                }
                updateVisc(&mesh, svw_m, pde, phys_unknowns, stages, 0); 
            }


            // for(int i = 0; i < npatches; i++)
            // {
            //     auto v = patches[i]->getFlowPtr();
            //     if(t == 0.0)
            //     { 
            //         filters[i].filter(v, filt_unknowns);                      
            //     }   
            //     else
            //     {
            //         filters[i + npatches].filter(v, filt_unknowns);  
            //     }   
            // }
            for(int i = 0; i < npatches; i++)
            {
                auto v = patches[i]->getFlowPtr();
                if(t == 0.0)
                { 
                    filters[i].filter(v, filt_unknowns, &fft_data, fft_locs[i]);                      
                }   
                else
                {
                    filters[i + npatches].filter(v, filt_unknowns, &fft_data,
                        fft_locs[i]);  
                }   
            }



            // Determination of the adaptive timestep
            if(adaptive)
            {
                dt = pde.getAdaptiveTimeStep(mesh, phys_unknowns, stages, CFL);
            }

            // std::cout << "dt = " << dt << std::endl;
            
            if(t + dt > T)
            {
                dt = T - t;
                t = T; 
            }
            else
            {
                t += dt;
            }           

            // Advancing
            for(int i = 0; i < npatches; i++)
            {
                auto v = patches[i]->getFlowPtr();
                // ts.advance(v, sp_diffs[i], pde, dt, true);
                ts.advance(v, sp_diffs[i], pde, dt, true, fft_data, fft_locs[i]);  
            }
 
            mesh.setPatches(patches);
            mesh.setIntraPatchBC(phys_unknowns);
            pde.getBC()(&mesh, t);
            
        }

        patches = mesh.getPatches();
        for(int i = 0; i < npatches; i++)
        {
            patches[i]->VectorFieldToNodes();
        }
        mesh.setPatches(patches);


        auto t2 = std::chrono::high_resolution_clock::now();

        // Getting number of milliseconds as an integer. //
        auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1); 
        // Getting number of milliseconds as a double. //
        std::chrono::duration<double, std::milli> ms_double = t2 - t1;

        std::cout << ms_int.count() << "ms\n" << std::endl;
        std::cout << ms_double.count() << "ms" << std::endl;  

    }  


    
};



#endif 