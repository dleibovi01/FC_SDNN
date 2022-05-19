/* PDE solver */

#ifndef SOLVER_H
#define SOLVER_H

#include "SpatDiffScheme.h"
#include "TimeStepScheme.h"
#include "Mesh.h"
#include "Mesh1D.h"
#include "PDE.h"
#include "printing.h"
#include "SDNN.h"
#include <chrono>
#include <array>
#include <string>
#include <omp.h>
#include "FC_2D.h"

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

    Solver(const Meshtype &_mesh, const PDE &_pde, TS &_ts, 
        const std::vector<Sp_diff> &_sp_diff) : mesh{_mesh}, pde{_pde}, ts{_ts},
        sp_diffs{_sp_diff} {}

    // const Meshtype& getMesh() const {return mesh;}
    const Meshtype& getMesh() const {return mesh;}

    // void solve(double dt, double CFL, bool adaptive)
    // {
    //     std::cout << "Solver" << std::endl;
    //     double t = 0.0;
    //     double T = pde.getT();
    //     auto patches = mesh.getPatches();
    //     int npatches = sp_diffs.size();
    //     int phys_unknowns = pde.getPhysUnknowns();
    //     int stages = ts.getStages();
    //     pde.getIC()(&mesh);

    //     auto t1 = std::chrono::high_resolution_clock::now();
    


    //     int ndata = 11;
    //     std::vector<VectorField> data_patch;
    //     std::vector<std::vector<VectorField> > data;
    //     for(int i = 0; i < npatches; i++)
    //     {
    //         for(int j = 0; j < ndata; j++)
    //         {
    //             data_patch.push_back(VectorField(phys_unknowns, 
    //                 patches[i]->getNnodes()));
    //         } 
    //         data.push_back(data_patch);
    //         data_patch.empty();
    //     }


    //     std::cout << " entering solver loop" << std::endl;
    //     // Time loop
    //     while (t < T)
    //     {
    //         patches = mesh.getPatches();

    //         // Determination of the adaptive timestep
    //         // if(adaptive)
    //         // {
    //         //     dt = pde.getAdaptiveTimeStep(mesh, phys_unknowns, stages, CFL);
    //         // }

    //         // std::cout << "dt = " << dt << std::endl;
            
    //         if(t + dt > T)
    //         {
    //             dt = T - t;
    //             t = T; 
    //         }
    //         else
    //         {
    //             t += dt;
    //         }      

    

    //         // Advancing
    //         for(int i = 0; i < npatches; i++)
    //         {
    //             ts.advance(patches[i], sp_diffs[i], pde, dt, t, &(data[i]));
    //         }

    //         mesh.setPatches(patches);
    //         mesh.setIntraPatchBC(phys_unknowns, 0);
    //         pde.getBC()(&mesh, t);
            
    //     }

    //     patches = mesh.getPatches();
    //     for(int i = 0; i < npatches; i++)
    //     {
    //         patches[i]->VectorFieldToNodes();
    //     }
    //     mesh.setPatches(patches);


    //     auto t2 = std::chrono::high_resolution_clock::now();

    //     // Getting number of milliseconds as an integer. //
    //     auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1); 
    //     // Getting number of milliseconds as a double. //
    //     std::chrono::duration<double, std::milli> ms_double = t2 - t1;

    //     std::cout << ms_int.count() << "ms\n" << std::endl;
    //     std::cout << ms_double.count() << "ms" << std::endl;  

    // }  






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
        int threads = 1;

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
                    sp_diffs[i]->getC()]);
                loc.push_back(i*phys_unknowns + j);
            }
            fft_locs.push_back(loc);
            loc.clear();
        }


        // Creating containers to overload advance_sdnn
        std::vector<double* > mux;
        for(int i = 0; i < npatches; i++)
        {
            mux.push_back(new double[patches[i]->getNnodes()]);
        }


        auto t1 = std::chrono::high_resolution_clock::now();
      
        std::vector<double> bc_l;
        std::vector<double> bc_r;
        double h = 0.0;
        
        // Time loop
        while (t < T)
        {
            patches = mesh.getPatches();
            // Smooth viscosity assignment
            if(visc)
            {
                // #pragma omp parallel for num_threads(threads)
                // #pragma omp parallel for
                for(int i = 0; i < npatches; i++)
                {
                    pde.getSDNN().updateTau(patches[i], sp_diffs[i], pde,
                        phys_unknowns, stages);
                }
                pde.getSDNN().updateVisc(&mesh, svw_m, pde, phys_unknowns,
                    stages, 0); 
            }

            // #pragma omp parallel for num_threads(threads)
            // #pragma omp parallel for

            // Print_VectorField(patches[0]->getFlow(), true);
            
            bc_l = pde.getBC().getBC_L(t);
            bc_r = pde.getBC().getBC_R(t);

            // std::cout << std::endl;

            for(int i = 0; i < npatches; i++)
            {
                h = patches[i]->getH();
                if(t == 0.0)
                { 
                    filters[i]->filter(patches[i]->getFlowPtr(), filt_unknowns,
                        &fft_data, fft_locs[i], h, bc_l, bc_r);                      
                }   
                else
                {
                    filters[i + npatches]->filter(patches[i]->getFlowPtr(),
                        filt_unknowns, &fft_data, fft_locs[i], h, bc_l, bc_r);  
                }   
            }

            mesh.setIntraPatchBC(phys_unknowns, 0);
            // pde.getBC()(&mesh, t);

            // Determination of the adaptive timestep
            if(adaptive)
            {
                dt = pde.getAdaptiveTimeStep(mesh, phys_unknowns, stages, CFL);
            }
            
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

            // #pragma omp parallel for num_threads(threads)
            // #pragma omp parallel for
            for(int i = 0; i < npatches; i++)
            {
                // auto v = patches[i]->getFlowPtr();
                ts.advance_sdnn(patches[i], sp_diffs[i], pde, dt, 
                    t, true, fft_data, fft_locs[i]);  
                // ts.advance_sdnn(v, sp_diffs[i], pde, dt, &mux);  
                // ts.advance_sdnn(&mesh, sp_diffs, pde, dt, &mux); 
            } 
            // ts.advance_sdnn(&mesh, sp_diffs, pde, dt, &mux); 
            // Print_VectorField(patches[0]->getFlow(), true, 17);
            // for(int i = 0; i < npatches; i++)
            // {
            //     patches[i]->VectorFieldToNodes();
            // }


            mesh.setPatches(patches);
            // mesh.setIntraPatchBC(phys_unknowns);
            pde.getBC()(&mesh, t);
            // Print_VectorField(mesh.getPatches()[0]->getFlow(), true, 7);
            // // Print_Mesh1D(mesh);

            std::cout << std::endl;
            
        }

        patches = mesh.getPatches();
        for(int i = 0; i < npatches; i++)
        {
            patches[i]->VectorFieldToNodes();
        }
        mesh.setPatches(patches);

        while(!mux.empty())
        {
            delete[] mux.back();
            mux.pop_back();
        }

        auto t2 = std::chrono::high_resolution_clock::now();

        // Getting number of milliseconds as an integer. //
        auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1); 
        // Getting number of milliseconds as a double. //
        std::chrono::duration<double, std::milli> ms_double = t2 - t1;

        std::cout << ms_int.count() << "ms\n" << std::endl;
        std::cout << ms_double.count() << "ms" << std::endl;  


    }  

    void solve2DSDNN(double dt, double CFL, bool adaptive, bool visc)
    {
        std::cout << "2D Solver" << std::endl;
        double t = 0.0;
        double T = pde.getT();
        auto patches = mesh.getPatches();
        int npatches = patches.size();
        int phys_unknowns = pde.getPhysUnknowns();
        int stages = ts.getStages();
        pde.getIC()(&mesh);
        int threads = 1;

        // Filter data      
        std::vector<int> filt_unknowns;
        for(int i = 0; i < phys_unknowns; i++)
        {
            filt_unknowns.push_back(i);
        }

        // Set viscosity normalization coefficients
        std::cout << "Setting SVW" << std::endl;
        SVW_mesh svw_m{mesh};


        // Creating containers to overload advance_sdnn
        std::vector<double* > mux;
        std::vector<double* > muy;
        std::cout << "Setting containers for mux and muy" << std::endl;
        for(int i = 0; i < npatches; i++)
        {
            std::cout << "i = " << i << std::endl;
            mux.push_back(new double[patches[i]->getNnodes()]);
            muy.push_back(new double[patches[i]->getNnodes()]);
        }


        auto t1 = std::chrono::high_resolution_clock::now();
      
        std::vector<std::vector<double* > > bc_l(npatches);
        std::vector<std::vector<double* > > bc_r(npatches);
        std::vector<std::vector<double* > > bc_d(npatches);
        std::vector<std::vector<double* > > bc_u(npatches);  
        std::cout << "Setting containers for bcs" << std::endl;    
        for(int i = 0; i < npatches; i++)
        {
            for(int j = 0; j < phys_unknowns; j++)
            {
                bc_l[i].push_back(new double[patches[i]->getNy()]);
                bc_r[i].push_back(new double[patches[i]->getNy()]);
                bc_d[i].push_back(new double[patches[i]->getNx()]);
                bc_u[i].push_back(new double[patches[i]->getNx()]);
            }
        }  

        std::cout << "Starting loop" << std::endl;
        double hx = 0.0;
        double hy = 0.0;
        
        // Time loop
        while (t < T)
        {
            // std::cout << " t = " << t << std::endl;
            patches = mesh.getPatches();
            // Smooth viscosity assignment
            if(visc)
            {
                for(int i = 0; i < npatches; i++)
                {
                    pde.getSDNN().updateTau(patches[i], sp_diffs[i], pde,
                        phys_unknowns, stages);
                }
                pde.getSDNN().updateVisc(&mesh, svw_m, pde, phys_unknowns,
                    stages, 0); 
            }       


            for(int i = 0; i < npatches; i++)
            {
                hx = patches[i]->getHx();
                hy = patches[i]->getHy();               
                pde.getBC().getBC_L(&bc_l[i], patches[i], t);
                pde.getBC().getBC_R(&bc_r[i], patches[i], t);
                pde.getBC().getBC_D(&bc_d[i], patches[i], t);
                pde.getBC().getBC_U(&bc_u[i], patches[i], t);

                if(t == 0.0)
                { 
                    filters[i]->filter(patches[i]->getFlowPtr(), filt_unknowns,
                        hx, hy, bc_l[i], bc_r[i], bc_d[i], bc_u[i]);                      
                }   
                else
                {
                    filters[i + npatches]->filter(patches[i]->getFlowPtr(),
                        filt_unknowns, hx, hy, bc_l[i], bc_r[i], bc_d[i],
                        bc_u[i]);  
                }   
            }

            mesh.setIntraPatchBC(phys_unknowns, 0);

            // Determination of the adaptive timestep
            if(adaptive)
            {
                dt = pde.getAdaptiveTimeStep(mesh, phys_unknowns, stages, CFL);
            }
            
            if(t + dt > T)
            {
                dt = T - t;
                t = T; 
            }
            else
            {
                t += dt;
            }         

            // std::cout << "dt = " << dt << std::endl;  

            // Advancing
            for(int i = 0; i < npatches; i++)
            {
                ts.advance2DSDNN(patches[i], sp_diffs[i], pde, dt, 
                    t, mux[i], muy[i], bc_l[i], bc_r[i], bc_d[i], bc_u[i]);  
            } 

            // for printing
            // patches[0]->VectorFieldToNodes();

            mesh.setPatches(patches);
            // mesh.setIntraPatchBC(phys_unknowns);
            pde.getBC()(&mesh, t);

            // Print_Mesh2DUniform(mesh, phys_unknowns);


            // std::cout << std::endl;
            
        }

        patches = mesh.getPatches();
        for(int i = 0; i < npatches; i++)
        {
            patches[i]->VectorFieldToNodes();
        }
        mesh.setPatches(patches);

        while(!mux.empty())
        {
            delete[] mux.back();
            mux.pop_back();
            delete[] muy.back();
            muy.pop_back();
        }

        auto t2 = std::chrono::high_resolution_clock::now();

        // Getting number of milliseconds as an integer. //
        auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1); 
        // Getting number of milliseconds as a double. //
        std::chrono::duration<double, std::milli> ms_double = t2 - t1;

        std::cout << ms_int.count() << "ms\n" << std::endl;
        std::cout << ms_double.count() << "ms" << std::endl;  
    }  














//    void solve_sdnn2(double dt, double CFL, bool adaptive, bool visc)
//     {
//         std::string filename = "viscosity.txt";
//         std::cout << "Solver" << std::endl;
//         double t = 0.0;
//         double T = pde.getT();
//         auto patches = mesh.getPatchesPtr();
//         int npatches = sp_diffs.size();
//         int phys_unknowns = pde.getPhysUnknowns();
//         int stages = ts.getStages();
//         pde.getIC()(&mesh);

//         // Filter data      
//         std::vector<int> filt_unknowns;
//         for(int i = 0; i < phys_unknowns; i++)
//         {
//             filt_unknowns.push_back(i);
//         }

//         // Set viscosity normalization coefficients
//         SVW_mesh svw_m{mesh};

//         // Pre-allocate a container to save the FFTS
//         std::vector<std::complex<double> * > fft_data;
//         std::vector<std::vector<int> > fft_locs;
//         std::vector<int> loc;
//         for(int i = 0; i < npatches; i++)
//         {
//             for(int j = 0; j < phys_unknowns; j++)
//             {
//                 fft_data.push_back(
//                     new std::complex<double> [mesh.getPatches()[i]->getNnodes() + 
//                     sp_diffs[i].getC()]);
//                 loc.push_back(i*phys_unknowns + j);
//             }
//             fft_locs.push_back(loc);
//             loc.clear();
//         }


//         // Creating containers to overload advance_sdnn
//         std::vector<double* > mux;
//         for(int i = 0; i < npatches; i++)
//         {
//             mux.push_back(new double[patches->at(i)->getNnodes()]);
//         }

//         std::vector<int> viscosity;
//         viscosity.push_back(stages*phys_unknowns + 1);
//         auto t1 = std::chrono::high_resolution_clock::now();
      
        
//         // Time loop
//         while (t < T)
//         {
//             patches = mesh.getPatchesPtr();
//             // Smooth viscosity assignment
//             if(visc)
//             {
//                 for(int i = 0; i < npatches; i++)
//                 {
//                     updateTau(patches->at(i), sp_diffs[i], pde, 
//                         sp_diffs[i].getShiftCoeffs(), phys_unknowns, stages);
//                 }
//                 updateVisc(&mesh, svw_m, pde, phys_unknowns, stages, 0); 
//             }

//             // Exchange viscosity at boundaries
//             mesh.setIntraPatchBC(1, phys_unknowns*stages + 1);

//             for(int i = 0; i < npatches; i++)
//             {
//                 auto v = patches->at(i)->getFlowPtr();
//                 if(t == 0.0)
//                 { 
//                     filters[i].filter(v, filt_unknowns, &fft_data, fft_locs[i]);                      
//                 }   
//                 // else
//                 {
//                     filters[i + npatches].filter(v, filt_unknowns, &fft_data,
//                         fft_locs[i]);  
//                 }   
//             }
//             // mesh.setPatches(patches);
//             mesh.setIntraPatchBC(phys_unknowns, 0);


//             pde.getBC()(&mesh, t);

//             // Determination of the adaptive timestep
//             if(adaptive)
//             {
//                 dt = pde.getAdaptiveTimeStep(mesh, phys_unknowns, stages, CFL);
//             }

//             // std::cout << "dt = " << dt << std::endl;
            
//             if(t + dt > T)
//             {
//                 dt = T - t;
//                 t = T; 
//             }
//             else
//             {
//                 t += dt;
//             }   


//             // ts.advance_sdnn(&mesh, sp_diffs, pde, dt, &mux, fft_data, fft_locs,
//             //     filters);
//             ts.advance_sdnn(&mesh, sp_diffs, pde, dt, t, &mux, fft_data, fft_locs);
//             // ts.advance_sdnn(&mesh, sp_diffs, pde, dt, &mux);
 
//             // std::cout << " patch 0" << std::endl;
//             // Print_VectorField(mesh.getPatches()[0]->getFlow().extract(viscosity), true);
//             // std::cout << std::endl;
//             // std::cout << std::endl;
//             // std::cout << " patch 1" << std::endl;
//             // Print_VectorField(mesh.getPatches()[1]->getFlow().extract(viscosity), true);
//             // std::cout << std::endl;
//             // mesh.setPatches(patches);
//             // mesh.setIntraPatchBC(phys_unknowns);
//             // pde.getBC()(&mesh, t);
//             // if(t > 0.13)
//             // {
//             //     Print_Mesh1D(mesh, phys_unknowns, 6, stages*phys_unknowns + 1,
//             //         filename);
//             //     std::cout << "t = " << t << std::endl;
//             // }

                    
//         }

//         patches = mesh.getPatchesPtr();
//         for(int i = 0; i < npatches; i++)
//         {
//             patches->at(i)->VectorFieldToNodes();
//         }
//         mesh.setPatches(*patches);

//         while(!mux.empty())
//         {
//             delete[] mux.back();
//             mux.pop_back();
//         }

//         auto t2 = std::chrono::high_resolution_clock::now();

//         // Getting number of milliseconds as an integer. //
//         auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1); 
//         // Getting number of milliseconds as a double. //
//         std::chrono::duration<double, std::milli> ms_double = t2 - t1;

//         std::cout << ms_int.count() << "ms\n" << std::endl;
//         std::cout << ms_double.count() << "ms" << std::endl;  

//         Print_Mesh1D(mesh, phys_unknowns, 6, stages*phys_unknowns + 1,
//             filename);

//     }  



    
};


void initialize1DFCSDNN(std::vector<FC_1D* > * diff_schemes, 
    std::vector<FC_1D* > * filters, const Mesh1DUniform &mesh, int npatches,
    std::string xx, int N, int d, int C, double delta, double alpha0,
    double p_0, double p);

void initialize2DFCSDNN(std::vector<FC_2D* > * diff_schemes, 
    std::vector<FC_2D* > * filters, const Mesh2DUniform &mesh, 
    int npatches_x, int npatches_y, std::string left, std::string right,
    std::string down, std::string up, int Nx, int Ny, int dx, int dy, int Cx,
    int Cy, double delta, double alpha0, double p_0, double p);

template<typename sp_diff>
void freeFCSDNN(std::vector<sp_diff* > * diff_schemes, 
    std::vector<sp_diff* > * filters)
{
    while(!diff_schemes->empty())
    {
        delete diff_schemes->back();
        diff_schemes->pop_back();
    }
    while(!filters->empty())
    {
        delete filters->back();
        filters->pop_back();
    }
}



#endif 