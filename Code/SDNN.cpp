/* Definition of SDNN routines. */

#include "SDNN.h"

SDNN::SDNN(const ANN &_ann, double _discard_noise, double _w1, double _w2,
    double _w3, double _w4) : discard_noise{_discard_noise},
    w1{_w1}, w2{_w2}, w3{_w3}, w4{_w4} 
{
    ann = new ANN(_ann.getAlpha());
}

SDNN::SDNN(double _discard_noise, double _w1, double _w2, double _w3,
    double _w4, double alpha) : discard_noise{_discard_noise}, w1{_w1}, w2{_w2},
    w3{_w3}, w4{_w4} 
{
    ann = new ANN(alpha);
}

SDNN::SDNN(const SDNN &sdnn) : discard_noise{sdnn.discard_noise}, w1{sdnn.w1},
    w2{sdnn.w2}, w3{sdnn.w3}, w4{sdnn.w4} 
{
    ann = new ANN(sdnn.ann->getAlpha());
}

SDNN & SDNN::operator= (const SDNN &sdnn)
{
    if(this == &sdnn)
    {
        return *this;
    }
    else
    {
        delete ann;  
        ann = new ANN(sdnn.ann->getAlpha()); 
        return *this;                    
    }   
}

void SDNN::getMaxedMWSB(Patch1D * patch, int s, const double * MWSB,
    double * MWSB_maxed) const
{
    double max_0;
    double max_end;
    int N = patch->getNnodes();
    max_0 = *(std::max_element(MWSB, MWSB + s + 1));
    max_end = *(std::max_element(MWSB + N - s, MWSB + N));
    MWSB_maxed[0] = max_0;
    MWSB_maxed[1] = max_0;
    MWSB_maxed[2] = max_0;
    MWSB_maxed[N - 3] = max_end;
    MWSB_maxed[N - 2] = max_end;
    MWSB_maxed[N - 1] = max_end;

    for(int i = 3; i < N - 3; i++)
    {
        MWSB_maxed[i] = *(std::max_element(MWSB + i - 3, MWSB + i + 4));
    }
}

void SDNN::getMaxedMWSB(Patch2D * patch, int s, const double * MWSB,
    double * MWSB_maxed) const
{
    int Nx = patch->getNx();
    int Ny = patch->getNy();
    double locStencil[s*s];
    for(int j = 0; j < Nx; j++)
    {
        for(int i = 0; i < Ny; i++)
        {
            form2DLocStencil(MWSB, locStencil, Nx, Ny, s, i, j);
            MWSB_maxed[j*Ny + i] = 
                *(std::max_element(locStencil, locStencil + s*s));
        }
    }
}

void SDNN::form2DLocStencil(const double * y, double * stencil, int Nx, int Ny,
    int s, int i, int j) const
{
    int i0;
    int j0;

    if(i < (s - 1)/2)
        i0 = 0;
    else if(i >= Ny - (s - 1)/2)
        i0 = Ny - s;
    else
        i0 = i - (s - 1)/2;

    if(j < (s - 1)/2)
        j0 = 0;
    else if(j >= Nx - (s - 1)/2)
        j0 = Nx - s;
    else
        j0 = j - (s - 1)/2;

    for(int l = 0; l < s; l++)
    {
        for(int k = 0; k < s; k++)
        {
            stencil[l*s + k] = y[(j0 + l)*Ny + i0 + k];
        }
    }
}

std::vector<double> SDNN::Visc_weights(int N, const double* tau) const
{
    std::vector<double> w(N, 0.0);
    #pragma omp simd
    for(int i = 0; i < N; i++)
    {
        w[i] = Q(tau[i]);
    }
    return w;
}

void SDNN::form_stencils(int N, int C, int s, const double* proxy_ext,
    double* stencils) const
{
    #pragma omp simd
    for(int i = 0; i < N; i++)
    {
        #pragma omp simd
        for (int j = 0; j < s; j++)
        {
            stencils[s*i+ j] = proxy_ext[(N + C - 3 + i + j) % (N + C)];
        }
    } 
}

void SDNN::preprocess_stencils(int N, double *stencils, bool* discard, int* tau)
    const
{
    double slope;
    double s0;
    double M = 0.0;
    double m = 0.0;

    #pragma omp simd
    // #pragma omp simd collapse(2)
    for(int i = 0; i < N; i++)
    {
        slope = double(stencils[(i + 1)*s - 1] - stencils[i*s]) / double(s - 1);
        s0 = stencils[i*s];
        #pragma unroll
        for(int j = 0; j < s; j++)
            stencils[i*s + j] = stencils[i*s + j] - (s0 + slope*double(j));
    }   

    // #pragma omp simd lastprivate(M, m)
    for(int i = 0; i < N; i++)
    {
        M = *(std::max_element(stencils + i*s, stencils + (i + 1)*s));
        m = *(std::min_element(stencils + i*s, stencils + (i + 1)*s));
        if(M - m > discard_noise)
        {
            discard[i] = false;
            // #pragma omp simd
            #pragma unroll
            for (int j = 0; j < s; j++)
            {
                stencils[i*s + j] = (2.0*stencils[i*s + j] - M - m) / (M - m);
            }
        }
        else
        {
            discard[i] = true;
            tau[i] = 4;
        }
    }  
}

void SDNN::formRegStencils(int N, int s, const double * stencils, 
    const bool * discard, double * regStencils, int * indices) const
{
    int current_regindex = 0;
    int current_index = 0;
    for(int i = 0; i < N; i++)
    {
        if(!discard[i])
        {
            std::copy(stencils + i*s, stencils + (i+1)*s, 
                regStencils + current_regindex);
            current_regindex += s;
            indices[current_index] = i;
            current_index++;
        }
    }
}

double SDNN::Q(double tau) const
{
    if(tau == 1.0)
        return w1;
    else if(tau == 2.0)
        return w2;
    else if (tau == 3.0)
        return w3;
    else
        return w4;
}

void SDNN::getRegularity(int * tau, bool * discard, double * stencils, int N,
    int s) const
{
    int M = 16;
    int output = 4;

    static int incx = 1;
    static int incy = 1;
    double a = 1.;
    double b = 0.;    

    int N0 = VectorSum(N, discard);
    int indices[N0];
    double regStencils[s*N0];
    int tau_ind[N0];
    formRegStencils(N, s, stencils, discard, regStencils, indices);
    ann->fwdClassif(tau_ind, regStencils, N0);
    for(int i = 0; i < N0; i++)
        tau[indices[i]] = tau_ind[i];
}