/* SDNN functions */

#ifndef SDNN_H
#define SDNN_H

#include <cmath>
#include "PDE.h"
#include <string>
#include "mkl_operations.h"
#include "FC.h"
#include "VectorField1D.h"
#include "Mesh.h"
#include "SDNN_data.h"
#include "SVW.h"
#include "MVOperations.h"
#include "VectorOperations.h"

constexpr int s = 7;
constexpr double discard_noise = 0.01;
// constexpr CBLAS_LAYOUT Layout = CblasColMajor;
constexpr CBLAS_LAYOUT Layout = CblasRowMajor;
constexpr CBLAS_TRANSPOSE TRANS = CblasNoTrans;

void readMatrix(double* A, std::string filename_A)
{
    int i = 0;
    double data;
    std::ifstream Adata (filename_A.c_str());
    if(Adata.is_open())
    {
        while(Adata >> data)
        {
            A[i] = data;
            i = i + 1;       
        }
    }        
}

inline void elu(const int N, double* input, const double alpha)
{
    #pragma omp simd
    for(int i = 0; i < N; i++)
    {
        if(input[i] <= 0)
        {
            input[i] =  alpha*(std::exp(input[i]) - 1.0);
        } 
    }
}



void softmax(int N, double* input)
{
    double sum_exp = 0;
    for(int i = 0; i < N; i++)
    {
        input[i] = std::exp(input[i]);
        sum_exp +=  input[i]; 
    }
    for(int i = 0; i < N; i++)
    {
        input[i] = input[i] / sum_exp;
    }    
}


template<typename Node>
double visc_window(const Node & node1, const Node & node2, double cutoff)
{
    double distance = node1.getDist(node2);
    const double pi = std::acos(-1);
    double H;
    if (distance < cutoff)
    {
        H = std::cos(pi*distance / 2.0 / cutoff);
        return H*H;
    }
    else
    {
        return 0.0;
    }
}

double Q(double tau)
{
    if(tau == 1.0)
    {
        return 2.0;
    }
    else if(tau == 2.0)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    } 
}


void getMaxedMWSB(int N, int s, const double * MWSB, double * MWSB_maxed)
{
    double max_0;
    double max_end;
    max_0 = *(std::max_element(MWSB, MWSB + s + 1));
    max_end = *(std::max_element(MWSB + N - s, MWSB + N));
    MWSB_maxed[0] = max_0;
    MWSB_maxed[1] = max_0;
    MWSB_maxed[2] = max_0;
    MWSB_maxed[N - 3] = max_end;
    MWSB_maxed[N - 2] = max_end;
    MWSB_maxed[N - 1] = max_end;

    #pragma omp simd
    for(int i = 3; i < N - 3; i++)
    {
        MWSB_maxed[i] = *(std::max_element(MWSB + i - 3, MWSB + i + 4));
    }
}


std::vector<double> Visc_weights(int N, const double* tau)
{
    std::vector<double> w(N, 0.0);
    #pragma omp simd
    for(int i = 0; i < N; i++)
    {
        w[i] = Q(tau[i]);
    }
    return w;
}

void form_stencils(int N, int C, const double* proxy_ext, double* stencils)
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


template<typename Mesh>
void Lambda(Mesh * mesh, int phys_unknowns, int stages, double cutoff, 
    bool initialization)
{
    auto patches = mesh->getPatches();
    auto patch_origin = patches[0];
    auto patch_target = patches[0];
    auto Node_origin = patch_origin->getNode(0);
    auto Node_target = patch_target->getNode(0);
    int npatches = patches.size();
    int N_origin;
    int N_target;
    auto v_origin = patch_origin->getFlow();
    auto v_target = patch_target->getFlow();
    double *tau;
    double *target_field;
    
    // Reset to zero the value of lambda on all nodes, all patches   
    for(int i = 0; i < npatches; i++)
    {
        patch_origin = patches[i];
        N_origin = patch_origin->getNnodes();
        v_origin = patch_origin->getFlow();
        for(int i = 0; i < N_origin; i++)
        {
            v_origin.setFieldValue(stages*phys_unknowns + 1, i, 0.0);
        }
        patch_origin->setField(v_origin);
    }
    
    for(int i = 0; i < npatches; i++)
    {
        patch_origin = patches[i];
        N_origin = patch_origin->getNnodes();
        v_origin = patch_origin->getFlow();
        tau = v_origin.getField(stages*phys_unknowns);
        for(int j = 0; j < N_origin; j++)
        {
            Node_origin = patch_origin->getNode(j);
            for(int k = 0; k < npatches; k++)
            {
                patch_target = patches[k];
                N_target = patch_target->getNnodes();
                v_target = patch_target->getFlow();                
                for(int l = 0; l < N_target; l++)
                {
                    Node_target = patch_target->getNode(l);
                    if(initialization)
                    {
                        v_target.setFieldValue(stages*phys_unknowns + 2, l, 
                            v_target.getFieldValue(stages*phys_unknowns + 2, l) 
                            + visc_window(*Node_origin, *Node_target, cutoff));
                    }
                    else
                    {
                        v_target.setFieldValue(stages*phys_unknowns + 1, l, 
                            v_target.getFieldValue(stages*phys_unknowns + 1, l) 
                            + Q(tau[j])*visc_window(*Node_origin, *Node_target, 
                            cutoff) / 
                            v_target.getFieldValue(stages*phys_unknowns + 2, l));                        
                    }
                }            
                patch_target->setField(v_target);
            }
        }
    }   
}




template<typename Patch, typename Mesh, typename SVW, typename PDE>
void updateViscPatch(Patch *patch, Mesh *mesh, const SVW &svw, const PDE &pde, 
    int phys_unknowns, int stages, int stage)
{
  
    auto v0 = patch->getFlowRef();
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
        v0 = patches[patchIds[i]]->getFlowRef();
        N0 = v0.getLength();
        weighted_tau = Visc_weights(N0, v0.getField(phys_unknowns*stages));
        W[i]->MV(alpha, weighted_tau.data(), beta, mu);
    }
 
    // Need to form stencils and get maximums of MWSB
    int s = 7;

    auto v = patch->getFlowPtr();
    pde.getMWSB(*v, MWSB);
    vdAbs(N, MWSB, MWSB);
    double MWSB_maxed[N];
    getMaxedMWSB(N, s, MWSB, MWSB_maxed);

    double y[N];
    for(int i = 0; i < N; i++)
    {
        y[i] = 0.0;
    }
   // Muiltiply by the wave speed
    VectorMul(N, MWSB_maxed, mu, mu);
    // Muiltiply by h
    cblas_daxpy(N, h, mu, 1, y, 1); 
    v->setField(phys_unknowns*stages + 1, N, y);

}


template<typename Mesh, typename PDE>
void updateVisc(Mesh *mesh, const SVW_mesh &svw_mesh, const PDE &pde,
    int phys_unknowns, int stages, int stage)
{
    auto patches = mesh->getPatches();
    int npatches = patches.size();
    auto svws = svw_mesh.getSVWs();
    for(int i = 0; i < npatches; i++)
    {
        // patch = patches[i];
        updateViscPatch(patches[i], mesh, svws[i], pde, phys_unknowns, stages,
            stage);
    }
}



void preprocess_stencils(int N, double *stencils, bool* discard, int* tau)
{
    double slope;
    double s0;
    double M = 0.0;
    double m = 0.0;

    #pragma omp simd
    for(int i = 0; i < N; i++)
    {
        slope = double(stencils[(i + 1)*s - 1] - stencils[i*s]) / double(s - 1);
        s0 = stencils[i*s];
        #pragma omp simd
        for(int j = 0; j < s; j++)
        {
            stencils[i*s + j] = stencils[i*s + j] - (s0 + slope*double(j));
        }
    }   

    // #pragma omp simd lastprivate(M, m)
    for(int i = 0; i < N; i++)
    {
        M = *(std::max_element(stencils + i*s, stencils + (i + 1)*s));
        m = *(std::min_element(stencils + i*s, stencils + (i + 1)*s));
        if(M - m > discard_noise)
        {
            discard[i] = false;
            #pragma omp simd
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


template<typename Patch, typename SpatDiffScheme, typename PDE>
void updateTau (Patch * patch, const SpatDiffScheme &sp, const PDE &pde,
    std::complex<double>* shift_coeffs, int unknowns, int stages)
{
    MKL_LONG status;
    int fourPts;
    // auto v = patch->getFlow();
    auto v = patch->getFlowPtr();
    
    int N = v->getLength();
    double proxy[N];
    pde.getProxy(*v, proxy);
    int tau[N];
    bool discard[N];

    int d = sp.getD();
    int C = sp.getC();
    double* AQ = sp.getAQ();
    double* FAQF = sp.getFAQF();
    DFTI_DESCRIPTOR_HANDLE desc_handle = sp.getDescHandle();

    // Form the continuation
    fourPts = N + C;
    double proxy_ext[N+C];   
    Fcont_shift(proxy, proxy_ext, shift_coeffs, N, d, C, double(fourPts), AQ, 
        FAQF, desc_handle);

    // Form the stencils 
    double stencils[s*N];
    form_stencils(N, C, proxy_ext, stencils);

    // Pre-process the stencils
    preprocess_stencils(N, stencils, discard, tau);

    // get classification (forward propagation)
    pde.getANN().getRegularity(tau, discard, stencils, N, s);

    double tau_dbl[N];  
    #pragma omp simd 
    for(int i = 0; i < N; i++)
    {
        tau_dbl[i] = double(tau[i]);
    }
    // v.setField(unknowns*stages, N, tau_dbl);
    // patch->setField(v);
    v->setField(unknowns*stages, N, tau_dbl);

}




// template< int M1, int N1, int M2, int N2, int M3, int N3, int M4, int N4>
class ANN {

double *W1;
double *B1;
double *W2;
double *B2;
double *W3;
double *B3;
double *W4;
double *B4;
// std::array< std::array<double, N1> M1> W1;

const double alpha;

public :
    
    ANN(std::string W1_filename, std::string W2_filename, 
        std::string W3_filename, std::string W4_filename, 
        std::string B1_filename, std::string B2_filename, 
        std::string B3_filename, std::string B4_filename, double _alpha) : 
        alpha{_alpha}
    {
        set_NN_weights(W1_filename, B1_filename, W2_filename, B2_filename, 
            W3_filename, B3_filename, W4_filename, B4_filename);
    }

    ANN(double _alpha) : alpha{_alpha}
    {
        set_NN_weights();
    }

    void set_NN_weights(std::string W1_file, std::string B1_file, 
        std::string W2_file, std::string B2_file, std::string W3_file, 
        std::string B3_file, std::string W4_file, std::string B4_file)
    { 
        W1 = new double[7*16];
        B1 = new double[16];
        W2 = new double[16*16];
        B2 = new double[16];
        W3 = new double[16*16];
        B3 = new double[16];
        W4 = new double[4*16];
        B4 = new double[4*16];    
        readMatrix(W1, W1_file);  
        readMatrix(W2, W2_file); 
        readMatrix(W3, W3_file); 
        readMatrix(W4, W4_file); 
        readMatrix(B1, B1_file); 
        readMatrix(B2, B2_file); 
        readMatrix(B3, B3_file); 
        readMatrix(B4, B4_file);                 
    }  

    void set_NN_weights()
    {
        W1 = new double[7*16];
        B1 = new double[16];
        W2 = new double[16*16];
        B2 = new double[16];
        W3 = new double[16*16];
        B3 = new double[16];
        W4 = new double[4*16];
        B4 = new double[4];   
        std::copy(W1_data.data(), W1_data.data() + 112, W1);
        std::copy(W2_data.data(), W2_data.data() + 256, W2);
        std::copy(W3_data.data(), W3_data.data() + 256, W3);
        std::copy(W4_data.data(), W4_data.data() + 64, W4);
        std::copy(B1_data.data(), B1_data.data() + 16, B1);
        std::copy(B2_data.data(), B2_data.data() + 16, B2);
        std::copy(B3_data.data(), B3_data.data() + 16, B3);
        std::copy(B4_data.data(), B4_data.data() + 4, B4);
    } 


    ANN(const ANN & ann) :  alpha{ann.alpha}
    {
        // deep copy of W1, B1, etc ...
        set_NN_weights();        
    }


    ANN & operator= (const ANN &ann)
    {
        if(this == &ann)
        {
            return *this;
        }
        else
        {
            delete[] W1;
            delete[] W2;
            delete[] W3;
            delete[] W4;
            delete[] B1;
            delete[] B2;
            delete[] B3;
            delete[] B4;   
            set_NN_weights();  
            return *this;                    
        }   
    }


    virtual ~ANN() {free_mem();}

    void free_mem()
    {
        delete[] W1;
        delete[] W2;
        delete[] W3;
        delete[] W4;
        delete[] B1;
        delete[] B2;
        delete[] B3;
        delete[] B4;       
    }


    void getRegularity(int * tau, bool * discard, double * stencils, int N,
        int s) const
    {
        int M = 16;
        int output = 4;
        // double res[M];
        // double res_end[output];


        static int incx = 1;
        static int incy = 1;
        double a = 1.;
        double b = 0.;    
        int lda;


        double input_1[M];
        double input_2[M];
        double input_3[M];
        double output_4[4];
        // for(int i = 0; i < N; i++)
        // {
        //     if(discard[i] == false)
        //     {
        //         std::copy(B1, B1 + M, input_1);
        //         std::copy(B2, B2 + M, input_2);
        //         std::copy(B3, B3 + M, input_3);
        //         std::copy(B4, B4 + 4, output_4);

        //         // lda = M;          
        //         lda = s; 
        //         cblas_dgemv (Layout, TRANS, M, s, a, W1, lda, stencils + i*s, 
        //             incx, 1.0, input_1, incy);   
        //         // MVMult_rowdom(M, s, W1, stencils + i*s, input_1);
        //         // VectorAdd(M, input_1, B1, input_1);     
        //         elu(M, input_1, alpha);

        //         lda = M; 
        //         cblas_dgemv (Layout, TRANS, M, M, a, W2, lda, input_1, incx, 1.0, 
        //             input_2, incy);
        //         // MVMult_rowdom(M, M, W2, input_1, input_2);
        //         // VectorAdd(M, input_2, B2, input_2);  
        //         elu(M, input_2, alpha); 

        //         lda = M;
        //         cblas_dgemv (Layout, TRANS, M, M, a, W3, lda, input_2, incx, 1.0, 
        //             input_3, incy);
        //         // MVMult_rowdom(M, M, W3, input_2, input_3);
        //         // VectorAdd(M, input_3, B3, input_3); 
        //         elu(M, input_3, alpha);           

        //         // lda = output;
        //         lda = M;
        //         cblas_dgemv (Layout, TRANS, output, M, a, W4, lda, input_3, 
        //             incx, 1.0, output_4, incy); 
        //         // MVMult_rowdom(output, M, W4, input_3, output_4);
        //         // VectorAdd(output, output_4, B4, output_4); 


        //         tau[i] = std::distance(output_4, 
        //             std::max_element(output_4, output_4 + output )) + 1;

        //     }
        // }

        int N0 = VectorSum(N, discard);
        int indices[N0];
        double regStencils[s*N0];
        formRegStencils(N, s, stencils, discard, regStencils, indices);
        if(N0 > 0)
        {
            // std::cout <<"indices[0] = " << indices[0] << std::endl;
            double input_1[M*N0];
            double input_2[M*N0];
            double input_3[M*N0];
            double output_4[output*N0];
            double ones[N0];

            #pragma omp simd
            for(int i = 0; i < M*N0; i++)
            {
                input_1[i] = 0.0;
                input_2[i] = 0.0;
                input_3[i] = 0.0;
            }
            #pragma omp simd
            for(int i = 0; i < output*N0; i++)
            {
                output_4[i] = 0.0;
            }
            #pragma omp simd
            for(int i = 0; i < N0; i++)
            {
                ones[i] = 1.0;
            }


            cblas_dger (CblasColMajor, M, N0, 1.0, B1, 1, ones, 1, input_1, M);
            cblas_dger (CblasColMajor, M, N0, 1.0, B2, 1, ones, 1, input_2, M);
            cblas_dger (CblasColMajor, M, N0, 1.0, B3, 1, ones, 1, input_3, M);
            cblas_dger (CblasColMajor, output, N0, 1.0, B4, 1, ones, 1,
                output_4, output);

            cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, M, N0, s,
                1.0, W1, M, regStencils, s, 1.0, input_1, M);
            elu(M*N0, input_1, alpha);

            cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, M, N0, M, 1.0,
                W2, M, input_1, M, 1.0, input_2, M);
            elu(M*N0, input_2, alpha);

            cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, M, N0, M, 1.0,
                W3, M, input_2, M, 1.0, input_3, M);
            elu(M*N0, input_3, alpha);

            cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, output, N0, M,
                1.0, W4, output, input_3, M, 1.0, output_4, output);

            #pragma omp simd
            for(int i = 0; i < N0; i++)
            {
                tau[indices[i]] = std::distance(output_4 + output*i, 
                    std::max_element(output_4 + output*i,
                        output_4 + output*(i + 1))) + 1;
            }
        }






    }

private :

    void formRegStencils(int N, int s, const double * stencils, 
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



};







#endif