/* 1D FC spatial differentiation schemes */

// #include "mkl_types.h"
#include <complex>
#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#include "FC_1D.h"


// FC_1D definitions

FC_1D::FC_1D(int _N, int _d, int _C, double _h, double delta)
{
    init(_N, _d, _C, _h);
    int k[fourPts];
    getK(k, fourPts);
    const double pi = std::acos(-1);
    const std::complex<double> I(0, 1);     
    for(int i = 0; i < fourPts; i++)
    {
        shift_coeffs[i] = std::exp(2.0*pi*I* double(k[i]) * delta / 
        double(fourPts));
    }
}

FC_1D::FC_1D(const FC_1D &slv)
{
    // patch = slv.getPatch();
    N = slv.getN();
    d = slv.getD();
    C = slv.getC();
    // h = slv.getH();
    h = slv.h;
    fourPts = slv.getFourPts();
    fourPts_dbl = slv.fourPts_dbl;
    prd = slv.getPrd();
    desc_handle = slv.getDescHandle();
    der_coeffs = new double [N+C];
    der_coeffs_2 = new double [N+C];
    filter_coeffs = new double [N+C]; 
    shift_coeffs = new std::complex<double> [N+C];
    AQ = new double [C*d];
    FAQF = new double [C*d];
    double * d_c = slv.getDerCoeffs();
    std::copy(d_c, d_c + N + C, der_coeffs);
    double * d_c_2 = slv.getDerCoeffs2();
    std::copy(d_c_2, d_c_2 + N + C, der_coeffs_2);
    double * f_c = slv.getFilterCoeffs();
    std::copy(f_c, f_c + N + C, filter_coeffs);  
    std::complex<double> * s_c = slv.getShiftCoeffs();
    std::copy(s_c, s_c + N + C, shift_coeffs);  
    double * aq = slv.getAQ();
    std::copy(aq, aq + C*d, AQ);    
    double * faqf = slv.getFAQF();
    std::copy(faqf, faqf + C*d, FAQF);      
}

FC_1D::~FC_1D()
{
    delete [] AQ;
    delete [] FAQF;
    delete [] filter_coeffs;
    delete [] der_coeffs;
    delete [] der_coeffs_2;
    delete [] shift_coeffs;
}

void FC_1D::init(int _N, int _d, int _C, double _h)
{
    N = _N;
    d = _d;
    C = _C;
    fourPts = N + C;
    fourPts_dbl = double(fourPts);
    h = _h;
    prd = fourPts_dbl*h;
    der_coeffs = new double [N+C];
    der_coeffs_2 = new double [N+C];
    filter_coeffs = new double [N+C];
    shift_coeffs = new std::complex<double> [N+C];
    double A[C*d];
    double Q[d*d];   
    AQ = new double [C*d];
    FAQF = new double [C*d];
    MKL_LONG status;
    const double pi = std::acos(-1);
    const std::complex<double> I(0, 1); 
    int k[N + C];
    getK(k, N + C);
    for(int j = 0; j < N + C; j++)
    {
        der_coeffs[j] = 2.0*pi/prd*double(k[j]);
        der_coeffs_2[j] = -4.0*pi*pi/prd/prd*double(k[j])*double(k[j]);
    }      
    for(int j = 0; j < N + C; j++)
    {
        filter_coeffs[j] = 1.0;
    }    
    status = DftiCreateDescriptor(&desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, 1,
        fourPts); 
    status = DftiCommitDescriptor(desc_handle);    
}

void FC_1D::filter(VectorField *v, const std::vector<int> &unknowns, 
    std::vector<std::complex<double> *> *ffts, 
    const std::vector<int> &fft_loc, double h, 
    const std::vector<double> &bc_l, const std::vector<double> &bc_r) const
{
    for(int i = 0; i < unknowns.size(); i++)
    {
        filter(v->getField(unknowns[i]), ffts->at(fft_loc[i]), h, bc_l[i], 
            bc_r[i]);   
    } 
}



////////////////////////////////////////////////////////////////////////////////
// FC_1D_DD definitions

FC_1D_DD::FC_1D_DD(int _N, int _d, int _C, double _h) : 
    FC_1D{_N, _d, _C, _h}
{
    double A[C*d];
    double Q[d*d]; 
    set_FC_Data(A, Q, d, C);
    build_Cont_Mat(A, Q, d, C, AQ, FAQF);  
}

FC_1D_DD::FC_1D_DD(int _N, int _d, int _C, double _h, int alpha, int p) : 
    FC_1D{_N, _d, _C, _h}
{
    double A[C*d];
    double Q[d*d]; 
    set_FC_Data(A, Q, d, C);
    build_Cont_Mat(A, Q, d, C, AQ, FAQF);  
    getFiltCoeffs(filter_coeffs, fourPts, alpha, p); 
}

FC_1D_DD::FC_1D_DD(int _N, int _d, int _C, double _h, double delta) : 
    FC_1D{_N, _d, _C, _h}
{
    double A[C*d];
    double Q[d*d]; 
    set_FC_Data(A, Q, d, C);
    build_Cont_Mat(A, Q, d, C, AQ, FAQF); 
    int k[fourPts];
    getK(k, fourPts);
    const double pi = std::acos(-1);
    const std::complex<double> I(0, 1);     
    for(int i = 0; i < fourPts; i++)
    {
        shift_coeffs[i] = std::exp(2.0*pi*I* double(k[i]) * delta / 
        double(fourPts));
    }
}

void FC_1D_DD::diff(const std::complex<double> *y_hat, double * y_der,
    double* y_der_2) const
{
    FC_Der(y_der, y_hat, der_coeffs, N, C, desc_handle);
    FC_Der(y_der_2, y_hat, der_coeffs_2, N, C, desc_handle, true);
}

void FC_1D_DD::diff(const double* y, double *y_der, double h, double bc_l,
    double bc_r) const
{
    std::complex<double> y_hat[N + C];
    Fcont_Gram_Blend_DD(y, y_hat, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle);
    FC_Der(y_der, y_hat, der_coeffs, N, C, desc_handle);
}

void FC_1D_DD::diff(const double* y, double * y_der, double* y_der_2, double h,
    double bc_l, double bc_r) const
{
    std::complex<double> y_hat[N + C];
    Fcont_Gram_Blend_DD(y, y_hat, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle);
    diff(y_hat, y_der, y_der_2);
}

void FC_1D_DD::filter(double* y, std::complex<double> *fft, double h,
    double bc_l, double bc_r) const
{
    std::complex<double> f_ext[N + C]; 
    Fcont_Gram_Blend_DD(y, fft, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle);
    VectorMulReCmp(N + C, filter_coeffs, fft, fft);
    std::copy(fft, fft + N + C, f_ext);
    int status = DftiComputeBackward(desc_handle, f_ext);
    for (int j = 0; j < N; j++)
    {
        y[j]= f_ext[j].real();
    }      
}

void FC_1D_DD::set_FC_Data(double* A, double* Q, int d, int C)
{
    if(C == 27)
    {
        if(d == 5)
        {
        std::copy(Ad5C27_data.data(), Ad5C27_data.data() + d*C, A);
        std::copy(Qd5C27_data.data(), Qd5C27_data.data() + d*d, Q);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// FC_1D_DN definitions

FC_1D_DN::FC_1D_DN(int _N, int _d, int _C, double _h) : 
    FC_1D{_N, _d, _C, _h}
{
    double A[C*d];
    double Q[d*d]; 
    double Q_tilde[d*d]; 
    set_FC_Data(A, Q, Q_tilde, d, C);
    build_Cont_Mat_DN(A, Q, Q_tilde, d, C, AQ, FAQF);  
}

FC_1D_DN::FC_1D_DN(int _N, int _d, int _C, double _h, int alpha, int p) : 
    FC_1D{_N, _d, _C, _h}
{
    double A[C*d];
    double Q[d*d]; 
    double Q_tilde[d*d]; 
    set_FC_Data(A, Q, Q_tilde, d, C);
    build_Cont_Mat_DN(A, Q, Q_tilde, d, C, AQ, FAQF);  
    getFiltCoeffs(filter_coeffs, fourPts, alpha, p); 
}

FC_1D_DN::FC_1D_DN(int _N, int _d, int _C, double _h, double delta) : 
    FC_1D{_N, _d, _C, _h}
{
    double A[C*d];
    double Q[d*d]; 
    double Q_tilde[d*d]; 
    set_FC_Data(A, Q, Q_tilde, d, C);
    build_Cont_Mat_DN(A, Q, Q_tilde, d, C, AQ, FAQF);  
    int k[fourPts];
    getK(k, fourPts);
    const double pi = std::acos(-1);
    const std::complex<double> I(0, 1);     
    for(int i = 0; i < fourPts; i++)
    {
        shift_coeffs[i] = std::exp(2.0*pi*I* double(k[i]) * delta / 
        double(fourPts));
    }
}

void FC_1D_DN::diff(const std::complex<double> *y_hat, double * y_der,
    double* y_der_2) const
{
    FC_Der(y_der, y_hat, der_coeffs, N, C, desc_handle);
    FC_Der(y_der_2, y_hat, der_coeffs_2, N, C, desc_handle, true);
}

void FC_1D_DN::diff(const double* y, double *y_der, double h, double bc_l,
    double bc_r) const
{
    std::complex<double> y_hat[N + C];
    Fcont_Gram_Blend_DN(y, y_hat, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle, 
        h, bc_r);
    FC_Der(y_der, y_hat, der_coeffs, N, C, desc_handle);
}

void FC_1D_DN::diff(const double* y, double * y_der, double* y_der_2, double h, 
    double bc_l, double bc_r) const
{
    std::complex<double> y_hat[N + C];
    Fcont_Gram_Blend_DN(y, y_hat, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle,
        h, bc_r);
    diff(y_hat, y_der, y_der_2);
}

void FC_1D_DN::filter(double* y, std::complex<double> *fft, double h,
    double bc_l, double bc_r) const
{
    std::complex<double> f_ext[N + C]; 
    Fcont_Gram_Blend_DN(y, fft, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle, 
        h, bc_r);
    VectorMulReCmp(N + C, filter_coeffs, fft, fft);
    std::copy(fft, fft + N + C, f_ext);
    int status = DftiComputeBackward(desc_handle, f_ext);
    for (int j = 0; j < N; j++)
    {
        y[j]= f_ext[j].real();
    }      
}

void FC_1D_DN::set_FC_Data(double* A, double* Q, double* Q_tilde, int d, int C)
{
    if(C == 27 && d == 5)
    {
        std::copy(Ad5C27_data.data(), Ad5C27_data.data() + d*C, A);
        std::copy(Qd5C27_data.data(), Qd5C27_data.data() + d*d, Q);
        std::copy(Qd5C27_tilde_data.data(), Qd5C27_tilde_data.data() + d*d,
        Q_tilde);
    }
}

////////////////////////////////////////////////////////////////////////////////
// FC_1D_ND definitions

FC_1D_ND::FC_1D_ND(int _N, int _d, int _C, double _h) : FC_1D{_N, _d, _C, _h}
{
    double A[C*d];
    double Q[d*d]; 
    double Q_tilde[d*d]; 
    set_FC_Data(A, Q, Q_tilde, d, C);
    build_Cont_Mat_ND(A, Q, Q_tilde, d, C, AQ, FAQF);  
}

FC_1D_ND::FC_1D_ND(int _N, int _d, int _C, double _h, int alpha, int p) : 
    FC_1D{_N, _d, _C, _h}
{
    double A[C*d];
    double Q[d*d]; 
    double Q_tilde[d*d]; 
    set_FC_Data(A, Q, Q_tilde, d, C); 
    build_Cont_Mat_ND(A, Q, Q_tilde, d, C, AQ, FAQF);  
    getFiltCoeffs(filter_coeffs, fourPts, alpha, p); 
}

FC_1D_ND::FC_1D_ND(int _N, int _d, int _C, double _h, double delta) : 
    FC_1D{_N, _d, _C, _h}
{
    double A[C*d];
    double Q[d*d]; 
    double Q_tilde[d*d]; 
    set_FC_Data(A, Q, Q_tilde, d, C);
    build_Cont_Mat_ND(A, Q, Q_tilde, d, C, AQ, FAQF);  
    int k[fourPts];
    getK(k, fourPts);
    const double pi = std::acos(-1);
    const std::complex<double> I(0, 1);     
    for(int i = 0; i < fourPts; i++)
    {
        shift_coeffs[i] = std::exp(2.0*pi*I* double(k[i]) * delta / 
        double(fourPts));
    }
}

void FC_1D_ND::diff(const std::complex<double> *y_hat, double * y_der,
    double* y_der_2) const
{
    FC_Der(y_der, y_hat, der_coeffs, N, C, desc_handle);
    FC_Der(y_der_2, y_hat, der_coeffs_2, N, C, desc_handle, true);
}

void FC_1D_ND::diff(const double* y, double *y_der, double h, double bc_l,
    double bc_r) const
{
    std::complex<double> y_hat[N + C];
    Fcont_Gram_Blend_ND(y, y_hat, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle,
        h, bc_l);
    FC_Der(y_der, y_hat, der_coeffs, N, C, desc_handle);
}

void FC_1D_ND::diff(const double* y, double * y_der, double* y_der_2, double h, 
    double bc_l, double bc_r) const
{
    std::complex<double> y_hat[N + C];
    Fcont_Gram_Blend_ND(y, y_hat, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle, 
        h, bc_l);
    diff(y_hat, y_der, y_der_2);
}

void FC_1D_ND::filter(double* y, std::complex<double> *fft, double h,
    double bc_l, double bc_r) const
{
    std::complex<double> f_ext[N + C]; 
    Fcont_Gram_Blend_ND(y, fft, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle, 
        h, bc_l);
    VectorMulReCmp(N + C, filter_coeffs, fft, fft);
    std::copy(fft, fft + N + C, f_ext);
    int status = DftiComputeBackward(desc_handle, f_ext);
    for (int j = 0; j < N; j++)
    {
        y[j]= f_ext[j].real();
    }      
}

void FC_1D_ND::set_FC_Data(double* A, double* Q, double* Q_tilde, int d, int C)
{
    if(C == 27 && d == 5)
    {
        std::copy(Ad5C27_data.data(), Ad5C27_data.data() + d*C, A);
        std::copy(Qd5C27_data.data(), Qd5C27_data.data() + d*d, Q);
        std::copy(Qd5C27_tilde_data.data(), Qd5C27_tilde_data.data() + d*d, 
        Q_tilde);
    }
}


////////////////////////////////////////////////////////////////////////////////
// FC_1D_NN definitions

FC_1D_NN::FC_1D_NN(int _N, int _d, int _C, double _h) : 
    FC_1D{_N, _d, _C, _h}
{
    double A[C*d];
    double Q[d*d]; 
    set_FC_Data(A, Q, d, C);
    build_Cont_Mat(A, Q, d, C, AQ, FAQF);  
}

FC_1D_NN::FC_1D_NN(int _N, int _d, int _C, double _h, int alpha, int p) : 
    FC_1D{_N, _d, _C, _h}
{
    double A[C*d];
    double Q[d*d]; 
    set_FC_Data(A, Q, d, C);
    build_Cont_Mat(A, Q, d, C, AQ, FAQF);  
    getFiltCoeffs(filter_coeffs, fourPts, alpha, p); 
}

FC_1D_NN::FC_1D_NN(int _N, int _d, int _C, double _h, double delta) : 
    FC_1D{_N, _d, _C, _h}
{
    double A[C*d];
    double Q[d*d]; 
    set_FC_Data(A, Q, d, C);
    build_Cont_Mat(A, Q, d, C, AQ, FAQF); 
    int k[fourPts];
    getK(k, fourPts);
    const double pi = std::acos(-1);
    const std::complex<double> I(0, 1);     
    for(int i = 0; i < fourPts; i++)
    {
        shift_coeffs[i] = std::exp(2.0*pi*I* double(k[i]) * delta / 
            double(fourPts));
    }
}

void FC_1D_NN::diff(const std::complex<double> *y_hat, double * y_der,
    double* y_der_2) const
{
    FC_Der(y_der, y_hat, der_coeffs, N, C, desc_handle);
    FC_Der(y_der_2, y_hat, der_coeffs_2, N, C, desc_handle, true);
}

void FC_1D_NN::diff(const double* y, double *y_der, double h, double bc_l,
    double bc_r) const
{
    std::complex<double> y_hat[N + C];
    Fcont_Gram_Blend_NN(y, y_hat, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle, 
        h, bc_l, bc_r);
    FC_Der(y_der, y_hat, der_coeffs, N, C, desc_handle);
}

void FC_1D_NN::diff(const double* y, double * y_der, double* y_der_2, double h, 
    double bc_l, double bc_r) const
{
    std::complex<double> y_hat[N + C];
    Fcont_Gram_Blend_NN(y, y_hat, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle, 
        h, bc_l, bc_r);
    diff(y_hat, y_der, y_der_2);
}

void FC_1D_NN::filter(double* y, std::complex<double> *fft, double h,
    double bc_l, double bc_r) const
{
    std::complex<double> f_ext[N + C]; 
    Fcont_Gram_Blend_NN(y, fft, N, d, C, fourPts_dbl, AQ, FAQF, desc_handle, 
        h, bc_l, bc_r);
    VectorMulReCmp(N + C, filter_coeffs, fft, fft);
    std::copy(fft, fft + N + C, f_ext);
    int status = DftiComputeBackward(desc_handle, f_ext);
    for (int j = 0; j < N; j++)
    {
        y[j]= f_ext[j].real();
    }      
}

void FC_1D_NN::set_FC_Data(double* A, double* Q, int d, int C)
{
    if(C == 27)
    {
        if(d == 5)
        {
            std::copy(Ad5C27_data.data(), Ad5C27_data.data() + d*C, A);
            std::copy(Qd5C27_tilde_data.data(), Qd5C27_tilde_data.data() + d*d,
                Q);
        }
    }
}

