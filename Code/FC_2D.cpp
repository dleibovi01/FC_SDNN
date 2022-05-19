/* 2D FC spatial differentiation schemes */

#include "FC_2D.h"



FC_2D::FC_2D(const std::string xx, const std::string yy, int _Nx, int dx,
    int Cx, double hx, int _Ny, int dy, int Cy, double hy, int alpha, int p) :
    Nx{_Nx}, Ny{_Ny}
{
    if(xx.compare("DD") == 0)
        FC_x = new FC_1D_DD(Nx, dx, Cx, hx, alpha, p);
    else if(xx.compare("DN") == 0)
        FC_x = new FC_1D_DN(Nx, dx, Cx, hx, alpha, p);
    else if(xx.compare("ND") == 0)
        FC_x = new FC_1D_ND(Nx, dx, Cx, hx, alpha, p);
    else if(xx.compare("NN") == 0)
        FC_x = new FC_1D_NN(Nx, dx, Cx, hx, alpha, p);

    if(yy.compare("DD") == 0)
        FC_y = new FC_1D_DD(Ny, dy, Cy, hy, alpha, p);
    else if(yy.compare("DN") == 0)
        FC_y = new FC_1D_DN(Ny, dy, Cy, hy, alpha, p);
    else if(yy.compare("ND") == 0)
        FC_y = new FC_1D_ND(Ny, dy, Cy, hy, alpha, p);
    else if(yy.compare("NN") == 0)
        FC_y = new FC_1D_NN(Ny, dy, Cy, hy, alpha, p);    
}

FC_2D::FC_2D(const std::string xx, const std::string yy, int _Nx, int dx,
    int Cx, double hx, int _Ny, int dy, int Cy, double hy, double delta) :
    Nx{_Nx}, Ny{_Ny}
{
    if(xx.compare("DD") == 0)
        FC_x = new FC_1D_DD(Nx, dx, Cx, hx, delta);
    else if(xx.compare("DN") == 0)
        FC_x = new FC_1D_DN(Nx, dx, Cx, hx, delta);
    else if(xx.compare("ND") == 0)
        FC_x = new FC_1D_ND(Nx, dx, Cx, hx, delta);
    else if(xx.compare("NN") == 0)
        FC_x = new FC_1D_NN(Nx, dx, Cx, hx, delta);

    if(yy.compare("DD") == 0)
        FC_y = new FC_1D_DD(Ny, dy, Cy, hy, delta);
    else if(yy.compare("DN") == 0)
        FC_y = new FC_1D_DN(Ny, dy, Cy, hy, delta);
    else if(yy.compare("ND") == 0)
        FC_y = new FC_1D_ND(Ny, dy, Cy, hy, delta);
    else if(yy.compare("NN") == 0)
        FC_y = new FC_1D_NN(Ny, dy, Cy, hy, delta);      
}

FC_2D::~FC_2D()
{
    delete FC_x;
    delete FC_y;
}

void FC_2D::diff_y(const std::complex<double> *y_hat, double * y_der,
    double* y_der_2) const
{
    for(int j = 0; j < Nx; j++)
    {
        FC_y->diff(y_hat + j*Ny, y_der + j*Ny, y_der_2 + j*Ny);
    }    
}

void FC_2D::diff_y(const double* y, double *y_der, double h, const double* bc_d,
    const double* bc_u) const
{
    for(int j = 0; j < Nx; j++)
    {
        FC_y->diff(y + j*Ny, y_der + j*Ny, h, bc_d[j], bc_u[j]);
    }
}

void FC_2D::diff_y(const double* y, double * y_der, double* y_der_2, double h, 
    const double* bc_d, const double* bc_u) const
{
    for(int j = 0; j < Nx; j++)
    {
        FC_y->diff(y + j*Ny, y_der + j*Ny, y_der_2 + j*Ny, h, bc_d[j], bc_u[j]);
    }    
}

void FC_2D::filter_y(VectorField *v, const std::vector<int> &unknowns, 
    std::vector<std::complex<double> *> *ffts, 
    const std::vector<int> &fft_loc, double h, 
    const std::vector<double* > &bc_d, const std::vector<double* > &bc_u) const
{
    for(int i = 0; i < unknowns.size(); i++)
    {
        for(int j = 0; j < Nx; j++)
        {
            FC_y->filter(v->getField(unknowns[i]) + j*Ny, 
                ffts->at(fft_loc[i]) + j*Ny, h, bc_d[i][j], bc_u[i][j]);   
        }
    } 
}

void FC_2D::filter_y(VectorField *v, const std::vector<int> &unknowns, double h, 
    const std::vector<double* > &bc_d, const std::vector<double* > &bc_u) const
{
    for(int i = 0; i < unknowns.size(); i++)
    {
        for(int j = 0; j < Nx; j++)
        { 
            FC_y->filter(v->getField(unknowns[i]) + j*Ny, h, bc_d[i][j],
                bc_u[i][j]);   
        }
    } 
}

void FC_2D::diff_x(const std::complex<double> *y_hat, double * y_der,
    double* y_der_2) const
{
    for(int j = 0; j < Ny; j++)
    {
        FC_x->diff(y_hat + j*Nx, y_der + j*Nx, y_der_2 + j*Nx);
    }   
    mkl_dimatcopy('C', 'T', Nx, Ny, 1.0, y_der, Nx, Ny);
    mkl_dimatcopy('C', 'T', Nx, Ny, 1.0, y_der_2, Nx, Ny); 
}

void FC_2D::diff_x(const double* y, double *y_der, double h, const double* bc_d,
    const double* bc_u) const
{
    double yy[Nx*Ny];
    std::copy(y, y + Nx*Ny, yy);
    mkl_dimatcopy('C', 'T', Ny, Nx, 1.0, yy, Ny, Nx);
    for(int j = 0; j < Ny; j++)
    {
        FC_x->diff(yy + j*Nx, y_der + j*Nx, h, bc_d[j], bc_u[j]);
    }
    mkl_dimatcopy('C', 'T', Nx, Ny, 1.0, y_der, Ny, Nx);
}

void FC_2D::diff_x(const double* y, double * y_der, double* y_der_2, double h, 
    const double* bc_d, const double* bc_u) const
{
    double yy[Nx*Ny];
    std::copy(y, y + Nx*Ny, yy);
    mkl_dimatcopy('C', 'T', Ny, Nx, 1.0, yy, Ny, Nx);
    for(int j = 0; j < Ny; j++)
    {
        FC_x->diff(yy + j*Nx, y_der + j*Nx, y_der_2 + j*Nx, h, bc_d[j], bc_u[j]);
    }    
    mkl_dimatcopy('C', 'T', Nx, Ny, 1.0, y_der, Nx, Ny);
    mkl_dimatcopy('C', 'T', Nx, Ny, 1.0, y_der_2, Nx, Ny); 
}

void FC_2D::filter_x(VectorField *v, const std::vector<int> &unknowns, 
    std::vector<std::complex<double> *> *ffts, 
    const std::vector<int> &fft_loc, double h, 
    const std::vector<double* > &bc_d, const std::vector<double* > &bc_u) const
{
    for(int i = 0; i < unknowns.size(); i++)
    {
        mkl_dimatcopy('C', 'T', Ny, Nx, 1.0, v->getField(unknowns[i]), Ny, Nx);
        for(int j = 0; j < Ny; j++)
        {
            FC_x->filter(v->getField(unknowns[i]) + j*Nx, 
                ffts->at(fft_loc[i]) + j*Nx, h, bc_d[i][j], bc_u[i][j]);   
        }
        mkl_dimatcopy('C', 'T', Nx, Ny, 1.0, v->getField(unknowns[i]), Nx, Ny);
    } 
}

void FC_2D::filter_x(VectorField *v, const std::vector<int> &unknowns, double h, 
    const std::vector<double* > &bc_d, const std::vector<double* > &bc_u) const
{
    for(int i = 0; i < unknowns.size(); i++)
    {
        mkl_dimatcopy('C', 'T', Ny, Nx, 1.0, v->getField(unknowns[i]), Ny, Nx);
        for(int j = 0; j < Ny; j++)
        {
            FC_x->filter(v->getField(unknowns[i]) + j*Nx, h, bc_d[i][j],
                bc_u[i][j]);   
        }
        mkl_dimatcopy('C', 'T', Nx, Ny, 1.0, v->getField(unknowns[i]), Nx, Ny);
    } 
}

void FC_2D::filter(VectorField *v, const std::vector<int> &unknowns, double hx,
    double hy, const std::vector<double* > &bc_l,
    const std::vector<double* > &bc_r, const std::vector<double* > &bc_d, 
    const std::vector<double* > &bc_u) const
{
    filter_x(v, unknowns, hx, bc_l, bc_r);  
    filter_y(v, unknowns, hy, bc_d, bc_u);
}
