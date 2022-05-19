/* Some solver routines. */

#include "Solver.h"
#include "Mesh.h"
#include <string.h>
#include <vector>
#include "FC_1D.h"
#include "FC_2D.h"


void initialize1DFCSDNN(std::vector<FC_1D* > * diff_schemes, 
    std::vector<FC_1D* > * filters, const Mesh1DUniform &mesh, int npatches,
    std::string xx, int N, int d, int C, double delta, double alpha0,
    double p_0, double p)
{
    if(xx.compare("DD") == 0)
    {
        for(int i = 0; i < npatches; i++)
        {
            diff_schemes->push_back(new FC_1D_DD(N, d, C,
                mesh.getPatches()[i]->getH(), delta));
            filters->push_back(new FC_1D_DD(N, d, C,
                mesh.getPatches()[i]->getH(), alpha0, p_0));
        }
        for(int i = 0; i < npatches; i++)
        {
            filters->push_back(new FC_1D_DD (N, d, C,
                mesh.getPatches()[i]->getH(), alpha0, p));
        }
    }
    else if(xx.compare("ND") == 0)
    {

    }
    else if(xx.compare("DN") == 0)
    {

    }
    else if(xx.compare("NN") == 0)
    {
        if (npatches == 1)
        {
            diff_schemes->push_back(new FC_1D_NN(N, d, C,
                mesh.getPatches()[0]->getH(), delta));
            filters->push_back(new FC_1D_NN(N, d, C,
                mesh.getPatches()[0]->getH(), alpha0, p_0));
            filters->push_back(new FC_1D_NN(N, d, C,
                mesh.getPatches()[0]->getH(), alpha0, p));
        }
        else
        {
            diff_schemes->push_back(new FC_1D_ND(N, d, C,
                mesh.getPatches()[0]->getH(), delta));
            filters->push_back(new FC_1D_ND(N, d, C,
                mesh.getPatches()[0]->getH(), alpha0, p_0));
            for(int i = 1; i < npatches - 1; i++)
            {
                diff_schemes->push_back(new FC_1D_DD(N, d, C,
                    mesh.getPatches()[i]->getH(), delta));
                filters->push_back(new FC_1D_DD(N, d, C,
                    mesh.getPatches()[i]->getH(), alpha0, p_0));
            }
            diff_schemes->push_back(new FC_1D_DN(N, d, C, 
                mesh.getPatches()[npatches - 1]->getH(), delta));
            filters->push_back(new FC_1D_DN(N, d, C,
                mesh.getPatches()[0]->getH(), alpha0, p_0));
            filters->push_back(new FC_1D_ND(N, d, C,
                mesh.getPatches()[0]->getH(), alpha0, p));
            for(int i = 1; i < npatches - 1; i++)
            {
                filters->push_back(new FC_1D_DD(N, d, C,
                    mesh.getPatches()[i]->getH(), alpha0, p));
            }
            filters->push_back(new FC_1D_DN(N, d, C,
                mesh.getPatches()[0]->getH(), alpha0, p));
        }
    }
}


void initialize2DFCSDNN(std::vector<FC_2D* > * diff_schemes, 
    std::vector<FC_2D* > * filters, const Mesh2DUniform &mesh, 
    int npatches_x, int npatches_y, std::string left, std::string right,
    std::string down, std::string up, int Nx, int Ny, int dx, int dy, int Cx,
    int Cy, double delta, double alpha0, double p_0, double p)
{
    if(npatches_y == 1)
    {
        if(npatches_x == 1)
        {
            diff_schemes->push_back(new FC_2D(left + right, down + up, Nx, dx,
                Cx, mesh.getPatches()[0]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[0]->getHy(), delta));
            filters->push_back(new FC_2D(left + right, down + up, Nx, dx, Cx,
                mesh.getPatches()[0]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[0]->getHy(), alpha0, p_0));
            filters->push_back(new FC_2D(left + right, down + up, Nx, dx, Cx,
                mesh.getPatches()[0]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[0]->getHy(), alpha0, p));
        }
        else
        {
            diff_schemes->push_back(new FC_2D(left + "D", down + up, Nx, dx, Cx,
                mesh.getPatches()[0]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[0]->getHy(), delta));
            filters->push_back(new FC_2D(left + "D", down + up, Nx, dx, Cx,
                mesh.getPatches()[0]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[0]->getHy(), alpha0, p_0));           
            for(int i = 1; i < npatches_x - 1; i++)
            {
                diff_schemes->push_back(new FC_2D("DD", down + up, Nx, dx, Cx,
                    mesh.getPatches()[i]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[i]->getHy(), delta));
                filters->push_back(new FC_2D("DD", down + up, Nx, dx, Cx,
                    mesh.getPatches()[i]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[i]->getHy(), alpha0, p_0));
            }  
            diff_schemes->push_back(new FC_2D("D" + right, down + up, Nx, dx,
                Cx, mesh.getPatches()[npatches_x - 1]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[npatches_x - 1]->getHy(), delta));
            filters->push_back(new FC_2D("D" + right, down + up, Nx, dx, Cx,
                mesh.getPatches()[npatches_x - 1]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[npatches_x - 1]->getHy(), alpha0, p_0));  

            filters->push_back(new FC_2D(left + "D", down + up, Nx, dx, Cx,
                mesh.getPatches()[0]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[0]->getHy(), alpha0, p));  
            for(int i = 1; i < npatches_x - 1; i++)
            {
                filters->push_back(new FC_2D("DD", down + up, Nx, dx,
                    Cx, mesh.getPatches()[i]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[i]->getHy(), alpha0, p));                  
            }   
            filters->push_back(new FC_2D("D" + right, down + up, Nx, dx, Cx,
                mesh.getPatches()[npatches_x - 1]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[npatches_x - 1]->getHy(), alpha0, p));                                                    
        }
    }
    else
    {
        if(npatches_x == 1)
        {
            diff_schemes->push_back(new FC_2D(left + right, down + "D", Nx, dx, Cx,
                mesh.getPatches()[0]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[0]->getHy(), delta));
            filters->push_back(new FC_2D(left + right, down + "D", Nx, dx, Cx,
                mesh.getPatches()[0]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[0]->getHy(), alpha0, p_0));           
            for(int i = 1; i < npatches_y - 1; i++)
            {
                diff_schemes->push_back(new FC_2D(left + right, "DD", Nx, dx, Cx,
                    mesh.getPatches()[i]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[i]->getHy(), delta));
                filters->push_back(new FC_2D(left + right, "DD", Nx, dx, Cx,
                    mesh.getPatches()[i]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[i]->getHy(), alpha0, p_0));
            }  
            diff_schemes->push_back(new FC_2D(left + right, "D" + up, Nx, dx,
                Cx, mesh.getPatches()[npatches_x - 1]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[npatches_x - 1]->getHy(), delta));
            filters->push_back(new FC_2D(left + right, "D" + up, Nx, dx, Cx,
                mesh.getPatches()[npatches_x - 1]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[npatches_x - 1]->getHy(), alpha0, p_0));  

            filters->push_back(new FC_2D(left + right, down + "D", Nx, dx, Cx,
                mesh.getPatches()[0]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[0]->getHy(), alpha0, p));  
            for(int i = 1; i < npatches_y - 1; i++)
            {
                filters->push_back(new FC_2D(left + right, "DD", Nx, dx,
                    Cx, mesh.getPatches()[i]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[i]->getHy(), alpha0, p));                  
            }   
            filters->push_back(new FC_2D(left + right, "D" + up, Nx, dx, Cx,
                mesh.getPatches()[npatches_x - 1]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[npatches_x - 1]->getHy(), alpha0, p));                                                    
        }
        else
        {
            // 1st column
            diff_schemes->push_back(new FC_2D(left + "D", down + "D", Nx, dx, Cx,
                mesh.getPatches()[0]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[0]->getHy(), delta));
            filters->push_back(new FC_2D(left + "D", down + "D", Nx, dx, Cx,
                mesh.getPatches()[0]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[0]->getHy(), alpha0, p_0));  
            for(int i = 1; i < npatches_y; i++)
            {
                diff_schemes->push_back(new FC_2D(left + "D", "DD", Nx, dx, Cx,
                    mesh.getPatches()[i]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[i]->getHy(), delta));
                filters->push_back(new FC_2D(left + "D", "DD", Nx, dx, Cx,
                    mesh.getPatches()[i]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[i]->getHy(), alpha0, p_0));
            }
            diff_schemes->push_back(new FC_2D(left + "D", "D" + up, Nx, dx,
                Cx, mesh.getPatches()[npatches_x - 1]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[npatches_x - 1]->getHy(), delta));
            filters->push_back(new FC_2D(left + "D", "D" + up, Nx, dx, Cx,
                mesh.getPatches()[npatches_x - 1]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[npatches_x - 1]->getHy(), alpha0, p_0));   

            // center columns
            for(int j = 1; j < npatches_x - 1; j++)
            {
                diff_schemes->push_back(new FC_2D("DD", down + "D", Nx, dx, Cx,
                    mesh.getPatches()[j*npatches_y]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[j*npatches_y]->getHy(), delta));
                filters->push_back(new FC_2D("DD", down + "D", Nx, dx, Cx,
                    mesh.getPatches()[j*npatches_y]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[j*npatches_y]->getHy(), alpha0, p_0));  
                for(int i = 1; i < npatches_y; i++)
                {
                    diff_schemes->push_back(new FC_2D("DD", "DD", Nx, dx, Cx,
                        mesh.getPatches()[j*npatches_y + i]->getHx(), Ny, dy, Cy, 
                        mesh.getPatches()[j*npatches_y + i]->getHy(), delta));
                    filters->push_back(new FC_2D("DD", "DD", Nx, dx, Cx,
                        mesh.getPatches()[j*npatches_y + i]->getHx(), Ny, dy, Cy, 
                        mesh.getPatches()[j*npatches_y + i]->getHy(), alpha0, p_0));
                }
                diff_schemes->push_back(new FC_2D("DD", "D" + up, Nx, dx,
                    Cx, mesh.getPatches()[j*npatches_y + npatches_x - 1]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[j*npatches_y + npatches_x - 1]->getHy(), delta));
                filters->push_back(new FC_2D("DD", "D" + up, Nx, dx, Cx,
                    mesh.getPatches()[j*npatches_y + npatches_x - 1]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[j*npatches_y + npatches_x - 1]->getHy(), alpha0, p_0));                 
            }      

            // Last column     
            diff_schemes->push_back(new FC_2D("D" + right, down + "D", Nx, dx, Cx,
                mesh.getPatches()[(npatches_x - 1)*npatches_y]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[(npatches_x - 1)*npatches_y]->getHy(), delta));
            filters->push_back(new FC_2D("D" + right, down + "D", Nx, dx, Cx,
                mesh.getPatches()[(npatches_x - 1)*npatches_y]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[(npatches_x - 1)*npatches_y]->getHy(), alpha0, p_0));  
            for(int i = 1; i < npatches_y; i++)
            {
                diff_schemes->push_back(new FC_2D("D" + right, "DD", Nx, dx, Cx,
                    mesh.getPatches()[(npatches_x - 1)*npatches_y + i]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[(npatches_x - 1)*npatches_y + i]->getHy(), delta));
                filters->push_back(new FC_2D("D" + right, "DD", Nx, dx, Cx,
                    mesh.getPatches()[(npatches_x - 1)*npatches_y + i]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[(npatches_x - 1)*npatches_y + i]->getHy(), alpha0, p_0));
            }
            diff_schemes->push_back(new FC_2D("D" + right, "D" + up, Nx, dx,
                Cx, mesh.getPatches()[npatches_y*npatches_x - 1]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[npatches_y*npatches_x - 1]->getHy(), delta));
            filters->push_back(new FC_2D("D" + right, "D" + up, Nx, dx, Cx,
                mesh.getPatches()[npatches_y*npatches_x - 1]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[npatches_y*npatches_x - 1]->getHy(), alpha0, p_0));      


            // And now the same thing, but for the filters t > 0    
            // 1st column
            filters->push_back(new FC_2D(left + "D", down + "D", Nx, dx, Cx,
                mesh.getPatches()[0]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[0]->getHy(), alpha0, p));  
            for(int i = 1; i < npatches_y; i++)
            {
                filters->push_back(new FC_2D(left + "D", "DD", Nx, dx, Cx,
                    mesh.getPatches()[i]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[i]->getHy(), alpha0, p));
            }
            filters->push_back(new FC_2D(left + "D", "D" + up, Nx, dx, Cx,
                mesh.getPatches()[npatches_x - 1]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[npatches_x - 1]->getHy(), alpha0, p));   

            // center columns
            for(int j = 1; j < npatches_x - 1; j++)
            {
                filters->push_back(new FC_2D("DD", down + "D", Nx, dx, Cx,
                    mesh.getPatches()[j*npatches_y]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[j*npatches_y]->getHy(), alpha0, p));  
                for(int i = 1; i < npatches_y; i++)
                {
                    filters->push_back(new FC_2D("DD", "DD", Nx, dx, Cx,
                        mesh.getPatches()[j*npatches_y + i]->getHx(), Ny, dy, Cy, 
                        mesh.getPatches()[j*npatches_y + i]->getHy(), alpha0, p));
                }
                filters->push_back(new FC_2D("DD", "D" + up, Nx, dx, Cx,
                    mesh.getPatches()[j*npatches_y + npatches_x - 1]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[j*npatches_y + npatches_x - 1]->getHy(), alpha0, p));                 
            }      

            // Last column     
            filters->push_back(new FC_2D("D" + right, down + "D", Nx, dx, Cx,
                mesh.getPatches()[(npatches_x - 1)*npatches_y]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[(npatches_x - 1)*npatches_y]->getHy(), alpha0, p));  
            for(int i = 1; i < npatches_y; i++)
            {
                filters->push_back(new FC_2D("D" + right, "DD", Nx, dx, Cx,
                    mesh.getPatches()[(npatches_x - 1)*npatches_y + i]->getHx(), Ny, dy, Cy, 
                    mesh.getPatches()[(npatches_x - 1)*npatches_y + i]->getHy(), alpha0, p));
            }
            filters->push_back(new FC_2D("D" + right, "D" + up, Nx, dx, Cx,
                mesh.getPatches()[npatches_y*npatches_x - 1]->getHx(), Ny, dy, Cy, 
                mesh.getPatches()[npatches_y*npatches_x - 1]->getHy(), alpha0, p));                   
        }

    }
}

